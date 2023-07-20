module temporal_filter
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
    include 'ctes3D'

    ! NOTE:
    !   FilterCoef in the filter setting file has to be a columns vector (NX1)!!!!!

    ! NOTE:
    !   The results read from the restart file will be designated as step 0 and time 0
    !   Loop ii will take results from step (ii-1), time (ii-1)DT and propagate it to step ii, time ii*DT
    !   The filtered result at step ii depends on the unfiltered result at step (ii-FilterOrder):ii,
    !   a total of FilterOrder+1 steps
    !
    !   For the old DNS code, the result saved at step N is the data at the end of step N-1, time (N-1)DT
    !   The time agrees with the new DNS data, but the step number is off by 1
    !
    !   If for some reason, no filtering is needed, then simply put:
    !      FilterOrder = 0, FilterCoef = 1
    !   in the filter settings h5 file
    !   The module will not increase computation cost and memory usage by much (if any)


    ! ---------------------------- Global variables ----------------------------
    ! For kz and y parallel assignments
    integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
    common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                    kbeg(0:numerop-1),kend(0:numerop-1), &
                    jb,je,kb,ke,mmy,mmz
    save /point/
    ! -------------- Filter Variables accessable in entire module --------------
    ! Filter order and coefficients
    integer :: FilterOrder
    real(kind=dp), dimension(:), allocatable :: FilterCoef
    ! Number of snapshots to store in buffer
    integer :: BufferDim4Size
    ! Buffer
    complex(kind=dp), dimension(:,:,:,:), allocatable :: vBuffer, vorBuffer
    real(kind=dp), dimension(:,:), allocatable :: u00Buffer, w00Buffer ! Note that 00 mode filtering is only done by master
    ! The step number of each element stored in dim 4 of the buffers (dim 2 in case of 00 mode buffers)
    integer, dimension(:), allocatable :: BufferStepNum
    ! Suppress this many snapshot saving due to uninitialized or partially initialized filter
    integer :: SuppressSavingNum
    ! ----------- File Saving Variables accessable in entire module -----------
    ! Snapshot FileNumber
    integer :: FileNumber
    ! File name
    character(:), allocatable :: BaseFileNameWithPath
    ! Sampling frequency
    integer :: SampleFreqInSteps
    ! -------------------- DNS Parameters to be saved in h5 --------------------
    real   (kind=sp) :: fundamentalkx, fundamentalkz
    real   (kind=sp) :: y(my)
    complex(kind=sp) :: kx(mx/2), kz(mz)
    real   (kind=sp) :: Re_ref
    real   (kind=dp) :: bulkVelocity

    ! processsor ID
    integer :: myid


    public :: initialize_temporal_filter_moudule, initialize_temporal_filter_DNS_parameters
    public :: cleanup_temporal_filter_module
    public :: update_filter, save_filter_restart_file


contains
! ****************** Initialization and Cleanup of the Module ******************
    ! function initialize_temporal_filter_moudule( myid_Input ) result( InitSuccess )
    ! Initializes the moudule, all processors call this during program initialization
    !
    ! CALLED IN: MAIN
    !
    ! Arguments
    !   myid_Input: [int, Input] Processor ID
    ! Output
    !   InitSuccess: [logical]
    !                .true. if sampling freq setting agree between filter settings file and DNS settings
    !                .false. otherwise
    function initialize_temporal_filter_moudule( myid_Input ) result( InitSuccess )
        ! -------------------------- Global variables --------------------------
        ! for file name and file number
        integer iinp,iout,id22,isn,ispf
        character*100 filinp,filout,filstt,filfilter
        common /ficheros/ iinp,iout,id22,isn,ispf,filinp,filout,filstt,filfilter
        save /ficheros/
        ! For DNS time
        real*4 Deltat,CFL,time,dtr,FixTimeStep
        common /tem/ Deltat,CFL,time,dtr,FixTimeStep
        save /tem/
        ! For step number stuff
        integer nimag,nstep,nhist,ihist,icfl,nsnapshot
        common /timacc/ nimag,nstep,nhist,ihist,icfl,nsnapshot
        save   /timacc/
        ! ----------------------------- Arguments -----------------------------
        ! processor ID
        integer, intent(in) :: myid_Input
        ! output
        logical :: InitSuccess
        ! ------------------------------- Others -------------------------------
        integer :: ii


        ! Default to successful initialization unles filter settings disagree
        InitSuccess = .true.
        ! Set myid as global variable in this module
        myid = myid_Input
        ! Set Sampling frequency
        SampleFreqInSteps = nsnapshot
        ! FileName
        BaseFileNameWithPath = trim(adjustl(filstt)) // '_flowfield.'
        ! Initialize Snapshot FileNumber to 1
        FileNumber = 1
        ! Write sampling information
        if (myid == 0) then
            write(*,'(A,I4)') "Sampling frequency of flowfield in timesteps: ", SampleFreqInSteps
        endif

        ! Read filter settings
        if (read_filter_settings( filfilter ) .eqv. .false.) then
            InitSuccess = .false.
            return
        endif
        ! Number of snapshots to keep in buffer
        ! EX: FilterOrder=1000,SampleFreqInSteps=50, BufferDim4Size=21
        !     At step 50, snapshot 1 (step 50) to snapshot 21 (step 1050) all depends on step 50
        BufferDim4Size = FLOOR( real(FilterOrder,dp)/real(SampleFreqInSteps,dp)) + 1
        ! Allocate filter buffer
        ALLOCATE(  vBuffer(my, mx/2, kb:ke, BufferDim4Size))
        ALLOCATE(vorBuffer(my, mx/2, kb:ke, BufferDim4Size))
        ALLOCATE( BufferStepNum( BufferDim4Size ) )
        if (myid .eq. 0) then ! 00 mode filtering is only done by master
            ALLOCATE(u00Buffer(my, BufferDim4Size))
            ALLOCATE(w00Buffer(my, BufferDim4Size))
        endif

        ! Compare current filter settings with settings from the restart file
        if (compare_filter_settings(filinp) .eqv. .true.) then
            ! if settings agree, then read the filterbuffer from the restart file
            call read_filter_restart(filinp)
        else
            ! Initialize buffer to 0
              vBuffer = 0.0_dp
            vorBuffer = 0.0_dp
            if (myid .eq. 0) then ! 00 mode filtering is only done by master
                u00Buffer = 0.0_dp
                w00Buffer = 0.0_dp
            endif
            ! Suppress snapshot saving until filter is initialized
            if ( FilterOrder .lt. SampleFreqInSteps ) then
                ! if SampleFreqInSteps - FilterOrder >= 1
                ! then the first snapshot will have all steps needed for the filtering
                ! Therefore no saving suppression
                SuppressSavingNum = 0
            else
                ! The last snapshot in the buffer will have its first dependence within this run,
                ! therefore it is the first snapshot that can be saved
                SuppressSavingNum = BufferDim4Size - 1
            endif
            ! The buffer will be organized such that Dim4 index N is the buffer for the N-th snapshot
            DO ii = 1, BufferDim4Size
                BufferStepNum(ii) = ii*SampleFreqInSteps
            ENDDO
            if (myid .eq. 0) then
                write(*,*)
                write(*,"(A)") "Initializing filter buffer with 0"
                write(*,"(A,I5,A)") "First ", SuppressSavingNum, " snapshots will not be saved."
            endif
        endif
    end function initialize_temporal_filter_moudule

    ! subroutine initialize_temporal_filter_DNS_parameters
    !   ( fundamentalkx_in, fundamentalkz_in, y_in, kx_in, kz_in, Re_ref_in, bulkVelocity_in )
    ! Set DNS parameters as global variables so that these parameters can be saved into h5 snapshot files during run time
    !
    ! CALLED IN: CROSS
    !
    ! Arguments
    !   fundamentalkx_in, fundamentalkz_in: [single real, Input]
    !                                       fundamental wavenumbers in the streamwise and spanwise directions
    !   y_in                              : [singele real, size my, Input]
    !                                       wall normal coordinate vector
    !   kx_in                             : [single complex, size mx/2, Input]
    !                                       streamwise wavenumber vector
    !   kz_in                             : [single complex, size mz, Input]
    !                                       spanwise wavenumber vector
    !   Re_ref_in                         : [real single, Input]
    !                                       Reference reynolds number
    !                                       (based on the reference velocity, which is approximately, but not equal to the centerline velocity)
    !   bulkVelocity_in                   : [real double, Input]
    !                                       bulkVelocity
    subroutine initialize_temporal_filter_DNS_parameters(  &
        fundamentalkx_in, fundamentalkz_in, y_in, kx_in, kz_in, Re_ref_in, bulkVelocity_in )
        real   (kind=sp), intent(in) :: fundamentalkx_in, fundamentalkz_in
        real   (kind=sp), intent(in) :: y_in(:)
        complex(kind=sp), intent(in) :: kx_in(:), kz_in(:)
        real   (kind=sp), intent(in) :: Re_ref_in
        real   (kind=dp), intent(in) :: bulkVelocity_in

        ! Set all as global variables for this module
        fundamentalkx = fundamentalkx_in
        fundamentalkz = fundamentalkz_in
        y = y_in
        kx = kx_in
        kz = kz_in
        Re_ref = Re_ref_in
        bulkVelocity = bulkVelocity_in
    end

    ! subroutine cleanup_temporal_filter_module( )
    ! Deallocate memory, should be called by all processors
    !
    ! CALLED IN: MAIN
    subroutine cleanup_temporal_filter_module( )
        DEALLOCATE(   vBuffer  )
        DEALLOCATE( vorBuffer  )
        DEALLOCATE( BufferStepNum )
        DEALLOCATE( FilterCoef )
        if (myid .eq. 0) then ! 00 mode filtering is only done by master
            DEALLOCATE( u00Buffer  )
            DEALLOCATE( w00Buffer  )
        endif
    end subroutine cleanup_temporal_filter_module

! *************************** Runtime Filter Update ***************************
    ! subroutine update_filter( v, vor, u00, w00, istep, time )
    ! Update the filter buffer with this current times step, and save snapshot if this is a step for saving
    !
    ! CALLED IN: CROSS
    !
    ! Arguements
    !   v, vor  : [single real, size (2*my,mx/2,kb:ke), Input]
    !             v, omega_2 fields at the end of the current time step, the first dim contains both real and imag
    !   u00, w00: [double, size my, Input]
    !             kx = kz = 0 mode of u and w
    !   istep   : [int, Input]
    !             step number of the time loop
    !   time    : [single, Input]
    !             time of the current time loop, should be istep*DeltaT
    subroutine update_filter( v, vor, u00, w00, istep, time )
        ! ----------------------------- Arguments -----------------------------
        ! v, omega_3
        ! Note: all real where the first dimension contains both the real and imag parts
        real(kind=sp), intent(in), dimension(2*my,mx/2,kb:ke) :: v, vor
        ! kx = kz = 0 modes for u and w
        real(kind=dp), intent(in), dimension(my) :: u00, w00
        ! current step number in the time loop
        integer      , intent(in) :: istep
        ! current time of the time loop, should be istep*DeltaT
        real(kind=sp), intent(in) :: time
        ! ------------------------------- Others -------------------------------
        ! Loop variable
        integer :: ii, FilterCoefIndex


        ! loop over all the snapshots stored in the buffer and update each
        DO ii = 1, BufferDim4Size
            ! when istep = BufferStepNum(ii)
            !   we should be using the last coef in FilterCoef -> FilterCoefIndex = FilterOrder+1
            ! when istep = BufferStepNum(ii) - FilterOrder (the first depending step)
            !   we should be using the first coef in FilterCoef -> FilterCoefIndex = 1
            FilterCoefIndex = istep - BufferStepNum(ii) + FilterOrder + 1

            if (FilterCoefIndex < 1) then
                ! if the FilterCoefIndex is less than 1, then it means istep < BufferStepNum(ii) - FilterOrder
                ! aka, the current step does not effect this snapshot,
                ! use the cycle statment to skip the remainder of the body for this current loop
                !   and procceed to the next iteration of the do loop
                cycle
            endif

            ! master updates the 00 mode filters
            if (myid .eq. 0) then
                u00Buffer(:,ii) = u00Buffer(:,ii) + FilterCoef( FilterCoefIndex ) * u00
                w00Buffer(:,ii) = w00Buffer(:,ii) + FilterCoef( FilterCoefIndex ) * w00
            endif

            ! Everyone updates the v and vor filter buffers, with conversion to double complex
              vBuffer(:,:,:,ii) =   vBuffer(:,:,:,ii) + FilterCoef( FilterCoefIndex ) &
                * CMPLX(   v(1:2*my:2,:,:),   v(2:2*my:2,:,:), dp)
            vorBuffer(:,:,:,ii) = vorBuffer(:,:,:,ii) + FilterCoef( FilterCoefIndex ) &
                * CMPLX( vor(1:2*my:2,:,:), vor(2:2*my:2,:,:), dp)

            ! if the current step matches the step number being buffered, then we save this snapshot
            if ( istep .eq. BufferStepNum(ii) ) then
                if (SuppressSavingNum .eq. 0) then
                    ! If no more SuppressSavingNum, then save the snapshot
                    call save_snapshot( ii, istep, time )
                else
                    ! Otherwise, SuppressSavingNum decreases by 1
                    SuppressSavingNum = SuppressSavingNum - 1
                    ! and the snapshot file number increases by 1
                    FileNumber = FileNumber + 1
                    if (myid .eq. 0) then
                        write(*,"(A,I5,A,(D14.6))") "    Suppressing Data Saving, Step ", istep, ", Time ", time
                    endif
                endif

                ! Clear the buffers at this index
                  vBuffer(:,:,:,ii) = 0.0_dp
                vorBuffer(:,:,:,ii) = 0.0_dp
                if (myid .eq. 0) then
                    u00Buffer(:,ii) = 0.0_dp
                    w00Buffer(:,ii) = 0.0_dp
                endif
                ! For the new Step numbers being stored in the buffer,
                ! The smallest should be BufferStepNum(ii) + SampleFreqInSteps (the next one to save)
                ! Therefore the biggest should be:
                BufferStepNum(ii) = BufferStepNum(ii) + BufferDim4Size*SampleFreqInSteps
            endif
        ENDDO
    end subroutine

! ******************************* Filter Restart *******************************
    ! subroutine save_filter_restart_file( RestartFileName )
    ! Save the necessary variables for filter restart to the restart file
    !
    ! CALLED IN: restart_file:write_restart_file
    !
    ! Arguements
    !   RestartFileName: [String, Input] restart file name
    subroutine save_filter_restart_file( RestartFileName )
        use h5save, only: h5save_R_dp, h5save_R1_dp, h5save_R2_dp, &
            h5save_CPartial_Init, h5save_C4Partial_SingleDim3_dp
        ! ----------------------------- Arguments -----------------------------
        character(len=*), intent(in) :: RestartFileName
        ! ------------------------------- Other -------------------------------
        integer :: ii
        ! MPI
        integer :: ierr


        if (myid == 0) then
            ! Master save filter data
            call h5save_R_dp( RestartFileName, "SampleFreqInSteps", real(SampleFreqInSteps, dp))
            call h5save_R_dp( RestartFileName, "SuppressSavingNum", real(SuppressSavingNum, dp))
            call h5save_R_dp( RestartFileName, "FilterOrder"      , real(FilterOrder      , dp))

            call h5save_R1_dp( RestartFileName, "FilterCoef"   , real(FilterCoef   , dp))
            call h5save_R1_dp( RestartFileName, "BufferStepNum", real(BufferStepNum, dp))

            call h5save_R2_dp( RestartFileName, "u00Buffer", u00Buffer )
            call h5save_R2_dp( RestartFileName, "w00Buffer", w00Buffer )
            ! Master initialize file for partial saving
            call h5save_CPartial_Init( RestartFileName,   "vBuffer", (/ my, mx/2, mz, BufferDim4Size /), dp )
            call h5save_CPartial_Init( RestartFileName, "vorBuffer", (/ my, mx/2, mz, BufferDim4Size /), dp )
        endif

        ! ---------------- Save the 4 dimensional filter buffer ----------------
        ! loop over all kz planes
        DO ii = 1, mz
            ! use MPI_BARRIER to ensure only one processor writes to the file at a single instance
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)

            ! if this plane is part of the data the current processor has
            if ( (ii .ge. kb) .and. (ii .le. ke) ) then
                ! write data for this y plane to the file
                call h5save_C4Partial_SingleDim3_dp( RestartFileName,   'vBuffer',   vBuffer(:,:,ii,:), ii )
                call h5save_C4Partial_SingleDim3_dp( RestartFileName, 'vorBuffer', vorBuffer(:,:,ii,:), ii )
            endif

            ! use MPI_BARRIER to ensure only one processor writes to the file at a single instance
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ENDDO
    end subroutine

    ! subroutine read_filter_restart(RestartFileName)
    ! Read the filter restart file
    !
    ! Arguements
    !   RestartFileName: [String, Input] restart file name
    subroutine read_filter_restart( RestartFileName )
        use h5load, only: h5load_R_dp, h5load_R1_dp, h5load_R2_dp, h5load_C4_partial_dp
        ! ----------------------------- Arguments -----------------------------
        character(len=*), intent(in) :: RestartFileName
        ! ------------------------------- Others -------------------------------
        ! Loop variable
        integer :: ii
        ! MPI
        integer :: ierr
        ! istep of the restart file
        integer :: istep_restart

        ! Loop over all processors to ensure only one processor reads the file
        DO ii = 0, numerop-1
            if (myid .eq. ii) then
                istep_restart = NINT( h5load_R_dp(RestartFileName, "step"))

                SuppressSavingNum = NINT(  h5load_R_dp(RestartFileName, "SuppressSavingNum") )
                BufferStepNum     = NINT( h5load_R1_dp(RestartFileName, "BufferStepNum") ) - istep_restart
                ! the restart file now becomes step 0,
                ! so all the buffered step numbers should subtract the step number of the restart file

                ! Each processor reads the data for its range of kz wavenumber kb to ke
                  vBuffer = h5load_C4_partial_dp( RestartFileName,   "vBuffer", (/0,0,kb-1,0/), (/my,mx/2,ke-kb+1,BufferDim4Size/) )
                vorBuffer = h5load_C4_partial_dp( RestartFileName, "vorBuffer", (/0,0,kb-1,0/), (/my,mx/2,ke-kb+1,BufferDim4Size/) )
            endif
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ENDDO

        ! Read buffer for the 00 modess
        if (myid .eq. 0) then ! 00 mode filtering is only done by master
            u00Buffer = h5load_R2_dp( RestartFileName, "u00Buffer" )
            w00Buffer = h5load_R2_dp( RestartFileName, "w00Buffer" )

            write(*,*)
            write(*,"(A,A)") "Initializing filter buffer from restart file: ", RestartFileName
            write(*,"(A,I5,A)") "First ", SuppressSavingNum, " snapshots will not be saved."
        endif
    end subroutine read_filter_restart

! ****************************** Filter Settings ******************************
    ! function read_filter_settings( FilterSettingFile) result( SettingsAgree )
    ! Read filter settings from the setting file, and compare the sample freq in
    !   filter setting with the DNS setting, return false is setting mismatch
    !
    ! Arguements
    !   FilterSettingFile: [String, Input] filename with path of filter setting h5 file
    !
    ! Return: .true. if sampling freq setting agree between filter settings file and DNS settings
    !         .false. otherwise
    function read_filter_settings( FilterSettingFile) result( SettingsAgree )
        use h5load, only: h5load_R_dp, h5load_R1_dp
        ! ----------------------------- Arguments -----------------------------
        character(len=*), intent(in) :: FilterSettingFile
        logical :: SettingsAgree
        ! ------------------------------- Others -------------------------------
        ! MPI
        integer :: istat(MPI_STATUS_SIZE), ierr
        ! The sampling frequency from the filter design
        integer :: FilterSampleFreqInSteps


        ! Filter order, master read and distribute to all slaves
        if (myid .eq. 0) then
            FilterOrder             = NINT( h5load_R_dp( FilterSettingFile, "FilterOrder"             ))
            FilterSampleFreqInSteps = NINT( h5load_R_dp( FilterSettingFile, "FilterSampleFreqInSteps" ))
            write(*,*)
            write(*,"(A,A)") "Reading temporal filter settings from file:  ", FilterSettingFile
            write(*,"(A,I5)") "Filter order: ", FilterOrder
            ! check if the filter designed sampling frequency agrees with the DNS run setting
            if (FilterSampleFreqInSteps .eq. SampleFreqInSteps) then
                write(*,"(A,I5)") "Filter sample fequency in steps agrees with DNS setting: ", SampleFreqInSteps
                SettingsAgree = .true.
            else
                ! if disagree, print error message
                write(*,"(A,I5,A,I5)") "Filter sample fequency in steps: ", FilterSampleFreqInSteps, &
                    " disagrees with DNS setting: ", SampleFreqInSteps
                write(*,"(A)") "ABORTING RUN"
                SettingsAgree = .false.
            endif
            write(*,*)
        endif

        ! Distribute the test result
        call MPI_BCAST( SettingsAgree, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )
        if (SettingsAgree .eqv. .false.) then
            return ! return if the settings do not match
        endif

        call MPI_BCAST( FilterOrder, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
        ! Filter coefs, master read and distribute to all slaves
        ALLOCATE( FilterCoef(FilterOrder+1) )
        if (myid .eq. 0) then
            FilterCoef = h5load_R1_dp( FilterSettingFile, "FilterCoef" )
        endif
        call MPI_BCAST( FilterCoef, FilterOrder+1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    end function read_filter_settings

    ! function compare_filter_settings(RestartFileName) result( SettingsAgree )
    ! Compare current filter settings and those from the restart file,
    !   if settings disagree then do not use the filter restart file
    !
    ! Arguements
    !   RestartFileName: [String, Input] filename with path of restart file
    !
    ! Return: .true. if current filter setting agree with those from restart file
    !         .false. otherwise
    function compare_filter_settings(RestartFileName) result( SettingsAgree )
        use h5load, only: h5load_check_variable_existence, h5load_R_dp, h5load_R1_dp
        ! -------------------------- Global variables --------------------------
        ! For DNS time
        real*4 Deltat,CFL,time,dtr,FixTimeStep
        common /tem/ Deltat,CFL,time,dtr,FixTimeStep
        save /tem/
        ! ----------------------------- Arguments -----------------------------
        character(len=*), intent(in) :: RestartFileName
        logical :: SettingsAgree
        ! ------------------------------- Others -------------------------------
        ! MPI
        integer :: istat(MPI_STATUS_SIZE), ierr
        ! temp variable to store filter coef from restart file
        real(kind=dp), dimension(:), allocatable :: temp


        ! Default to settings agree
        ! changed to disagree if any conditions are not met
        SettingsAgree = .true.

        ! Master does all the comparision
        if (myid .eq. 0 ) then
            if ( INDEX(RestartFileName, ".h5 ") .eq. 0) then
                ! If not an h5 file (INDEX return 0 if substring is not found)
                ! Then this is old version of restart file, no filter settings
                write(*,"(A)") "Old version of restart file, no filter settings."
                SettingsAgree = .false.
            else
                ! Check if filter settings exist in the restart file
                if ( .not. h5load_check_variable_existence(RestartFileName,"FilterOrder") ) then
                    write(*,"(A)") "Filter Order not saved in restart file."
                    SettingsAgree = .false.
                endif
                if ( .not. h5load_check_variable_existence(RestartFileName,"FilterCoef" ) ) then
                    write(*,"(A)") "Filter Coefficients not saved in restart file."
                    SettingsAgree = .false.
                endif
            endif

            if ( SettingsAgree .eqv. .true. ) then
                ! Compare DNS settings ( SampleFreqInSteps, FixTimeStep, mx, my, mz )
                if (SampleFreqInSteps .ne. NINT( h5load_R_dp(RestartFileName, "SampleFreqInSteps") )) then
                    write(*,"(A)") "Restart file SampleFreqInSteps disagree."
                    SettingsAgree = .false.
                endif
                if (FixTimeStep .ne. real(h5load_R_dp(RestartFileName, "FixTimeStep"),sp) ) then
                    write(*,"(A)") "Restart file FixTimeStep disagree."
                    SettingsAgree = .false.
                endif
                if (mx .ne. NINT( h5load_R_dp(RestartFileName, "mx") ) ) then
                    write(*,"(A)") "Restart file mx disagree."
                    SettingsAgree = .false.
                endif
                if (my .ne. NINT( h5load_R_dp(RestartFileName, "my") ) ) then
                    write(*,"(A)") "Restart file my disagree."
                    SettingsAgree = .false.
                endif
                if (mz .ne. NINT( h5load_R_dp(RestartFileName, "mz") ) ) then
                    write(*,"(A)") "Restart file mz disagree."
                    SettingsAgree = .false.
                endif
            endif

            if ( SettingsAgree .eqv. .true. ) then
                ! Compare filter settings (FilterOrder and FilterCoef )
                if (FilterOrder .ne. NINT( h5load_R_dp(RestartFileName, "FilterOrder") ) ) then
                    write(*,"(A)") "Restart file FilterOrder disagree."
                    SettingsAgree = .false.
                else
                    ALLOCATE(temp(FilterOrder+1))
                    temp = h5load_R1_dp(RestartFileName, "FilterCoef" )
                    if (any( FilterCoef .ne. temp )) then
                        write(*,"(A)") "Restart file FilterCoef disagree."
                        SettingsAgree = .false.
                    endif
                    DEALLOCATE(temp)
                endif
            endif
        endif

        ! Master send the results to all other processors
        call MPI_BCAST( SettingsAgree, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )

    end function compare_filter_settings

! ****************************** Snapshot Saving ******************************
    ! subroutine save_snapshot( BufferIndex, istep, time )
    ! Save a snapshot that just finished buffering into a h5 file
    !
    ! Arguments
    !   BufferIndex: [int, Input] the dim4 index in the buffer to save
    !   istep      : [int, Input] the current step number
    !   time       : [single real, Input] the current DNS time
    subroutine save_snapshot( BufferIndex, istep, time )
        ! ----------------------------- Arguments -----------------------------
        ! the dim 4 index in buffer to save
        integer      , intent(in) :: BufferIndex
        ! current step number in the time loop
        integer      , intent(in) :: istep
        ! current time of the time loop, should be istep*DeltaT
        real(kind=sp), intent(in) :: time
        ! ----------------------------- File Name -----------------------------
        character(:), allocatable :: filename
        character(len=4) :: FileNumberString


        ! Generate filename
        write(FileNumberString, '(i4.4)') FileNumber
        filename = BaseFileNameWithPath // FileNumberString // '.h5'
        ! Master print message
        if (myid == 0) then
            write(*,"(A,I5,A,(D14.6),2A)") &
                '    Saving Snapshot, Step: ', istep, &
                ', Time: ', time, ', File: ', filename
        endif

        ! Initialize file (initialize the flow field matrices for later partial saving)
        call initialize_h5_file(filename)
        ! Write flowfield data to the h5 file (in serial)
        call write_flowfields_to_h5(filename, BufferIndex)
        ! Write simulation paramters to h5 file
        call write_parameters_to_h5(filename, istep, time )

        ! increment fileNumber for the module
        FileNumber = FileNumber + 1
    end subroutine

    ! subroutine initialize_h5_file(filename)
    ! This file initializes the flowfield matrices in the h5 file for later partial saving
    !
    ! Parameters:
    !   filename: [string, Input] h5 filename with path
    !
    ! Note:
    !   Only mpi master performs the data matrix initialization, but all
    !   processes should call this function to synchronize progress
    subroutine initialize_h5_file(filename)
        use h5save, only: h5save_CPartial_Init
        ! ----------------------------- Arguments -----------------------------
        character(len=*), intent(in) :: filename
        ! ------------------------------- Others -------------------------------
        ! MPI
        integer :: ierr

        ! Only master initialize the vairable for the h5 files
        if (myid == 0) then
            ! Initialize variables in the h5 for partial saving
            call h5save_CPartial_Init( filename,  'v', (/ mx/2, mz, my /), sp )
            call h5save_CPartial_Init( filename, 'o2', (/ mx/2, mz, my /), sp )
        endif

        ! Ensure variables are initialized before proceding
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end subroutine

    ! subroutine write_flowfields_to_h5(filename, BufferIndex)
    ! This file writes the flowfield data to the h5 file
    ! Write is done in serial, one kz plane at a time
    !
    ! Parameters:
    !   filename   : [string, Input] h5 filename with path
    !   BufferIndex: [int, Input] the dim 4 index in buffer to save
    subroutine write_flowfields_to_h5(filename, BufferIndex)
        use h5save, only: h5save_C3Partial_SingleDim2_sp, h5save_R1_dp
        ! ----------------------------- Arguments -----------------------------
        character(len=*), intent(in) :: filename
        ! the dim 4 index in buffer to save
        integer      , intent(in) :: BufferIndex
        ! ------------------------------- Others -------------------------------
        ! MPI
        integer :: ierr
        ! Loop Index
        integer :: ii
        ! Temp variable
        complex(kind=sp), dimension(mx/2,my) :: temp


        ! loop over all kx planes
        DO ii = 1, mz
            ! use MPI_BARRIER to ensure only one processor writes to the file at a single instance
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)

            ! if this plane is part of the data the current processor has
            if ( (ii .ge. kb) .and. (ii .le. ke) ) then
                ! write data for this y plane to the file
                temp = TRANSPOSE( CMPLX(  vBuffer(:,:,ii,BufferIndex), KIND=sp )) ! Convert to single and transpose
                call h5save_C3Partial_SingleDim2_sp( filename,  'v', temp, ii )
                temp = TRANSPOSE( CMPLX(vorBuffer(:,:,ii,BufferIndex), KIND=sp )) ! Convert to single and transpose
                call h5save_C3Partial_SingleDim2_sp( filename, 'o2', temp, ii )
            endif

            ! use MPI_BARRIER to ensure only one processor writes to the file at a single instance
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ENDDO

        ! Master save the two 00 modes
        if (myid .eq. 0) then
            call h5save_R1_dp( filename, 'u00', u00Buffer(:,BufferIndex) )
            call h5save_R1_dp( filename, 'w00', w00Buffer(:,BufferIndex) )
        endif

    end subroutine

    ! subroutine write_parameters_to_h5(filename, istep, time )
    ! This function writes the simulation parameters to the h5 file
    !
    ! Parameters:
    !   filename: [string, Input] h5 filename with path
    !   istep   : [int, Input] current step number in the time loop
    !   time    : [single, Input] current time of the time loop, should be istep*DeltaT
    !
    ! Note:
    !   Only mpi master performs the data matrix initialization, but all
    !   processes should call this function to synchronize progress
    subroutine write_parameters_to_h5(filename, istep, time )
        use h5save, only: h5save_R_dp, h5save_R1_dp
        ! ----------------------------- Arguments -----------------------------
        character(len=*), intent(in) :: filename
        ! current step number in the time loop
        integer      , intent(in) :: istep
        ! current time of the time loop, should be istep*DeltaT
        real(kind=sp), intent(in) :: time
        ! ---------------------- Physical x and z vectors ----------------------
        real(kind=dp) :: boxLengthX, deltaX, boxLengthZ, deltaZ
        real(kind=dp) :: x(mgalx), z(mgalz)
        integer :: i
        real(kind=dp), parameter :: pi = acos(-1.0_dp)

        ! MPI
        integer :: ierr


        ! All single precision variables are converted to double precision for saving purpose
        ! Only master saves these variables
        if (myid == 0) then
            ! Compute x and z phyiscal space coordinate vectors
            boxLengthX = 2.0_dp * pi / real(fundamentalkx, dp)
            deltaX = boxLengthX / mgalx
            do i = 1, mgalx
                x(i) = (i-1) * deltaX
            enddo

            boxLengthZ = 2.0_dp * pi / real(fundamentalkz, dp)
            deltaZ = boxLengthZ / mgalz
            do i = 1, mgalz
                z(i) = (i-1) * deltaZ
            enddo

            ! Save physical space coordinates to h5 file
            call h5save_R1_dp( filename, 'x', real( x, dp) )
            call h5save_R1_dp( filename, 'y', real( y, dp) )
            call h5save_R1_dp( filename, 'z', real( z, dp) )

            ! Save wavenumbers to h5 file
            ! the DNS code stores the wavenumbers as purely imaginary vector,
            ! so we will extract the imaginary part
            call h5save_R1_dp( filename, 'kx', real( aimag(kx), dp) )
            call h5save_R1_dp( filename, 'kz', real( aimag(kz), dp) )

            ! Save current time, step number, reference Reynolds number and bulk velocity
            call h5save_R_dp( filename, 'time', real( time, dp) )
            call h5save_R_dp( filename, 'step', real( istep, dp))
            call h5save_R_dp( filename, 'Re_ref', real( Re_ref, dp) )
            call h5save_R_dp( filename, 'bulkVelocity', real( bulkVelocity, dp) )
        endif

        ! Ensure all saving are completed before proceding
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end subroutine

end module