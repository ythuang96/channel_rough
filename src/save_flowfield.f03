module save_flowfield
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
    include 'ctes3D'


    ! Global variables, distribution pointers
    integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
    common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                   kbeg(0:numerop-1),kend(0:numerop-1), &
                   jb,je,kb,ke,mmy,mmz
    save   /point/

    ! Buffer for Fourier space velocity and vorticity
    complex(kind=sp), dimension(:,:,:), allocatable :: uBuffer, vBuffer, wBuffer, o1Buffer, o2Buffer, o3Buffer

    ! whether or not to collect flowfield
    logical :: collectFlowfield
    ! File number and file path
    integer :: fileNumber
    character(:), allocatable :: baseFileNameWithPath
    ! Sampling Fequencies
    integer :: SampleFreqInSteps


    ! Public Variables
    public :: collectFlowfield, SampleFreqInSteps
    ! Public Subroutines
    public :: initialize_save_flowfield_module, cleanup_save_flowfield_module
    public :: assess_whether_to_collect_flowfield_step
    public :: save_plane_data_to_buffer, write_h5

contains
    ! subroutine initialize_save_flowfield_module(runNameWithPath)
    ! This function initializes this module
    ! Including setting the sampling frequencies and initializing the buffers
    !
    ! Parameters
    !   runNameWithPath: [string, Input] directory for file saving appended with run name'
    !
    ! Note: this subrouting only needs to be run once at the very beginning of
    !       the code by all processors
    subroutine initialize_save_flowfield_module(runNameWithPath)
        character(len=*), intent(in) :: runNameWithPath
        integer :: mpiRank, mpiError
        integer, parameter :: mpiMaster = 0

        ! Set Sampling parameters
        SampleFreqInSteps = 50

        collectFlowfield = .false.
        fileNumber = 1
        baseFileNameWithPath = trim(adjustl(runNameWithPath)) // '_flowfield.'

        ! Initialize Fourier Space Buffers
        allocate( uBuffer(mx/2, mz, jb:je))
        allocate( vBuffer(mx/2, mz, jb:je))
        allocate( wBuffer(mx/2, mz, jb:je))
        allocate(o1Buffer(mx/2, mz, jb:je))
        allocate(o2Buffer(mx/2, mz, jb:je))
        allocate(o3Buffer(mx/2, mz, jb:je))
        ! Initialize Fourier Buffers
        uBuffer = (0.0_sp, 0.0_sp)
        vBuffer = (0.0_sp, 0.0_sp)
        wBuffer = (0.0_sp, 0.0_sp)
        o1Buffer = (0.0_sp, 0.0_sp)
        o2Buffer = (0.0_sp, 0.0_sp)
        o3Buffer = (0.0_sp, 0.0_sp)

        ! Write sampling information
        call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiError)
        if (mpiRank == mpiMaster) then
            write(*,'(a,i4)') 'sampling frequency of flowfield in timesteps: ', SampleFreqInSteps
        endif
    end subroutine


    ! subroutine cleanup_save_flowfield_module
    ! This function deallocates the buffers for this module
    !
    ! Note: this subrouting only needs to be run once at the very end of
    !       the code by all processors
    subroutine cleanup_save_flowfield_module
        ! Deallocate Fourier Space Buffers
        deallocate( uBuffer)
        deallocate( vBuffer)
        deallocate( wBuffer)
        deallocate(o1Buffer)
        deallocate(o2Buffer)
        deallocate(o3Buffer)
    end subroutine


    ! subroutine assess_whether_to_collect_flowfield_step(timeStepNr)
    ! This funciton access whether or not to collect flowfield based on step
    ! number, and set the public variable "collectFlowfield" to true or false
    !
    ! Parameters
    !   timeStepNr: [integer, Input] current time step number
    subroutine assess_whether_to_collect_flowfield_step(timeStepNr)
        integer, intent(in) :: timeStepNr
        if ( (mod(timeStepNr  ,SampleFreqInSteps) == 0) .or.  &
            ((mod(timeStepNr-1,SampleFreqInSteps) == 0) .and. (timeStepNr.ne.1)) ) then
            ! "Double-Pulse" Case
            collectFlowfield = .true.
        else
            collectFlowfield = .false.
        endif
    end subroutine


    ! subroutine save_plane_data_to_buffer( u, v, w, o1, o2, o3, yplane )
    ! This function saves data at a sinlge y plane to the buffers
    !
    ! Parameters:
    !    u,  v,  w: [complex, size (mx/2, mz), Input]  velocity fields at a single y plane
    !   o1, o2, o3: [complex, size (mx/2, mz), Input] vorticity fields at a single y plane
    !   yplane    : [integer, Input] the y plane index of the data
    subroutine save_plane_data_to_buffer( u, v, w, o1, o2, o3, yplane )
        complex(kind=sp), intent(in), dimension(:,:) :: u, v, w, o1, o2, o3
        integer, intent(in) :: yplane

        uBuffer(:,:,yplane) = u
        vBuffer(:,:,yplane) = v
        wBuffer(:,:,yplane) = w
        o1Buffer(:,:,yplane) = o1
        o2Buffer(:,:,yplane) = o2
        o3Buffer(:,:,yplane) = o3
    end subroutine


    ! subroutine write_h5(fundamentalkx, fundamentalkz, y, kx, kz, time, Re_ref, bulkVelocity)
    ! This function writes the data to a h5 file
    !
    ! Parameters:
    !   fundamentalkx, fundamentalkz: [real, Input]
    !                                 fundamental wavenumbers in the streamwise and spanwise directions
    !   y                           : [real, size my, Input]
    !                                 wall normal coordinate vector
    !   kx                          : [compex, size mx/2, Input]
    !                                 streamwise wavenumber vector
    !   kz                          : [compex, size mz, Input]
    !                                 spanwise wavenumber vector
    !   time                        : [real, Input]
    !                                 current simulation time
    !   Re_ref                      : [real, Input]
    !                                 Reference reynolds number
    !                                 (based on the reference velocity, which is approximately, but not equal to the centerline velocity)]
    !   bulkVelocity                : [real, Input]
    !                                 bulkVelocity
    subroutine write_h5(fundamentalkx, fundamentalkz, y, kx, kz, time, Re_ref, bulkVelocity)
        real(kind=sp), intent(in) :: fundamentalkx, fundamentalkz
        real(kind=sp), intent(in) :: y(:)
        complex(kind=sp), intent(in) :: kx(:), kz(:)
        real(kind=sp), intent(in) :: time, Re_ref
        real(kind=dp), intent(in) :: bulkVelocity

        character(:), allocatable :: filename
        character(len=4) :: fileNumberString

        ! Generate and update the fileName
        write(fileNumberString, '(i4.4)') fileNumber
        filename = baseFileNameWithPath // fileNumberString // '.h5'

        ! Initialize file (initialize the flow field matrices for later partial saving)
        call initialize_h5_file(filename)
        ! Write flowfield data to the h5 file (in serial)
        call write_flowfields_to_h5(filename)
        ! Write simulation paramters to h5 file
        call write_parameters_to_h5(filename, fundamentalkx, fundamentalkz, y, kx, kz, time, Re_ref, bulkVelocity)

        ! update variables for the module by:
        ! increment fileNumber
        fileNumber = fileNumber + 1
        ! setting collectFlowfield to false
        collectFlowfield = .false.
    end subroutine


    ! subroutine initialize_h5_file(filename)
    ! This file initializes the flowfield matrices in the h5 file for later
    ! partial saving
    !
    ! Parameters:
    !   filename: [string, Input] h5 filename with path
    !
    ! Note:
    !   Only mpi master performs the data matrix initialization, but all
    !   processes should call this function to synchronize progress
    subroutine initialize_h5_file(filename)
        use h5save, only: h5save_CPartial_Init_sp
        ! filename
        character(len=*), intent(in) :: filename
        ! MPI Variables
        integer :: mpiRank, mpiError
        integer, parameter :: mpiMaster = 0

        ! Get the rank of current processor
        call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiError)

        ! Only master initialize the vairable for the h5 files
        if (mpiRank == mpiMaster) then
            ! Initialize variables in the h5 for the velocities for partial saving
            call h5save_CPartial_Init_sp( filename,  'u', (/ mx/2, mz, my /) )
            call h5save_CPartial_Init_sp( filename,  'v', (/ mx/2, mz, my /) )
            call h5save_CPartial_Init_sp( filename,  'w', (/ mx/2, mz, my /) )
            ! Initialize variables in the h5 for the vorticities for partial saving
            call h5save_CPartial_Init_sp( filename, 'o1', (/ mx/2, mz, my /) )
            call h5save_CPartial_Init_sp( filename, 'o2', (/ mx/2, mz, my /) )
            call h5save_CPartial_Init_sp( filename, 'o3', (/ mx/2, mz, my /) )
        endif

        ! Ensure variables are initialized before proceding
        call MPI_BARRIER(MPI_COMM_WORLD, mpiError)
    end subroutine


    ! subroutine write_flowfields_to_h5(filename)
    ! This file writes the flowfield data to the h5 file
    ! Write is done in serial, one y plane at a time
    !
    ! Parameters:
    !   filename: [string, Input] h5 filename with path
    subroutine write_flowfields_to_h5(filename)
        use h5save, only: h5save_C3Partial_SingleDim3_sp
        character(len=*), intent(in) :: filename
        integer :: ii
        integer :: mpiError

        ! loop over all y planes
        DO ii = 1, my
            ! use MPI_BARRIER to ensure only one processor writes to the file at a single instance
            call MPI_BARRIER(MPI_COMM_WORLD, mpiError)

            ! if this plane is part of the data the current processor has
            if ( (ii .ge. jb) .and. (ii .le. je) ) then
                ! write data for this y plane to the file
                call h5save_C3Partial_SingleDim3_sp( fileName,  'u',  uBuffer(:,:,ii), ii )
                call h5save_C3Partial_SingleDim3_sp( fileName,  'v',  vBuffer(:,:,ii), ii )
                call h5save_C3Partial_SingleDim3_sp( fileName,  'w',  wBuffer(:,:,ii), ii )
                call h5save_C3Partial_SingleDim3_sp( fileName, 'o1', o1Buffer(:,:,ii), ii )
                call h5save_C3Partial_SingleDim3_sp( fileName, 'o2', o2Buffer(:,:,ii), ii )
                call h5save_C3Partial_SingleDim3_sp( fileName, 'o3', o3Buffer(:,:,ii), ii )
            endif

            ! use MPI_BARRIER to ensure only one processor writes to the file at a single instance
            call MPI_BARRIER(MPI_COMM_WORLD, mpiError)
        ENDDO
    end subroutine


    ! subroutine write_parameters_to_h5(filename, fundamentalkx, fundamentalkz, y, kx, kz, time, Re_ref, bulkVelocity)
    ! This function writes the simulation parameters to the h5 file
    !
    ! Parameters:
    !   filename                    : [string, Input]
    !                                 h5 filename with path
    !   fundamentalkx, fundamentalkz: [real, Input]
    !                                 fundamental wavenumbers in the streamwise and spanwise directions
    !   y                           : [real, size my, Input]
    !                                 wall normal coordinate vector
    !   kx                          : [compex, size mx/2, Input]
    !                                 streamwise wavenumber vector
    !   kz                          : [compex, size mz, Input]
    !                                 spanwise wavenumber vector
    !   time                        : [real, Input]
    !                                 current simulation time
    !   Re_ref                      : [real, Input]
    !                                 Reference reynolds number
    !                                 (based on the reference velocity, which is approximately, but not equal to the centerline velocity)]
    !   bulkVelocity                : [real, Input]
    !                                 bulkVelocity
    !
    ! Note:
    !   Only mpi master performs the data matrix initialization, but all
    !   processes should call this function to synchronize progress
    subroutine write_parameters_to_h5(filename, fundamentalkx, fundamentalkz, y, kx, kz, time, Re_ref, bulkVelocity)
        use h5save, only: h5save_R_dp, h5save_R1_dp
        character(len=*), intent(in) :: filename
        real(kind=sp), intent(in) :: fundamentalkx, fundamentalkz
        real(kind=sp), intent(in) :: y(:)
        complex(kind=sp), intent(in) :: kx(:), kz(:)
        real(kind=sp), intent(in) :: time, Re_ref
        real(kind=dp), intent(in) :: bulkVelocity

        real(kind=dp) :: boxLengthX, deltaX, boxLengthZ, deltaZ
        real(kind=dp) :: x(mgalx), z(mgalz)
        integer :: i

        real(kind=dp), parameter :: pi = acos(-1.0_dp)

        ! MPI Variables
        integer :: mpiRank, mpiError
        integer, parameter :: mpiMaster = 0

        ! All single precision variables are converted to double precision for saving purpose


        ! Get the rank of current processor
        call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiError)

        ! Only master saves these variables
        if (mpiRank == mpiMaster) then
            ! Compute x and z phyiscal space coordinate vectors
            boxLengthX = 2.0_dp * pi / real(fundamentalkx, dp)
            deltaX = boxLengthX / mgalx
            do i = 1, mgalx
                x(i) = (i-1) * deltaX
            enddo

            boxLengthZ = 2.0_sp * pi / real(fundamentalkz, dp)
            deltaZ = boxLengthZ / mgalz
            do i = 1, mgalz
                z(i) = (i-1) * deltaZ
            enddo

            ! Save physical space coordinates to h5 file
            call h5save_R1_dp( fileName, 'x', real( x, dp) )
            call h5save_R1_dp( fileName, 'y', real( y, dp) )
            call h5save_R1_dp( fileName, 'z', real( z, dp) )

            ! Save wavenumbers to h5 file
            ! the DNS code stores the wavenumbers as purely imaginary vector,
            ! so we will extract the imaginary part
            call h5save_R1_dp( fileName, 'kx', real( aimag(kx), dp) )
            call h5save_R1_dp( fileName, 'kz', real( aimag(kz), dp) )

            ! Save current time, reference Reynolds number and bulk velocity
            call h5save_R_dp( fileName, 'time', real( time, dp) )
            call h5save_R_dp( fileName, 'Re_ref', real( Re_ref, dp) )
            call h5save_R_dp( fileName, 'bulkVelocity', real( bulkVelocity, dp) )
        endif

        ! Ensure all saving are completed before proceding
        call MPI_BARRIER(MPI_COMM_WORLD, mpiError)
    end subroutine


end module save_flowfield


