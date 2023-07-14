module restart_file
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
    include 'ctes3D'

    public :: read_restart_file_old, read_restart_file_new, write_restart_file

contains
    ! read_restart_file_old(myid, vor,phi,u00,w00,rf0u,rf0w,hg,hv,u00wk,w00wk,phiwk,vorwk)
    ! Read the restart file in the original DNS format
    !
    ! Arguments:
    !   myid     : [int] processor ID
    !   vor, phi : [single, size 2*my,mx1+1,kb:*, Input/Output]
    !              wall normal vorticity and phi for kx-y planes assigned to the current processor at a few kz
    !              variables are allocated outside, and filled with data from the restart file
    !   u00, w00 : [double size my, Input/Ouput]
    !              u and w 00 modes, variables allocated outside and filled with data from the restart file
    !   rf0u,rf0w: [double size my, Input/Ouput]
    !              y derivatives of u and w 00 modes, variables allocated outside and filled with data from the restart file
    !   hg, hv   : [single, size(0:2*my-1,0:mx1,kb:ke), Input/Output)
    !              v and dvdy, variables allocated outside and filled with data computed from the restart file
    !   u00wk = u00
    !   w00wk = w00
    !   phiwk = phi
    !   vorwk = vor
    subroutine read_restart_file_old(myid, vor,phi,u00,w00,rf0u,rf0w,hg,hv,u00wk,w00wk,phiwk,vorwk)
        implicit none
        ! -------------------------- Global variables --------------------------
        ! For DNS time
        real*4 Deltat,CFL,time,dtr,FixTimeStep
        common /tem/ Deltat,CFL,time,dtr,FixTimeStep
        save /tem/
        ! For kz and y parallel assignments
        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save /point/
        ! for file name and file number
        integer iinp,iout,id22,isn,ispf
        character*100 filinp,filout,filstt
        common /ficheros/ iinp,iout,id22,isn,ispf,filinp,filout,filstt
        save /ficheros/

        ! ----------------------------- Arguments -----------------------------
        real*4 vor(mx,mz,*),phi(mx,mz,*)
        real(kind=sp), dimension(0:2*my-1,0:mx1,kb:ke), intent(inout) :: &
            hv,hg,phiwk,vorwk
        real(kind=dp), dimension(my), intent(inout) :: u00,w00,rf0u,rf0w,u00wk,w00wk
        integer :: myid

        ! ------------------------------- Others -------------------------------
        ! MPI
        integer :: istat(MPI_STATUS_SIZE),ierr
        ! Work variable
        real(kind=sp), dimension(:), allocatable :: wk1

        integer :: i,pe,k,j,mxe,mye,mze,ntotr,iproc,mym,mmy2,nbuffsize
        real(kind=sp) :: a0e,Ree,alpe,bete,ce,xx


        ! zero everything
        vor(1:mx,1:mz,1:mmy) = 0.0_sp
        phi(1:mx,1:mz,1:mmy) = 0.0_sp
        u00(1:my) = 0.0_dp
        w00(1:my) = 0.0_dp

        ! Allocate work variable
        nbuffsize = mx*max(mmy*mz,mmz*my)
        ALLOCATE(wk1(nbuffsize))

        ! the time could be read from the restart file below, but here we
        ! set it explicitly to zero, so that the final time contains
        ! information about the averaging period.
        time = 0.0

        ! Read from restart file and distribute
        if (myid.eq.0) then
            ! Master reads from the file and send data to slaves
            open(iinp,file=filinp,status='unknown',form='unformatted',access='direct',recl=8*iwd)
            write(*,*) 'infile opened'
            read(iinp,rec=1) xx,Ree,alpe,bete,a0e,mxe,mye,mze
            ! mxe, mye, mze are the dimensions of the restart file

            write(*,*)
            write(*,*) 'reading input file ...'
            write(*,*) 'time=',time,Ree,alpe,bete,a0e,mxe,mye,mze
            write(*,*) 'mesh:',ce,pe
            write(*,*)
            ntotr=2*mxe*mze
            close(iinp)

            open(iinp,file=filinp,status='unknown',form='unformatted',access='direct',recl=ntotr*iwd)
            ! Read 00 mode from restart file
            read(iinp,rec=1) xx,xx,xx,xx,xx,i,i,i,i,xx,i,i,xx,xx,xx,xx,(wk1(j),j=1,2*mye)

            ! Organize 00 modes
            mym = min(my,mye)
            u00(1:mym) = wk1(1:2*mym:2)
            w00(1:mym) = wk1(2:2*mym:2)

            ! send restart file dimensions and 00 modes to all slaves
            do iproc=1,numerop-1
                call MPI_SEND(mxe, 1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
                call MPI_SEND(mye, 1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
                call MPI_SEND(mze, 1,MPI_INTEGER,iproc,iproc,MPI_COMM_WORLD,ierr)
                call MPI_SEND(u00,my,MPI_REAL8  ,iproc,iproc,MPI_COMM_WORLD,ierr)
                call MPI_SEND(w00,my,MPI_REAL8  ,iproc,iproc,MPI_COMM_WORLD,ierr)
            enddo

            ! Read and organize vor and phi modes, for master node
            mmy2 = min(je,mye)-jb+1
            do j=1,mmy2
                read(iinp,rec=j+jb) (wk1(i),i=1,ntotr)
                call assign(wk1,vor(1,1,j),phi(1,1,j),mx,mz,mxe,mze)
            enddo

            ! Read and organize vor and phi modes, and send to all slaves
            do iproc=1,numerop-1
                do j=jbeg(iproc),min(jend(iproc),mye)
                    read(iinp,rec=j+1) (wk1(i),i=1,ntotr)
                    call MPI_SEND(wk1,ntotr,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,ierr)
                enddo
            enddo
            close(iinp)

        else  ! Slaves receive data from master
            ! Receive 00 mode
            call MPI_RECV(mxe, 1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(mye, 1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(mze, 1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(u00,my,MPI_REAL8  ,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(w00,my,MPI_REAL8  ,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

            ! Receive vor and phi modes and organize
            ntotr=2*mxe*mze
            mmy2 = min(je,mye)-jb+1
            do j=1,mmy2
                call MPI_RECV(wk1,ntotr,MPI_REAL,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
                call assign(wk1,vor(1,1,j),phi(1,1,j),mx,mz,mxe,mze)
            enddo

        endif ! end master and slave separation if

        ! so far, each processor has data in a few y planes
        ! before going any further, exchange cuts in y with cuts in z!
        call chikj2jik(vor,vor,wk1,wk1,myid)
        call chikj2jik(phi,phi,wk1,wk1,myid)

        ! Organize the data, mainly solve for v from phi and compute y derivatives
        call orgainize_data_from_restart_file(vor,phi,u00,w00,rf0u,rf0w,hg,hv,u00wk,w00wk,phiwk,vorwk)

        DEALLOCATE(wk1)

    end subroutine


    ! Subroutine to perform the special first time step in the original code,
    ! Mainly solving for v from phi and compute y derivatives of the 00 modes.
    ! This is only used for reading the old version restart file
    subroutine orgainize_data_from_restart_file(vor,phi,u00,w00,rf0u,rf0w,hg,hv,u00wk,w00wk,phiwk,vorwk)
        use wall_roughness, only: set_wall_roughness, &
            uWallBottom, uWallTop, vWallBottom, vWallTop, wWallBottom, wWallTop
        implicit none
        ! -------------------------- Global variables --------------------------
        ! For kz and y parallel assignments
        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save   /point/
        ! For kx and kz wavenumbers
        complex*8     xalp, xbet
        real*4        alp2,bet2
        integer       iax,icx
        common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),iax(mx),icx(0:mz1)
        save   /wave/

        ! ----------------------------- Arguments -----------------------------
        real(kind=sp), dimension(0:2*my-1,0:mx1,kb:ke), intent(inout) :: &
            vor,phi,hv,hg,phiwk,vorwk
        real(kind=dp), dimension(my), intent(inout) :: u00,w00,rf0u,rf0w,u00wk,w00wk

        ! ------------------------------- Others -------------------------------
        integer i,k,j,k1
        real*8 rk
        complex*8 bcb,bct


        ! Set boundary conditions, ignoring the BC from restart file
        call set_wall_roughness( )

        ! Calculate v from phi
        do k=kb,ke
            k1 = icx(k-1)
            do i=0,mx1
                rk = bet2(k1)+alp2(i)
                ! note on indices:
                ! kb and ke are in the range [1, mz]
                ! the indices of the boundary planes are in the range
                ! [0, mz1], where mz1 = mz - 1
                ! therefore we need to subtract 1 from the loop index
                ! k (which is defined in terms of kb, ke) when
                ! accessing the boundary planes
                bcb = vWallBottom(i,k-1)
                bct = vWallTop(i,k-1)
                ! debugging
                if(abs(bcb) .gt. 1e-4) then
                write(*,*) 'initial boundary conditions'
                write(*,*) 'i = ', i
                write(*,*) 'abs(k) = ', k1
                write(*,*) 'vwallBottom = ', bcb
                write(*,*) 'vWallTop = ', bct
                endif
                ! solve nabla^2 v = phi
                ! phi is the  input
                call Lapvdv(phi(0,i,k),hg(0,i,k),hv(0,i,k),rk,bcb,bct)
                ! hg  is the output: v
                ! hv  is the output: dv/dy
            enddo
        enddo

        ! Compute 00 modes of omega1 and omega3
        ! Compute y derivative of u and w
        u00wk=u00
        w00wk=w00
        call deryr(u00wk,rf0u,my)
        call deryr(w00wk,rf0w,my)
        ! rf0u = du/dy(kx=kz=0) = - o3(kx=kz=0)
        ! rf0w = dw/dy(kx=kz=0) = + o1(kx=kz=0)

        ! Copy data to vorwk and phiwk
        vorwk=vor
        phiwk=phi
        ! vorwk and phi wk are size 0:2*my-1, 0:mx1, kb:ke
        ! the first index includes both real and imaginary part
        ! vorwk: omega2
        ! phiwk: phi = nabla^2 v

    end subroutine


    ! assign(work,vor,phi,mx,mz,mxe,mze)
    ! Organize the data read from the restart file and puts
    ! into vor and phi for each y plane
    !
    ! Arguments:
    !   work     : [single size (2,mxe,mze), Input]
    !              data read from the restart file, containing vor and phi for a single y plane
    !   vor, phi : [single size (mx, mz), Input/Ouput]
    !              vor and phi allocated before calling, getting filled with data for a single y plane
    !   mx, mz   : [int] x and z size for this run
    !   mxe, mze : [int] x and z size for the restart file
    !
    ! Note that mx is twice the size of kx, containing both the real and imaginary parts
    ! and dimension 1 is organized as: real(kx=0), imag(kx=0), real(kx=0.5), image(kx=0.5) ...
    subroutine assign(work,vor,phi,mx,mz,mxe,mze)
        implicit none
        integer   mxm,mzm,klen,kini1,kini2,i,k,k1,k2,mx,mz,mxe,mze
        real*4    vor(mx,mz),phi(mx,mz)
        real*4    work(2,mxe,mze)

        mxm=min(mx,mxe)
        mzm=min(mz,mze)

        klen=(mzm+1)/2
        kini1=mze - klen + 1
        kini2=mz  - klen + 1

        ! loop over the non negative kz
        ! Note: for kx, since the real and imag are next to each other, we only need to loop to the smaller mx size,
        ! the effect will be truncation of high kx wavenumbers (if restart file is larger)
        ! or 0 padding high kx wavenumbers (if restart file is smaller)
        do k=1,klen
            do i=1,mxm
                vor(i,k)=work(1,i,k)
                phi(i,k)=work(2,i,k)
            enddo
        enddo

        ! loop over the negative kz
        do k=1,klen-1
            k1 = k + kini1
            k2 = k + kini2
            do i=1,mxm
                vor(i,k2)=work(1,i,k1)
                phi(i,k2)=work(2,i,k1)
            enddo
         enddo
    end subroutine


    ! Write a new version H5 format restart file
    ! Arguments
    !   write_time       : [double, Input/Output]
    !                      timer for data writing, gets updated here
    !   istep            : [int, Input]
    !                      current DNS step number
    !   phi, vor, v, dvdy: [single, size (2*my,mx/2,kb:ke), Input]
    !                      nabla^2 v, omega_2, v and dv/dy
    !   u00, w00         : [double, size my, Input]
    !                      kx = kz = 0 mode of u and w
    !   myid             : [int, Input]
    !                      processor ID number
    subroutine write_restart_file( write_time, istep, phi, vor, v, dvdy, u00, w00, myid)
        use h5save, only: check_filename, h5save_R_dp, h5save_R1_dp, &
            h5save_CPartial_Init, h5save_C3Partial_SingleDim3_sp
        use save_flowfield, only: SampleFreqInSteps
        implicit none
        ! -------------------------- Global variables --------------------------
        ! For DNS time
        real*4 Deltat,CFL,time,dtr,FixTimeStep
        common /tem/ Deltat,CFL,time,dtr,FixTimeStep
        save /tem/
        ! for file name and file number
        integer iinp,iout,id22,isn,ispf
        character*100 filinp,filout,filstt
        common /ficheros/ iinp,iout,id22,isn,ispf,filinp,filout,filstt
        save /ficheros/
        ! For kz and y parallel assignments
        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save /point/

        ! ----------------------------- Arguments -----------------------------
        ! timer for file writing
        real(kind=dp), intent(inout) :: write_time
        ! current step number in the time loop
        integer      , intent(in   ) :: istep
        ! phi, omega_3, v and dv/dy
        ! Note , all real where the first dimension contains both the real and imag parts
        real(kind=sp), intent(in   ), dimension(2*my,mx/2,kb:ke) :: phi, vor, v, dvdy
        ! kx = kz = 0 modes for u and w
        real(kind=dp), intent(in   ), dimension(my) :: u00, w00
        ! processor ID
        integer      , intent(in   ) :: myid

        ! ------------------------------- Others -------------------------------
        ! MPI
        integer :: istat(MPI_STATUS_SIZE), ierr
        ! h5 filename
        character(len=200) :: FileOut
        ! Loop variables
        integer :: ii, jj


        ! ----------------------- Timer for file writing -----------------------
        if (myid.eq.0) then
            write_time = -MPI_WTIME()
        endif

        ! ---------------------------- H5 filename ----------------------------
        ! create and check filename
        write(FileOut,"(2A,I3.3,A)") filout(1:index(filout,' ')-1), "_restart_", id22, ".h5"
        if (myid .eq. 0) then
            write(*,*) " "
            write(*,"(A,I5)") "  Saving restart file, Step ", istep
        endif
        call check_filename( FileOut, myid )

        ! -------------------------- Master Processor --------------------------
        if (myid .eq. 0) then
            ! Initialize variable for partial saving
            call h5save_CPartial_Init( FileOut, "phi" , (/ my, mx/2, mz/), sp )
            call h5save_CPartial_Init( FileOut, "v"   , (/ my, mx/2, mz/), sp )
            call h5save_CPartial_Init( FileOut, "dvdy", (/ my, mx/2, mz/), sp )
            call h5save_CPartial_Init( FileOut, "vor" , (/ my, mx/2, mz/), sp )
            ! Save time and step number
            call h5save_R_dp( FileOut, "time" , real(time ,dp) )
            call h5save_R_dp( FileOut, "istep", real(istep,dp) )
            ! Save 00 modes for u and w
            call h5save_R1_dp( FileOut, "u00", u00)
            call h5save_R1_dp( FileOut, "w00", w00)
            ! Save grid size
            call h5save_R_dp( FileOut, "mx", real(mx,dp) ) ! Note this is twice the size of kx
            call h5save_R_dp( FileOut, "my", real(my,dp) )
            call h5save_R_dp( FileOut, "mz", real(mz,dp) )
            ! Save DNS run settings
            call h5save_R_dp( FileOut, "FixTimeStep", real(FixTimeStep,dp) )
            call h5save_R_dp( FileOut, "SampleFreqInSteps", real(SampleFreqInSteps,dp) )
        endif

        ! -------------------------- Save data to h5 --------------------------
        DO jj = 0, numerop-1 ! Loop over all processor
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            if (myid .eq. jj) then
                DO ii = kb, ke ! Loop over all planes of this processor and save the planes
                    call h5save_C3Partial_SingleDim3_sp( FileOut, "phi" , r2c_array( phi(:,:,ii)), ii )
                    call h5save_C3Partial_SingleDim3_sp( FileOut, "v"   , r2c_array(   v(:,:,ii)), ii )
                    call h5save_C3Partial_SingleDim3_sp( FileOut, "dvdy", r2c_array(dvdy(:,:,ii)), ii )
                    call h5save_C3Partial_SingleDim3_sp( FileOut, "vor" , r2c_array( vor(:,:,ii)), ii )
                ENDDO
            endif
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ENDDO

        ! ----------------------- Increment File Number -----------------------
        id22 = id22+1
    end subroutine


    ! Convert a real array with twice the size for dimension 1 to a complex array
    ! Used in writing the new version restart file
    function r2c_array( array_in ) result( array_out )
        implicit none

        real(kind=sp), dimension(2*my, mx/2), intent(in) :: array_in
        complex(kind=sp), dimension(my, mx/2):: array_out

        array_out = CMPLX( array_in(1:2*my:2,:), array_in(2:2*my:2,:), sp)

    end function


    ! read_restart_file_new(myid, vor,phi,u00,w00,dudy00,dwdy00,v,dvdy, u00wk,w00wk,phiwk,vorwk)
    ! Read the restart file in the new format
    !
    ! Arguments:
    !   myid         : [int] processor ID
    !   vor, phi     : [single, size 2*my,mx1+1,kb:*, Input/Output]
    !                  wall normal vorticity and phi for kx-y planes assigned to the current processor at a few kz
    !                  variables are allocated outside, and filled with data from the restart file
    !   u00, w00     : [double size my, Input/Ouput]
    !                  u and w 00 modes, variables allocated outside and filled with data from the restart file
    !   dudy00,dwdy00: [double size my, Input/Ouput]
    !                  y derivatives of u and w 00 modes, variables allocated outside and filled with data from the restart file
    !   v, dvdy      : [single, size(0:2*my-1,0:mx1,kb:ke), Input/Output)
    !                  v and dvdy, variables allocated outside and filled with data computed from the restart file
    !   u00wk = u00
    !   w00wk = w00
    !   phiwk = phi
    !   vorwk = vor
    subroutine read_restart_file_new(myid, vor,phi,u00,w00,dudy00,dwdy00,v,dvdy, u00wk,w00wk,phiwk,vorwk)
        use h5load, only: h5load_R_dp
        implicit none
        ! -------------------------- Global variables --------------------------
        ! for file name and file number
        integer iinp,iout,id22,isn,ispf
        character*100 filinp,filout,filstt
        common /ficheros/ iinp,iout,id22,isn,ispf,filinp,filout,filstt
        save /ficheros/
        ! For kz and y parallel assignments
        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save /point/
        ! For DNS time
        real*4 Deltat,CFL,time,dtr,FixTimeStep
        common /tem/ Deltat,CFL,time,dtr,FixTimeStep
        save /tem/
        ! ----------------------------- Arguments -----------------------------
        ! phi, omega_3, v and dv/dy
        ! Note , all real where the first dimension contains both the real and imag parts
        real(kind=sp), intent(inout), dimension(2*my,mx/2,kb:ke) :: phi, vor, v, dvdy, phiwk,vorwk
        ! kx = kz = 0 modes for u and w
        real(kind=dp), intent(inout), dimension(my) :: u00, w00, dudy00, dwdy00, u00wk,w00wk
        ! processor ID
        integer      , intent(in   ) :: myid
        ! ------------------------------- Others -------------------------------
        ! MPI
        integer :: istat(MPI_STATUS_SIZE), ierr
        ! Restart file grid size
        integer :: mxr, myr, mzr


        ! the time could be read from the restart file below, but here we
        ! set it explicitly to zero, so that the final time contains
        ! information about the averaging period.
        time = 0.0

        ! Read restart file grid dimensions and distribute
        if (myid .eq. 0) then
            ! (NINT rounds real to int)
            mxr = NINT( h5load_R_dp(filinp, "mx") )
            myr = NINT( h5load_R_dp(filinp, "my") )
            mzr = NINT( h5load_R_dp(filinp, "mz") )
        endif
        ! Send to all slave processors
        call MPI_BCAST( mxr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
        call MPI_BCAST( myr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
        call MPI_BCAST( mzr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

        ! Read 00 modes
        if (myid .eq. 0) then
            u00 = read_and_assign_00mode( filinp, "u00", myr )
            w00 = read_and_assign_00mode( filinp, "w00", myr )
        endif
        ! Send to all slave processors
        call MPI_BCAST( u00, my, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
        call MPI_BCAST( w00, my, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
        ! compute y derivatives
        call deryr(u00,dudy00,my)
        call deryr(w00,dwdy00,my)

        ! Master reads phi, vor, v, dvdy and distribute
        phi  = read_and_assign_matrix( filinp, "phi" , mxr, myr, mzr, myid)
        vor  = read_and_assign_matrix( filinp, "vor" , mxr, myr, mzr, myid)
        v    = read_and_assign_matrix( filinp, "v"   , mxr, myr, mzr, myid)
        dvdy = read_and_assign_matrix( filinp, "dvdy", mxr, myr, mzr, myid)

        ! assign work variables
        u00wk = u00
        w00wk = w00
        phiwk = phi
        vorwk = vor

    end subroutine


    ! Read 00 mode from restart file and correct for potential y grid size differenc
    ! Arguements
    !   FileName: [string]
    !   VarName : [string] either "u00" or "w00"
    !   myr     : [int] y grid size of restart file
    !
    ! Output
    !   u00 or w00 [double real, size my]
    function read_and_assign_00mode( FileName, VarName, myr ) result( mode00 )
        use h5load, only: h5load_R1_dp
        ! ----------------------------- Arguments -----------------------------
        character(len=*), intent(in) :: FileName, VarName
        ! y grid size of restart file
        integer, intent(in) :: myr
        ! output
        real(kind=dp), dimension(my) :: mode00
        ! ------------------------------- Others -------------------------------
        ! temp file for reading from restart file
        real(kind=dp), dimension(:), allocatable :: temp
        ! min of y grid size (smaller of restart file y size and current y size)
        integer :: my_min


        ! min of y grid size (smaller of restart file y size and current y size)
        my_min = min(my, myr)

        ! Intialize
        if ( VarName .eq. "u00" ) then
            mode00(1:my) = 1.0_dp ! u is initialized to 1 (u centerline is close to 1)
        elseif ( VarName .eq. "w00" ) then
            mode00(1:my) = 0.0_dp ! w is initialized to 0 (w centerline is close to 0)
        endif

        ! allocate and read from restart file
        ALLOCATE(temp(myr))
        temp = h5load_R1_dp(FileName, VarName)

        ! assign to 00 mode, taking into account that y gride size might be different
        ! for different y grid size, will assign data starting from the top/bottom wall
        mode00( 1:my_min/2 ) = temp( 1:my_min/2 ) ! bottom half data
        mode00( my-my_min/2+1:my ) = temp( myr-my_min/2+1:myr ) ! top half data

        DEALLOCATE(temp)

    end function read_and_assign_00mode


    ! Read phi, vor, v or dvdy from restart file and correct for potential grid size differenc
    ! Arguements
    !   FileName     : [string]
    !   VarName      : [string] either "phi", "vor", "v" or "dvdy"
    !   mxr, myr, mzr: [int] x, y,z grid size of restart file
    !   myid         : [int] processor ID
    !
    ! Output
    !   phi, v, vor or dvdy :[single real, size size (2*my,mx/2,kb:ke)]
    !                        The first dimension contains both real and imag parts
    function read_and_assign_matrix( FileName, VarName, mxr, myr, mzr, myid) result(matrix)
        use h5load, only: h5load_C3_sp
        ! -------------------------- Global variables --------------------------
        ! For kz and y parallel assignments
        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save /point/
        ! ----------------------------- Arguments -----------------------------
        character(len=*), intent(in) :: filename, varname
        integer, intent(in) :: mxr, myr, mzr
        integer, intent(in) :: myid
        real(kind=sp), dimension(2*my,mx/2,kb:ke) :: matrix
        ! ------------------------------- Others -------------------------------
        ! MPI
        integer :: istat(MPI_STATUS_SIZE), ierr
        ! matrices for restart file reading
        complex(kind=sp), dimension(:,:,:), allocatable :: mat_restart, mat_full
        ! buffer for master mpi send
        real(kind=sp), dimension(:,:,:), allocatable :: buffer
        ! mpi send and recv size
        integer :: nsend, nrecv
        ! min of x,y,z grid size (smaller of restart file and current size)
        integer :: mx_min, my_min, mz_min
        ! Loop index
        integer :: kk


        ! min of x,y,z grid size (smaller of restart file and current size)
        mx_min = min(mx,mxr)
        my_min = min(my,myr)
        mz_min = min(mz,mzr)

        if (myid .eq. 0) then
            ! Master read and distibute data
            ALLOCATE(mat_restart(myr, mxr/2, mzr)) ! restart file size
            ALLOCATE(mat_full(my, mx/2, mz)) ! current run grid size

            ! Initialize to 0
            mat_full = 0.0_sp
            ! Read from restart file
            mat_restart = h5load_C3_sp(FileName, VarName )

            ! Deal with potential grid mismatch
            ! for y  match data starting from the two walls, the center part maybe discarded (if restart is bigger)
            !        or left to initial value of 0 (if restart is smaller)
            ! for kx match data for low kx, and discard or left as 0 for high kx
            ! for kz match data from the two ends (low kz), discard or left as 0 for high kz

            ! --------------- non-negtive kz (including kz = 0) ---------------
            ! Bottom half of channel
               mat_full(1:my_min/2, 1:mx_min/2, 1:(mz_min+1)/2) = &
            mat_restart(1:my_min/2, 1:mx_min/2, 1:(mz_min+1)/2)
            ! top half of channel
               mat_full(my -my_min/2+1:my , 1:mx_min/2, 1:(mz_min+1)/2) =  &
            mat_restart(myr-my_min/2+1:myr, 1:mx_min/2, 1:(mz_min+1)/2)
            ! -------------------------- negative kz --------------------------
            ! Bottom half of channel
               mat_full(1:my_min/2, 1:mx_min/2, mz -(mz_min-1)/2+1:mz) = &
            mat_restart(1:my_min/2, 1:mx_min/2, mzr-(mz_min-1)/2+1:mz)
            ! top half of channel
               mat_full(my -my_min/2+1:my , 1:mx_min/2, mz -(mz_min-1)/2+1:mz) = &
            mat_restart(myr-my_min/2+1:myr, 1:mx_min/2, mzr-(mz_min-1)/2+1:mz)

            DEALLOCATE(mat_restart)

            ! save data for master itself
            matrix(1:2*my:2,1:mx/2,kb:ke) =  real( mat_full(:,:,kb:ke), sp)
            matrix(2:2*my:2,1:mx/2,kb:ke) = aimag( mat_full(:,:,kb:ke)    )

            ! loop over all slave processors and send data
            DO kk = 1, numerop-1
                ! allocate send buffer
                ALLOCATE(buffer(1:2*my,1:mx/2,kbeg(kk):kend(kk)))
                ! fill it with real and imag parts, for the kz planes for this target slave processor
                buffer(1:2*my:2,:,kbeg(kk):kend(kk)) =  real( mat_full(:,:,kbeg(kk):kend(kk)), sp)
                buffer(2:2*my:2,:,kbeg(kk):kend(kk)) = aimag( mat_full(:,:,kbeg(kk):kend(kk))    )
                ! send buffer size
                nsend = 2*my*(mx/2)*(kend(kk)-kbeg(kk)+1)
                ! send data
                call MPI_SEND( buffer, nsend, MPI_REAL, kk, 0, MPI_COMM_WORLD, ierr )
                DEALLOCATE(buffer)
            ENDDO

            DEALLOCATE(mat_full)

        else
            ! all slaves recieve data
            nrecv = 2*my*(mx/2)*(ke-kb+1)
            call MPI_RECV( matrix, nrecv, MPI_REAL, 0, 0, MPI_COMM_WORLD, istat, ierr )

        endif
    END function read_and_assign_matrix

end module restart_file