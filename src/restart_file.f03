module restart_file
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
    include 'ctes3D'

    public :: read_restart_file_old

contains
    ! read_restart_file_old(vor,phi,u00,w00,wk1,myid)
    ! Read the restart file in the original DNS format
    !
    ! Arguments:
    !   vor, phi : [single, size 2*my,mx1+1,kb:*, Input/Output]
    !              wall normal vorticity and phi for kx-y planes assigned to the current processor at a few kz
    !              variables are allocated outside, and filled with data from the restart file
    !   u00, w00 : [single size my, Input/Ouput]
    !              u and w 00 modes, variables allocated outside and filled with data from the restart file
    !   wk1      : [single, work variable]
    !   myid     : [int] processor ID
    subroutine read_restart_file_old(vor,phi,u00,w00,wk1,myid)
        use boundary_planes, only: set_bc_from_restart_file
        implicit none

        integer istat(MPI_STATUS_SIZE),ierr

        integer myid,master,i,pe,k
        real*4  a0e

        real*4 vor(mx,mz,*),phi(mx,mz,*),wk1(*)
        real*8 u00(*),w00(*)

        real*4  Ree,alpe,bete,ce,xx

        real*4 Deltat,CFL,time,dtr,FixTimeStep
        common /tem/ Deltat,CFL,time,dtr,FixTimeStep
        save /tem/

        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save /point/

        integer iinp,iout,id22,isn,ispf
        character*100 filinp,filout,filstt
        common /ficheros/ iinp,iout,id22,isn,ispf,filinp,filout,filstt
        save /ficheros/

        integer j,mxe,mye,mze,ntotr,iproc,mym,mmy2

        real*4, allocatable :: uWallBottom(:,:), uWallTop(:,:)
        real*4, allocatable :: vWallBottom(:,:), vWallTop(:,:)
        real*4, allocatable :: wWallBottom(:,:), wWallTop(:,:)


        master = 0

        ! zero everything
        do k=1,mz
            do j=1,mmy
                do i=1,mx
                    vor(i,k,j) = 0.
                    phi(i,k,j) = 0.
                enddo
            enddo
        enddo

        do j=1,my
            u00(j) = 0d0
            w00(j) = 0d0
        enddo

        ! the time could be read from the restart file below, but here we
        ! set it explicitly to zero, so that the final time contains
        ! information about the averaging period.
        time = 0.0

        ! Read from restart file and distribute
        if (myid.eq.master) then
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
            do j=1,mym
                u00(j) = wk1(2*j-1)
                w00(j) = wk1(2*j)
            enddo

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

            ! velocity boundary conditions: master process
            ! allocate temporary memory to read data from restart file
            allocate(uWallBottom(mxe,mze))
            allocate(uWallTop(mxe,mze))
            allocate(vWallBottom(mxe,mze))
            allocate(vWallTop(mxe,mze))
            allocate(wWallBottom(mxe,mze))
            allocate(wWallTop(mxe,mze))
            ! read values from file
            read(iinp, rec=mye+2)((uWallBottom(i,k), uWallTop(i,k), i=1,mxe), k=1,mze)
            read(iinp, rec=mye+3)((vWallBottom(i,k), vWallTop(i,k), i=1,mxe), k=1,mze)
            read(iinp, rec=mye+4)((wWallBottom(i,k), wWallTop(i,k), i=1,mxe), k=1,mze)
            close(iinp)

        else  ! Slaves receive data from master
            ! Receive 00 mode
            call MPI_RECV(mxe, 1,MPI_INTEGER,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(mye, 1,MPI_INTEGER,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(mze, 1,MPI_INTEGER,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(u00,my,MPI_REAL8  ,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(w00,my,MPI_REAL8  ,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

            ! Receive vor and phi modes and organize
            ntotr=2*mxe*mze
            mmy2 = min(je,mye)-jb+1
            do j=1,mmy2
                call MPI_RECV(wk1,ntotr,MPI_REAL,master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
                call assign(wk1,vor(1,1,j),phi(1,1,j),mx,mz,mxe,mze)
            enddo

            ! velocity boundary conditions: slave processes
            ! allocate temporary memory
            allocate(uWallBottom(mxe,mze))
            allocate(uWallTop(mxe,mze))
            allocate(vWallBottom(mxe,mze))
            allocate(vWallTop(mxe,mze))
            allocate(wWallBottom(mxe,mze))
            allocate(wWallTop(mxe,mze))

        endif ! end master and slave separation if

        ! all processes set velocity boundary condition
        call set_bc_from_restart_file(uWallBottom, uWallTop,vWallBottom, vWallTop,wWallBottom, wWallTop)
        ! and dealocate temporary buffers
        deallocate(uWallBottom)
        deallocate(uWallTop)
        deallocate(vWallBottom)
        deallocate(vWallTop)
        deallocate(wWallBottom)
        deallocate(wWallTop)

        ! so far, each processor has data in a few y planes
        ! before going any further, exchange cuts in y with cuts in z!
        call chikj2jik(vor,vor,wk1,wk1,myid)
        call chikj2jik(phi,phi,wk1,wk1,myid)

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


end module restart_file