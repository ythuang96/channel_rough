module restart_file
    use types, only: sp, dp
    use mpi
    use hdf5
    implicit none
    private
    include 'ctes3D'

    public :: read_restart_file_old, write_restart_file_old

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
    subroutine read_restart_file_old(vor,phi,u00,w00,wk1,myid, rf0u,rf0w,u00wk,w00wk,hv,hg,phiwk,vorwk)
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

        real(kind=sp), dimension(0:2*my-1,0:mx1,kb:ke), intent(inout) :: &
            hv,hg,phiwk,vorwk
        real(kind=dp), dimension(my), intent(inout) :: rf0u,rf0w,u00wk,w00wk

        master = 0

        ! zero everything
        do k=1,mz
            do j=1,mmy
                do i=1,mx
                    vor(i,k,j) = 0.0_sp
                    phi(i,k,j) = 0.0_sp
                enddo
            enddo
        enddo

        do j=1,my
            u00(j) = 0.0_sp
            w00(j) = 0.0_sp
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

        endif ! end master and slave separation if

        ! so far, each processor has data in a few y planes
        ! before going any further, exchange cuts in y with cuts in z!
        call chikj2jik(vor,vor,wk1,wk1,myid)
        call chikj2jik(phi,phi,wk1,wk1,myid)

        ! Organize the data, mainly solve for v from phi and compute y derivatives
        call orgainize_data_from_restart_file(vor,phi,u00,w00,rf0u,rf0w,u00wk,w00wk,hv,hg,phiwk,vorwk)

    end subroutine


    ! Subroutine to perform the special first time step in the original code,
    ! Mainly solving for v from phi and compute y derivatives of the 00 modes.
    ! This is only used for reading the old version restart file
    subroutine orgainize_data_from_restart_file(vor,phi,u00,w00,rf0u,rf0w,u00wk,w00wk,hv,hg,phiwk,vorwk)
        use wall_roughness, only: set_wall_roughness, &
            uWallBottom, uWallTop, vWallBottom, vWallTop, wWallBottom, wWallTop
        implicit none

        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save   /point/

        complex*8     xalp, xbet
        real*4        alp2,bet2
        integer       iax,icx
        common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),iax(mx),icx(0:mz1)
        save   /wave/


        integer i,k,j,k1
        real*8 rk
        complex*8 bcb,bct

        real(kind=sp), dimension(0:2*my-1,0:mx1,kb:ke), intent(inout) :: &
            vor,phi,hv,hg,phiwk,vorwk
        real(kind=dp), dimension(my), intent(inout) :: u00,w00,rf0u,rf0w,u00wk,w00wk


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


    ! Write a old version restart file
    subroutine write_restart_file_old(write_time, phi, vor, dvordy, chwk, u00, w00, massu0, myid)
        use wall_roughness, only: uWallBottom, uWallTop, vWallBottom, vWallTop, wWallBottom, wWallTop
        implicit none

        real*4 Deltat,CFL,time,dtr,FixTimeStep
        common /tem/ Deltat,CFL,time,dtr,FixTimeStep
        save   /tem/

        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save   /point/

        real*4 uampl,vampl,wampl,vspeed
        integer mxwall,mzwall
        common /boundary/ uampl,vampl,wampl,vspeed,mxwall,mzwall
        save /boundary/

        integer iinp,iout,id22,isn,ispf
        character*100 filinp,filout,filstt
        common /ficheros/ iinp,iout,id22,isn,ispf,filinp,filout,filstt
        save /ficheros/

        integer istat(MPI_STATUS_SIZE),ierr

        integer myid,iproc,leng,leng1,leng2,istep
        real*4 phi(0:2*my-1,0:mx1,kb:ke),vor(0:2*my-1,0:mx1,kb:ke)
        real*4 dvordy(0:2*my-1,0:mx1,kb:ke)
        real*4 chwk(*),work(20*my)

        real*8 u00(0:*),w00(0:*)
        real*8 massu0

        real*8 write_time

        real(kind=sp), dimension(1:mx, 1:mz) :: &
            uWallBottom_temp, uWallTop_temp, &
            vWallBottom_temp, vWallTop_temp, &
            wWallBottom_temp, wWallTop_temp

        integer j

        ! covert complex array to real array with twice the dimension for kx
        call c2r_array( uWallBottom, uWallBottom_temp)
        call c2r_array( vWallBottom, vWallBottom_temp)
        call c2r_array( wWallBottom, wWallBottom_temp)
        call c2r_array( uWallTop   , uWallTop_temp   )
        call c2r_array( vWallTop   , vWallTop_temp   )
        call c2r_array( wWallTop   , wWallTop_temp   )

        if (myid.eq.0) then
            write_time = -MPI_WTIME()
            write(*,*) time, vWallBottom(mxwall,mzwall),vWallTop(mxwall, mzwall)
        endif

        ! Change from slices of kx-y planes to slices of kx-kz planes
        call chjik2ikj(phi,phi,dvordy,dvordy,myid)
        call chjik2ikj(vor,vor,dvordy,dvordy,myid)

        do j=0,mx*max(mmy*mz,mmz*my)
            chwk(1+j)      = vor(j,0,kb)
            dvordy(j,0,kb) = phi(j,0,kb)
        enddo

        if(myid.ne.0) then
            ! everybody sends data to the master
            call MPI_SEND(  chwk,mx*mz*mmy,MPI_REAL,0,myid,MPI_COMM_WORLD,ierr)
            call MPI_SEND(dvordy,mx*mz*mmy,MPI_REAL,0,myid,MPI_COMM_WORLD,ierr)
        else
            ! the master first writes its stuff
            call escru(chwk,dvordy,u00,w00,jb,je,0, &
                uWallBottom_temp,uWallTop_temp, &
                vWallBottom_temp,vWallTop_temp, &
                wWallBottom_temp,wWallTop_temp,massu0)

            ! then receives everything from everybody
            do iproc=1,numerop-1

                leng=(jend(iproc)-jbeg(iproc)+1)*mx*mz
                call MPI_RECV(  chwk,leng,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                call MPI_RECV(dvordy,leng,MPI_REAL,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                ! and writes it to file
                call escru(chwk,dvordy,u00,w00,jbeg(iproc), &
                    jend(iproc),iproc, &
                    uWallBottom_temp,uWallTop_temp, &
                    vWallBottom_temp,vWallTop_temp, &
                    wWallBottom_temp,wWallTop_temp,massu0)
            enddo
        endif

        ! Change back to slices of kx-y planes
        call chikj2jik(phi,phi,dvordy,dvordy,myid)
        call chikj2jik(vor,vor,dvordy,dvordy,myid)

        ! Display write timer
        if (myid.eq.0) then
            write(*,*) 'time write:',MPI_WTIME()+write_time
            id22 = id22+1
        endif
    end subroutine


    ! Convert a complex array to a real array with twice the size for dimension 1
    ! Used in writing the old version restart file
    subroutine c2r_array( array_in, array_out )
        implicit none

        complex(kind=sp), dimension(mx/2, mz), intent(in) :: array_in
        real(kind=sp), dimension(mx, mz), intent(out) :: array_out

        array_out(1:mx-1:2,:) =  real( array_in, sp )
        array_out(2:mx:2  ,:) = aimag( array_in)

    end subroutine


    ! Write to the restart file
    subroutine escru(vor,phi,u00,w00,j1,j2,iproc,uwallb,uwallt,vwallb,vwallt,wwallb,wwallt,uBulk)
        implicit none
        integer iproc,j1,j2
        real*4 vor(mx,mz,j1:j2),phi(mx,mz,j1:j2)
        real*8 u00(*),w00(*)
        real*4 uwallb(mx,mz),uwallt(mx,mz)
        real*4 vwallb(mx,mz),vwallt(mx,mz)
        real*4 wwallb(mx,mz),wwallt(mx,mz)
        real*8 uBulk

        integer i,j,k

        character*3 ext1
        character*4 ext2
        character*104 fnameima

        real*4 gamma
        integer imesh
        common /mesh/ gamma,imesh
        save /mesh/

        real*4 Deltat,CFL,time,dtr,FixTimeStep
        common /tem/ Deltat,CFL,time,dtr,FixTimeStep
        save /tem/

        real*4 uampl,vampl,wampl,vspeed
        integer mxwall,mzwall
        common /boundary/ uampl,vampl,wampl,vspeed,mxwall,mzwall
        save /boundary/

        real*8 fmap,y2
        real*4  Re,alp,bet,a0,y,hy
        common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
        save   /fis/

        integer iinp,iout,id22,isn,ispf
        character*100 filinp,filout,filstt
        common /ficheros/ iinp,iout,id22,isn,ispf,filinp,filout,filstt
        save /ficheros/

        integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
        common /point /jbeg(0:numerop-1),jend(0:numerop-1), &
                       kbeg(0:numerop-1),kend(0:numerop-1), &
                       jb,je,kb,ke,mmy,mmz
        save /point/


        if(iproc.eq.0) then    !!!!! coming from node 0
            ! start writing image
            if(id22.gt.999.or.id22.lt.0) then
                write(*,*) 'number of images out of range'
                stop
            endif

            write(ext1,'(i3.3)') id22
            fnameima=filout(1:index(filout,' ')-1)//'.'//ext1

            open (iout,file=fnameima,status='unknown', &
                form='unformatted',access='direct',recl=mx*mz*2*iwd)

            write(*,*) j1,j2, 'in escru image',fnameima
            write(iout,rec=1) time,Re,alp,bet,a0,mx,my,mz,imesh,gamma, &
                mxwall,mzwall,uampl,vampl,wampl,vspeed, &
                (real(u00(j)),real(w00(j)),j=1,my)

            do j=j1,j2
                write(iout,rec=j+1)((vor(i,k,j),phi(i,k,j),i=1,mx),k=1,mz)
            enddo


        else
            write(*,*) j1,j2, 'in escru image'

            do j=j1,j2
                write(iout,rec=j+1)((vor(i,k,j),phi(i,k,j),i=1,mx),k=1,mz)
            enddo

        endif

        if(iproc.eq.numerop-1) then
            ! last node writes boundaries
            write(iout,rec=my+2)((uwallb(i,k),uwallt(i,k),i=1,mx),k=1,mz)
            write(iout,rec=my+3)((vwallb(i,k),vwallt(i,k),i=1,mx),k=1,mz)
            write(iout,rec=my+4)((wwallb(i,k),wwallt(i,k),i=1,mx),k=1,mz)

         write(*,*) 'closing files'
         close(iout)
      endif

    end subroutine


end module restart_file