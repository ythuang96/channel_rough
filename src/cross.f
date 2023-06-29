c/*............................................c/*
c/*  MODULE FOR MPI SP2                        c/*
c/*............................................c/*
c/*............................................c/*

      subroutine cross1(vor,
     .     phi,
     .     u00,w00,
     .     rf0u,rf0w,u00wk,w00wk,
     .     hv,
     .     hg,
     .     phiwk,
     .     vorwk,
     .     dvordy,
     .     chwk,
     .     myid)
      use boundary_planes, only: uWallBottom, uWallTop, vWallBottom,
     &  vWallTop, wWallBottom, wWallTop
      use wall_roughness, only: set_wall_roughness
      use save_flowfield, only: collectFlowfield, write_h5,
     &  assess_whether_to_collect_flowfield_step,
     &  assess_whether_to_collect_flowfield_time
      implicit none
      include "mpif.h"
      include "ctes3D"

      real*4  c1,c2,c3
      parameter(c1=  1./3. ,c2=  1./2. ,c3=  1.  )

      integer myid,iproc,leng,leng1,leng2,istep
      integer irun,rkstep,i,ii,k,kk,j,jj,i1,k1,k2,
     .     ipo1,ipo2,ipo3,ipo4,ipo
      real*4  r1, dtr1, du0,duL,tmp
      real*8 rk,rk2,bcbr,bctr
      complex*8 bcb,bct,bcbdv,bctdv,bcbo,bcto
      real*4  reynota
      real*8  H00u,H00w
      real*8  massu0,massw0,massu,massw
      real*8  massu1,massu2,massw1,massw2,cteu,ctew

      integer istat(MPI_STATUS_SIZE),ierr

      real*4  Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/

      real*4 Deltat,CFL,time,dtr,FixTimeStep
      common /tem/ Deltat,CFL,time,dtr,FixTimeStep
      save   /tem/

      integer nimag,nstep,nhist,ihist,icfl
      common /timacc/ nimag,nstep,nhist,ihist,icfl
      save   /timacc/

      real*8  fmap,y2
      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/

      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .     kbeg(0:numerop-1),kend(0:numerop-1),
     .     jb,je,kb,ke,mmy,mmz
      save   /point/

      integer iinp,iout,id22,isn,ispf
      character*100 filinp,filout,filstt
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .     filinp,filout,filstt
      save /ficheros/
      character*104 fname


      real*4 hg (0:2*my-1,0:mx1,kb:ke),vorwk(0:2*my-1,0:mx1,kb:ke),
     .     hv (0:2*my-1,0:mx1,kb:ke),phiwk(0:2*my-1,0:mx1,kb:ke),
     .     phi(0:2*my-1,0:mx1,kb:ke),vor  (0:2*my-1,0:mx1,kb:ke),
     .     dvordy(0:2*my-1,0:mx1,kb:ke)
      real*4 chwk(*),work(20*my)

      real*8 u00(0:*),w00(0:*)
      real*8 rf0u(0:*),u00wk(0:*),
     .     rf0w(0:*),w00wk(0:*)

      complex*8     xalp, xbet
      real*4        alp2,bet2
      integer       iax,icx
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .     iax(mx),icx(0:mz1)
      save   /wave/
      real*4 dkx2,dkz2
      complex*8 dkx,dkz,dkk

      real*4 c(3)
      data c/c1,c2,c3/
      save c

      real*4 uampl,vampl,wampl,vspeed
      integer mxwall,mzwall
      common /boundary/ uampl,vampl,wampl,vspeed,mxwall,mzwall
      save /boundary/

      real*4 u1r,u2r,u3r,o1r,o2r,o3r
      common /planes/
     .     u1r(mgalx+2,mgalz),u2r(mgalx+2,mgalz),u3r(mgalx+2,mgalz),
     .     o1r(mgalx+2,mgalz),o2r(mgalx+2,mgalz),o3r(mgalx+2,mgalz)
      save/planes/

      real*4 trp
      real*8 trp2
      common /mass/ trp(0:my1),trp2(0:my1)
      save /mass/


      character*3 ext,ext1

      real*8 commtimer,transtimer,totaltimer
      common/timers/ commtimer, transtimer, totaltimer
      save/timers/
      real*8 iter_time,write_time,laps_time,comm_time

      comm_time = 0D0
      commtimer=0.0D0
      transtimer=0.0D0
      totaltimer=0.0D0

      ihist    = 0

      irun   = 0                ! first time step is special in tim3rkp
      icfl   = 1                ! first time step always needs a step size

      if(myid.eq.0) then
         write(ext1,'(i3.3)') id22
         fname=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.cf'
         write (*,*) fname
         open(39,file=fname,status='unknown')
      endif


c/********************************************************************/
c/*     THIS IS THE TIME LOOP                                        */
c/********************************************************************/
      do istep=1,nstep

         if (myid.eq.0) then
            totaltimer = totaltimer-MPI_WTIME()
            iter_time=-MPI_WTIME()
         endif

         if (mod(istep-1,nhist) .eq. 0) then
            ihist=1
            icfl= 1
         endif

c     -----------------------------------------------------------------
c     Assess whether it to save flow field or not
         ! call assess_whether_to_collect_flowfield_time(time)
         call assess_whether_to_collect_flowfield_step(istep)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     write restart file
         if (mod(istep-1,nimag) .eq. 0 .and. istep.ne.1) then

            if (myid.eq.0) then
               write_time = -MPI_WTIME()
               write(*,*) time, vWallBottom(mxwall,mzwall),
     &             vWallTop(mxwall, mzwall)
            endif

            ! procedure to receive all the slices
            call chjik2ikj(phi,phi,dvordy,dvordy,myid)
            call chjik2ikj(vor,vor,dvordy,dvordy,myid)

            do j=0,mx*max(mmy*mz,mmz*my)
               chwk(1+j)        = vor(j,0,kb)
               dvordy(j,0,kb) = phi(j,0,kb)
            enddo

            if(myid.ne.0) then
               ! everybody sends data to the master
               call MPI_SEND(chwk,mx*mz*mmy,MPI_REAL,
     &              0,myid,MPI_COMM_WORLD,ierr)
               call MPI_SEND(dvordy,mx*mz*mmy,MPI_REAL,
     &              0,myid,MPI_COMM_WORLD,ierr)
            else
               ! the master first writes its stuff
               call escru(chwk,dvordy,u00,w00,jb,je,0,
     .              uWallBottom,uWallTop,vWallBottom,vWallTop,
     .              wWallBottom,wWallTop,massu0)

               ! then receives everything from everybody
               do iproc=1,numerop-1

                  leng=(jend(iproc)-jbeg(iproc)+1)*mx*mz
                  call MPI_RECV(chwk,leng,MPI_REAL,
     &                 iproc,iproc,MPI_COMM_WORLD,istat,ierr)

                  call MPI_RECV(dvordy,leng,MPI_REAL,
     &                 iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                  ! and writes it
                  call escru(chwk,dvordy,u00,w00,jbeg(iproc),
     &                 jend(iproc),iproc,
     &                 uWallBottom,uWallTop,vWallBottom,vWallTop,
     &                 wWallBottom,wWallTop,massu0)
               enddo
            endif


            call chikj2jik(phi,phi,dvordy,dvordy,myid)
            call chikj2jik(vor,vor,dvordy,dvordy,myid)

            if (myid.eq.0) then
               write(*,*) 'time write:',MPI_WTIME()+write_time
               id22 = id22+1
            endif

         endif
c     finished writing image
c     -----------------------------------------------------------------


c/********************************************************************/
c/*     Special treatment for first time step                        */
c/********************************************************************/
         if(irun.eq.0) then

            irun = 1

c     -----------------------------------------------------------------
           ! all processors have explicitly set the time to zero, so
           ! this step can be skipped
           ! call  MPI_BCAST(time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
           ! first step: set boundary conditions for v
           ! vWall is read from restart file, but let's set it
           ! explicitly here, just in case we want to do something
           ! else than in the restart file
           call set_wall_roughness(uWallBottom, uWallTop, vWallBottom,
     &       vWallTop)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Calculate the v from phi
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
                  call Lapvdv(phi(0,i,k),hg(0,i,k),hv(0,i,k),
     .                 rk,bcb,bct)
                  ! hg  is the output: v
                  ! hv  is the output: dv/dy
               enddo
            enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute 00 modes of omega1 and omega3
            ! Compute y derivative of u and w
            do j=0,my1
               u00wk(j)=u00(j)
               w00wk(j)=w00(j)
            enddo
            call deryr(u00wk,rf0u,my)
            call deryr(w00wk,rf0w,my)
            ! rf0u = du/dy(kx=kz=0) = - o3(kx=kz=0)
            ! rf0w = dw/dy(kx=kz=0) = + o1(kx=kz=0)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Copy data to vorwk and phiwk
            do  k=kb,ke
               do i=0,mx1
                  do j=0,2*my-1
                     vorwk(j,i,k)=vor(j,i,k)
                     phiwk(j,i,k)=phi(j,i,k)
                  enddo
               enddo
            enddo
            ! vorwk and phi wk are size 0:2*my-1, 0:mx1, kb:ke
            ! the first index includes both real and imaginary part
            ! vorwk: omega2
            ! phiwk: phi = nabla^2 v
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     compute u mass flux
            ! Integrate the channel
            massu0=0d0
            massw0=0d0
            do j=0,my1
               massu0 = massu0 + trp2(j)*u00(j)
            enddo

            ! target mass flux, historical Madrid value
            massu0 = .8987636566162d0
            if(myid.eq.0) then
               write (*,*) 'mass flux:',massu0,massw0
            endif
c     -----------------------------------------------------------------


         endif
c/********************************************************************/
c/*     End special first step                                       */
c/********************************************************************/


c/********************************************************************/
c/*     Runge-Kutta third order                                      */
c/********************************************************************/
         do rkstep=1,3


c     -----------------------------------------------------------------
c     Compute y derivative of omega2
            call deryr2(vorwk,dvordy,(mx1+1)*mmz,my,chwk)
            ! dvordy = d omega2 / dy
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     All arrays changed from kx-y planes to kx-kz planes
            call chjik2ikj(phiwk ,phiwk ,chwk,chwk,myid)
            call chjik2ikj(hv    ,hv    ,chwk,chwk,myid)
            call chjik2ikj(hg    ,hg    ,chwk,chwk,myid)
            call chjik2ikj(dvordy,dvordy,chwk,chwk,myid)
            call chjik2ikj(vorwk ,vorwk ,chwk,chwk,myid)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Storaging 0 modes in work
            ipo1 = 1    + my
            ipo2 = ipo1 + my
            ipo3 = ipo2 + my
            ipo4 = ipo3 + my

            ! work is a column vector of the kx=kz=0 modes
            ! 1st one is: du/dy(00) = - omega3(00)
            ! 2nd one is: dw/dy(00) = + omega1(00)
            ! 3rd one is:     u(00)
            ! 4th one is:     w(00)
            do j=0,my1
               work(   1+j)=rf0u (j)
               work(ipo1+j)=rf0w (j)
               work(ipo2+j)=u00wk(j)
               work(ipo3+j)=w00wk(j)
            enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute non-linear terms
            ! phiwk, vorwk, hv, hg, dvordy, work are inputs
            ! u1r, u2r ... are all work variables
            call hvhg(phiwk,vorwk,hv,hg,rf0u,
     .           rf0w,dvordy,work,myid,rkstep,
     .           u1r,u2r,u3r,o1r,o2r,o3r,
     .           u1r,u2r,u3r,o1r,o2r,o3r)
            ! Ouputs:
            ! hv    : (d^2/dx^2 + d^2/dz^2) H2
            ! hg    : - d H3/dx + d H1/dz
            ! rf0u  : H1(kx=kz=0)
            ! rf0w  : H3(kx=kz=0)
            ! dvordy: + d H1/dx + d H3/dz
            !
            ! The 3 matrices are kx-kz-y
            ! with each processor containing a few kx-kz planes
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Mean subtraction for rf0u and rf0w
            ! Integrate in y
            H00u=0d0
            H00w=0d0
            do j=0,my1
               H00u = H00u + trp2(j)*rf0u(j)
               H00w = H00w + trp2(j)*rf0w(j)
            enddo
            ! Mean subtract
            do j=0,my1
               rf0u(j) = rf0u(j)-H00u
               rf0w(j) = rf0w(j)-H00w
            enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     All arrays changed from kx-kz planes to kx-y planes
            call chikj2jik(dvordy,dvordy,chwk,chwk,myid)
            call chikj2jik(hv    ,hv    ,chwk,chwk,myid)
            call chikj2jik(hg    ,hg    ,chwk,chwk,myid)
            ! hv    : (d/dx^2 + d/dz^2) H2
            ! hg    : - d H3/dx + d H1/dz
            ! dvordy: + d H1/dx + d H3/dz
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute y derivative of dvordy
            ! before: dvordy = + d H1/dx + d H3/dz
            call deryr2(dvordy,dvordy,(mx1+1)*mmz,my,chwk)
            ! now:    dvordy = d/dy ( + d H1/dx + d H3/dz )
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute v velocity forcing hv
            ! before: hv = (d/dx^2 + d/dz^2) H2
            do k=kb,ke
               do i=0,mx1
                  do j=0,2*my-1
                     hv(j,i,k) = hv(j,i,k) - dvordy(j,i,k)
                  enddo
               enddo
            enddo
            ! now hv = (d/dx^2 + d/dz^2) H2 - d(dH1/dx + dH3/dz)/dy
            ! this is the complete v velocity forcing
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute the RHS forcing for u00, w00, omega3, phi
            ! gamma_s * Deltat
            r1=c(rkstep)*Deltat
            ! Compute RHS forcing for u00 and w00
            ! Note: they are missing a scalar multiplier
            do j=0,my1
               u00wk(j)=u00(j)+r1*rf0u(j)
               w00wk(j)=w00(j)+r1*rf0w(j)
            enddo
            ! Compute RHS forcing for omega2 and phi
            ! Note they are missing a scalar multiplier
            do k=kb,ke
               do i=0,mx1
                  do j=0,2*my-1
                     vorwk(j,i,k)=vor(j,i,k)+r1*hg(j,i,k)
                     phiwk(j,i,k)=phi(j,i,k)+r1*hv(j,i,k)
                  enddo
               enddo
            enddo

            ! Re/ ( gamma_s * Deltat )
            dtr1=dtr/c(rkstep)
            ! Apply the scalar multiplier
            do j=0,my1
               rf0u(j)=-dtr1*u00wk(j)
               rf0w(j)=-dtr1*w00wk(j)
            enddo
            ! Apply the scalar multiplier
            do k=kb,ke
               do i=0,mx1
                  do j=0,2*my-1
                     hg(j,i,k)=-dtr1*vorwk(j,i,k)
                     hv(j,i,k)=-dtr1*phiwk(j,i,k)
                  enddo
               enddo
            enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Update boundary condition for the viscous time step
            call set_wall_roughness(uWallBottom, uWallTop, vWallBottom,
     &          vWallTop)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Advance phi, v, dv/dy, omega2
            do k=kb,ke ! Loop over the assigned kz of this processor
               k1 = icx(k-1)
               do i=0,mx1 ! Loop over all kx

                  rk = bet2(k1)+alp2(i) ! kx^2 + kz^2
                  rk2 = dtr1+rk ! kx^2 + kz^2 + Re/( gamma_s * Deltat )

                  ! boundary condition for rough wall
                  ! BC for omega2 from BC of u and w
                  bcbo = -xalp(i)*wWallBottom(i,k-1) +
     &              xbet(k-1)*uWallBottom(i,k-1)
                  bcto = -xalp(i)*wWallTop(i,k-1) +
     &              xbet(k-1)*uWallTop(i,k-1)
                  ! BC for v
                  bcb = vWallBottom(i,k-1)
                  bct = vWallTop(i,k-1)
                  ! BC for dv/dy from u and w using continuity
                  bcbdv = -xalp(i)*uWallBottom(i,k-1) -
     &              xbet(k-1)*wWallBottom(i,k-1)
                  bctdv = -xalp(i)*uWallTop(i,k-1) -
     &              xbet(k-1)*wWallTop(i,k-1)

                  ! temporary: no slip bc
                  !bcbo  = (0.0, 0.0)
                  !bcto  = (0.0, 0.0)
                  !bcb   = (0.0, 0.0)
                  !bct   = (0.0, 0.0)
                  !bcbdv = (0.0, 0.0)
                  !bctdv = (0.0, 0.0)

                  ! Solve for phi, v, dv/dy, omega2 at the next RK step
                  ! for each kx, kz wavenumbers
                  call lapsov(phiwk(0,i,k),hg(0,i,k),hv(0,i,k),
     &                 hv(0,i,k),vorwk(0,i,k),hg(0,i,k),
     &                 rk2,rk,bcb,bct,bcbdv,bctdv,bcbo,bcto)
                  ! Outputs:
                  ! phiwk: phi = nabla^2 v
                  ! hg   : v
                  ! hv   : dv/dy
                  ! vorwk: omega2
                  ! All 4 matrices are size 0:2*my-1, 0:mx1, kb:ke
               enddo ! end loop over kx
            enddo ! end loop over kz

c     Advance phi, v, dv/dy, omega2 complete
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Advance kx = kz = 0 modes of u and w
            ! Velocity ​​without pressure-mass correction
            ! Solve [d/dt - 1/Re * d^2/dy^2] u = H1(00)
            bcbr = 0d0
            bctr = 0d0
            rk = dtr1
            call Lapv1(rf0u,u00wk,rk,bcbr,bctr)
            call Lapv1(rf0w,w00wk,rk,bcbr,bctr)
            ! The solution is u00wk and w00wk

            ! Velocity with constant RHS forcing
            ! Solve [d/dt - 1/Re * d^2/dy^2] u = constant
            bcbr = 0d0
            bctr = 0d0
            do j=0,my1
               rf0u(j) = -dtr1
               rf0w(j) = -dtr1
            enddo
            rk = dtr1
            call Lapv1(rf0u,rf0u,rk,bcbr,bctr)
            call Lapv1(rf0w,rf0w,rk,bcbr,bctr)
            ! The solution is rf0u and rf0w

            ! Compute the mass flux for the two solutions
            ! integrate in y
            massu1=0d0
            massu2=0d0
            massw1=0d0
            massw2=0d0
            do j=0,my1
               massu1 = massu1 + trp2(j)*u00wk(j)
               massw1 = massw1 + trp2(j)*w00wk(j)
               massu2 = massu2 + trp2(j)*rf0u(j)
               massw2 = massw2 + trp2(j)*rf0w(j)
            enddo

            ! Linearly conbine the two solutions
            ! to get constant mass flux
            ! massu0, massw0 is the target mass flux
            cteu =(massu0 - massu1)/massu2
            ctew =(massw0 - massw1)/massw2
            do j=0,my1
               u00wk(j) = u00wk(j) + cteu*rf0u(j)
               w00wk(j) = w00wk(j) + ctew*rf0w(j)
            enddo

            ! Compute the derivatives of u00 and w00
            call deryr(u00wk,rf0u,my)
            call deryr(w00wk,rf0w,my)

            ! u00wk: the kx = kz = 0 u mode
            ! w00wk: the kx = kz = 0 w mode
            ! rf0u : d/dy of u00
            ! rf0w : d/dy of w00
c     Advance kx = kz = 0 modes of u and w complete
c     -----------------------------------------------------------------


         ENDDO
c/********************************************************************/
c/*     Runge-Kutta third order Complete                             */
c/********************************************************************/


c     -----------------------------------------------------------------
c     Write Flowfield from buffer to h5 file
         if (collectFlowfield) then
           ! data collected at the previous timestep is written at the
           ! subsequent timestep. time-Deltat corresponds to the time
           ! at which the data was collected
           if (myid .eq. 0) then
326          format(a22,i5,a7,(d14.6))
             write(*,326) '    Saving Data, Step ', istep, ', Time ',
     &           time
           endif
           call write_h5(alp, bet, y, xalp, xbet, time, Re, massu0)
         endif
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Copy variables for the new time step
         ! 00 modes for u and w
         do j=0,my1
            u00(j)=u00wk(j)
            w00(j)=w00wk(j)
         enddo
         ! phi and omega3 dat
         do k=kb,ke
            do i=0,mx1
               do j= 0,2*my-1
                  vor(j,i,k)=vorwk(j,i,k)
                  phi(j,i,k)=phiwk(j,i,k)
               enddo
            enddo
         enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     WxL, WzL, uvL
         ! WxL: plane averaged omega1 at top wall, also = + d/dy w00
         ! WzL: plane averaged omega3 at top wall, also = - d/dy u00
         ! uvL: plane averaged Re stress at top wall
         if (ihist.eq.1) then
            ! Last processor send WxL, WzL, uvL to master
            if (myid.eq.numerop-1) then
               chwk (1) = WxL
               chwk (2) = WzL
               chwk (3) = uvL
               call MPI_SEND(chwk,3,MPI_REAL,
     &              0,0,MPI_COMM_WORLD,ierr)
            endif
            ! Master recieves WxL, WzL, uvL
            if (myid.eq.0) then
               call MPI_RECV(chwk,3,MPI_REAL,
     &              numerop-1,0,MPI_COMM_WORLD,istat,ierr)
               WxL =  chwk(1)
               WzL =  chwk(2)
               uvL =  chwk(3)
            endif
         endif
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Master writes history record
         if(myid.eq.0) then

            reynota=0.5*re*re*(abs(wz0/re+uv0)
     .           +abs(wzl/re+uvL))

            ! mass flux in u and w
            massu = massu1 +cteu*massu2
            massw = massw1 +ctew*massw2

            ! write histroy record
            if (ihist.eq.1) then
               tmp=my1/2

               ! write to the output txt file
 325           format(a10,i5,9(d14.6))
               write(*,325) 'Time Step ',
     .              istep,time,-1.*Wz0,WzL,sqrt(reynota),Deltat,
     .              u00(floor(tmp))*Re/sqrt(reynota),massw,uv0,uvL
               ! write to cf file
 324           format(i5,9(d22.14))
               write(39,324)
     .              istep,time,-1.*Wz0,WzL,sqrt(reynota),Deltat,
     .              u00(floor(tmp))*Re/sqrt(reynota),massw,uv0,uvL
               call flush(39)
            endif

            ! switch to new .cf file
            if (mod(istep-1,nimag).eq.0 .and.
     &           istep.ne.1 .and. istep.ne.nstep) then

               close(39)
               write(ext1,'(i3.3)') id22
               fname=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.cf'
               write (*,*) fname
               open(39,file=fname,status='unknown')
            endif

         endif
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Increment time
         time=time+Deltat
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Reset icfl, ihist
         if(icfl.eq.1)   icfl   = 0
         if(ihist.eq.1)  ihist  = 0
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Master write step time
         if (myid.eq.0) then
            totaltimer=totaltimer+MPI_WTIME()
c            write(*,'(i7,3f20.5)') istep,
c     >           MPI_WTIME()+iter_time-commtimer+comm_time,
c     >           commtimer-comm_time,MPI_WTIME()+iter_time
            comm_time = commtimer
         end if
c     -----------------------------------------------------------------

      ENDDO
c/********************************************************************/
c/*     TIME LOOP COMPLETE                                           */
c/********************************************************************/


c     -----------------------------------------------------------------
c     Master writes data about the total computation time
      if (myid.eq.0) then
         print *,"Total time: ",totaltimer
         print *,"Trans. time: ",transtimer
         print *,"Comm. time: ",commtimer
         print *,"Comm/Total: ",commtimer/totaltimer
         print *,"Trans/Total: ",transtimer/totaltimer
         print *,"Aver time per step: ",totaltimer/nstep
      end if
c     -----------------------------------------------------------------

      end



c/********************************************************************/
c/*                                                                  */
c/*         computes the forcing terms (convective terms)            */
c/*                                                                  */
c/*    input:                                                        */
c/*      phic : phi = nabla^2 v  (F-F-P)                             */
c/*      ome2c: omega2           (F-F-P)                             */
c/*      rhvc : dv/dy            (F-F-P)                             */
c/*      rhgc : v                (F-F-P)                             */
c/*      ome1c: d omega2/dy      (F-F-P)                             */
c/*      work2: kx = kz = 0 modes of omega3, omega1, u, w            */
c/*             stacked in a column vector                           */
c/*                                                                  */
c/*   output:                                                        */
c/*     rhvc : d^2H2/dx^2 + d^2H2/dz^2  (F-F-P)                      */
c/*     rhgc : -dH3/dx + dH1/dz         (F-F-P)                      */
c/*     rf0u : kx = kz = 0 mode of H1   (column vector)              */
c/*     rf0w : kx = kz = 0 mode of H3   (column vector)              */
c/*     ome1c:  dH1/dx + dH3/dz         (F-F-P)                      */
c/*                                                                  */
c/*                                                                  */
c/*      (F-F-P) indicates the matrix dimensions are                 */
c/*              Fourier x - Fourier z - Physical y                  */
c/*              and each processor contains data from a few         */
c/*              kx-kz planes.                                       */
c/*                                                                  */
c/********************************************************************/
      subroutine hvhg(phic,ome2c,rhvc,rhgc,
     .     rf0u,rf0w,ome1c,
     .     work2,myid,rkstep,
     .     u1r,u2r,u3r,o1r,o2r,o3r,
     .     u1c,u2c,u3c,o1c,o2c,o3c )
      use save_flowfield, only: collectFlowfield,
     &  save_plane_data_to_buffer
      implicit none
      include "ctes3D"
      include "mpif.h"

      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .     kbeg(0:numerop-1),kend(0:numerop-1),
     .     jb,je,kb,ke,mmy,mmz
      save   /point/

      real*4 Deltat,CFL,time,dtr,FixTimeStep
      common /tem/ Deltat,CFL,time,dtr,FixTimeStep
      save /tem/

      real*8  fmap,y2
      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/

      integer nimag,nstep,nhist,ihist,icfl
      common /timacc/ nimag,nstep,nhist,ihist,icfl
      save   /timacc/

      integer iax,icx
      real alp2,bet2
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     .     iax(mx),icx(0:mz1)
      save /wave/

      real*4  Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/

      complex*8 phic (0:mx1,0:mz1,jb:*),
     .     ome1c(0:mx1,0:mz1,jb:*),
     .     ome2c(0:mx1,0:mz1,jb:*),
     .     rhgc (0:mx1,0:mz1,jb:*),
     .     rhvc (0:mx1,0:mz1,jb:*)

      real*8 rf0u(*),rf0w(*)
      real*4 work2(*)

c---------------6 * (mgalx+2)  * mgalz planes


      real*4
     &     u1r(mgalx+2,mgalz),u2r(mgalx+2,mgalz),u3r(mgalx+2,mgalz),
     &     o1r(mgalx+2,mgalz),o2r(mgalx+2,mgalz),o3r(mgalx+2,mgalz)

      complex*8
     &     u1c(0:mx1,0:mz1),u2c(0:mx1,0:mz1),u3c(0:mx1,0:mz1),
     &     o1c(0:mx1,0:mz1),o2c(0:mx1,0:mz1),o3c(0:mx1,0:mz1)

      complex*8 dk
      real*4 dk2
      integer myid,rkstep
      integer ipo1,ipo2,ipo3
      integer i,j,k,jj,kk,jndex
      integer istat(MPI_STATUS_SIZE),ierr
      integer mmyr
      integer iproc,pproc
      integer pnodes

      real*4 cflx,cfly,cflz,hxalp,hzalp,hyy,cfl0,reigmx1
      integer uggg
      real*8 aa
      complex*16 cc

      real*4 temp

      integer ilocalu,ilocalw,icount

      ipo1 = 1    + my
      ipo2 = ipo1 + my
      ipo3 = ipo2 + my

c     -----------------------------------------------------------------
c     at this point:
c     rhv   is dv / dy
c     rhg   is v
c     phic  is nabla^2(v)
c     ome1c is d(omega_2)/dy
c     ome2  is omega_2
c     All Fourier x -Fourier z -Physical y
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
      cflx = 0.
      cfly = 0.
      cflz = 0.
      cfl0 = 0.

      hxalp=alp*mx*0.5
      hzalp=bet*mz*0.5
c     -----------------------------------------------------------------


c     Loop over y planes
      DO J = JB,JE

c     -----------------------------------------------------------------
c     computes omega_1 (o1c) and omega_3 (o3c)
c     except kx = 0 and kz = 0 modes
         ! kx = 0 modes
         do k=1,mz1
            dk  = 1.0/xbet(k)
            o3c(0,k) = -ome1c(0,k,j) * dk
            o1c(0,k) = - phic(0,k,j) * dk
         enddo
         ! All other modes
         do k=0,mz1
            dk  = xbet(k)
            dk2 = bet2(k)
            do i=1,mx1

               temp = 1.0/(alp2(i) + dk2)

               o3c(i,k) = ( ome1c(i,k,j)*dk-phic(i,k,j)*xalp(i))
     &              *temp
               o1c(i,k) = ( ome1c(i,k,j)*xalp(i) + phic(i,k,j)*dk)
     &              *temp
            enddo
         enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     computes u1c (u) and u3c (w)
c     except kx = 0 and kz = 0 modes
         ! kx = 0 modes
         do k=1,mz1
            dk  = 1.0/xbet(k)
            u3c(0,k) = -rhvc (0,k,j) * dk
            u1c(0,k) = ome2c (0,k,j) * dk
         enddo
         ! All other modes
         do k=0,mz1
            dk  = xbet(k)
            dk2 = bet2(k)
            do i=1,mx1

               temp = 1.0/(alp2(i) + dk2)

               u3c(i,k) = ( rhvc(i,k,j)*dk + ome2c(i,k,j)*xalp(i) )
     &              *temp
               u1c(i,k) = ( rhvc(i,k,j)*xalp(i) - ome2c(i,k,j)*dk )
     &              *temp
            enddo
         enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     00 modes of -omega3,omega1,u,w unpacked from input work2
c     They are stacked in a column vector
         jj = j-1
         o3c(0,0) = -work2(      j)
         o1c(0,0) =  work2(ipo1+jj)
         u1c(0,0) =  work2(ipo2+jj)
         u3c(0,0) =  work2(ipo3+jj)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     copy v and omega2 from the input
         do k=0,mz1
            do i=0,mx1
               u2c(i,k) = rhgc(i,k,j)
               o2c(i,k) = ome2c(i,k,j) !!! warning, ome2c(0,0,j) must be zero
            enddo
         enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     at this point:
c     u1c is u
c     u2c is v
c     u3c is w
c     o1c is omega_1
c     o2c is omega_2
c     o3c is omega_3
c     All variables in Fourier x - Fourier z for one x-z plane
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     compute vorticity & Re stress at walls
         if (rkstep.eq.1) then
            ! Bottom Wall
            if (j.eq.1) then
               Wx0 = o1c(0,0)
               Wz0 = o3c(0,0)
            endif
            ! Top Wall
            if (j.eq.my) then
               WxL = o1c(0,0)
               WzL = o3c(0,0)
            endif
         endif
         ! Wx0, Wz0, WxL, WzL are global variables
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     save velocity and vorticity fields if required
         if (collectFlowfield .and. rkstep==1) then
            call save_plane_data_to_buffer(
     &           u1c, u2c, u3c, o1c, o2c, o3c, j )
         endif
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Move everything to Physical x - Physical z

         ! substract umax/2 to u00 to increase dt !!!!
         u1c(0,0) = u1c(0,0) - a0
         ! IFFT2
         call fourxz(u1c,u1r,1,1) ! u
         call fourxz(u2c,u2r,1,1) ! v
         call fourxz(u3c,u3r,1,1) ! w
         call fourxz(o1c,o1r,1,1) !omega_1
         call fourxz(o2c,o2r,1,1) !omega_2
         call fourxz(o3c,o3r,1,1) !omega_3
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
         if (rkstep.eq.1) then
            ! Compute uv Stress at the bottom wall
            if (j.eq.1) then
               uv0=0.
               do k=1,mgalz
                  do i=1,mgalx
                     uv0 = uv0 + u1r(i,k)*u2r(i,k)
                  enddo
               enddo
               uv0= uv0/(mgalz*mgalx)
            endif
            ! Compute uv Stress at the top wall
            if (j.eq.my) then
               uvL=0.
               do k=1,mgalz
                  do i=1,mgalx
                     uvL = uvL + u1r(i,k)*u2r(i,k)
                  enddo
               enddo
               uvL= uvL/(mgalz*mgalx)
            endif
         endif
         ! uv0 and uvL are global variables
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     CFL Condition
         uggg=0

         if (icfl.eq.1.and.rkstep.eq.1) then
            ! hyy is y grid spacing at the current x-z plane
            hyy = hy(j)
            do k=1,mgalz
               do i=1,mgalx
                  cflx = max(cflx,abs(u1r(i,k))*hxalp )
                  cfly = max(cfly,abs(u2r(i,k))/hyy )
                  cflz = max(cflz,abs(u3r(i,k))*hzalp )
               enddo
            enddo

         endif
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     computes H = u X omega
c     u2r = H1 = v.omega_3 - w.omega_2
c     u3r = H2 = w.omega_1 - u.omega_3
c     u1r = H3 = u.omega_2 - v.omega_1
c     All 3 variables are physical x - physical z for one x-z plane
         do k=1,mgalz
            do i=1,mgalx
               aa = u2r(i,k)
               u2r(i,k) = u2r(i,k)*o3r(i,k)-u3r(i,k)*
     &              o2r(i,k)
               u3r(i,k) = u3r(i,k)*o1r(i,k)-u1r(i,k)*
     &              o3r(i,k)
               u1r(i,k) = u1r(i,k)*o2r(i,k)-aa*
     &              o1r(i,k)
            enddo
         enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     FFT2 of H to Fourier x - Fourier z
         ! H3
         call fourxz(u1c,u1r,-1,1)
         ! H1
         call fourxz(u2c,u2r,-1,1)
         ! H2
         call fourxz(u3c,u3r,-1,1)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     saves coefficients for kx = kz = 0 modes of H1 and H3
         ! H1
         rf0u(j)=real(u2c(0,0))
         ! H3
         rf0w(j)=real(u1c(0,0))
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute derivatives of H1, H2, H3
         ! o1 =   dH1/dx + dH3/dz
         ! u2 = - dH3/dx + dH1/dz
         do k=0,mz1
            do i=0,mx1
               o1c(i,k) =  xalp(i)*u2c(i,k)+xbet(k)*u1c(i,k)
               u2c(i,k) = -xalp(i)*u1c(i,k)+xbet(k)*u2c(i,k)
            enddo
         enddo
         ! u1 = d^2 H2/dx^2 + d^H2/dz^2
         do k=0,mz1
            do i=0,mx1
               u1c(i,k) = - u3c(i,k)*(alp2(i)+bet2(k))
            enddo
         enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Copy planes into outputs
         ! rhvc  = d^2 H2/dx^2 + d^H2/dz^2
         ! rhgc  = - dH3/dx + dH1/dz
         ! ome1c =   dH1/dx + dH3/dz
         do k=0,mz1
            do i=0,mx1
               rhvc(i,k,j)=u1c(i,k)
               rhgc(i,k,j)=u2c(i,k)
               ome1c(i,k,j)=o1c(i,k)
            enddo
         enddo
c     -----------------------------------------------------------------

c     Finosh loop over y planes
      ENDDO


c     -----------------------------------------------------------------
c     Things that require info from all y planes
c     Exchange data between all processors
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Computes Deltat
      if (icfl.eq.1.and.rkstep.eq.1) then

         cfl0 = max(cflx,max(cfly,cflz))

         call MPI_ALLREDUCE(cfl0,reigmx1,1,MPI_REAL,MPI_MAX,
     .        MPI_COMM_WORLD,ierr)

         if (reigmx1.lt.1e-1)  uggg=1

         if (CFL.ne.0 .and. FixTimeStep.eq.0) then
            ! Case of adaptive time stepping
            Deltat=CFL/reigmx1
         else if (CFL.eq.0 .and. FixTimeStep.ne.0) then
            ! Case of fixed time stepping
            Deltat=FixTimeStep
         endif
         dtr=Re/Deltat

         if (uggg.ne.0) then
            write(*,*) 'UGGG', uggg,myid,ihist,jb,je,kb,ke
         endif

      endif
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Saves coefficients for kx = kz = 0

c     MODULE FOR SP2 MPI uses SENDRECV
c     sends a block (jb:je) to everybody and receive from everybody

      ! Check for NAN
      ilocalu =0
      ilocalw =0
      icount  =0
      do j=1,my
         if (rf0u(j).ne.rf0u(j)) then
            icount  =icount +1
            ilocalu  =1
         endif
         if (rf0w(j).ne.rf0w(j)) then
            icount  =icount +1
            ilocalw  =1
         endif
      enddo
      if (ilocalu.eq.1 .OR. ilocalw.eq.1) then
         write(*,*) 'NaN before SENDRECV in hvhg',ilocalu,ilocalw,
     .        'myid=',myid,'   count:',icount
         call system("hostname")
      endif

      ! All processors exchange data for 00 modes
      do iproc=0,numerop-1

         if (iproc.ne.myid) then
            mmyr=jend(iproc)-jbeg(iproc)+1
            call MPI_SENDRECV(rf0u(jb),mmy,MPI_REAL8,
     .           iproc,0,rf0u(jbeg(iproc)),
     .           mmyr,MPI_REAL8,
     .           iproc,0,MPI_COMM_WORLD,istat,ierr)

            call MPI_SENDRECV(rf0w(jb),mmy,MPI_REAL8,
     .           iproc,0,rf0w(jbeg(iproc)),
     .           mmyr,MPI_REAL8,
     .           iproc,0,MPI_COMM_WORLD,istat,ierr)
         endif
      enddo

      ! Check for NAN
      ilocalu =0
      ilocalw =0
      icount  =0
      do j=1,my
         if (rf0u(j).ne.rf0u(j)) then
            icount  =icount +1
            ilocalu  =1
         endif
         if (rf0w(j).ne.rf0w(j)) then
            icount  =icount +1
            ilocalw  =1
         endif
      enddo
      if (ilocalu.eq.1 .OR. ilocalw.eq.1) then
         write(*,*) 'NaN after  SENDRECV in hvhg',ilocalu,ilocalw,
     ,        'myid=',myid,'   count:',icount
         call system("hostname")
      endif
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     at this point: ome1c, rhvc and rhgc are the outputs
c     they are all Fourier x - Fourier z - Physical y
c     each processor has data in a few kx-kz planes
c     ome1c:  dH1/dx + dH3/dz
c     rhgc : -dH3/dx + dH1/dz
c     rhvc : d^2H2/dx^2 + d^2H2/dz^2
c     -----------------------------------------------------------------

      end

