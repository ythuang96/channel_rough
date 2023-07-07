c/********************************************************************/
c/*                                                                  */
c/*      Resuelve las ecuaciones de N-S incompresibles para una      */
c/*    canal tridimensional con una formulacion basada en            */
c/*    la vorticidad y la laplaciana de la velocidad en la dire-     */
c/*    ccion normal a la pared, Kim, Moin y Moser (1987).            */
c/*                                                                  */
c/*       tocado j.jimenez (1/92)                                    */
c/*       hundido j.c.a. (1/01)
c/*                                                                  */
c/*     generates new output files (with dimensions )                */
c/*       head =(time,Re,alp,bet,a0,mx,my,mz)
c/*                                                                  */
c/*     REFERENCIAS                                                  */
c/*       - Kim,J., Moin,P. and Moser,R. (1987) Turbulence Statis    */
c/*         tics in fully developped channel flow at low Reynolds    */
c/*         numbers, J. Fluid Mech., 177, 133-166                    */
c/*                                                                  */
c/*    This version uses CFD in y (o. flores 06)                     */
ccccc VERSION for  MPI  !!!!                                    ccccccc
c/********************************************************************/
      use save_flowfield, only: initialize_save_flowfield_module,
     &  cleanup_save_flowfield_module
      use wall_roughness, only: initialize_wall_roughness
      use restart_file, only: read_restart_file_old
      implicit none
      include "mpif.h"
      include "ctes3D"

      integer nbuffsize,j
      real*4, allocatable:: vor(:),phi(:)
      real*4, allocatable:: hv(:),hg(:)
      real*4, allocatable:: phiwk(:),vorwk(:)
      real*4, allocatable:: dvordy(:),chwk(:)
      real*8  u00(my),w00(my)
      real*8  rf0u(my),rf0w(my),u00wk(my),w00wk(my)

      integer istat(MPI_STATUS_SIZE),ierr
      integer myid,numprocs

      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .               kbeg(0:numerop-1),kend(0:numerop-1),
     .               jb,je,kb,ke,mmy,mmz
      save /point/

      real*8 fmap, y2
      real*4  Re, alp, bet, a0, y, hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/

      integer iinp,iout,id22,isn,ispf
      character*100 filinp,filout,filstt
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                              filinp,filout,filstt
      save /ficheros/

      integer iax, icx
      real alp2, bet2
      complex*8 xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     >     iax(mx),icx(0:mz1)
      save /wave/


c     ! ------------------------ Initializes everything ------------------------
c     ! Initializes MPI
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

      if (numprocs.ne.numerop) then
         write(*,*) 'wrong number of processors',numprocs
         write(*,*) 'compiled with', numerop
         stop
      endif

c     ! Initializes commons and things
      call initcr(myid)

c     ! Allocates buffers
      nbuffsize = mx*max(mmy*mz,mmz*my)

      allocate(vor(nbuffsize))
      allocate(phi(nbuffsize))
      allocate(hv(nbuffsize))
      allocate(hg(nbuffsize))
      allocate(phiwk(nbuffsize))
      allocate(vorwk(nbuffsize))
      allocate(dvordy(nbuffsize))
      allocate(chwk(nbuffsize))

c     ! Initialize custom modules
      call initialize_save_flowfield_module(filstt)
      call initialize_wall_roughness(xalp, xbet)


c     ! Read data from restart file
      call read_restart_file_old(myid, vor,phi,
     .     u00,w00,rf0u,rf0w,hv,hg,
     .     u00wk,w00wk, phiwk,vorwk)
c     ! ------------------------ Start time advancement ------------------------
      call cross1(vor,phi,u00,w00,
     .     rf0u,rf0w,u00wk,w00wk,
     .     hv,hg,
     .     phiwk,vorwk,dvordy,
     .     chwk,
     .     myid)

c     ! -------------------------- Finalize procedure --------------------------
c     ! Clean up the save_flowfield module (deallocate variables)
      call cleanup_save_flowfield_module

      DEALLOCATE(vor   )
      DEALLOCATE(phi   )
      DEALLOCATE(hv    )
      DEALLOCATE(hg    )
      DEALLOCATE(phiwk )
      DEALLOCATE(vorwk )
      DEALLOCATE(dvordy)
      DEALLOCATE(chwk  )

c     ! Finalize MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_FINALIZE(ierr)

      end









c/********************************************************************/
c/*
c/     reads input parameters and initializes things
c/*
c/********************************************************************/
      subroutine initcr(myid)
      use matrices

      implicit none
      include "mpif.h"
      include "ctes3D"

      integer myid

      integer istat(MPI_STATUS_SIZE),ierr

      real*4 trp
      real*8 trp2
      common /mass/ trp(my),trp2(my)
      save /mass/

      integer nimag,nstep,nhist,ihist,icfl
      common /timacc/ nimag,nstep,nhist,ihist,icfl
      save   /timacc/

      real*8  fmap,y2
      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/

      real*4 Deltat,CFL,time,dtr,FixTimeStep
      common /tem/ Deltat,CFL,time,dtr,FixTimeStep
      save /tem/

      integer iinp,iout,id22,isn,ispf
      character*100 filinp,filout,filstt
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                              filinp,filout,filstt
      save /ficheros/

      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .                kbeg(0:numerop-1),kend(0:numerop-1),
     .               jb,je,kb,ke,mmy,mmz
      save /point/

      integer myslice
      common /MPI_datatype/myslice(0:numerop-1)
      save /MPI_datatype/

      real*4 gamma
      integer imesh
      common /mesh/ gamma,imesh
      save /mesh/

      real*4 uampl,vampl,wampl,vspeed
      integer mxwall,mzwall
      common /boundary/ uampl,vampl,wampl,vspeed,mxwall,mzwall
      save /boundary/

      real*4  Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/

      integer idat(30)
      real*4 dat(20),pi
      character*100 text

      integer mybu,mzbu,i,j,iproc,jj,k

      real*8 d11(my,5),d12(my,5),d21(my,5),d22(my,5),aaa(my)

c     PARAMETROS QUE SON CALZABLES:


      Pi=(4.*atan(1.))

c                           /* reads in data                       */
      if(myid.eq.0) then
         open(19,file='hre.dat',status='old')

965      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 965
         read(text,*) (dat(j),j=1,4)

966      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 966
         read(text,*) (idat(j),j=5,7)

964      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 964
         read(text,*) idat(1),dat(7)

967      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 967
         read(text,*) dat(6), dat(12)

968      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 968
         read(text,*) idat(12)

969      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 969
         read(text,*) idat(2),idat(3),dat(8),dat(9),dat(10),dat(11)

65       read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 65
         read(text,'(a100)') filout

166      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 166
         read(text,'(a100)') filinp

167      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 167
         read(text,'(a100)') filstt

         close(19)

         do iproc=1,numerop-1
            call MPI_SEND(dat,12,MPI_REAL,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(idat,19,MPI_INTEGER,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(filout,len(filout),MPI_CHARACTER,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(filinp,len(filinp),MPI_CHARACTER,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(filstt,len(filstt),MPI_CHARACTER,iproc,
     &                 iproc,MPI_COMM_WORLD,ierr)
         enddo

      else

         call MPI_RECV(dat,12,MPI_REAL,0,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(idat,19,MPI_INTEGER,0,
     &                 MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(filout,len(filout),MPI_CHARACTER,0,
     &                 MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(filinp,len(filinp),MPI_CHARACTER,0,
     &                 MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(filstt,len(filstt),MPI_CHARACTER,0,
     &                 MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

      endif


      Re    = dat(1)
      alp   = dat(2)
      bet   = dat(3)
      a0    = dat(4)
      cfl   = dat(6)
      gamma = dat(7)
      uampl = dat(8)
      vampl = dat(9)
      wampl = dat(10)
      vspeed= dat(11) -a0      ! To compensate the movement of the wall
      FixTimeStep = dat(12)

      imesh  = idat(1)
      mxwall = idat(2)
      mzwall = idat(3)
      nstep  = idat(5)
      nimag  = idat(6)
      nhist  = idat(7)
      id22   = idat(12)

      iinp=32                       !! file numbers
      ispf=35
      isn=33
      iout=31

c ---------------  compute y coordinates, pointers and
c ---------------  modes values for FOURIER
      call malla(my,imesh,gamma,y,fmap,y2)
      call genexp

      call pointers(jbeg,jend,kbeg,kend)
c     ===============================================
c     This function divides the y and z grid for all
c     the processors
c     jbeg, kbeg are arrays with the y,z start index
c     for each processor
c     jend, kend are the end index
c     ===============================================
      jb=jbeg(myid) ! y start index for this processor
      je=jend(myid) ! y end index for this processor
      kb=kbeg(myid) ! z start index for this processor
      ke=kend(myid) ! y end index for this processor

      mmz = ke-kb+1 ! number of y and z points for each processor
      mmy = je-jb+1


c    ------------  initializes fast fourier transforms and CFDiff ----
      call cfti(mgalz)
      call rfti(mgalx)
      call prederiv1m(my,d11(1,1),d11(1,2),d11(1,3),d11(1,4),d11(1,5),
     .                d12(1,1),d12(1,2),d12(1,3),d12(1,4),d12(1,5))
      call prederiv2 (my,d21(1,1),d21(1,2),d21(1,3),d21(1,4),d21(1,5),
     .               d22(1,1),d22(1,2),d22(1,3),d22(1,4),d22(1,5),y2(1))

      do j=1,my
         do i=1,5
           prem1(i,j)= d11(j,i)
           dt12(i,j) = d12(j,i)
           dt21(i,j) = d21(j,i)
           dt22(i,j) = d22(j,i)
         enddo
      enddo

      call bandec5(prem1,my)

      !!!! recomputes fmap numerically

      aaa(1) =dt12(3,1)*y2(1)+dt12(4,1)*y2(2) + dt12(5,1)*y2(3)
      aaa(2) =dt12(2,2)*y2(1)+dt12(3,2)*y2(2)+
     .        dt12(4,2)*y2(3)+dt12(5,2)*y2(4)
      do j=3,my-2
      aaa(j)=dt12(1,j)*y2(1+j-3)
      do i=2,5
            aaa(j) = aaa(j) + dt12(i,j)*y2(i+j-3)
      enddo
      enddo
      aaa(my-1)=dt12(1,my-1)*y2(my-3)+dt12(2,my-1)*y2(my-2)+
     .          dt12(3,my-1)*y2(my-1)+dt12(4,my-1)*y2(my)
      aaa(my)=  dt12(1,my)*y2(my-2)  +dt12(2,my)*y2(my-1)+
     .          dt12(3,my)*y2(my)

      call banbks5(prem1,my,aaa)
      do j=1,my
         fmap(j)=1d0/aaa(j)
      enddo

!!!!!!!!  end recompute

      trp(1) = 0.25*(y(2)-y(1))
      trp2(1) = 25d-2*(y2(2)-y2(1))
      do j=2,my-1
         trp(j) = 0.25*(y(j+1)-y(j-1))
         trp2(j) = 25d-2*(y2(j+1)-y2(j-1))
      enddo
      trp(my) = 0.25*(y(my)-y(my-1))
      trp2(my)= 25d-2*(y2(my)-y2(my-1))



      do j=2,my-1
         hy(j) = (y(j+1)-y(j-1))/(2.*2.5)
      enddo
      hy(1)  = (y(2)-y(1))/2.5
      hy(my) = (y(my)-y(my-1))/2.5


c ------------------ Re/dt --------------------------------
      dtr  = Re     !!!! initialize, just in case


c ------------------ MPI Datatypes ------------------------

      do iproc=0,numerop-1

         if (iproc.ne.myid) then

             mzbu = kend(iproc)-kbeg(iproc)+1
             call MPI_TYPE_VECTOR(mmy,mx*mzbu,
     .                            mx*mz,MPI_REAL,myslice(iproc),
     .                            ierr)

             call MPI_TYPE_COMMIT(myslice(iproc),ierr)

         endif
      enddo

c --------------  write header for output -------------
      if(myid.eq.0) then
         write(*,'(a7,f8.2,a8,f6.3,a8,f6.3)')
     .                    '  Re =',Re,'channel'
         write(*,'(a7,f8.3,a8,f6.3,a8,f6.3)')
     .                    'alp =',alp,'  bet =',bet
         write(*,*)

         write(*,'(a8,i5,a8,i5,a8,i5)')
     .                    'mgalx =',mgalx,'mgalz =',mgalz,'my =',my
         write(*,'(a8,i5,a8,i5,a8,i5)')
     .                    'mx =',mx,'mz =',mz
         if (imesh .eq. 1) then
            write(*,*) 'generating tanh mesh'
         else if (imesh .eq. 2) then
            write(*,*) 'generating sen mesh'
         endif
         write(*,*)

         write(*,'(a10,i6,a9,i6,a9,i5)')
     .     'nstep =',nstep,'nimag =',nimag,'nhist =',nhist

         if (CFL.ne.0 .and. FixTimeStep.eq.0) then
            ! Case of adaptive time stepping
            write(*,*)
            write(*,'(a23)') 'Adaptive time stepping:'
            write(*,'(a16,e10.4)') '    CFL = ',CFL
         else if (CFL.eq.0 .and. FixTimeStep.ne.0) then
            ! Case of fixed time stepping
            write(*,*)
            write(*,'(a23)') 'Constant time stepping:'
            write(*,'(a16,e12.7)') '  FixTimeStep = ',FixTimeStep
         else
            ! Case of incorrect input
            ! print error message and exit
            write(*,*)
            write(*,'(a44)')
     &         'Incorrect CFL or FixTimeStep Input. Exiting.'
            stop
         endif

         write(*,*)
         write(*,'(a6,3f8.3,a7,f8.3)')
     .   'vwall',uampl,vampl,wampl,'vspeed',vspeed+a0
         write(*,'(a7,i4,a7,i4)') 'mxwall',mxwall,'mzwall',mzwall
         write(*,*)
         write(*,'(a,a)') 'reading from:  ',filinp
         write(*,*)
         write(*,'(a,a)') '  write in :  ',filout
      endif

      end






      subroutine genexp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     computes wavenumbers and fourier indices
c
c     updated j.j.s.     22/12/00
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include "ctes3D"

      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my)
      save   /fis/

      integer iax,icx
      real alp2,bet2
      complex*8    xalp, xbet
      common /wave/ xalp(0:mx1),xbet(0:mz1),alp2(0:mx1),bet2(0:mz1),
     >     iax(mx),icx(0:mz1)
      save /wave/

      real zero
      integer i,j,k

      zero = 0e0

      do 10 k=0,nz1
         xbet(k) = cmplx(zero,bet*k)
         icx(k) = k
 10   continue

      do 20 k=nz1+1,mz1
         xbet(k) = cmplx(zero ,-bet*(mz1+1-k))
 20   continue

      do 30 k=1,nz1
         icx(mz-k) = k
 30   continue

      do 40 i=0,mx1
         iax(2*i+1) = i
         iax(2*i+2) = i
 40   continue

      do i=0,mx1
         xalp(i) = cmplx(zero ,alp*i)
      enddo

      do i=0,mx1
         alp2(i) = -xalp(i)**2
      enddo

      do j=0,mz1
         bet2(j) = -xbet(j)**2
      enddo

      end


