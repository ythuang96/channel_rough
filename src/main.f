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
      use save_flowfield, only: initialize_save_flowfield
      use wall_roughness, only: initialize_wall_roughness
      implicit none 
      include "mpif.h"
      include "ctes3D"

      integer nbuffsize,nwkasize,j
      real*4, allocatable:: vor(:),phi(:)
      real*8  u00(my),w00(my)
      real*4, allocatable::  wk(:)

      integer istat(MPI_STATUS_SIZE),ierr
      integer master,myid,iproc,numprocs
      integer ihv,ihg,iphiwk,ivorwk,irf0u,irf0w,iu00wk,iw00wk,idvordy,
     .        ichwk,nstr
      integer itags,newtag,imess,i
      real*4 val 

      integer nacumsp,jsp,jsptot,jspiproc,jspend,jspbeg,jspb,jspe,jspbb
       
      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .               kbeg(0:numerop-1),kend(0:numerop-1),
     .               jb,je,kb,ke,mmy,mmz
      save /point/

      common/spectra/   nacumsp,jsp(my),
     .                  jsptot(2*nspec+1),jspiproc(2*nspec+1),
     .                  jspbeg(0:numerop-1), jspend(0:numerop-1),
     .                  jspb,jspe,jspbb
      save/spectra/
      integer nspsize
      real*4, allocatable:: sp(:)

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


c                              /*   initializes everything    */
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

      if (numprocs.ne.numerop) then
         write(*,*) 'wrong number of processors',numprocs
         write(*,*) 'compiled with', numerop
         stop
      endif

      write(*,*) 'MPI enviroment initialized ..',myid

c--------------- initializes commons and things 

      call initcr(myid)
 
c--------------  allocates spectra
    
      jspbb = jspb 
      if (jspe.lt.jspb) jspbb=jspe
      nspsize = (mx1+1)*(nz1+1)*7*(jspe-jspbb+1)

      allocate(sp(nspsize))

c--------------- allocates buffers

      nbuffsize = mx*max(mmy*mz,mmz*my) 
      nwkasize  = 6*nbuffsize + 2*4*my 

c ----------  storage for the pressure
      allocate(vor(nbuffsize))
      allocate(phi(nbuffsize))
      allocate(wk(nwkasize))

      nacumsp = 0
      do i=1,nspsize
         sp(i)=0.
      enddo

c               /*   read input data           */

      write(*,*) 'going to read file ...'

      call getfil(vor,phi,u00,w00,wk,myid)

      write(*,*) '... file read'

      ! initialize custom modules
      call initialize_save_flowfield(filstt)
      call initialize_wall_roughness(xalp, xbet)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c                    /* start time advancement  */

 
      irf0u   = 1
      irf0w   = irf0u   + 2*my 
      iu00wk  = irf0w   + 2*my
      iw00wk  = iu00wk  + 2*my
      ihv     = iw00wk  + 2*my
      ihg     = ihv     + nbuffsize
      iphiwk  = ihg     + nbuffsize
      ivorwk  = iphiwk  + nbuffsize
      idvordy = ivorwk  + nbuffsize
      ichwk   = idvordy + nbuffsize

      call cross1(vor,
     .            phi,
     .            u00,w00,
     .            wk(irf0u),wk(irf0w),wk(iu00wk),wk(iw00wk),
     .            wk(ihv),
     .            wk(ihg),
     .            wk(iphiwk),
     .            wk(idvordy),
     .            wk(ivorwk),
     .            wk(idvordy),
     .            wk(ichwk),
     .            sp,myid)



c                    /* finalize procedure      */

      master=0
      if (myid.ne.0) then

         itags=myid
         call MPI_SEND(1.,1,MPI_REAL,master,
     &              itags,MPI_COMM_WORLD,ierr)

      else

         do iproc=1,numprocs-1

            call MPI_RECV(val,1,MPI_REAL,MPI_ANY_SOURCE,
     &              MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

            imess=istat(MPI_TAG)
            write(*,*) 'process',imess,'over'
            newtag=100
            call MPI_SEND(1.,1,MPI_REAL,imess,
     &                newtag,MPI_COMM_WORLD,ierr)
         enddo

      endif

      if (myid.ne.0) then

         call MPI_RECV(val,1,MPI_REAL,master,
     &              MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         if(istat(MPI_TAG).eq.100) goto 200

      else

         write(*,*) 'process',master,'over'

      endif

200   call MPI_FINALIZE(ierr)

      end

c ====================================================================
c
c                   get data field from file
c                                             jjs  27/1/01
c ====================================================================

      subroutine getfil(vor,phi,u00,w00,wk1,myid)
      use boundary_planes, only: set_bc_from_restart_file
      implicit none
      include "ctes3D"
      include "mpif.h"
      integer istat(MPI_STATUS_SIZE),ierr

      integer myid,master,i,pe,k
      real*4  a0e

      real*4 vor(mx,mz,*),phi(mx,mz,*),wk1(*)
      real*8 u00(*),w00(*)

      real*4  Ree,alpe,bete,ce,xx

      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
      save /tem/

      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .               kbeg(0:numerop-1),kend(0:numerop-1),
     .               jb,je,kb,ke,mmy,mmz
      save /point/

      integer iinp,iout,id22,isn,ispf
      character*100 filinp,filout,filstt
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .                              filinp,filout,filstt
      save /ficheros/

      integer j,mxe,mye,mze,ntotr,iproc,mym,mmy2

      real*4, allocatable ::  uWallBottom(:,:), uWallTop(:,:),
     &    vWallBottom(:,:), vWallTop(:,:), wWallBottom(:,:),
     &    wWallTop(:,:)


      master = 0

c      ---------------  zero everything, fitf
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

c      ---------------  read from file and distribute

      if (myid.eq.master) then     ! ----- this is done by the master

         open(iinp,file=filinp,status='unknown',form='unformatted',
     &        access='direct',recl=8*iwd)
         write(*,*) 'infile opened'
         ! don't read time from restart file
         !read(iinp,rec=1) time,Ree,alpe,bete,a0e,mxe,mye,mze
         read(iinp,rec=1) xx,Ree,alpe,bete,a0e,mxe,mye,mze

         write(*,*) 
         write(*,*) 'reading input file ...' 
         write(*,*) 'time=',time,Ree,alpe,bete,a0e,mxe,mye,mze
         write(*,*) 'mesh:',ce,pe
         write(*,*)
         ntotr=2*mxe*mze
         close(iinp)  
         
         !!!!!! time = 0.   !!!!!!!!!!  ACHTUNG 
       
         open(iinp,file=filinp,status='unknown',form='unformatted',
     &        access='direct',recl=ntotr*iwd)
         read(iinp,rec=1) xx,xx,xx,xx,xx,i,i,i,i,xx,
     &        i,i,xx,xx,xx,xx,(wk1(j),j=1,2*mye) 
      


c     ------------- 00 modes ------------------
         mym = min(my,mye)
         do j=1,mym
            u00(j) = wk1(2*j-1)
            w00(j) = wk1(2*j)
         enddo

         do iproc=1,numerop-1
            call MPI_SEND(mxe,1,MPI_INTEGER,iproc,
     &                iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(mye,1,MPI_INTEGER,iproc,
     &                iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(mze,1,MPI_INTEGER,iproc,
     &                iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(u00,my,MPI_REAL8,iproc,
     &                iproc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(w00,my,MPI_REAL8,iproc,
     &                iproc,MPI_COMM_WORLD,ierr)
         enddo

c         ------------  other modes,  master node -----
         mmy2 = min(je,mye)-jb+1
         write(*,*) 'master reads its data' 
         do j=1,mmy2
            read(iinp,rec=j+jb) (wk1(i),i=1,ntotr)
            call assign(wk1,vor(1,1,j),phi(1,1,j),mx,mz,mxe,mze)
         enddo

c         ------------  other modes,  distribute to slaves --
         do iproc=1,numerop-1

            write(*,*)'master reads proc no',iproc,' data and send them' 

            do j=jbeg(iproc),min(jend(iproc),mye)
               read(iinp,rec=j+1) (wk1(i),i=1,ntotr)
               call MPI_SEND(wk1,ntotr,MPI_REAL,
     &                  iproc,iproc,MPI_COMM_WORLD,ierr)
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
         read(iinp, rec=mye+2)
     &       ((uWallBottom(i,k), uWallTop(i,k), i=1,mxe), k=1,mze)
         read(iinp, rec=mye+3)
     &       ((vWallBottom(i,k), vWallTop(i,k), i=1,mxe), k=1,mze)
         read(iinp, rec=mye+4)
     &       ((wWallBottom(i,k), wWallTop(i,k), i=1,mxe), k=1,mze)

         close(iinp)

      else          ! -------- this is done by all slaves

c         ------------  receive 00 mode ----------
         call MPI_RECV(mxe,1,MPI_INTEGER,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(mye,1,MPI_INTEGER,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(mze,1,MPI_INTEGER,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(u00,my,MPI_REAL8,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(w00,my,MPI_REAL8,master,
     &                MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)

         ntotr=2*mxe*mze

c         ------------  receive other modes ----------
         mmy2 = min(je,mye)-jb+1
         do j=1,mmy2
            call MPI_RECV(wk1,ntotr,MPI_REAL,master,
     &                    MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call assign(wk1,vor(1,1,j),phi(1,1,j),mx,mz,mxe,mze)

         enddo

         write(*,*) 'proc no.',myid,'receives data from master'

         ! velocity boundary conditions: slave processes
         ! allocate temporary memory
         allocate(uWallBottom(mxe,mze))
         allocate(uWallTop(mxe,mze))
         allocate(vWallBottom(mxe,mze))
         allocate(vWallTop(mxe,mze))
         allocate(wWallBottom(mxe,mze))
         allocate(wWallTop(mxe,mze))

      endif

      ! all processes set velocity boundary condition
      call set_bc_from_restart_file(uWallBottom, uWallTop,
     &                              vWallBottom, vWallTop,
     &                              wWallBottom, wWallTop)
      ! and dealocate temporary buffers
      deallocate(uWallBottom)
      deallocate(uWallTop)
      deallocate(vWallBottom)
      deallocate(vWallTop)
      deallocate(wWallBottom)
      deallocate(wWallTop)

c ----- before going any further, exchange cuts in y with cuts in z!

      call chikj2jik(vor,vor,wk1,wk1,myid)
      call chikj2jik(phi,phi,wk1,wk1,myid)

      end


c/********************************************************************/
c/*                                                                  */
c/*               writes a intermediate solution                     */
c/*                                                                  */
c/*    single  jjs  4/01/01                                          */
c/********************************************************************/
      subroutine escru(vor,phi,u00,w00,sp,j1,j2,iproc,iopt,jjsp,
     .     uwallb,uwallt,vwallb,vwallt,wwallb,wwallt,uBulk)
      
      implicit none
      
      include "ctes3D"
      integer iproc,iopt,j1,j2
      real*4 vor(mx,mz,j1:j2),phi(mx,mz,j1:j2)
      real*8 u00(*),w00(*)
      real*4 uwallb(mx,mz),uwallt(mx,mz)
      real*4 vwallb(mx,mz),vwallt(mx,mz)
      real*4 wwallb(mx,mz),wwallt(mx,mz)
      real*8 uBulk
      
      integer i,j,k,jjsp
      
      character*3 ext1
      character*4 ext2
      character*104 fnameima,fnamespe,fnamesta
      
      real*4 timed,Reed,alped,beted,a0ed
      
      real*8  fac
      real*8  um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .     Wx0a,Wz0a,ep,uuv,wwv,vvv
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .     up(my), vp(my), wp(my),
     .     w1m(my),w2m(my),w3m(my),
     .     w1p(my),w2p(my),w3p(my),
     .     uvr(my),uwr(my),vwr(my),
     .     ep(my),uuv(my),wwv(my),vvv(my),
     .     Wx0a,Wz0a,
     .     istati,ntimes,nacum,nstart
      save /statis/
      
      real*4 gamma
      integer imesh
      common /mesh/ gamma,imesh
      save /mesh/
      
      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
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
      common /ficheros/ iinp,iout,id22,isn,ispf,
     .     filinp,filout,filstt
      save /ficheros/
      
c----------------now
      integer nacumsp,jsptot,jsp,jspend,jspbeg,jspiproc,jspb,jspe,jspbb
      
      common/spectra/   nacumsp,jsp(my),
     .     jsptot(2*nspec+1),jspiproc(2*nspec+1),
     .     jspbeg(0:numerop-1), jspend(0:numerop-1),
     .     jspb,jspe,jspbb
      save/spectra/
      real*4 sp(0:mx1,0:nz1,7)
      
      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .     kbeg(0:numerop-1),kend(0:numerop-1),
     .     jb,je,kb,ke,mmy,mmz
      save /point/

      integer numberOfStatisticalQuantities
      
      
      if (iopt.eq.1) then
c--------phi, vor
         
         if(iproc.eq.0) then    !!!!! coming from node 0
c            /* start writing image */
            if(id22.gt.999.or.id22.lt.0) then
               write(*,*) 'number of images out of range'
               stop
            endif
            
            write(ext1,'(i3.3)') id22
            fnameima=filout(1:index(filout,' ')-1)//'.'//ext1
            
            open (iout,file=fnameima,status='unknown',
     &           form='unformatted',access='direct',recl=mx*mz*2*iwd)
            
                                !! rewind(iout)
            write(*,*) j1,j2, 'in escru image',fnameima
            write(iout,rec=1) time,Re,alp,bet,a0,mx,my,mz,imesh,gamma,
     &           mxwall,mzwall,uampl,vampl,wampl,vspeed,
     &           (real(u00(j)),real(w00(j)),j=1,my)
            
            do j=j1,j2
               write(iout,rec=j+1) 
     &              ((vor(i,k,j),phi(i,k,j),i=1,mx),k=1,mz)
            enddo
            
            
c     /*       write statistics       */
            if (nstart.ne.0.and.nacum.ne.0) then
               
               write(*,*) 'stat esc', nstart,nacum
               timed = time
               Reed  = Re
               alped = alp
               beted = bet
               a0ed  = a0
               
               fac = 1./dble(mgalx*mgalz)

               ! number of statistical quantities written to record (r1)
               ! below
               numberOfStatisticalQuantities = 19
               
               write(ext1,'(i3.3)') id22
               fnamesta=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.sta'
               open (isn,file=fnamesta,status='unknown',
     &              form='unformatted')
               rewind(isn)
               write(isn) nacum,Wx0a/nacum,Wz0a/nacum
               write(isn) my,timed,Reed,alped,beted,a0ed
               write(isn) numberOfStatisticalQuantities, uBulk
               ! record (r1). Change numberOfStatisticalQuantities above
               ! if you add/remove quantities here
               write(isn) (um(j),vm(j),wm(j),up(j),vp(j),wp(j),
     &              w1m(j),w2m(j),w3m(j),w1p(j),w2p(j),w3p(j),
     &              uvr(j),uwr(j),vwr(j),ep(j),uuv(j)*fac,
     &              wwv(j)*fac,vvv(j)*fac,j=1,my)
               write(isn) (y(j),j=1,my)
               write(isn) (fmap(j),j=1,my)
               
            endif
            
         else
            write(*,*) j1,j2, 'in escru image'
            
            do j=j1,j2
               write(iout,rec=j+1) 
     &              ((vor(i,k,j),phi(i,k,j),i=1,mx),k=1,mz)
            enddo
            
         endif
         
         if(iproc.eq.numerop-1) then
c     last node writes boundaries!!
            write(iout,rec=my+2) 
     &           ((uwallb(i,k),uwallt(i,k),i=1,mx),k=1,mz)
            write(iout,rec=my+3) 
     &           ((vwallb(i,k),vwallt(i,k),i=1,mx),k=1,mz)
            write(iout,rec=my+4) 
     &           ((wwallb(i,k),wwallt(i,k),i=1,mx),k=1,mz)
            
            
            write(*,*) 'closing files'
            close(iout)
            close(isn )
         endif
         
      endif
      
      if (iopt.eq.2) then
         
         if (nstart.ne.0.and.nacumsp.ne.0) then
c     /*       write spectra        */
            if (jjsp.eq.1) then
               
               write(ext1,'(i3.3)') id22
               fnamespe=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.spe'
               open (ispf,file=fnamespe,status='unknown',
     &              form='unformatted')
               rewind(ispf)
               write(ispf) time,Re, alp, bet, mx,my,mz,nspec+1,
     &              nacumsp
               write(ispf) (jsptot(j), j=1,nspec+1)
               
               write(*,*) 'escribe espectro, proc ',iproc,
     &              'jj= ',jjsp
               
               do k=1,7
                  write(ispf) (sp(i,0,k),
     &                 i=0,(mx1+1)*(nz1+1) - 1)
               enddo        
               
            else
               
               write(*,*) 'escribe espectro, proc ',iproc,
     &              'jj= ',jjsp
               do k=1,7
                  write(ispf) (sp(i,0,k),
     &                 i=0,(mx1+1)*(nz1+1) - 1)
               enddo        
               
            endif        
            
            if (jjsp.eq.nspec+1) then
               write(*,*) 'closing files'
               close(ispf)
            endif
            
         endif
         
      endif
      

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

      real*8  um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .        ep,uuv,wwv,vvv,Wx0a,Wz0a
      integer istati,ntimes,nacum,nstart
      common /statis/   um(my), vm(my), wm(my),
     .                  up(my), vp(my), wp(my),
     .                  w1m(my),w2m(my),w3m(my),
     .                  w1p(my),w2p(my),w3p(my),
     .                  uvr(my),uwr(my),vwr(my),
     .                  ep(my),uuv(my),wwv(my),vvv(my),
     .                  Wx0a,Wz0a,
     .                  istati,ntimes,nacum,nstart
      save /statis/

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

      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
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

      integer nacumsp,jsp,jsptot,jspiproc,jspend,jspbeg,jspb,jspe,jspbb
       
      common/spectra/   nacumsp,jsp(my),
     .                  jsptot(2*nspec+1),jspiproc(2*nspec+1),
     .                  jspbeg(0:numerop-1), jspend(0:numerop-1),
     .                  jspb,jspe,jspbb
      save/spectra/

      real*4 gamma
      integer imesh
      common /mesh/ gamma,imesh
      save /mesh/

      real*4 uampl,vampl,wampl,vspeed
      integer mxwall,mzwall
      common /boundary/ uampl,vampl,wampl,vspeed,mxwall,mzwall
      save /boundary/

      real*4  ener,Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
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
         read(text,*) dat(6) 

968      read(19,'(a)') text
         if(text(1:2).eq.'CC') goto 968
         read(text,*) (idat(j), j=12,14)

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
            call MPI_SEND(dat,11,MPI_REAL,iproc,
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

         call MPI_RECV(dat,11,MPI_REAL,0,
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
      
       


      imesh  = idat(1)
      mxwall = idat(2)
      mzwall = idat(3)
      nstep  = idat(5)
      nimag  = idat(6)
      nhist  = idat(7)
      id22   = idat(12)
      nstart = idat(13)
      ntimes = idat(14)
      
      iinp=32                       !! file numbers
      ispf=35
      istati=0
      isn=33
      iout=31
 
c ---------------  compute y coordinates, pointers and
c ---------------  modes values for FOURIER
      call malla(my,imesh,gamma,y,fmap,y2)
      call genexp

      call pointers(jbeg,jend,kbeg,kend)
      jb=jbeg(myid)
      je=jend(myid)
      kb=kbeg(myid)
      ke=kend(myid)

      write(*,*) 'pointers: ',myid,' jb,je,my=',jb,je,my
      write(*,*) 'pointers: ',myid,' kb,ke,mz=',kb,ke,mz

      mmz = ke-kb+1
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


c --------------  prepare spectra -----
      do j=1,nspec
         jsptot(j)           =jspecy(j)
         jsptot(2*nspec+2-j) = my-jspecy(j)+1
      enddo
      jsptot(nspec+1) = (my+1)/2

      do i=0,numerop-1
        jspbeg(i)=2*nspec+2
         do j=2*nspec+1,1,-1
            if (jsptot(j).ge.jbeg(i)) jspbeg(i) = j
         enddo

         jspend(i)=0
         do j=1,2*nspec+1
            if (jsptot(j).le.jend(i)) jspend(i) = j
         enddo
      enddo

      do i=0,numerop-1 
         do j=1,2*nspec+1
            do jj=jspbeg(i),jspend(i)
               if (j.eq.jj) jspiproc(j) = i
            enddo
         enddo
      enddo 

      do j=1,my
         jsp(j) = 0    
         do jj = 1,2*nspec+1
            if (jsptot(jj).eq.j) jsp(j) = 1
         enddo
      enddo
 
      jspb = jspbeg(myid)
      jspe = jspend(myid)

      write(*,*) 'jspb',jspb,'jspe',jspe,'myid',myid


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


c --------------  initialize stats -------------

      do j=1,my
           um(j)  = 0.
           vm(j)  = 0.
           wm(j)  = 0.
           up(j)  = 0.
           vp(j)  = 0.
           wp(j)  = 0.
           uvr(j) = 0.
           uwr(j) = 0.
           vwr(j) = 0.
           w1m(j) = 0.
           w2m(j) = 0.
           w3m(j) = 0.
           w1p(j) = 0.
           w2p(j) = 0.
           w3p(j) = 0.
           ep(j) = 0.
           uuv(j) = 0.
           wwv(j) = 0.
           vvv(j) = 0.

        enddo

        nacum = 0

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
         write(*,*)

         write(*,'(a10,i6,a9,i6,a9,i5)')
     .     'nstep =',nstep,'nimag =',nimag,'nhist =',nhist
         write(*,'(a10,e10.4,a8,f5.2)') '  CFL =',CFL
         write(*,*)
         write(*,'(a6,3f8.3,a7,f8.3)')
     .   'vwall',uampl,vampl,wampl,'vspeed',vspeed+a0
         write(*,'(a7,i4,a7,i4)') 'mxwall',mxwall,'mzwall',mzwall
         write(*,*)
         write(*,'(a,a)')
     .     'reading from:  ',filinp
         write(*,*)
         write(*,'(a,a)')
     .     '  write in :  ',filout
      endif

      end

c ==============================================================
      subroutine assign(work,vor,phi,mx,mz,mxe,mze)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          single jjs 4/01/01, rewritten jjs 28/01/01
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer   mxm,mzm,klen,kini1,kini2,i,k,k1,k2,mx,mz,mxe,mze
      real*4    vor(mx,mz),phi(mx,mz)
      real*4    work(2,mxe,mze)

      mxm=min(mx,mxe)
      mzm=min(mz,mze)

      klen=(mzm+1)/2
      kini1=mze - klen + 1
      kini2=mz - klen + 1

      do k=1,klen
         do i=1,mxm
            vor(i,k)=work(1,i,k)
            phi(i,k)=work(2,i,k)
         enddo
      enddo
      
      
      do k=1,klen-1
         k1 = k + kini1
         k2 = k + kini2
         do i=1,mxm
            vor(i,k2)=work(1,i,k1)
            phi(i,k2)=work(2,i,k1)
         enddo
      enddo
      
      
      end
      
      subroutine sendsta(myid)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          single jjs 4/01/01
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include "mpif.h"
      include "ctes3D"
      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .     kbeg(0:numerop-1),kend(0:numerop-1),
     .     jb,je,kb,ke,mmy,mmz
      save /point/
      real*8 um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,
     .     uvr,uwr,vwr,Wx0a,Wz0a,ep,uuv,wwv,vvv
      integer istat(MPI_STATUS_SIZE)
      common /statis/   um(my), vm(my), wm(my),
     .     up(my), vp(my), wp(my),
     .     w1m(my),w2m(my),w3m(my),
     .     w1p(my),w2p(my),w3p(my),
     .     uvr(my),uwr(my),vwr(my),
     .     ep(my),uuv(my),wwv(my),vvv(my),
     .     Wx0a,Wz0a,
     .     istati,ntimes,nacum,nstart
      save/statis/
      
      if(myid.ne.0) then
         
         ipo=jb
         leng=mmy
         
         call MPI_SEND(um(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(vm(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(wm(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(up(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(vp(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(wp(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(w1m(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(w2m(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(w3m(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(w1p(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(w2p(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(w3p(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(uvr(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(uwr(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(vwr(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(ep(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(uuv(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(wwv(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
         call MPI_SEND(vvv(ipo),leng,MPI_DOUBLE_PRECISION,0,
     &        myid,MPI_COMM_WORLD,ierr)
      else
         
         do iproc=1,numerop-1
            ipo=jbeg(iproc)
            leng=jend(iproc)-jbeg(iproc)+1
            
            call MPI_RECV(um(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(vm(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(wm(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(up(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(vp(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(wp(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(w1m(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(w2m(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(w3m(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(w1p(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(w2p(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(w3p(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(uvr(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(uwr(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(vwr(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(ep(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(uuv(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(wwv(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
            call MPI_RECV(vvv(ipo),leng,MPI_DOUBLE_PRECISION,iproc,
     &           iproc,MPI_COMM_WORLD,istat,ierr)
         enddo
         
      endif
      
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      
      
