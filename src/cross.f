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
     .     spwk,
     .     vorwk,
     .     dvordy,
     .     chwk,
     .     sp,myid)
      use boundary_planes, only: uWallBottom, uWallTop, vWallBottom,
     &  vWallTop, wWallBottom, wWallTop
      use wall_roughness, only: set_wall_roughness
      use velocity_gradient_tensor, only: collectData,
     &  allocate_data_buffers, deallocate_data_buffers,
     &  write_velocity_gradient_tensor, write_velocity_fields,
     &  collectionFrequencyInTimesteps, fileNrVisualization
      use save_flowfield, only: assess_whether_to_collect_flowfield_dt,
     &  assess_whether_to_collect_flowfield,
     &  write_flowfield_to_hdf_file, collectFlowfield,
     &  save_velocity_forcing_to_buffer, collectWallVelocity,
     &  assess_whether_to_collect_wall_velocity,
     &  save_velocity_bottom_wall_to_buffer,
     &  save_velocity_top_wall_to_buffer,
     &  compute_acceleration_bottom_wall,
     &  compute_acceleration_top_wall

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
      
      real*4  ener,Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/
      
      real*8  um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .     ep,uuv,wwv,vvv,Wx0a,Wz0a
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
      
      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
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
      
      integer nacumsp,jsp,jsptot,jspiproc,jspend,jspbeg,jspb,jspe,jspbb
      common/spectra/   nacumsp,jsp(my),
     .     jsptot(2*nspec+1),jspiproc(2*nspec+1),
     .     jspbeg(0:numerop-1), jspend(0:numerop-1),
     .     jspb,jspe,jspbb
      save/spectra/
      
      real*4 sp  (0:mx1,0:nz1,7,jspbb:jspe),
     .     spwk(0:mx1,0:nz1,7,1:*)
      
      
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
      istati   = 0
      
      irun   = 0                ! first time step is special in tim3rkp
      icfl   = 1                ! first time step always needs a step size
      
      if(myid.eq.0) then
         
         write(ext1,'(i3.3)') id22
         fname=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.cf'
	 write (*,*) fname
         open(39,file=fname,status='unknown')
      endif
c========================================================================
c     THIS IS THE TIME LOOP
c========================================================================
      
      do 30 istep=1,nstep
         
      	 if (myid.eq.0) then
            totaltimer = totaltimer-MPI_WTIME()
            iter_time=-MPI_WTIME()
      	 endif
         
         if (mod(istep-1,nhist) .eq. 0) then
            ihist=1
            icfl= 1
         endif
         
         if (mod(istep-1,ntimes) .eq.0 .and.nstart.ne.0) then
            istati=1
         endif

         ! --- velocity gradient tensor and visualization--- !
         ! if collected, write velocity gradient tensor and
         ! visualization data to file
         if(collectData) then
             ! careful: you have to call write_velocity_fields before
             ! write_velocity_gradient_tensor. The latter has a nasty
             ! side effect (reallocates the velocity arrays) that has
             ! to be fixed
             call write_velocity_fields(myid, jbeg, jend, filstt, 
     &           fileNrVisualization, alp, bet, y, time)
             call write_velocity_gradient_tensor(myid, filstt, 
     &           fileNrVisualization, alp, bet, y, jbeg, jend)
             call deallocate_data_buffers(myid)
             collectData = .false.
             fileNrVisualization = fileNrVisualization + 1
         endif
         ! velocity gradient tensor and visualization data is collected
         ! and written whenever we write a restart file
         if(mod(istep, collectionFrequencyInTimesteps) == 0) then
             call allocate_data_buffers(jb, je, myid)
             collectData = .true.
         endif
         ! --- end velocity gradient tensor and visualization --- !

         ! --- save flowfields
         if (collectFlowfield) then
           ! data collected at the previous timestep is written at the
           ! subsequent timestep. time-Deltat corresponds to the time
           ! at which the data was collected
           call write_flowfield_to_hdf_file(alp, bet, y, jb, je, 
     &         xalp, xbet, time-Deltat, Re, massu0)
         endif
         ! call assess_whether_to_collect_flowfield_dt(time)
         call assess_whether_to_collect_flowfield(istep-1)
         call assess_whether_to_collect_wall_velocity(istep)
         ! --- end save flowfields

c      ==================================   write image to a file    */
         if (mod(istep-1,nimag) .eq. 0 .and. istep.ne.1) then
            
            if (myid.eq.0) then
               write_time = -MPI_WTIME()
               write(*,*) time, vWallBottom(mxwall,mzwall),
     &             vWallTop(mxwall, mzwall)
            endif
            
            
c     /* procedure to receive all the slices and stats  */
            call sendsta(myid)
            call chjik2ikj(phi,phi,dvordy,dvordy,myid)
            call chjik2ikj(vor,vor,dvordy,dvordy,myid)
            
            do j=0,mx*max(mmy*mz,mmz*my)
               chwk(1+j)        = vor(j,0,kb)
               dvordy(j,0,kb) = phi(j,0,kb)
            enddo
            
c-------everybody sends data to the master
            
            if(myid.ne.0) then
               
               call MPI_SEND(chwk,mx*mz*mmy,MPI_REAL,
     &              0,myid,MPI_COMM_WORLD,ierr)
               call MPI_SEND(dvordy,mx*mz*mmy,MPI_REAL,
     &              0,myid,MPI_COMM_WORLD,ierr)
               
            else
               
c------the master first writes its stuff
               
               call escru(chwk,dvordy,u00,w00,spwk,jb,je,0,1,1,
     .              uWallBottom,uWallTop,vWallBottom,vWallTop,
     .              wWallBottom,wWallTop,massu0)
               
c------then receives everything from everybody
               
               do iproc=1,numerop-1
                  
                  leng=(jend(iproc)-jbeg(iproc)+1)*mx*mz
                  call MPI_RECV(chwk,leng,MPI_REAL,
     &                 iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                  
                  call MPI_RECV(dvordy,leng,MPI_REAL,
     &                 iproc,iproc,MPI_COMM_WORLD,istat,ierr)
                  
                  
c-------and writes it
                  
                  call escru(chwk,dvordy,u00,w00,spwk,jbeg(iproc),
     &                 jend(iproc),iproc,1,1,
     &                 uWallBottom,uWallTop,vWallBottom,vWallTop,
     &                 wWallBottom,wWallTop,massu0)
                  
               enddo
               
            endif 
            
c-------upper half sends spectra
            
            do j=2*nspec+1,nspec+1,-1
               
               jj = 2*nspec+2 - j
               
               if (myid.eq.jspiproc(j).and.
     &              myid.ne.jspiproc(jj)) then
                  
                  leng = (mx1+1)*(nz1+1)*7
                  
                  call MPI_SEND(sp(0,0,1,j),leng,MPI_REAL,
     &                 jspiproc(jj),0,MPI_COMM_WORLD,ierr)
                  
               endif
               
c-------lower half receives upper half spectra and computes average
               
               if (myid.eq.jspiproc(jj)) then
                  
                  leng = (mx1+1)*(nz1+1)*7
                  leng1 = (mx1+1)*(nz1+1)*3
                  leng2 = (mx1+1)*(nz1+1)*4
                  
                  if (jspiproc(j).ne.myid) then
                     
                     call MPI_RECV(spwk(0,0,1,1),leng,MPI_REAL,
     &                    jspiproc(j),0,MPI_COMM_WORLD,istat,ierr)


                  else
                     
                     do i=0,leng-1
                        spwk(i,0,1,1) = sp (i,0,1,j)
                     enddo
                     
                  endif
                  
c-------------velocity spectra are symmetric
                  
                  do i = 0,leng1-1
                     spwk(i,0,1,1) =.5*(spwk(i,0,1,1) + sp(i,0,1,jj))
                  enddo
                  
c-------------velocity cospectra are skew-symmetric
                  
                  do i = leng1,leng2-1
                     spwk(i,0,1,1) =.5*(-spwk(i,0,1,1) + sp(i,0,1,jj))
                  enddo
                  
c-------------vorticity spectra are symmetric
                  
                  do i = leng2,leng-1
                     spwk(i,0,1,1) =.5*(spwk(i,0,1,1) + sp(i,0,1,jj))
                  enddo
                  
c-------everybody sends data to the master
                  
                  if(myid.ne.0) then
                     
c-------only lower half sends averaged spectra to the master
                     
                     leng = (mx1+1)*(nz1+1)*7
                     call MPI_SEND(spwk,leng,MPI_REAL,
     &                    0,myid,MPI_COMM_WORLD,ierr)
                     
                  else
c------the master first writes its stuff
                     
                     call escru(chwk,dvordy,u00,w00,spwk,jb,je,0,2,jj,
     .                    uWallBottom,uWallTop,vWallBottom,vWallTop,
     .                    wWallBottom,wWallTop,massu0)
                     
                  endif             
                  
               endif 
               
               if (myid.eq.0.and.jspiproc(jj).ne.myid) then
                  
c------then receives everything from everybody
                  
                  leng = (mx1+1)*(nz1+1)*7
                  call MPI_RECV(spwk,leng,MPI_REAL,jspiproc(jj),
     &                 jspiproc(jj),MPI_COMM_WORLD,istat,ierr)
                  
c-------and writes it
                  
                  call escru(chwk,dvordy,u00,w00,spwk,jbeg(iproc),
     &                 jend(iproc),iproc,2,jj,
     &                 uWallBottom,uWallTop,vWallBottom,vWallTop,
     &                 wWallBottom,wWallTop,massu0)
                  
               endif
               
            enddo
            
            
            call chikj2jik(phi,phi,dvordy,dvordy,myid)
            call chikj2jik(vor,vor,dvordy,dvordy,myid)
            
            if (myid.eq.0) then
               write(*,*) 'time write:',MPI_WTIME()+write_time
               id22 = id22+1
            endif
            
            
         endif
c     ==================================  finished writing image */

         
c     /********************************************************************/
c     /*      time stepper                                                */
c     /*      da un paso en el tiempo. Runge - Kutta  para terminos       */
c     /*      convectivos euler inverso para la parte viscosa.            */
c     /*                                                                  */
c     /*       Resuelve:    Gt  + Hg = 1/Re G"                            */
c     /*                    V"t + Hv = 1/Re V""                           */
c     /*                                                                  */
c     /*   input:                                                         */
c     /*     vor: vorticidad segun la direccion y (n) 		     */
c     /*     phi: laplaciana velocidad segun la direccion y  (n)          */
c     /*     vorwk: copia de vor para el calc. de los term. no lineales   */
c     /*     phiwk: copia de phi para el calc. de los term. no lineales   */
c     /*     hg: Hg                                                       */
c     /*     hv: Hv                                                       */
c     /*     dvordy: area de trabajo de dimension mx*nxymax               */
c     /*     chwk: area de trabajo para los chz2y                         */
c     /*                                                                  */
c     /*  output:                                                         */
c     /*     vor: vorticidad segun la direccion y (n+1)                   */
c     /*     phi: laplaciana velocidad segun la direccion y  (n+1)        */
c     /*      hg: contiene  v (velocidad segun y)                         */
c     /*      hv: contiene  dv/dy                                         */
c     /*..................................................................*/
c     /*  MODULE FOR MPI SP2                                              */
c     /*..................................................................*/
c     /*                                                                  */
c     /*   updated jjs 07/01/01                                           */
c     /*   in jik form jca                                                */
c     /*   to CFdiff by of						  */
c     /*   cleaned by jj 05/08 for ctr, without major changes             */
c     */                                                                  */
c     */      ACHTUNG:  This runge kutta is ok but could be better        */
c     */                                                                  */
c     /********************************************************************/
         



c     /********************************************************************/
c     /*      time stepper                                                */
c     /*      take a step in time. Runge - Kutta for terms       */
c     /*      inverse euler convective for the viscous part.            */
c     /*                                                                  */
c     /*       Solves  :    Gt  + Hg = 1/Re G"                            */
c     /*                    V"t + Hv = 1/Re V""                           */
c     /*                                                                  */
c     /*   input:                                                         */
c     /*     vor: vorticity according to the y (n) direction          */
c     /*     phi: laplacian velocity according to the direction y (n)     */
c     /*     vorwk: vor copy for the calc. of the non-linear terms.    */
c     /*     phiwk: copy of phi for the calc. of the non-linear terms.    */
c     /*     hg: Hg                                                       */
c     /*     hv: Hv                                                       */
c     /*     dvordy: working area of ​​dimension mx * nxymax               */
c     /*     chwk: work area for chz2y                         */
c     /*                                                                  */
c     /*  output:                                                         */
c     /*     vor: vorticity according to the y direction (n + 1)                   */
c     /*     phi: laplacian velocity according to the y direction (n + 1)        */
c     /*      hg: contains v (velocity according to y)                         */
c     /*      hv: contains dv / dy                                         */
c     /********************************************************************/


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this is done only for the first step
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
         if(irun.eq.0) then
            
            irun = 1
            
c     /*   Todos necesitan el tiempo!!!!
           ! all processors have explicitly set the time to zero, so
           ! this step can be skipped
           ! call  MPI_BCAST(time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
            
           ! first step: set boundary conditions for v
           ! vWall is read from restart file, but let's set it
           ! explicitly here, just in case we want to do something
           ! else than in the restart file
           call set_wall_roughness(uWallBottom, uWallTop, vWallBottom,
     &       vWallTop)

c     /*   calcula la v a partir de phi */
c     /*   calculate the v from phi */
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
                  call Lapvdv(phi(0,i,k),hg(0,i,k),hv(0,i,k),
     .                 rk,bcb,bct)
               enddo
            enddo

c     /*     prepara phiwk,vorwk,u00wk,w00wk */
            
            do j=0,my1
               u00wk(j)=u00(j)
               w00wk(j)=w00(j)
            enddo
            
            call deryr(u00wk,rf0u,my)
            call deryr(w00wk,rf0w,my)
            
            do  k=kb,ke
               do i=0,mx1
                  do j=0,2*my-1
                     vorwk(j,i,k)=vor(j,i,k)
                     phiwk(j,i,k)=phi(j,i,k)
                  enddo
               enddo
            enddo
c     computes mass 
            massu0=0d0
            massw0=0d0
            do j=0,my1
               massu0 = massu0 + trp2(j)*u00(j)
            enddo

c     !!! target mass flux.  historical Madrid value 
            massu0 = .8987636566162d0 
            if(myid.eq.0) then
               write (*,*) 'mass flux:',massu0,massw0
            endif
            
         endif
cccccccccc--end special first step ccccccccccccccccccccc
         
         
         
c     /*  Runge-Kutta third order  */
c-----
         
         do 10 rkstep=1,3
            
c-----------------------computes d (ome_2) / dy
c-----------------------dvordy      : d (vorwk) / dy -- F-F
            
            call deryr2(vorwk,dvordy,(mx1+1)*mmz,my,chwk)
            
c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
c     all arrays changed from z slices to y slices
c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
            
            call chjik2ikj(phiwk,phiwk,chwk,chwk,myid)
            call chjik2ikj(hv,hv,chwk,chwk,myid)
            call chjik2ikj(hg,hg,chwk,chwk,myid)
            call chjik2ikj(dvordy,dvordy,chwk,chwk,myid)
            call chjik2ikj(vorwk,vorwk,chwk,chwk,myid)
            
c-----------------------storage for 0 modes
            
            ipo1 = 1    + my
            ipo2 = ipo1 + my
            ipo3 = ipo2 + my
            ipo4 = ipo3 + my
            
c-----------------------all nodes compute 0  mode of vorticity
c     work(1):du00;  work(ipo1):dw00;  work(ipo2):u00;  work(ipo3):w00;
            
            do j=0,my1
               work(1+j)   =rf0u(j)           
               work(ipo1+j)=rf0w(j)
               work(ipo2+j)=u00wk(j)
               work(ipo3+j)=w00wk(j)
               
            enddo
            
            
            call hvhg(phiwk,vorwk,hv,hg,rf0u,
     .           rf0w,dvordy,work,sp,myid,rkstep, 
     .           u1r,u2r,u3r,o1r,o2r,o3r, 
     .           u1r,u2r,u3r,o1r,o2r,o3r)
            
            
            H00u=0d0
            H00w=0d0
            do j=0,my1
               H00u = trp2(j)*rf0u(j)
               H00w = trp2(j)*rf0w(j)
            enddo
            
            do j=0,my1
               rf0u(j) = rf0u(j)-H00u
               rf0w(j) = rf0w(j)-H00w
            enddo

c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
c     at this point: dvordy, rhv and rhg are the outs
c     they must be trasformed from y slices to z slices before completion
c     dvordy: dH1/dx + dH3/dz
c     hg: -dH3/dx + dH1/dz
c     hv: d^2H2/dx^2 + d^2H2/dz^2
c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
            
            call chikj2jik(dvordy,dvordy,chwk,chwk,myid)
            call chikj2jik(hv,hv,chwk,chwk,myid)
            call chikj2jik(hg,hg,chwk,chwk,myid)
            
c------------------------computes dvordy = d (dH1/dx + dH3/dz) / dy
            call deryr2(dvordy,dvordy,(mx1+1)*mmz,my,chwk)
            
c---  computes  rhv =d^2H2/dx^2 + d^2H2/dz^2 - d(dH1/dx + dH3/dz)/dy
            
            do k=kb,ke
               do i=0,mx1
                  do j=0,2*my-1
                     hv(j,i,k) = hv(j,i,k) - dvordy(j,i,k)
                  enddo
               enddo
            enddo
            
c-----------

            ! save velocity forcing if required
            ! note: the velocity, vorticity and vorticity forcing
            ! fields are saved in subroutine hvhg, see
            ! save_velocity_at_wallparallel_plane_to_buffer
            ! save_vorticity_at_wallparallel_plane_to_buffer
            ! save_vorticity_forcing_at_plane_to_buffer
            if (collectFlowfield .and. rkstep==1) then
              call save_velocity_forcing_to_buffer(hv, jb, je)
            endif
            
c-----------advances in time : nonlinear terms explicitly
            
            r1=c(rkstep)*Deltat
            
            
            do j=0,my1
               u00wk(j)=u00(j)+r1*rf0u(j)
               w00wk(j)=w00(j)+r1*rf0w(j)
            enddo
            
            do k=kb,ke
               do i=0,mx1
                  do j=0,2*my-1
                     vorwk(j,i,k)=vor(j,i,k)+r1*hg(j,i,k)
                     phiwk(j,i,k)=phi(j,i,k)+r1*hv(j,i,k)
                  enddo
               enddo
            enddo
            
            dtr1=dtr/c(rkstep)
            
            
            do j=0,my1
               rf0u(j)=-dtr1*u00wk(j)
               rf0w(j)=-dtr1*w00wk(j)
            enddo
            
            do k=kb,ke
               do i=0,mx1
                  do j=0,2*my-1
                     hg(j,i,k)=-dtr1*vorwk(j,i,k)
                     hv(j,i,k)=-dtr1*phiwk(j,i,k)
                  enddo
               enddo
            enddo
            

c------------ updating boundary condition for the viscous time step 
c                    ----    BCS-HERE ----
            ! subsequent timesteps: set velocity boundary conditions
            ! uWall = wWall = 0, this is built in implicitly below
            ! in bcbo, bcbdv etc.
            ! v is set according to opposition control
            call set_wall_roughness(uWallBottom, uWallTop, vWallBottom,
     &          vWallTop)

c      ----   END OF BCS-HERE FOR THE BOTTOM WALL ---- 

            ! save wall acceleration
            ! note: the velocity, vorticity and vorticity forcing
            ! fields are saved in subroutine hvhg, see:
            ! save_velocity_at_wallparallel_plane_to_buffer
            ! save_vorticity_at_wallparallel_plane_to_buffer
            ! save_vorticity_forcing_at_plane_to_buffer
            ! the velocity forcing is saved above, see
            ! save_velocity_forcing_to_buffer
            ! time step before data collection: save previous velocity
            ! to compute acceleration
            if(collectWallVelocity .and. (.not. collectFlowfield)) then
              if(rkstep == 1) then
                call save_velocity_bottom_wall_to_buffer(uWallBottom,
     &            vWallBottom, wWallBottom)
                call save_velocity_top_wall_to_buffer(uWallTop,
     &            vWallTop, wWallTop)
              endif
            endif
            ! data collection timestep: compute acceleration from
            ! previous velocity and current velocity
            if(collectWallVelocity .and. collectFlowfield) then
              if(rkstep == 1) then
                call compute_acceleration_bottom_wall(uWallBottom,
     &            vWallBottom, wWallBottom, Deltat)
                call compute_acceleration_top_wall(uWallTop,
     &            vWallTop, wWallTop, Deltat)
              endif
            endif
            
            do k=kb,ke
               k1 = icx(k-1)
               do i=0,mx1
                  
                  rk = bet2(k1)+alp2(i)
                  rk2 = dtr1+rk

                  ! boundary condition for rough wall
                  bcbo = -xalp(i)*wWallBottom(i,k-1) +
     &              xbet(k-1)*uWallBottom(i,k-1)
                  bcto = -xalp(i)*wWallTop(i,k-1) +
     &              xbet(k-1)*uWallTop(i,k-1)

                  bcb = vWallBottom(i,k-1)
                  bct = vWallTop(i,k-1)

                  bcbdv = -xalp(i)*uWallBottom(i,k-1) -
     &              xbet(k-1)*wWallBottom(i,k-1)
                  bctdv = -xalp(i)*uWallTop(i,k-1) -
     &              xbet(k-1)*wWallTop(i,k-1)

                  ! temporary: no slip bc
                  !bcbo = (0.0, 0.0)
                  !bcto = (0.0, 0.0)

                  !bcb = (0.0, 0.0)
                  !bct = (0.0, 0.0)

                  !bcbdv = (0.0, 0.0)
                  !bctdv = (0.0, 0.0)
                  
                  call lapsov(phiwk(0,i,k),hg(0,i,k),hv(0,i,k),
     &                 hv(0,i,k),vorwk(0,i,k),hg(0,i,k),
     &                 rk2,rk,bcb,bct,bcbdv,bctdv,bcbo,bcto)
                  
               enddo
            enddo
c                    ----  END  BCS-HERE ---- 
            
            

c     /* Kx = Kz = 0 modes         */
c     Se calculan en dos fases:
            
            
c     Velocidad sin correccion de presion-masa
            bcbr = 0d0          !! Estas no tienen por que ser nulas!!
            bctr = 0d0
            
            rk = dtr1 
            call Lapv1(rf0u,u00wk,rk,bcbr,bctr)
            call Lapv1(rf0w,w00wk,rk,bcbr,bctr)
            
c     Velocidad correccion presion-masa, cc homogeneas:
            bcbr = 0d0
            bctr = 0d0
            
            do j=0,my1
               rf0u(j) = -dtr1
               rf0w(j) = -dtr1
            enddo 
            
            rk = dtr1
            call Lapv1(rf0u,rf0u,rk,bcbr,bctr)
            call Lapv1(rf0w,rf0w,rk,bcbr,bctr)
            
            
c     Calculo masa de cada flujo:
            
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
            
            cteu =(massu0 - massu1)/massu2
            ctew =(massw0 - massw1)/massw2
            
            do j=0,my1
               u00wk(j) = u00wk(j) + cteu*rf0u(j)
               w00wk(j) = w00wk(j) + ctew*rf0w(j)
            enddo
            
            call deryr(u00wk,rf0u,my)
            call deryr(w00wk,rf0w,my)
            
c     u00wk,w00wk ;   y ademas rf0$ = d$00        
            
            
 10      continue               !!! finalizado subpaso del RK3
         
         do j=0,my1
            u00(j)=u00wk(j)
            w00(j)=w00wk(j)
         enddo

         do k=kb,ke
            do i=0,mx1
               do j= 0,2*my-1
                  vor(j,i,k)=vorwk(j,i,k)
                  phi(j,i,k)=phiwk(j,i,k)
               enddo
            enddo
         enddo
         
c--------------------send WxL, WzL
         if (ihist.eq.1) then
            
            if (myid.eq.numerop-1) then
               
               chwk (1) = WxL
               chwk (2) = WzL
               chwk (3) = uvL 
               
               call MPI_SEND(chwk,3,MPI_REAL,
     &              0,0,MPI_COMM_WORLD,ierr)
               
            endif
            
            if (myid.eq.0) then
               
               call MPI_RECV(chwk,3,MPI_REAL,
     &              numerop-1,0,MPI_COMM_WORLD,istat,ierr)
               
               WxL =  chwk(1)
               WzL =  chwk(2)
               uvL =  chwk(3)
            endif
         endif
         
         
c     /*     write history record     */
         
         if(myid.eq.0) then
            
            
            reynota=0.5*re*re*(abs(wz0/re+uv0)
     .           +abs(wzl/re+uvL))
            
            massu = massu1 +cteu*massu2
            massw = massw1 +ctew*massw2
            
            
            if (ihist.eq.1) then
               
               tmp=my1/2
               
 325           format(i5,9(d14.6))
               write(*,325) istep,time,-1.*Wz0,WzL,sqrt(reynota),Deltat,
     .              u00(floor(tmp))*Re/sqrt(reynota),
     .              massw,uv0,uvL
               
 324           format(17(d22.14))
               write(39,324) time,-1.*Wz0,WzL,sqrt(reynota),ener,
     .              u00(floor(tmp))*Re/sqrt(reynota),massw
               call flush(39)
               
            endif
            
            if (mod(istep-1,nimag).eq.0 .and. 
     &           istep.ne.1 .and. istep.ne.nstep) then
               
               close(39)
               write(ext1,'(i3.3)') id22
               fname=filstt(1:index(filstt,' ')-1)//'.'//ext1//'.cf'
               write (*,*) fname
               open(39,file=fname,status='unknown')
               
            endif
            
         endif
         
c     /* time:
         time=time+Deltat
         
         
         if(istati.eq.1) istati = 0
         if(icfl.eq.1)   icfl   = 0
         if(ihist.eq.1)  ihist  = 0
         
         if (myid.eq.0) then
            totaltimer=totaltimer+MPI_WTIME()
            write(*,'(i7,3f20.5)') istep,
     >           MPI_WTIME()+iter_time-commtimer+comm_time,
     >           commtimer-comm_time,MPI_WTIME()+iter_time
            comm_time = commtimer
         end if
         
 30   continue
      
      if (myid.eq.0) then
         print *,"Total time: ",totaltimer
         print *,"Trans. time: ",transtimer
         print *,"Comm. time: ",commtimer
         print *,"Comm/Total: ",commtimer/totaltimer
         print *,"Trans/Total: ",transtimer/totaltimer
         print *,"Aver time per step: ",totaltimer/nstep
         
      end if
      
      end
      
c     /********************************************************************
c     /*                                                                  */
c     /*         computes the forcing terms (convective terms)            */
c     /*                                                                  */
c     /*    input:                                                        */
c     /*      phi: Delt(v)   (F-F-Tch)                                    */
c     /*      vor: vorticity (F-F-Tch)                                    */
c     /*      rhg: velocity along y axis (F-F-phys)                       */
c     /*      rhv: dv/dy                 (F-F-phys)                       */
c     /*      work: area de trabajo de dimension al menos  max(mgalx,4*my) */
c     /*                                                                  */
c     /*   output:                                                        i/
c     /*     rhv: nonlinear term for the phi equation    (F-F-Tch)       */
c     /*     rhg: non linear term for vorticity equation (F-F-Tch)       */
c     /*     rf0u: non linear term for the evolution of Kx=Kz=0  u        */
c     /*     rf0w: non linear term for the evolution of Kx=Kz=0  w        */
c     /*                                                                  */
c     /*..................................................................*/
c     /*  MODULE FOR MPI SP2                                              */
c     /*..................................................................*
c     /*                                                                  */
c     /*    updated jjs 22/12/00     (INCOMPLETE)                         */
c     /*    single  jjs  4/01/01                                          */
c     /*    low storage : 24/01/01 (incomplete)
c     /********************************************************************/
      subroutine hvhg(phic,ome2c,rhvc,rhgc,
     .     rf0u,rf0w,ome1c,
     .     work2,sp,myid,rkstep, 
     .     u1r,u2r,u3r,o1r,o2r,o3r, 
     .     u1c,u2c,u3c,o1c,o2c,o3c )
      use velocity_gradient_tensor, only: collectData,
     &  compute_velocity_gradient_tensor, save_velocity_fields
      use save_flowfield, only: collectFlowfield,
     &  save_velocity_at_wallparallel_plane_to_buffer,
     &  save_vorticity_at_wallparallel_plane_to_buffer,
     &  save_vorticity_forcing_at_plane_to_buffer
      implicit none
      include "ctes3D"
      include "mpif.h"
      
      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .     kbeg(0:numerop-1),kend(0:numerop-1),
     .     jb,je,kb,ke,mmy,mmz
      save   /point/
      
      integer nacumsp,jsp,jsptot,jspiproc,jspend,jspbeg,jspb,jspe,jspbb
      
      common/spectra/   nacumsp,jsp(my),
     .     jsptot(2*nspec+1),jspiproc(2*nspec+1),
     .     jspbeg(0:numerop-1), jspend(0:numerop-1),
     .     jspb,jspe,jspbb
      save/spectra/
      real*4 sp(0:mx1,0:nz1,7,jspbb:jspe)
      
      
      real*4 Deltat,CFL,time,dtr
      common /tem/ Deltat,CFL,time,dtr
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
      
      real*4 uner(9)
      real*4  ener,Wx0,Wz0,WxL,WzL,uv0,uvL
      common /diag/ ener(9),Wx0,Wz0,WxL,WzL,uv0,uvL
      save   /diag/
      
      real*8 um,vm,wm,up,vp,wp,w1m,w2m,w3m,w1p,w2p,w3p,uvr,uwr,vwr,
     .     ep,uuv,wwv,vvv,Wx0a,Wz0a
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
      integer i,j,k,jj,iy,kk,jndex
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
      
c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
c     at this point:
c     rhv is dv / dy -- F-F-P
c     rhg is v -- F-F-P
c     phi is nabla^2(v) -- F-F-P
c     ome1 is d(omega_2)/dy -- F-F-P
c     ome2 is omega_2 --F-F-P
c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
      
c----------------------start operating by planes
      
c----------------Initialize variables out of the y loop
      
      do kk=1,9
         uner(kk) = 0.
      enddo
      
      iy = jspb
      
      cflx = 0.
      cfly = 0.
      cflz = 0.
      cfl0 = 0.
      
      hxalp=alp*mx*0.5
      hzalp=bet*mz*0.5
      
      
      DO J = JB,JE
         
c----------------------computes non 0 modes of ome1, ome3
         
         do k=1,mz1
            dk  = 1.0/xbet(k)
            o3c(0,k) = -ome1c(0,k,j) * dk
            o1c(0,k) = - phic(0,k,j) * dk
         enddo
         
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
         
c----------------------computes non 0 modes of u,w
         do k=1,mz1
            dk  = 1.0/xbet(k)
            u3c(0,k) = -rhvc (0,k,j) * dk
            u1c(0,k) = ome2c (0,k,j) * dk
         enddo
         
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
         
c----------------------all nodes, 0 modes, u,w,ome1,ome3
         
         jj = j-1
         o3c(0,0) = -work2(      j)
         o1c(0,0) =  work2(ipo1+jj)
         u1c(0,0) =  work2(ipo2+jj) 
         u3c(0,0) =  work2(ipo3+jj)
         
c     -------------- copy v and omega_2 into their planes
         do k=0,mz1
            do i=0,mx1
               u2c(i,k) = rhgc(i,k,j)
               o2c(i,k) = ome2c(i,k,j) !!! warning, ome2c(0,0,j) must be zero
            enddo
         enddo

c     
c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
c     at this point:
c     3-D arrays    --------------
c     rhv is dv / dy -- F-F-P
c     rhg is v -- F-F-P
c     phi is nabla^2(v) -- F-F-P
c     ome1 is d(omega_2)/dy -- F-F-P
c     ome2 is omega_2 --F-F-P
c     3-D arrays    --------------
c     
c     2-D arrays    --------------
c     u1 is u
c     u2 is v
c     u3 is w
c     o1 is omega_1
c     o2 is omega_2
c     o3 is omega_3
c     all variables in Fourierx -- Fourier z -- Physical y
c     2-D arrays    --------------
c     everything in ONE x-z plane
c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
         
         
c======
c     /*     statistics of  v, omega_1, omega_3      */
c     /*     statistics of  u,w, omega_2             */
         
         if (istati.eq.1 .and. rkstep.eq.1) then
            
c     -----------------  spectra
c     ---------------- only if my plane contains spectra information
            
            if (jsp(j).eq.1) then
               
               do kk = 0,mz1
                  k = icx(kk)
                  sp(0,k,1,iy) = sp(0,k,1,iy)+u1c(0,kk)*
     &                 conjg(u1c(0,kk))
                  sp(0,k,2,iy) = sp(0,k,2,iy)+u2c(0,kk)*
     &                 conjg(u2c(0,kk))
                  sp(0,k,3,iy) = sp(0,k,3,iy)+u3c(0,kk)*
     &                 conjg(u3c(0,kk))
                  sp(0,k,4,iy) = sp(0,k,4,iy)+real(u1c(0,kk)*
     &                 conjg(u2c(0,kk)))
                  sp(0,k,5,iy) = sp(0,k,5,iy)+o1c(0,kk)*
     &                 conjg(o1c(0,kk))
                  sp(0,k,6,iy) = sp(0,k,6,iy)+o2c(0,kk)*
     &                 conjg(o2c(0,kk))
                  sp(0,k,7,iy) = sp(0,k,7,iy)+o3c(0,kk)*
     &                 conjg(o3c(0,kk))
                  do i = 1,mx1
                     sp(i,k,1,iy) = sp(i,k,1,iy)+2.*u1c(i,kk)*
     &                    conjg(u1c(i,kk))
                     sp(i,k,2,iy) = sp(i,k,2,iy)+2.*u2c(i,kk)*
     &                    conjg(u2c(i,kk))
                     sp(i,k,3,iy) = sp(i,k,3,iy)+2.*u3c(i,kk)*
     &                    conjg(u3c(i,kk))
                     sp(i,k,4,iy) =sp(i,k,4,iy)+2.*real(u1c(i,kk)*
     &                    conjg(u2c(i,kk)))
                     sp(i,k,5,iy) = sp(i,k,5,iy)+2.*o1c(i,kk)*
     &                    conjg(o1c(i,kk))
                     sp(i,k,6,iy) = sp(i,k,6,iy)+2.*o2c(i,kk)*
     &                    conjg(o2c(i,kk))
                     sp(i,k,7,iy) = sp(i,k,7,iy)+2.*o3c(i,kk)*
     &                    conjg(o3c(i,kk))
                  enddo
               enddo
               
c     ----------------  update spectra y index
               
               iy = iy + 1
               
            endif
            
c----------just add 1 to nacum once !!!!
            
            if (j.eq.je)  nacumsp = nacumsp +1
            
c     ----- intensities ------------
c     ----- & dissipation ------------
            
            do kk=1,9
               ener(kk) = 0.
            enddo
            
            do k=0,mz1
               
c     intensities ----------------
               
               aa = u1c(0,k)*conjg(u1c(0,k))
               up(j) = up(j) + aa
               ener(1)=ener(1) + aa
               
               aa= u2c(0,k)*conjg(u2c(0,k))
               vp(j) = vp(j) + aa
               ener(2)=ener(2) + aa
               
               aa= u3c(0,k)*conjg(u3c(0,k))
               wp(j) = wp(j) + aa
               ener(3)=ener(3) + aa
               
               aa= real(u1c(0,k)*conjg(u2c(0,k)))
               uvr(j)= uvr(j) + aa
               ener(4)=ener(4) + aa
               
               aa= real(u1c(0,k)*conjg(u3c(0,k)))
               uwr(j)= uwr(j) + aa
               ener(5)=ener(5) + aa
               
               aa= real(u3c(0,k)*conjg(u2c(0,k)))
               vwr(j)= vwr(j) + aa
               ener(6)=ener(6) + aa
               
               aa= o1c(0,k)*conjg(o1c(0,k))
               w1p(j)= w1p(j) + aa
               ener(7)=ener(7) + aa
               
               aa = o2c(0,k)*conjg(o2c(0,k))
               w2p(j)= w2p(j) + aa
               ener(8)=ener(8) + aa
               
               aa = o3c(0,k)*conjg(o3c(0,k))
               w3p(j)= w3p(j) + aa
               ener(9)=ener(9) + aa
               
c     dissipation  ----------------
               
               aa =  bet2(k) * 
     &              ( u1c(0,k)*conjg(u1c(0,k)) +
     &              u2c(0,k)*conjg(u2c(0,k)) + 
     &              u3c(0,k)*conjg(u3c(0,k)) ) +
     &              rhvc(0,k,j)*conjg(rhvc(0,k,j))
               
               cc = o1c(0,k) + xbet(k)*u2c(0,k)
               aa = aa + cc*conjg(cc) 
               cc = o3c(0,k) 
               aa = aa + cc*conjg(cc) 
               
               ep(j) = ep(j) + aa
               
               do i=1,mx1
                  aa = 2.*u1c(i,k)*conjg(u1c(i,k))
                  up(j) = up(j) + aa
                  ener(1)=ener(1) + aa
                  
                  aa= 2.*u2c(i,k)*conjg(u2c(i,k))
                  vp(j) = vp(j) + aa
                  ener(2)=ener(2) + aa
                  
                  aa= 2.*u3c(i,k)*conjg(u3c(i,k))
                  wp(j) = wp(j) + aa
                  ener(3)=ener(3) + aa
                  
                  aa= 2.*real(u1c(i,k)*conjg(u2c(i,k)))
                  uvr(j)= uvr(j) + aa
                  ener(4)=ener(4) + abs(aa)
                  
                  aa= 2.*real(u1c(i,k)*conjg(u3c(i,k)))
                  uwr(j)= uwr(j) + aa
                  ener(5)=ener(5) + abs(aa)
                  
                  aa= 2.*real(u3c(i,k)*conjg(u2c(i,k)))
                  vwr(j)= vwr(j) + aa
                  ener(6)=ener(6) + abs(aa)
                  
                  aa= 2.*o1c(i,k)*conjg(o1c(i,k))
                  w1p(j)= w1p(j) + aa
                  ener(7)=ener(7) + aa
                  
                  aa = 2.*o2c(i,k)*conjg(o2c(i,k))
                  w2p(j)= w2p(j) + aa
                  ener(8)=ener(8) + aa
                  
                  aa = 2.*o3c(i,k)*conjg(o3c(i,k))
                  w3p(j)= w3p(j) + aa
                  ener(9)=ener(9) + aa
                  
c              dissipation  ----------------
                  
                  aa = ( alp2(i) + bet2(k) ) * 
     &                 ( u1c(i,k)*conjg(u1c(i,k)) +
     &                 u2c(i,k)*conjg(u2c(i,k)) + 
     &                 u3c(i,k)*conjg(u3c(i,k)) ) +
     &                 rhvc(i,k,j)*conjg(rhvc(i,k,j) ) 
                  
                  cc = o1c(i,k) + xbet(k)*u2c(i,k)
                  aa = aa + cc*conjg(cc) 
                  cc = o3c(i,k) - xalp(i)*u2c(i,k) 
                  aa = aa + cc*conjg(cc)
                  
                  ep(j) = ep(j) + 2.*aa 
                  
               enddo
            enddo
            
c     c --------------- add this plane energy
            
            hyy = hy(j)
            do kk = 1,9
               uner(kk) = uner(kk) + ener(kk)*hyy
            enddo
            
c     ------------ means
            
            um(j) = um(j)+u1c(0,0)
            vm(j) = vm(j)+u2c(0,0)
            wm(j) = wm(j)+u3c(0,0)
            w1m(j)= w1m(j)+o1c(0,0)
            w2m(j)= w2m(j)+o2c(0,0)
            w3m(j)= w3m(j)+o3c(0,0)
            
c--------------update nacum just once !!!
            
            if (j.eq.je)   nacum = nacum+1
            
            if (myid.eq.0) then
               Wx0a=Wx0a+o1c(0,0)
               Wz0a=Wz0a+o3c(0,0)
            endif
         endif
         
         if (rkstep.eq.1) then
            
c-------------compute vorticity & Re stress  at walls
            if (j.eq.1) then
               Wx0 = o1c(0,0)
               Wz0 = o3c(0,0)
            endif
            if (j.eq.my) then
               WxL = o1c(0,0)
               WzL = o3c(0,0)
            endif
            
         endif

         ! compute velocity gradient tensor
         if(collectData .and. rkstep==1) then
             call compute_velocity_gradient_tensor(u1c, u2c, u3c,
     &           o1c, o3c, xalp, xbet, j, jb)
         endif

         ! save velocity and vorticity fields if required
         ! note: the forcing fields are saved elsewhere, see
         ! save_vorticity_forcing_at_plane_to_buffer (in hvhg)
         ! save_velocity_forcing_to_buffer (in cross1)
         if (collectFlowfield .and. rkstep==1) then
           call save_velocity_at_wallparallel_plane_to_buffer(u1c,
     &       u2c, u3c, j, jb)
           call save_vorticity_at_wallparallel_plane_to_buffer(o1c,
     &       o2c, o3c, j, jb)
         endif
         
c     /*      Move everything to PPP          */
         
c     ------- substract umax/2 to u00 to increase dt !!!!
         u1c(0,0) = u1c(0,0) - a0 
         
         call fourxz(u1c,u1r,1,1) ! u
         call fourxz(u2c,u2r,1,1) ! v
         call fourxz(u3c,u3r,1,1) ! w
         call fourxz(o1c,o1r,1,1) !omega_1
         call fourxz(o2c,o2r,1,1) !omega_2
         call fourxz(o3c,o3r,1,1) !omega_3

         ! collect visualization data
         if (collectData .and. rkstep==1) then
             call save_velocity_fields(j, jb, u2r)
         endif
         
         
c     ----------- triple products 
         
         if (rkstep.eq.1) then
            
            
            if (istati.eq.1) then
               do k = 1,mgalz
                  do i=1,mgalx 
                     aa = u2r(i,k)
                     uuv(j) = uuv(j) +aa*u1r(i,k)**2  
                     wwv(j) = wwv(j) +aa*u3r(i,k)**2  
                     vvv(j) = vvv(j) +aa**3
                  enddo 
               enddo 
            endif
            
            if (j.eq.1) then
               uv0=0.
               do k=1,mgalz
                  do i=1,mgalx
                     uv0 = uv0 + u1r(i,k)*u2r(i,k)
                  enddo
               enddo
               uv0= uv0/(mgalz*mgalx)
            endif
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
         
         uggg=0
         
         if (icfl.eq.1.and.rkstep.eq.1) then
            
c     /*      estimates maximum time step     */
c     /********************************************************************/
c     /*                                                                  */
c     /*    estimates spatial eigenvalues and computes maximum            */
c     /*    time step.                                                    */
c     /*                                                                  */
c     /*  input :                                                         */
c     /*    rhv,rhg,phi : velocities in (phys-phys-phys) plane            */
c     /*                                                                  */
c     /********************************************************************/
            hyy = hy(j)
            do k=1,mgalz
               do i=1,mgalx
                  cflx = max(cflx,abs(u1r(i,k)) )
                  cfly = max(cfly,abs(u2r(i,k))/hyy )
                  cflz = max(cflz,abs(u3r(i,k)) )
               enddo
            enddo
            
         endif
         
c     /* rhg= H1 = v.omega_3 - w.omega_2 (F-F-P)  */
c     /* phi= H2 = w.omega_1 - u.omega_3 (F-F-P)  */
c     /* rhv= H3 = u.omega_2 - v.omega_1 (F-F-P)  */
         
c     /********************************************************************/
c     /*                                                                  */
c     /*         computes u X  omega                                      */
c     /*                                                                  */
c     /*                                                                  */
c     /********************************************************************/
         
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
         
c---------------------at this point
c---------------------u1 : H3
c---------------------u3 : H2
c---------------------u2 : H1
         
c---------------------back to F-F-T
         
         call fourxz(u1c,u1r,-1,1)
         call fourxz(u2c,u2r,-1,1)
         call fourxz(u3c,u3r,-1,1)
         
         
c     /* saves coefficients for Kx = Kz = 0    */
         
         rf0u(j)=real(u2c(0,0))
         rf0w(j)=real(u1c(0,0))
         
c     =====
c     /*   u2 = - dH3/dx + dH1/dz      */
c     /*   o1 = dH1/dx + dH3/dz      */
         
         do k=0,mz1
            do i=0,mx1
               o1c(i,k) = xalp(i)*u2c(i,k)+xbet(k)*u1c(i,k)
               u2c(i,k) = -xalp(i)*u1c(i,k)+xbet(k)*u2c(i,k)
            enddo
         enddo
         
c     /*   rhv = d^2 H2/dx^2 + d^H2/dz^2      */
         
         do k=0,mz1
            do i=0,mx1
               u1c(i,k) = - u3c(i,k)*(alp2(i)+bet2(k))
            enddo
         enddo
         
c======
c     --------------------- copy planes into outputs
         do k=0,mz1
            do i=0,mx1
               rhvc(i,k,j)=u1c(i,k)
               rhgc(i,k,j)=u2c(i,k)
               ome1c(i,k,j)=o1c(i,k)
            enddo
         enddo

         ! save vorticity forcing if required
         ! note: the velocity, vorticity and velocity forcing fields
         ! are saved elsewhere, see
         ! save_velocity_at_wallparallel_plane_to_buffer (in hvhg)
         ! save_vorticity_at_wallparallel_plane_to_buffer (in hvhg)
         ! save_velocity_forcing_to_buffer (in cross1)
         if (collectFlowfield .and. rkstep==1) then
           call save_vorticity_forcing_at_plane_to_buffer(rhgc(:,:,j),
     &       j, jb)
         endif
         
c     -------------------------- finishes the y loop
      ENDDO
      
c     ----------------------- some things have to be done after the y loop
      
c     ---------------------- adds up the total energy
      
c     /********************************************************************/
c     /*        computes the energy
c     /********************************************************************/
      if (istati.eq.1.and.rkstep.eq.1.and.ihist.eq.1) then
         call MPI_ALLREDUCE(uner,ener,9,MPI_REAL,MPI_SUM,
     .        MPI_COMM_WORLD,ierr)
         
         do i=1,9
            ener(i)=sqrt(abs(ener(i)))
         enddo
         
         
      endif
      
c--------------------computes Deltat
      
      if (icfl.eq.1.and.rkstep.eq.1) then
         
         
         cflx = cflx*hxalp
         cflz = cflz*hzalp
         cfl0 = max(cflx,max(cfly,cflz))
         
         call MPI_ALLREDUCE(cfl0,reigmx1,1,MPI_REAL,MPI_MAX,
     .        MPI_COMM_WORLD,ierr)
         
         if (reigmx1.lt.1e-1)  uggg=1
         
         
         Deltat=CFL/reigmx1
         dtr=Re/Deltat
         
         if (uggg.ne.0) then
            write(*,*) 'UGGG', uggg,myid,ihist,jb,je,kb,ke
         endif
         
      endif
      
c                          /* saves coefficients for Kx = Kz = 0    */
c     /cccccccccccccccccccccccccccccccccccccccccccccccc/
c     MODULE FOR SP2 MPI uses SENDRECV               c/
c     sends a block (jb:je) everybody and receive    c/
c     from everybody                                 c/
c     /cccccccccccccccccccccccccccccccccccccccccccccccc/

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

c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
c     at this point: ome1, rhv and rhg are the outs
c     they must be trasformed from xy-xz before completion
c     o1: dH1/dx + dH3/dz
c     u2: -dH3/dx + dH1/dz
c     u1: d^2H2/dx^2 + d^2H2/dz^2
c     c c c c c c c c c c cc c c c c c c c c c c cc c c c c c c c c c c c
      
      end
      
