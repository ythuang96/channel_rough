      subroutine chikj2jik(xy,xz,wk1,wk2,myid)
c     /cccccccccccccccccccccccccccccccccccccccccccccccc/
c     MODULE FOR SP2 MPI uses SENDRECV               c/
c     /cccccccccccccccccccccccccccccccccccccccccccccccc/
      implicit none
      include 'mpif.h'
      include "ctes3D"
      integer myid,icontrol,iresult,ilocal,icount

      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .     kbeg(0:numerop-1),kend(0:numerop-1),
     .     jb,je,kb,ke,mmy,mmz
      save /point/

      real*4 xy(mx,mz,*),xz(2*my,mx1+1,kb:*),wk1(*),
     &     wk2(mx,kb:ke,*)

      integer myslice
      common /MPI_datatype/myslice(0:numerop-1)
      save /MPI_datatype/

      integer istat(MPI_STATUS_SIZE),ierr
      integer iproc,nsetotr,ipoxz,mzbu,ntoto
      integer i,j,k,mmz2

      real*8 commtimer,transtimer,totaltimer
      common/timers/ commtimer, transtimer, totaltimer
      save/timers/

c     /cccccccccccccccccccccccccccccccccccccccccccccccc/
c     NOTE:
c     c  wk2 dimensioned at least (mx*mzp*myp)
c     /cccccccccccccccccccccccccccccccccccccccccccccccc/

c--------------------------checks NANS 1/2 -----------------------
      icontrol=0
      icount = 0
      ilocal = 0

      do j=1,je-jb+1
         do k=1,mz
            do i=1,mx
               if (xy(i,k,j).ne.xy(i,k,j)) then
                  icontrol=1
                  icount = icount + 1
                  ilocal = 1
               endif
            enddo
         enddo
      enddo

      if (ilocal.eq.1) then
         write(*,*) 'NaN in chikj2jik before MPI_SENDRECV,myid=,',
     .        myid,'   count=',icount
         call system("hostname")
      endif
c--------------------------end checks NANS 1/2 -----------------------


      mmz2=mmz*mx

      if (myid.eq.0) then
         commtimer = commtimer-MPI_WTIME()
      endif


c-----------------------------------------------------------------------
c     Exchange data between processors
c     Each porcessor will get a wk1 vector that contains all kx, all y
c     but only the kz that is assigned to this processor
c-----------------------------------------------------------------------
      do iproc=0,numerop-1

         if(iproc.ne.myid)then ! exchange data with different processors

            nsetotr=mmz2*(jend(iproc)-jbeg(iproc)+1)
            ipoxz=1+mmz2*(jbeg(iproc)-1)

            call MPI_SENDRECV(xy(1,kbeg(iproc),1),1,myslice(iproc),
     .           iproc,0,wk1(ipoxz),nsetotr,
     .           MPI_REAL,iproc,0,
     .           MPI_COMM_WORLD,istat,ierr)

         else

            mzbu=kend(iproc)-kbeg(iproc)+1
            ipoxz=1+mmz2*(jb-1)
            call pacy2z(wk1(ipoxz),xy,mzbu,kbeg(iproc),iproc)

         endif

      enddo

      if (myid.eq.0) then
         commtimer = commtimer+MPI_WTIME()
      endif


      if (myid.eq.0) then
         transtimer = transtimer-MPI_WTIME()
      endif

      do k=kb,ke
         do j=1,my
            do i=1,mx1+1
               xz(2*j-1,i,k) = wk2(2*i-1,k,j)
               xz(2*j  ,i,k) = wk2(2*i  ,k,j)
            enddo
         enddo
      enddo


      if (myid.eq.0) then
         transtimer = transtimer+MPI_WTIME()
      endif

c-----------------------------checks NANS 2/2-----------------------
      icount=0
      ilocal=0
      do k=kb,ke
         do i=1,mx1+1
            do j=1,2*my
               if (xz(j,i,k).ne.xz(j,i,k)) then
                  icontrol=1
                  ilocal=1
                  icount = icount + 1
               endif
            enddo
         enddo
      enddo

      if (ilocal.eq.1) then
         write(*,*) 'NaN in chikj2jik after MPI_SENDRECV,myid=,',
     .        myid,'   count:',icount
         call system("hostname")
      endif

      call MPI_ALLREDUCE(icontrol,iresult,1,MPI_INTEGER,
     .     MPI_MAX,MPI_COMM_WORLD,ierr)

      if (iresult.eq.1) then
         write(*,*) 'Finishing process ... ',myid
         stop
      endif
c--------------------------end checks NANS 2/2 -----------------------

      end


      subroutine chjik2ikj(xz,xy,wk1,wk2,myid)
c     /cccccccccccccccccccccccccccccccccccccccccccccccc/
c     MODULE FOR SP2 MPI uses SENDRECV               c/
c     sends a block (mx,mz,jb:je) to a (mx,kb:ke,my) c/
c     /cccccccccccccccccccccccccccccccccccccccccccccccc/
      implicit none
      include "mpif.h"
      include "ctes3D"

      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .     kbeg(0:numerop-1),kend(0:numerop-1),
     .     jb,je,kb,ke,mmy,mmz
      save /point/
      real*4 xz(2*my,mx1+1,kb:*),xy(mx,mz,*),wk1(*),wk2(mx,kb:ke,*)
      integer myid,icontrol,ilocal,icount,iresult


      integer myslice
      common /MPI_datatype/myslice(0:numerop-1)
      save /MPI_datatype/

      integer istat(MPI_STATUS_SIZE),ierr
      integer iproc,nsetots,ipoxz,mzbu,ntoto
      integer i,j,k,mmy2,mmz2

      real*8 commtimer,transtimer,totaltimer
      common/timers/ commtimer, transtimer, totaltimer
      save/timers/

      mmz2=mmz*mx
      mmy2=mmy*mx

c--------------------------checks NANS 1/2 -----------------------
      icontrol=0
      icount = 0
      ilocal = 0

      do k=kb,ke
         do i=1,mx1+1
            do j=1,2*my
               if (xz(j,i,k).ne.xz(j,i,k)) then
                  icontrol=1
                  icount = icount + 1
                  ilocal = 1
               endif
            enddo
         enddo
      enddo

      if (ilocal.eq.1) then
         write(*,*) 'NaN in chjik2ikj before MPI_SENDRECV,myid=,',
     .        myid,'   count=',icount
         call system("hostname")
      endif
c--------------------------end checks NANS 1/2 -----------------------

      if (myid.eq.0) then
         transtimer = transtimer-MPI_WTIME()
      endif


      do k=kb,ke
         do j=1,my
            do i=1,mx1+1
               wk2(2*i-1,k,j) = xz(2*j-1,i,k)
               wk2(2*i  ,k,j) = xz(2*j  ,i,k)
            enddo
         enddo
      enddo

      if (myid.eq.0) then
         transtimer = transtimer+MPI_WTIME()
         commtimer = commtimer-MPI_WTIME()
      endif



      do iproc=0,numerop-1

         if(iproc.ne.myid)then

            nsetots=(jend(iproc)-jbeg(iproc)+1)*mmz2
            ipoxz=1+mmz2*(jbeg(iproc)-1)

            call MPI_SENDRECV(wk1(ipoxz),nsetots,
     .           MPI_REAL,iproc,0,xy(1,kbeg(iproc),1),
     .           1,myslice(iproc),
     .           iproc,0,MPI_COMM_WORLD,istat,ierr)

         else

            mzbu=kend(iproc)-kbeg(iproc)+1
            ipoxz=1+mmz2*(jb-1)
            call unpacz2y(wk1(ipoxz),xy,mzbu,kbeg(iproc),iproc)

         endif

      enddo

      if (myid.eq.0) then
         commtimer = commtimer+MPI_WTIME()
      endif

c--------------------------checks NANS 2/2 -----------------------
      icount = 0
      ilocal = 0

      do j=1,je-jb+1
         do k=1,mz
            do i=1,mx
               if (xy(i,k,j).ne.xy(i,k,j)) then
                  icontrol=1
                  icount = icount + 1
                  ilocal = 1
               endif
            enddo
         enddo
      enddo

      if (ilocal.eq.1) then
         write(*,*) 'NaN in chjik2ikj after MPI_SENDRECV,myid=,',
     .        myid,'   count=',icount
         call system("hostname")
      endif

      call MPI_ALLREDUCE(icontrol,iresult,1,MPI_INTEGER,
     .     MPI_MAX,MPI_COMM_WORLD,ierr)

      if (iresult.eq.1) then
         write(*,*) 'Finishing process ... ',myid
         stop
      endif
c--------------------------end checks NANS 2/2 -----------------------

      end


      subroutine unpacz2y(xyi,xyo,mzbu,kb1,iproc)
      implicit none
      include "ctes3D"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This is part of the transform from slices in z to slices in y
c     (reverse of pacy2z)
c     Input: xyi, this contains all kx, wavenumbers, and all y grid, 
c                 but only a few slices in kz
c     Output: xyo, this contains all kx, all kz, but only a few slicess
c                  in y that are assigned to proccesor #iproc
c     In otherwords, this function extract the few slices in y that are
c                    assigned to processor #iproc from xyi (which 
c                    contains all y)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .     kbeg(0:numerop-1),kend(0:numerop-1),
     .     jb,je,kb,ke,mmy,mmz
      save /point/
      integer mzbu,kb1,iproc
      real*4 xyi(mx,kb1:kb1+mzbu-1,*),xyo(mx,mz,*)
      integer  i,j,k

      do j=1,mmy
         do k=kbeg(iproc),kend(iproc)
            do i=1,mx

               xyo(i,k,j)=xyi(i,k,j)

            enddo
         enddo
      enddo

      end


      subroutine pacy2z(xyo,xyi,mzbu,kb1,iproc)
      implicit real*4(a-h,o-z)
      include "ctes3D"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This is part of the transform from slices in y to slices in z
c     Input: xyi, this contains all kx, kz wavenumbers, but only a few 
c                 slices in y
c     Output: xyo, this contains all kx, by only a few slices in kz 
c                  (the slices of kz that is assigned to processor #iproc)
c                  it contains all y point that is contained in xyi
c     In otherwords, this function extract the few slices in kz that are
c                    assigned to processor #iproc from xyi (which 
c                    contains all kz)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      integer jbeg,jend,kbeg,kend,jb,je,kb,ke,mmy,mmz
      common /point /jbeg(0:numerop-1),jend(0:numerop-1),
     .     kbeg(0:numerop-1),kend(0:numerop-1),
     .     jb,je,kb,ke,mmy,mmz

      save /point/
      integer mzbu,kb1,iproc
      dimension xyi(mx,mz,*),xyo(mx,kb1:kb1+mzbu-1,*)
      integer i,j,k

      do j=1,mmy
         do k=kbeg(iproc),kend(iproc)
            do i=1,mx

               xyo(i,k,j)=xyi(i,k,j)

            enddo
         enddo
      enddo

      end

      subroutine pointers(jb,je,kb,ke)
c     ===============================================
c     jjs,  aug/2001  (bug in old version)
c     ===============================================
c     ===============================================
c     This function divides the y and z grid for all
c     the processors
c     If the number of y, z grid are not divisible by
c     the number of processors, some processors will
c     have an extra grid point
c     jb, kb are arrays with the y,z start index for
c     each processor
c     je, ke are the end index
c     ===============================================
      implicit none

      include "ctes3D"
      integer jb(numerop),je(numerop),kb(numerop),ke(numerop),n,n1,n2

      n1=my/numerop
      n2=my-numerop*n1

      jb(1)=1
      do n=1,n2
         je(n)  = jb(n)+n1
         jb(n+1)= je(n)+1
      enddo
      do n=n2+1,numerop-1
         je(n)=jb(n)+n1-1
         jb(n+1)= je(n)+1
      enddo
      je(numerop)=jb(numerop)+n1-1


      n1=mz/numerop
      n2=mz-numerop*n1

      kb(1)=1
      do n=1,n2
         ke(n)  = kb(n)+n1
         kb(n+1)= ke(n)+1
      enddo
      do n=n2+1,numerop-1
         ke(n)=kb(n)+n1-1
         kb(n+1)= ke(n)+1
      enddo
      ke(numerop)=kb(numerop)+n1-1

      end
