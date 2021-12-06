c/********************************************************************/
c/*                                                                  */
c/*     Solve the problem  v'' - rK v = phi                          */
c/*             with v(-1) = bcb  and v(+1) = bct                    */
c/*                                                                  */
c/*  Input:                                                          */
c/*     phi: forcing vector phi                                      */
c/*     bcb: bottom boundary condition                               */
c/*     bct: top boundary condition                                  */
c/*                                                                  */
c/* Output:                                                          */
c/*       v: solution of the problem                                 */
c/*    dvdy: dv/dy                                                   */
c/*                                                                  */
c/* The input and output matrices are size 1:2,1:my,                 */
c/* first index of 1 indicates real part and                         */
c/* first index of 2 indicates imaginary part                        */
c/********************************************************************/
      subroutine Lapvdv(phi,v,dvdy,rK,bcb,bct)
      use matrices,only: dt21,dt22,dt12,prem1
      
      implicit none
      include "ctes3D"
      integer n,i,j
      real*4 phi(2,my),v(2,my),dvdy(2,my)
      real*8 rK
      
      real*4 zero
      complex*8 bcb,bct
      real*8 wk1(5,my)
      real*8 fwk1(my),fwk2(my),fwk(my)
      
      real*8  fmap,y2
      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/
      
      zero = 0.
      
      if(rK.eq.0) then
         do j=1,my
            v   (1,j) =  zero
            v   (2,j) =  zero
            dvdy(1,j) =  zero
            dvdy(2,j) =  zero
         enddo
      else


c     -----------------------------------------------------------------
c     Prepare wk1 
         wk1(1,1) = 0d0
         wk1(2,1) = 0d0
         wk1(3,1) = 1d0
         wk1(4,1) = 0d0
         wk1(5,1) = 0d0
         do j=2,my-1
            do i=1,5
               wk1(i,j)=dt22(i,j)-rk*dt21(i,j)
            enddo
         enddo
         wk1(1,my) = 0d0
         wk1(2,my) = 0d0
         wk1(3,my) = 1d0
         wk1(4,my) = 0d0
         wk1(5,my) = 0d0
         
         call bandec5(wk1,my)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve v'' - rk*v = phi  with  v(-1) = bcb, v(+1) = bct
         ! Real part
         do j=1,my
            fwk(j) = phi(1,j)
         enddo
         
         fwk1(1)=real(bcb)
         fwk1(2)=dt21(2,2)*fwk(1)+dt21(3,2)*fwk(2)+dt21(4,2)*fwk(3)+
     .        dt21(5,2)*fwk(4)
         do j=3,my-2
            fwk1(j)=dt21(1,j)*fwk(j-2)
            do i=2,5
               fwk1(j)=fwk1(j)+dt21(i,j)*fwk(j-3+i)
            enddo
         enddo
         fwk1(my-1)=dt21(1,my-1)*fwk(my-3)+dt21(2,my-1)*fwk(my-2)+
     .        dt21(3,my-1)*fwk(my-1)+dt21(4,my-1)*fwk(my) 
         fwk1(my)=real(bct)
         
         call banbks5(wk1,my,fwk1)
         
         ! Imaginary part
         do j=1,my
            fwk(j) = phi(2,j)
         enddo
         fwk2(1)=aimag(bcb)
         fwk2(2)=dt21(2,2)*fwk(1)+dt21(3,2)*fwk(2)+dt21(4,2)*fwk(3)+
     .        dt21(5,2)*fwk(4)
         do j=3,my-2
            fwk2(j)=dt21(1,j)*fwk(j-2)
            do i=2,5
               fwk2(j)=fwk2(j)+dt21(i,j)*fwk(j-3+i)
            enddo
         enddo
         fwk2(my-1)=dt21(1,my-1)*fwk(my-3)+dt21(2,my-1)*phi(2,my-2)+
     .        dt21(3,my-1)*fwk(my-1)+dt21(4,my-1)*phi(2,my)
         fwk2(my)=aimag(bct)
         
         call banbks5(wk1,my,fwk2)
         
         ! Copy data to v
         do j=1,my
            v(1,j) = fwk1(j)
            v(2,j) = fwk2(j)
         enddo
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute dv/dy
         ! Real Part:
         fwk(1)= dt12(3,1)*fwk1(1)+dt12(4,1)*fwk1(2) + 
     .       dt12(5,1)*fwk1(3)         
         fwk(2)= dt12(2,2)*fwk1(1) +
     .       dt12(3,2)*fwk1(2)+dt12(4,2)*fwk1(3) +
     .       dt12(5,2)*fwk1(4)
         do j=3,my-2
            fwk(j)=dt12(1,j)*fwk1(j-2)
            do i=2,5
               fwk(j)=fwk(j) + dt12(i,j)*fwk1(i+j-3)
            enddo
         enddo
         fwk(my-1)=dt12(1,my-1)*fwk1(my-3)+dt12(2,my-1)*fwk1(my-2)+
     .        dt12(3,my-1)*fwk1(my-1)+dt12(4,my-1)*fwk1(my)
         fwk(my)=  dt12(1,my)*fwk1(my-2)+dt12(2,my)*fwk1(my-1)  +
     .        dt12(3,my)*fwk1(my) 
        
         call banbks5(prem1,my,fwk)
        
         do j=1,my
            dvdy(1,j)=fwk(j)*fmap(j)
         enddo

         ! Imaginary Part
         fwk(1)= dt12(3,1)*fwk2(1)+dt12(4,1)*fwk2(2) + 
     .        dt12(5,1)*fwk2(3)         
         fwk(2)= dt12(2,2)*fwk2(1) +
     .        dt12(3,2)*fwk2(2)+dt12(4,2)*fwk2(3) +
     .        dt12(5,2)*fwk2(4)
         do j=3,my-2
            fwk(j)=dt12(1,j)*fwk2(j-2)
            do i=2,5
               fwk(j)=fwk(j) + dt12(i,j)*fwk2(i+j-3)
            enddo
         enddo
         fwk(my-1)=dt12(1,my-1)*fwk2(my-3)+dt12(2,my-1)*fwk2(my-2)+
     .        dt12(3,my-1)*fwk2(my-1)+dt12(4,my-1)*fwk2(my)
         fwk(my)=  dt12(1,my)*fwk2(my-2)+dt12(2,my)*fwk2(my-1)+
     .        dt12(3,my)*fwk2(my)
        
         call banbks5(prem1,my,fwk)
        
         do j=1,my
            dvdy(2,j)=fwk(j)*fmap(j)
         enddo
c     -----------------------------------------------------------------
      endif
      
      end
      
      
      
c/********************************************************************/
c/*                                                                  */
c/*     Solve the problem  v'' - rK v = phi                          */
c/*             with v(-1) = bcb  and v(+1) = bct                    */
c/*                                                                  */
c/*  Input:                                                          */
c/*     phi: forcing vector phi                                      */
c/*     bcb: bottom boundary condition                               */
c/*     bct: top boundary condition                                  */
c/*                                                                  */
c/* Output:                                                          */
c/*       v: solution of the problem                                 */
c/*                                                                  */
c/* The input and output matrices are size 1:2,1:my,                 */
c/* first index of 1 indicates real part and                         */
c/* first index of 2 indicates imaginary part                        */
c/*                                                                  */
c/* This is the same as Lapvdv,                                      */
c/* without computing the first derivative                           */
c/********************************************************************/
      subroutine Lapv(phi,v,rk,bcb,bct)
      use matrices,only: dt22,dt21
      
      implicit none
      include "ctes3D"
      integer n,i,j
      real*4 phi(2,my),v(2,my),dvdy(2,my)
      real*8 rk
      
      complex*8 bcb,bct
      real*8 wk1(5,my)
      real*8 fwk1(my),fwk2(my),fwk(my)
      
      if(rK.eq.0) then
         do j=1,my
            v   (1,j) = 0. 
            v   (2,j) = 0. 
         enddo
      else
c     -----------------------------------------------------------------
c     Prepare wk1
         wk1(1,1) = 0d0
         wk1(2,1) = 0d0
         wk1(3,1) = 1d0
         wk1(4,1) = 0d0
         wk1(5,1) = 0d0
         do j=2,my-1
            do i=1,5
               wk1(i,j)=dt22(i,j)-rk*dt21(i,j)
            enddo
         enddo
         wk1(1,my) = 0d0
         wk1(2,my) = 0d0
         wk1(3,my) = 1d0
         wk1(4,my) = 0d0
         wk1(5,my) = 0d0

         call bandec5(wk1,my)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve v'' - rk*v = phi  with  v(-1) = bcb, v(+1) = bct
         ! Real part
         do j=1,my
            fwk(j) = phi(1,j)
         enddo
         
         fwk1(1)=real(bcb)
         fwk1(2)=dt21(2,2)*fwk(1)+dt21(3,2)*fwk(2)+dt21(4,2)*fwk(3)+
     .        dt21(5,2)*fwk(4)
         do j=3,my-2
            fwk1(j)=dt21(1,j)*fwk(j-2)
            do i=2,5
               fwk1(j)=fwk1(j)+dt21(i,j)*fwk(j-3+i)
            enddo
         enddo
         fwk1(my-1)=dt21(1,my-1)*fwk(my-3)+dt21(2,my-1)*fwk(my-2)+
     .        dt21(3,my-1)*fwk(my-1)+dt21(4,my-1)*fwk(my) 
         fwk1(my)=real(bct)
         
         call banbks5(wk1,my,fwk1)
         
         ! Imaginary part
         do j=1,my
            fwk(j) = phi(2,j)
         enddo
         fwk2(1)=aimag(bcb)
         fwk2(2)=dt21(2,2)*fwk(1)+dt21(3,2)*fwk(2)+dt21(4,2)*fwk(3)+
     .        dt21(5,2)*fwk(4)
         do j=3,my-2
            fwk2(j)=dt21(1,j)*fwk(j-2)
            do i=2,5
               fwk2(j)=fwk2(j)+dt21(i,j)*fwk(j-3+i)
            enddo
         enddo
         fwk2(my-1)=dt21(1,my-1)*fwk(my-3)+dt21(2,my-1)*phi(2,my-2)+
     .        dt21(3,my-1)*fwk(my-1)+dt21(4,my-1)*phi(2,my)
         fwk2(my)=aimag(bct)
         
         call banbks5(wk1,my,fwk2)

         ! Copy data to v
         do j=1,my
            v(1,j) = fwk1(j)
            v(2,j) = fwk2(j)
         enddo
c     -----------------------------------------------------------------
      endif
      
      end
      
      
c/********************************************************************/
c/*  Solve for v and phi                                             */
c/*    for the next rk step at a particular kx kz pair               */
c/*                                                                  */
c/*  For phi:                                                        */
c/*    phi'' - rk1 phi = f                                           */
c/*      with boundary conditions determined by the dv/dy condition  */
c/*                                                                  */
c/*  For v:                                                          */
c/*    v'' -  rk2 v = phi                                            */
c/*      with     v(-1) = bcbv ,     v(+1) = bctv                    */
c/*      and  dv/dy(-1) = bcbdv, dv/dy(+1) = bctdv                   */
c/*                                                                  */
c/*  To implement the boundary conditions for dv/dy                  */
c/*  v and phi are decomposed into particular and homogeous solution */
c/*  then linearly conbined to give final answer                     */
c/*                                                                  */
c/* Input:                                                           */
c/*    f          : forcing for phi                                  */
c/*    rk1        : constant, should be kx^2+kz^2+Re/c/Deltat        */
c/*    rk2        : constant, should be kx^2+kz^2                    */
c/*    bcbv ,bctv : boundary conditions for v                        */
c/*    bcbdv,bctdv: boundary conditions for dv/dy                    */
c/* Output:                                                          */
c/*    phi, v, dvdy,                                                 */
c/*    solutions are size 2,my first index indicates real/imag part  */
c/*                                                                  */
c/* This function is the same as Lapsov,                             */
c/* without solving for omega                                        */
c/********************************************************************/
      subroutine Bilap(phi,v,dvdy,f,rk1,rk2,bcbv,bctv,bcbdv,bctdv)
      use matrices,only:dt21,dt22,dt12,prem1
      
      implicit none
      include "ctes3D"
      integer i,j,k,jj
      
      real*4 phi(2,my),v(2,my),dvdy(2,my),f(2,my)
      real*8 rk1,rk2
      real*8 phipr(my),phipi(my),vpr(my),vpi(my),dvpr(my),dvpi(my),
     .     phi1(my),v1(my),dv1(my)
      real*8 zero
      complex*8 bcbv,bctv,bcbdv,bctdv,zeroc
      real*8 det,Ar,Ai,Br,Bi
      real*8 wk1(5,my)
      
      real*8  fmap,y2
      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/
      
      
      zero = 0d0
      zeroc = (0.,0.)
      
c     -----------------------------------------------------------------
c     Prepare the wk1 matrix with rk1, 
c     which is the same for all phi
      wk1(1,1) = 0d0
      wk1(2,1) = 0d0
      wk1(3,1) = 1d0
      wk1(4,1) = 0d0
      wk1(5,1) = 0d0
      do j=2,my-1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rk1*dt21(i,j)
         enddo
      enddo
      wk1(1,my) = 0d0
      wk1(2,my) = 0d0
      wk1(3,my) = 1d0
      wk1(4,my) = 0d0
      wk1(5,my) = 0d0

      call bandec5(wk1,my)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve: phip'' - rk1*phip = f  with  phip(1) = phip(-1) = 0 
      ! Real part
      do j=1,my
         v1(j) = f(1,j) 
      enddo 
      phipr(1)=0d0
      phipr(2)=dt21(2,2)*v1(1)+dt21(3,2)*v1(2)+dt21(4,2)*v1(3)+
     .         dt21(5,2)*v1(4)
      do j=3,my-2
         phipr(j)=dt21(1,j)*v1(j-2)
         do i=2,5
            phipr(j)=phipr(j)+dt21(i,j)*v1(j-3+i)
         enddo
      enddo
      phipr(my-1)=dt21(1,my-1)*v1(my-3)+dt21(2,my-1)*v1(my-2)+
     .     dt21(3,my-1)*v1(my-1)+dt21(4,my-1)*v1(my)
      phipr(my)=0d0
      call banbks5(wk1,my,phipr)

      ! Imaginary part
      do j=1,my
         v1(j) = f(2,j) 
      enddo 
      phipi(1)=0d0
      phipi(2)=dt21(2,2)*v1(1)+dt21(3,2)*v1(2)+dt21(4,2)*v1(3)+
     .     dt21(5,2)*v1(4)
      do j=3,my-2
         phipi(j)=dt21(1,j)*v1(j-2)
         do i=2,5
            phipi(j)=phipi(j)+dt21(i,j)*v1(j-3+i)
         enddo
      enddo
      phipi(my-1)=dt21(1,my-1)*v1(my-3)+dt21(2,my-1)*v1(my-2)+
     .     dt21(3,my-1)*v1(my-1)+dt21(4,my-1)*v1(my)
      phipi(my)=0d0
      call banbks5(wk1,my,phipi)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve phi1'' - rk1*phi1 = 0  with  phi1(1) =0  phi1(-1)=1 
      phi1(1)=1d0
      do j=2,my
         phi1(j)=0d0
      enddo
      call banbks5(wk1,my,phi1)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Prepare the wk1 matrix with rk2, 
c     which is the same for all velocities
      wk1(1,1) = 0d0
      wk1(2,1) = 0d0
      wk1(3,1) = 1d0
      wk1(4,1) = 0d0
      wk1(5,1) = 0d0
      do j=2,my-1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rk2*dt21(i,j)
         enddo
      enddo
      wk1(1,my) = 0d0
      wk1(2,my) = 0d0
      wk1(3,my) = 1d0
      wk1(4,my) = 0d0
      wk1(5,my) = 0d0
      
      call bandec5(wk1,my)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve vp'' - rk2*vp = phip  with  vp(1) = bctv, vp(-1) = bcbv
      ! Real Part
      vpr(1)=real(bcbv)
      vpr(2) = dt21(2,2)*phipr(1)+dt21(3,2)*phipr(2)+
     .     dt21(4,2)*phipr(3)+dt21(5,2)*phipr(4)
      do j=3,my-2
         vpr(j)=dt21(1,j)*phipr(j-2)
         do i=2,5
            vpr(j)=vpr(j)+dt21(i,j)*phipr(j-3+i)
         enddo
      enddo
      vpr(my-1)=dt21(1,my-1)*phipr(my-3)+dt21(2,my-1)*phipr(my-2)+
     .     dt21(3,my-1)*phipr(my-1)+dt21(4,my-1)*phipr(my)
      vpr(my)=real(bctv)
      call banbks5(wk1,my,vpr)

      ! Imaginary Part
      vpi(1)=aimag(bcbv)
      vpi(2) = dt21(2,2)*phipi(1)+dt21(3,2)*phipi(2)+
     .     dt21(4,2)*phipi(3)+dt21(5,2)*phipi(4)
      do j=3,my-2
         vpi(j)=dt21(1,j)*phipi(j-2)
         do i=2,5
            vpi(j)=vpi(j)+dt21(i,j)*phipi(j-3+i)
         enddo
      enddo
      vpi(my-1)=dt21(1,my-1)*phipi(my-3)+dt21(2,my-1)*phipi(my-2)+
     .     dt21(3,my-1)*phipi(my-1)+dt21(4,my-1)*phipi(my)
      vpi(my)=aimag(bctv)
      call banbks5(wk1,my,vpi)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve v1'' - rk2*v1 = phi1  with  v1(1) = v1(-1) =0
      v1(1)= 0d0
      v1(2)=dt21(2,2)*phi1(1)+dt21(3,2)*phi1(2)+
     .     dt21(4,2)*phi1(3)+dt21(5,2)*phi1(4)
      do j=3,my-2
         v1(j)=dt21(1,j)*phi1(j-2)
         do i=2,5
            v1(j)=v1(j)+dt21(i,j)*phi1(j-3+i)
         enddo
      enddo
      v1(my-1)=dt21(1,my-1)*phi1(my-3)+dt21(2,my-1)*phi1(my-2)+
     .     dt21(3,my-1)*phi1(my-1)+dt21(4,my-1)*phi1(my)
      v1(my)=0d0
      call banbks5(wk1,my,v1)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute derivative of vp -> dvp
c     Note: fmap is not applied!!
         ! Real Part
      dvpr(1)=dt12(3,1)*vpr(1)+dt12(4,1)*vpr(2) + 
     .     dt12(5,1)*vpr(3)        
      dvpr(2)=dt12(2,2)*vpr(1) +
     .     dt12(3,2)*vpr(2)+dt12(4,2)*vpr(3) +
     .     dt12(5,2)*vpr(4)
      do j=3,my-2
         dvpr(j)=dt12(1,j)*vpr(j-2)
         do i=2,5
            dvpr(j)=dvpr(j) + dt12(i,j)*vpr(i+j-3)
         enddo
      enddo
      dvpr(my-1)=dt12(1,my-1)*vpr(my-3)+dt12(2,my-1)*vpr(my-2)+
     .     dt12(3,my-1)*vpr(my-1)+dt12(4,my-1)*vpr(my)
      dvpr(my)=  dt12(1,my)*vpr(my-2)  +dt12(2,my)*vpr(my-1)+
     .     dt12(3,my)*vpr(my)
      call banbks5(prem1,my,dvpr)
      
      ! Imaginary Part
      dvpi(1)=dt12(3,1)*vpi(1)+dt12(4,1)*vpi(2) + 
     .     dt12(5,1)*vpi(3)        
      dvpi(2)=dt12(2,2)*vpi(1) +
     .     dt12(3,2)*vpi(2)+dt12(4,2)*vpi(3) +
     .     dt12(5,2)*vpi(4)
      do j=3,my-2
         dvpi(j)=dt12(1,j)*vpi(j-2)
         do i=2,5
            dvpi(j)=dvpi(j) + dt12(i,j)*vpi(i+j-3)
         enddo
      enddo
      dvpi(my-1)=dt12(1,my-1)*vpi(my-3)+dt12(2,my-1)*vpi(my-2)+
     .     dt12(3,my-1)*vpi(my-1)+dt12(4,my-1)*vpi(my)
      dvpi(my)=  dt12(1,my)*vpi(my-2)  +dt12(2,my)*vpi(my-1)+
     .     dt12(3,my)*vpi(my)
      call banbks5(prem1,my,dvpi)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute derivative of v1 -> dv1
c     Note: fmap is not applied!!
      dv1(1) = dt12(3,1)*v1(1)+dt12(4,1)*v1(2)+
     .     dt12(5,1)*v1(3)        
      dv1(2) = dt12(2,2)*v1(1)+
     .     dt12(3,2)*v1(2)+dt12(4,2)*v1(3)+
     .     dt12(5,2)*v1(4)
      do j=3,my-2
         dv1(j)=dt12(1,j)*v1(j-2)
         do i=1,5
            dv1(j) = dv1(j) + dt12(i,j)*v1(i+j-3)
         enddo
      enddo
      dv1(my-1)= dt12(1,my-1)*v1(my-3)+dt12(2,my-1)*v1(my-2)+
     .     dt12(3,my-1)*v1(my-1)+dt12(4,my-1)*v1(my)  
      dv1(my)=   dt12(1,my)*v1(my-2)  +dt12(2,my)*v1(my-1)+
     .     dt12(3,my)*v1(my)   
      call banbks5(prem1,my,dv1)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute the complex coefficients A and B 
c     so that the solution satisfies dvdy(1)=bctdv, dvdy(-1)=bcbdv
c     Note that fmap is applied here
      ! Determinant
      det=-dv1(1)*dv1(1) + dv1(my)*dv1(my)
      ! Real part
      Ar=(-dv1(my)*(dvpr(my)-real(bctdv)/fmap(my)) 
     .     + dv1(1)*(dvpr(1)-real(bcbdv)/fmap(1)))/det
      Br=(-dv1(1)*(dvpr(my)-real(bctdv)/fmap(my))
     .     +  dv1(my)*(dvpr(1)-real(bcbdv)/fmap(1)))/det
      ! Imag part
      Ai=(-dv1(my)*(dvpi(my)-aimag(bctdv)/fmap(my))
     .     + dv1(1)*(dvpi(1)-aimag(bcbdv)/fmap(1)))/det
      Bi=(-dv1(1)*(dvpi(my)-aimag(bctdv)/fmap(my))
     .     + dv1(my)*(dvpi(1)-aimag(bcbdv)/fmap(1)))/det
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Combine the particular solution with the homogeneous solutions
c     Note that fmap is applied here
      if (rk2.ne.0.) then
         do j=1,my
            phi(1,j) = phipr(j) +Ar*phi1(j)+Br*phi1(my-j+1)
            phi(2,j) = phipi(j) +Ai*phi1(j)+Bi*phi1(my-j+1)
            v(1,j)   = vpr(j) +Ar*v1(j)+Br*v1(my-j+1)
            v(2,j)   = vpi(j) +Ai*v1(j)+Bi*v1(my-j+1)
            dvdy(1,j)= fmap(j)*(dvpr(j) +Ar*dv1(j)-Br*dv1(my-j+1))
            dvdy(2,j)= fmap(j)*(dvpi(j) +Ai*dv1(j)-Bi*dv1(my-j+1))
         enddo
      else
         do j=1,my
            phi(1,j)  = phipr(j) +Ar*phi1(j)+Br*phi1(my-j+1)
            phi(2,j)  = phipi(j) +Ai*phi1(j)+Bi*phi1(my-j+1)
            v(1,j)    = zero
            v(2,j)    = zero
            dvdy(1,j) = zero
            dvdy(2,j) = zero
         enddo
      endif
c     -----------------------------------------------------------------
      end
      


c/********************************************************************/
c/*        Solve the problem   fwk'' - rK fwk = f2                   */
c/*                                                                  */
c/*                                with fwk(-1)= bcb                 */
c/*                                     fwk(+1)= bct                 */
c/*                                                                  */
c/* Input:                                                           */
c/*      f2      : forcing                                           */
c/*      rK      : independent constant.                             */
c/*      bct, bcb: top and bottom boundary conditions                */
c/* Output:                                                          */
c/*      fwk     : solution                                          */
c/********************************************************************/
      subroutine Lapv1(f2,fwk,rK,bcb,bct)
      use matrices,only: dt22,dt21
      
      implicit none
      include "ctes3D"
      integer i,j
      real*8 rK

      real*8 bcb,bct
      real*8 wk1(5,my) !!!,wk2(my,2)
      real*8 fwk(my),f2(my)

      real*8 fmap,y2
      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/

      ! Compute wk1 matrix: dt22 - rk*dt21 with boundary conditions
      wk1(1,1)=0d0
      wk1(2,1)=0d0
      wk1(3,1)=1d0
      wk1(4,1)=0d0
      wk1(5,1)=0d0
      do j=2,my-1
         do i=1,5
            wk1(i,j)=dt22(i,j)-rk*dt21(i,j)
         enddo
      enddo
      wk1(1,my)=0d0
      wk1(2,my)=0d0
      wk1(3,my)=1d0
      wk1(4,my)=0d0
      wk1(5,my)=0d0
      call bandec5(wk1,my)
      
      ! multiply forcing with dt21 and apply BC
      fwk(1)=bcb
      fwk(2)=dt21(2,2)*f2(1)+dt21(3,2)*f2(2)+dt21(4,2)*f2(3)+
     .     dt21(5,2)*f2(4)
      do j=3,my-2
         fwk(j)=dt21(1,j)*f2(j-2)
         do i=2,5
            fwk(j)=fwk(j)+dt21(i,j)*f2(j-3+i)
         enddo
      enddo
      fwk(my-1)=dt21(1,my-1)*f2(my-3)+dt21(2,my-1)*f2(my-2)+
     .     dt21(3,my-1)*f2(my-1)+dt21(4,my-1)*f2(my)
      fwk(my)=bct
      
      ! solve linear system
      call banbks5(wk1,my,fwk)
      
      end


c/********************************************************************/
c/*  Solve for v, omega and phi                                      */
c/*    for the next rk step at a particular kx kz pair               */
c/*                                                                  */
c/*  For omega:                                                      */
c/*    omega'' - rk1 omega = g                                       */
c/*      with omega(-1) = bcbo, omega(+1) = bcto                     */
c/*                                                                  */
c/*  For phi:                                                        */
c/*    phi'' - rk1 phi = f                                           */
c/*      with boundary conditions determined by the dv/dy condition  */
c/*                                                                  */
c/*  For v:                                                          */
c/*    v'' -  rk2 v = phi                                            */
c/*      with     v(-1) = bcbv ,     v(+1) = bctv                    */
c/*      and  dv/dy(-1) = bcbdv, dv/dy(+1) = bctdv                   */
c/*                                                                  */
c/*  To implement the boundary conditions for dv/dy                  */
c/*  v and phi are decomposed into particular and homogeous solution */
c/*  then linearly conbined to give final answer                     */
c/*                                                                  */
c/* Input:                                                           */
c/*    f, g       : forcing for phi and omega                        */
c/*    rk1        : constant, should be kx^2+kz^2+Re/c/Deltat        */
c/*    rk2        : constant, should be kx^2+kz^2                    */
c/*    bcbv ,bctv : boundary conditions for v                        */
c/*    bcbdv,bctdv: boundary conditions for dv/dy                    */
c/*    bcbo ,bcto : boundary conditions for omega                    */
c/* Output:                                                          */
c/*    phi, v, dvdy, ome                                             */
c/*    solutions are size 2,my first index indicates real/imag part  */
c/********************************************************************/
      subroutine lapsov(phi,v,dvdy,f,ome,g,
     .     rk1,rk2,bcbv,bctv,bcbdv,bctdv,bcbo,bcto)
      use matrices,only:dt21,dt22,dt12,prem1
      
      implicit none
      include "ctes3D"
      integer i,j,k,jj
      
      real*4 phi(2,my),v(2,my),dvdy(2,my),f(2,my)
      real*4 ome(2,my),g(2,my)
      real*8 rk1,rk2
      real*8 phipr(my),phipi(my),vpr(my),vpi(my),dvpr(my),dvpi(my),
     .       phi1(my),v1(my),dv1(my)
      real*8 omer(my),omei(my)
      real*8 zero
      complex*8 bcbv,bctv,bcbdv,bctdv,zeroc,bcbo,bcto
      real*8 det,Ar,Ai,Br,Bi
      real*8 wk1(5,my)

      real*8  fmap,y2
      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/


      zero = 0d0
      zeroc = (0.,0.)

      if (rk2.eq.0) then 
         do j=1,my
            phi(1,j)  = zero
            phi(2,j)  = zero
            ome(1,j)  = zero
            ome(2,j)  = zero
            v(1,j)    = zero
            v(2,j)    = zero
            dvdy(1,j) = zero
            dvdy(2,j) = zero
         enddo
      else


c     First solve phi and omega


c     -----------------------------------------------------------------
c     Prepare the wk1 matrix with rk1, 
c     which is the same for all phi
         wk1(1,1) = 0d0
         wk1(2,1) = 0d0
         wk1(3,1) = 1d0
         wk1(4,1) = 0d0
         wk1(5,1) = 0d0
         do j=2,my-1
            do i=1,5
               wk1(i,j)=dt22(i,j)-rk1*dt21(i,j)
            enddo
         enddo
         wk1(1,my) = 0d0
         wk1(2,my) = 0d0
         wk1(3,my) = 1d0
         wk1(4,my) = 0d0
         wk1(5,my) = 0d0
         call bandec5(wk1,my)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve: phip'' - rk1*phip = f  with  phip(1) = phip(-1) = 0 
         ! Real part
         do j=1,my
            v1(j) = f(1,j) 
         enddo 
         phipr(1)=0d0
         phipr(2)=dt21(2,2)*v1(1)+dt21(3,2)*v1(2)+dt21(4,2)*v1(3)+
     .        dt21(5,2)*v1(4)
         do j=3,my-2
            phipr(j)=dt21(1,j)*v1(j-2)
            do i=2,5
               phipr(j)=phipr(j)+dt21(i,j)*v1(j-3+i)
            enddo
         enddo
         phipr(my-1)=dt21(1,my-1)*v1(my-3)+dt21(2,my-1)*v1(my-2)+
     .        dt21(3,my-1)*v1(my-1)+dt21(4,my-1)*v1(my)
         phipr(my)=0d0
         call banbks5(wk1,my,phipr)
         
         ! Imaginary part
         do j=1,my
            v1(j) = f(2,j) 
         enddo 
         phipi(1)=0d0
         phipi(2)=dt21(2,2)*v1(1)+dt21(3,2)*v1(2)+dt21(4,2)*v1(3)+
     .        dt21(5,2)*v1(4)
         do j=3,my-2
            phipi(j)=dt21(1,j)*v1(j-2)
            do i=2,5
               phipi(j)=phipi(j)+dt21(i,j)*v1(j-3+i)
            enddo
         enddo
         phipi(my-1)=dt21(1,my-1)*v1(my-3)+dt21(2,my-1)*v1(my-2)+
     .        dt21(3,my-1)*v1(my-1)+dt21(4,my-1)*v1(my)
         phipi(my)=0d0
         call banbks5(wk1,my,phipi)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve phi1'' - rk1*phi1 = 0  with  phi1(1) =0  phi1(-1)=1 
         phi1(1)=1d0
         do j=2,my
            phi1(j)=0d0
         enddo
         call banbks5(wk1,my,phi1)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve omer'' - rk1*omer = g  with  ome(-1)=bcbo,ome(1)=bcto
         ! Real part
         do j=1,my
            v1(j) = g(1,j) 
         enddo 
         omer(1)=real(bcbo)
         omer(2)=dt21(2,2)*v1(1)+dt21(3,2)*v1(2)+dt21(4,2)*v1(3)+
     .        dt21(5,2)*v1(4)
         do j=3,my-2
            omer(j)=dt21(1,j)*v1(j-2)
            do i=2,5
               omer(j)=omer(j)+dt21(i,j)*v1(j-3+i)
            enddo
         enddo
         omer(my-1)=dt21(1,my-1)*v1(my-3)+dt21(2,my-1)*v1(my-2)+
     .        dt21(3,my-1)*v1(my-1)+dt21(4,my-1)*v1(my)
         omer(my)=real(bcto)
         call banbks5(wk1,my,omer)
         
         ! Imaginary Part
         do j=1,my
            v1(j) = g(2,j) 
         enddo 
         omei(1)=aimag(bcbo)
         omei(2)=dt21(2,2)*v1(1)+dt21(3,2)*v1(2)+dt21(4,2)*v1(3)+
     .        dt21(5,2)*v1(4)
         do j=3,my-2
            omei(j)=dt21(1,j)*v1(j-2)
            do i=2,5
               omei(j)=omei(j)+dt21(i,j)*v1(j-3+i)
            enddo
         enddo
         omei(my-1)=dt21(1,my-1)*v1(my-3)+dt21(2,my-1)*v1(my-2)+
     .        dt21(3,my-1)*v1(my-1)+dt21(4,my-1)*v1(my)
         omei(my)=aimag(bcto)
         call banbks5(wk1,my,omei)
c     -----------------------------------------------------------------


c     Now solve for v


c     -----------------------------------------------------------------
c     Prepare the wk1 matrix with rk2, 
c     which is the same for all velocities
         wk1(1,1) = 0d0
         wk1(2,1) = 0d0
         wk1(3,1) = 1d0
         wk1(4,1) = 0d0
         wk1(5,1) = 0d0
         do j=2,my-1
            do i=1,5
               wk1(i,j)=dt22(i,j)-rk2*dt21(i,j)
            enddo
         enddo
         wk1(1,my) = 0d0
         wk1(2,my) = 0d0
         wk1(3,my) = 1d0
         wk1(4,my) = 0d0
         wk1(5,my) = 0d0
         call bandec5(wk1,my)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve vp'' - rk2*vp = phip  with vp(1) = bctv, vp(-1) = bcbv
         ! Real Part
         vpr(1)=real(bcbv)
         vpr(2) = dt21(2,2)*phipr(1)+dt21(3,2)*phipr(2)+
     .        dt21(4,2)*phipr(3)+dt21(5,2)*phipr(4)
         do j=3,my-2
            vpr(j)=dt21(1,j)*phipr(j-2)
            do i=2,5
               vpr(j)=vpr(j)+dt21(i,j)*phipr(j-3+i)
            enddo
         enddo
         vpr(my-1)=dt21(1,my-1)*phipr(my-3)+dt21(2,my-1)*phipr(my-2)+
     .        dt21(3,my-1)*phipr(my-1)+dt21(4,my-1)*phipr(my)
         vpr(my)=real(bctv)
         call banbks5(wk1,my,vpr)
         
         ! Imaginary Part
         vpi(1)=aimag(bcbv)
         vpi(2) = dt21(2,2)*phipi(1)+dt21(3,2)*phipi(2)+
     .        dt21(4,2)*phipi(3)+dt21(5,2)*phipi(4)
         do j=3,my-2
            vpi(j)=dt21(1,j)*phipi(j-2)
            do i=2,5
               vpi(j)=vpi(j)+dt21(i,j)*phipi(j-3+i)
            enddo
         enddo
         vpi(my-1)=dt21(1,my-1)*phipi(my-3)+dt21(2,my-1)*phipi(my-2)+
     .        dt21(3,my-1)*phipi(my-1)+dt21(4,my-1)*phipi(my)
         vpi(my)=aimag(bctv)
         call banbks5(wk1,my,vpi)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Solve v1'' - rk2*v1 = phi1  with  v1(1) = v1(-1) =0
         v1(1)= 0d0
         v1(2)=dt21(2,2)*phi1(1)+dt21(3,2)*phi1(2)+
     .        dt21(4,2)*phi1(3)+dt21(5,2)*phi1(4)
         do j=3,my-2
            v1(j)=dt21(1,j)*phi1(j-2)
            do i=2,5
               v1(j)=v1(j)+dt21(i,j)*phi1(j-3+i)
            enddo
         enddo
         v1(my-1)=dt21(1,my-1)*phi1(my-3)+dt21(2,my-1)*phi1(my-2)+
     .        dt21(3,my-1)*phi1(my-1)+dt21(4,my-1)*phi1(my)
         v1(my)=0d0
         call banbks5(wk1,my,v1)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute derivative of vp -> dvp
c     Note: fmap is not applied!!
         ! Real Part
         dvpr(1)=dt12(3,1)*vpr(1)+dt12(4,1)*vpr(2) + 
     .        dt12(5,1)*vpr(3)        
         dvpr(2)=dt12(2,2)*vpr(1) +
     .        dt12(3,2)*vpr(2)+dt12(4,2)*vpr(3) +
     .        dt12(5,2)*vpr(4)
         do j=3,my-2
            dvpr(j)=dt12(1,j)*vpr(j-2)
            do i=2,5
               dvpr(j)=dvpr(j) + dt12(i,j)*vpr(i+j-3)
            enddo
         enddo
         dvpr(my-1)=dt12(1,my-1)*vpr(my-3)+dt12(2,my-1)*vpr(my-2)+
     .        dt12(3,my-1)*vpr(my-1)+dt12(4,my-1)*vpr(my)
         dvpr(my)=  dt12(1,my)*vpr(my-2)  +dt12(2,my)*vpr(my-1)+
     .        dt12(3,my)*vpr(my)
         call banbks5(prem1,my,dvpr)

         ! Imaginary Part
         dvpi(1)=dt12(3,1)*vpi(1)+dt12(4,1)*vpi(2) + 
     .        dt12(5,1)*vpi(3)        
         dvpi(2)=dt12(2,2)*vpi(1) +
     .        dt12(3,2)*vpi(2)+dt12(4,2)*vpi(3) +
     .        dt12(5,2)*vpi(4)
         do j=3,my-2
            dvpi(j)=dt12(1,j)*vpi(j-2)
            do i=2,5
               dvpi(j)=dvpi(j) + dt12(i,j)*vpi(i+j-3)
            enddo
         enddo
         dvpi(my-1)=dt12(1,my-1)*vpi(my-3)+dt12(2,my-1)*vpi(my-2)+
     .        dt12(3,my-1)*vpi(my-1)+dt12(4,my-1)*vpi(my)
         dvpi(my)=  dt12(1,my)*vpi(my-2)  +dt12(2,my)*vpi(my-1)+
     .        dt12(3,my)*vpi(my)
         call banbks5(prem1,my,dvpi)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute derivative of v1 -> dv1
c     Note: fmap is not applied!!
         dv1(1) = dt12(3,1)*v1(1)+dt12(4,1)*v1(2)+
     .        dt12(5,1)*v1(3)        
         dv1(2) = dt12(2,2)*v1(1)+
     .        dt12(3,2)*v1(2)+dt12(4,2)*v1(3)+
     .        dt12(5,2)*v1(4)
         do j=3,my-2
            dv1(j)=dt12(1,j)*v1(j-2)
            do i=1,5
               dv1(j) = dv1(j) + dt12(i,j)*v1(i+j-3)
            enddo
         enddo
         dv1(my-1)= dt12(1,my-1)*v1(my-3)+dt12(2,my-1)*v1(my-2)+
     .        dt12(3,my-1)*v1(my-1)+dt12(4,my-1)*v1(my)  
         dv1(my)=   dt12(1,my)*v1(my-2)  +dt12(2,my)*v1(my-1)+
     .        dt12(3,my)*v1(my)   
         call banbks5(prem1,my,dv1)
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Compute the complex coefficients A and B 
c     so that the solution satisfies dvdy(1)=bctdv, dvdy(-1)=bcbdv
c     Note that fmap is applied here
         ! Determinant
         det=-dv1(1)*dv1(1) + dv1(my)*dv1(my)
         ! Real part
         Ar=(-dv1(my)*(dvpr(my)-real(bctdv)/fmap(my)) 
     .        + dv1(1)*(dvpr(1)-real(bcbdv)/fmap(1)))/det
         Br=(-dv1(1)*(dvpr(my)-real(bctdv)/fmap(my))
     .        +  dv1(my)*(dvpr(1)-real(bcbdv)/fmap(1)))/det
         ! Imag part
         Ai=(-dv1(my)*(dvpi(my)-aimag(bctdv)/fmap(my))
     .        + dv1(1)*(dvpi(1)-aimag(bcbdv)/fmap(1)))/det
         Bi=(-dv1(1)*(dvpi(my)-aimag(bctdv)/fmap(my))
     .        + dv1(my)*(dvpi(1)-aimag(bcbdv)/fmap(1)))/det
c     -----------------------------------------------------------------


c     -----------------------------------------------------------------
c     Combine the particular solution with the homogeneous solutions
c     Note that fmap is applied here
         do j=1,my
            ome(1,j) = omer(j)
            ome(2,j) = omei(j)
            phi(1,j) = phipr(j) +Ar*phi1(j)+Br*phi1(my-j+1)
            phi(2,j) = phipi(j) +Ai*phi1(j)+Bi*phi1(my-j+1)
            v(1,j)   = vpr(j) +Ar*v1(j)+Br*v1(my-j+1)
            v(2,j)   = vpi(j) +Ai*v1(j)+Bi*v1(my-j+1)
            dvdy(1,j)= fmap(j)*(dvpr(j) +Ar*dv1(j)-Br*dv1(my-j+1))
            dvdy(2,j)= fmap(j)*(dvpi(j) +Ai*dv1(j)-Bi*dv1(my-j+1))
         enddo
c     -----------------------------------------------------------------
         
      endif
      
      
      end

