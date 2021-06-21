***************************************************************************
c	PAQUETE DE SUBRUTINAS PARA DIFERENCIAS FINITAS COMPACTAS
c	O.F, nov 2002
c
c	updated dic 2002
c	touched ene 2003  ---------> this subroutines don't pass worksp
c       d11,d12,d21,d22 & solver real*8
c***************************************************************************


c************************************************
c*      h:      coordinates                      *
c*      n:      number of points                *
c*      p:      flag:  0->unif      1->no-unif  *
c*		p=1 -> tanh			*
c*		p=2 -> sin			*
c*      c:      narrowing parameter 	*
c*      g:	mapping function                *
c************************************************
      subroutine malla(n,p,c,h,g,h2)

      implicit none

      integer j,n,p
      real*4 h(n),pi,c
      real*8 h2(n),pi2,g(n)
      pi=acos(-1.)
      pi2=dacos(-1d0)


      if (p.eq.0) then
         do j=1,n
            h(j)=2.*(j-1)/(n-1)
            h(j)=h(j)-1.
            h2(j)=2.*float(j-1)/float(n-1)
            h2(j)=h2(j)-1d0
            g(j)=1d0
         enddo

         write(*,*) 'uniform mesh'
      endif

      if (p.eq.1) then          !!! malla tangente hiperbolica
         write(*,*) 'generating tanh mesh'
         do j=1,n
            h(j)=tanh(c*(2d0*(j-1)/(n-1)-1d0))/tanh(c)
            h2(j)=tanh(c*(2d0*float(j-1)/float(n-1)-1d0))/tanh(c)
            g(j)=(tanh(c)/c)*1d0/(1-(tanh(c*(2d0*(j-1)/(n-1)-1d0))**2))
         enddo
      endif

      if (p.eq.2) then          !!! malla tipo seno
         write(*,*) 'generating sen mesh'
         do j=1,n
            h(j) = sin(c*(2d0*(j-1)/(n-1)-1d0)*pi/2d0)/sin(c*pi/2d0)
            h2(j)= sin(c*(2d0*float(j-1)/float(n-1)-1d0)*pi2/2d0)/
     &           sin(c*pi2/2d0)
            g(j) = 2d0*sin(c*pi/2d0)/(cos(c*(2d0*(j-1)/
     &           (n-1)-1d0)*pi/2d0)*c*pi)

         enddo
      endif


      end




c***************************************************************************
c*	prederiv1:
c*
c*	Genera los coeficientes alpha,beta,a0,a1 y a2 para el esquema de
c*	diferencias finitas compactas con malla uniforme:
c*
c*	beta*(du(j+2)+du(j-2))+alpha*(du(j+1)+du(j-1))+du(j) =
c*	 = a1*(u(j+1)-u(j-1)) + a2(u(j+2)-u(j-2))
c*
c*		En la frontera del dominio:
c*		alpha*du(2) + du(1) =
c*		 = a1*u(1) + a2*u(2) + a3*u(3) + a4*u(4) + a5*u(5)
c*
c*		alpha*(du(3)+du(1)) + du(2) =
c*		 = a1*(u(3)-u(1))
c*
c***************************************************************************
      subroutine prederiv1m(my,ben,aln,aaa,alp,bep,a2n,a1n,a0,a1p,a2p)

      implicit none

      integer my
      real*8 ben(my),aln(my),aaa(my),alp(my),bep(my)
      real*8 a2n(my),a1n(my),a0(my),a1p(my),a2p(my)
      real*8 c(9,9),b(9),d(9,9)
      real*8 w1,w2,h,pi,omega1,omega2
      integer j

      pi = 4d0*atan(1d0)

      do j=1,my
      aaa(j)=1d0
      enddo

      h=2d0/(my-1)

c-------------------------centro del dominio-------------------------------
c	IMPONGO LOS COEFICIENTES A PARTIR DE CALCULOS EXTERNOS
      b(1) = 40941015625d-12
      b(2) = 4905703125d-10
      b(3) = (4d0-16d0*b(1)+2d0*b(2))/(6d0*h)
      b(4) = (22d0*b(1) + 4d0*b(2) - 1d0)/(12d0*h)

      do j=3,my-2
         ben(j)=b(1)
         aln(j)=b(2)
         alp(j)=b(2)
         bep(j)=b(1)
         a0(j)=0d0
         a1p(j)=b(3)
         a1n(j)=-a1p(j)
         a2p(j)=b(4)
         a2n(j)=-a2p(j)
      enddo

c---------------------fronteras del dominio--------------------------------

c     ESTANDAR:SOLO IMPONGO ECUACIONES DE ORDEN -> orden3

      j=1
      a2p(j) = 1d0/(2d0*h)
      a1p(j) = 2d0/h
      a0(j) = -15d0/(6d0*h)
      a2n(j) = 0d0
      a1n(j) = 0d0
      alp(j) = 2d0
      bep(j) = 0d0
      aln(j) = 0d0
      ben(j) = 0d0

      j=my
      a2n(j) = -1d0/(2d0*h)
      a1n(j) = -2d0/h
      a0(j)  = 15d0/(6d0*h)
      a2p(j) = 0d0
      a1p(j) = 0d0
      aln(j) = 2d0
      ben(j) = 0d0
      alp(j) = 0d0
      bep(j) = 0d0

c		PADE
      j=2
      a2p(j) = 0d0
      a1p(j) = 3d0/(4d0*h)
      a0(j)  = 0d0
      a1n(j) = -3d0/(4d0*h)
      a2n(j) = 0d0

      ben(j) = 0d0
      aln(j) = 25d-2
      alp(j) = 25d-2
      bep(j) = 0d0

      j=my-1
      a2n(j) = 0d0
      a1n(j) = -3d0/(4d0*h)
      a0(j)  = 0d0
      a1p(j) = 3d0/(4d0*h)
      a2p(j) = 0d0

      ben(j) = 0d0
      aln(j) = 25d-2
      alp(j) = 25d-2
      bep(j) = 0d0

      return
      end

c***************************************************************************
c*	prederiv1:
c*
c*	Genera los coeficientes alpha,beta,a0,a1 y a2 para el esquema de
c*	diferencias finitas compactas:
c*
c*	beta*(d2u(j+2)+d2u(j-2))+alpha*(d2u(j+1)+d2u(j-1))+d2u(j) =
c*	 = a0*u(j)+a1*(u(j+1)+u(j-1)) + a2(u(j+2)+u(j-2))
c*
c*
c*
c***************************************************************************
      subroutine prederiv2(my,ben,aln,aaa,alp,bep,a2n,a1n,a0,a1p,a2p,x)

      implicit none

      integer my,n
      real*8 x(my)
      real*8 ben(my),aln(my),aaa(my),alp(my),bep(my)
      real*8 a2n(my),a1n(my),a0(my),a1p(my),a2p(my)
      real*8 c(9,9),b(9)

      integer j

      do j=1,my
         aaa(j)=1d0
      enddo

c     Ahora calculo los coeficientes del centro del dominio

      do j=3,my-2
         c(1,1)=0d0
         c(1,2)=0d0
         c(1,3)=0d0
         c(1,4)=0d0
         c(1,5)=1d0
         c(1,6)=1d0
         c(1,7)=1d0
         c(1,8)=1d0
         c(1,9)=1d0
         b(1)=0d0

         c(2,1)=0d0
         c(2,2)=0d0
         c(2,3)=0d0
         c(2,4)=0d0
         c(2,5)=(x(j-2)-x(j))
         c(2,6)=(x(j-1)-x(j))
         c(2,7)=0d0
         c(2,8)=(x(j+1)-x(j))
         c(2,9)=(x(j+2)-x(j))
         b(2)=0d0

         c(3,1)=2d0
         c(3,2)=2d0
         c(3,3)=2d0
         c(3,4)=2d0
         c(3,5)=-(x(j-2)-x(j))**2
         c(3,6)=-(x(j-1)-x(j))**2
         c(3,7)=0d0
         c(3,8)=-(x(j+1)-x(j))**2
         c(3,9)=-(x(j+2)-x(j))**2
         b(3)=-2d0

         do n=3,8
            c(n+1,1)=n*(n-1)*(x(j-2)-x(j))**(n-2)
            c(n+1,2)=n*(n-1)*(x(j-1)-x(j))**(n-2)
            c(n+1,3)=n*(n-1)*(x(j+1)-x(j))**(n-2)
            c(n+1,4)=n*(n-1)*(x(j+2)-x(j))**(n-2)
            c(n+1,5)=-(x(j-2)-x(j))**n
            c(n+1,6)=-(x(j-1)-x(j))**n
            c(n+1,7)=0d0
            c(n+1,8)=-(x(j+1)-x(j))**n
            c(n+1,9)=-(x(j+2)-x(j))**n
            b(n+1)=0d0
         enddo

         call gaussj(c,9,9,b,1,1)


         ben(j)=b(1)
         aln(j)=b(2)
         alp(j)=b(3)
         bep(j)=b(4)
         a2n(j)=b(5)
         a1n(j)=b(6)
         a0(j)=b(7)
         a1p(j)=b(8)
         a2p(j)=b(9)

      enddo


c			fronteras del dominio

c*************************************************************************
c*************************************************************************
      j=1
      c(1,1)=0d0
      c(1,2)=1d0
      c(1,3)=1d0
      c(1,4)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j+1)-x(j))
      c(2,4)=(x(j+2)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=0d0
      c(3,3)=-(x(j+1)-x(j))**2
      c(3,4)=-(x(j+2)-x(j))**2
      b(3)=-2d0

      do n=3,3
         c(n+1,1)=n*(n-1)*(x(j+1)-x(j))**(n-2)
         c(n+1,2)=0d0
         c(n+1,3)=-(x(j+1)-x(j))**n
         c(n+1,4)=-(x(j+2)-x(j))**n
         b(n+1)=0d0
      enddo

      call gaussj(c,4,9,b,1,1)

      ben(j)=0d0
      aln(j)=0d0
      alp(j)=b(1)
      bep(j)=0d0
      a2n(j)=0d0
      a1n(j)=0d0
      a0(j)=b(2)
      a1p(j)=b(3)
      a2p(j)=b(4)

C************************************************************************
c************************************************************************
      j=2
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=1d0
      c(1,4)=1d0
      c(1,5)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j-1)-x(j))
      c(2,4)=0d0
      c(2,5)=(x(j+1)-x(j))
      b(2)=0d0
      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=-(x(j-1)-x(j))**2
      c(3,4)=0d0
      c(3,5)=-(x(j+1)-x(j))**2
      b(3)=-2d0

      do n=3,4
      c(n+1,1)=n*(n-1)*(x(j-1)-x(j))**(n-2)
      c(n+1,2)=n*(n-1)*(x(j+1)-x(j))**(n-2)
      c(n+1,3)=-(x(j-1)-x(j))**n
      c(n+1,4)=0d0
      c(n+1,5)=-(x(j+1)-x(j))**n
      b(n+1)=0d0
      enddo

      call gaussj(c,5,9,b,1,1)

      ben(j)=0d0
      aln(j)=b(1)
      alp(j)=b(2)
      bep(j)=0d0
      a2n(j)=0d0
      a1n(j)=b(3)
      a0(j)=b(4)
      a1p(j)=b(5)
      a2p(j)=0d0

c*************************************************************************
c*************************************************************************
      j=my-1
      c(1,1)=0d0
      c(1,2)=0d0
      c(1,3)=1d0
      c(1,4)=1d0
      c(1,5)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j-1)-x(j))
      c(2,4)=0d0
      c(2,5)=(x(j+1)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=2d0
      c(3,3)=-(x(j-1)-x(j))**2
      c(3,4)=0d0
      c(3,5)=-(x(j+1)-x(j))**2
      b(3)=-2d0

      do n=3,4
         c(n+1,1)=n*(n-1)*(x(j-1)-x(j))**(n-2)
         c(n+1,2)=n*(n-1)*(x(j+1)-x(j))**(n-2)
         c(n+1,3)=-(x(j-1)-x(j))**n
         c(n+1,4)=0d0
         c(n+1,5)=-(x(j+1)-x(j))**n
         b(n+1)=0d0
      enddo

      call gaussj(c,5,9,b,1,1)

      ben(j)=0d0
      aln(j)=b(1)
      alp(j)=b(2)
      bep(j)=0d0
      a2n(j)=0d0
      a1n(j)=b(3)
      a0(j)=b(4)
      a1p(j)=b(5)
      a2p(j)=0d0

c*************************************************************************
c*************************************************************************
      j=my
      c(1,1)=0d0
      c(1,2)=1d0
      c(1,3)=1d0
      c(1,4)=1d0
      b(1)=0d0

      c(2,1)=0d0
      c(2,2)=0d0
      c(2,3)=(x(j-1)-x(j))
      c(2,4)=(x(j-2)-x(j))
      b(2)=0d0

      c(3,1)=2d0
      c(3,2)=0d0
      c(3,3)=-(x(j-1)-x(j))**2
      c(3,4)=-(x(j-2)-x(j))**2
      b(3)=-2d0

      do n=3,3
         c(n+1,1)=n*(n-1)*(x(j-1)-x(j))**(n-2)
         c(n+1,2)=0d0
         c(n+1,3)=-(x(j-1)-x(j))**n
         c(n+1,4)=-(x(j-2)-x(j))**n
         b(n+1)=0d0
      enddo

      call gaussj(c,4,9,b,1,1)

      ben(j)=0d0
      aln(j)=b(1)
      alp(j)=0d0
      bep(j)=0d0
      a2n(j)=b(4)
      a1n(j)=b(3)
      a0(j)=b(2)
      a1p(j)=0d0
      a2p(j)=0d0

      return
      end


c***************************************************************************
c*	deryr2: complex*8
c* 	deryr:  real*4
c*
c* Find the first derivative
c*	u,du: vector dato y salida, (du puede sobrescribir u)
c*	m: number of vectors
c*	n: vector size
c*
c***************************************************************************
      subroutine deryr2(u,du,m,n,wk4)
      use matrices,only: dt12,prem1

      implicit none
      include "ctes3D"
      integer n,m,i,j,k
      real*4 u(2,n,m),du(2,n,m),wk4(*)

      real*8 fmap,y2
      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/

      real*8 wk1(my),wk2(my)
      real*8 uwk(my,2)

c				  !!!  calculo el termino independiente
      do k=1,m

         do j=1,my
            uwk(j,1) = u(1,j,k)
            uwk(j,2) = u(2,j,k)
         enddo

         wk1(1)=dt12(3,1)*uwk(1,1)+dt12(4,1)*uwk(2,1)+
     &        dt12(5,1)*uwk(3,1)
         wk1(2)=dt12(2,2)*uwk(1,1)+dt12(3,2)*uwk(2,1)+
     &        dt12(4,2)*uwk(3,1)+dt12(5,2)*uwk(4,1)
         do j=3,my-2
            wk1(j)= dt12(1,j)*uwk(1+j-3,1)
            do i=2,5
               wk1(j)=wk1(j) + dt12(i,j)*uwk(i+j-3,1)
            enddo
         enddo
         wk1(my-1)=dt12(1,my-1)*uwk(my-3,1)+dt12(2,my-1)*uwk(my-2,1)+
     &        dt12(3,my-1)*uwk(my-1,1)+dt12(4,my-1)*uwk(my,1)
         wk1(my)  =dt12(1,my)*uwk(my-2,1)+dt12(2,my)*uwk(my-1,1)+
     &        dt12(3,my)*uwk(my,1)

         wk2(1)=dt12(3,1)*uwk(1,2)+
     &        dt12(4,1)*uwk(2,2)+dt12(5,1)*uwk(3,2)
         wk2(2)=dt12(2,2)*uwk(1,2)+dt12(3,2)*uwk(2,2)+
     &        dt12(4,2)*uwk(3,2)+dt12(5,2)*uwk(4,2)
         do j=3,my-2
            wk2(j)=dt12(1,j)*uwk(1+j-3,2)
            do i=2,5
               wk2(j)=wk2(j) + dt12(i,j)*uwk(i+j-3,2)
            enddo
         enddo
         wk2(my-1)=dt12(1,my-1)*uwk(my-3,2)+dt12(2,my-1)*uwk(my-2,2)+
     .        dt12(3,my-1)*uwk(my-1,2)+dt12(4,my-1)*uwk(my,2)
         wk2(my)  =dt12(1,my)*uwk(my-2,2)+dt12(2,my)*uwk(my-1,2)+
     .        dt12(3,my)*uwk(my,2)

c     Real:
c-------PREPARO EL TERMINO INDEPENDIENTE [d2]*{u} !!!d2 casi-pentadiagonal!!!
c			backsubstitution, storage and MAPING

         call banbks5(prem1,my,wk1)
         call banbks5(prem1,my,wk2)

         do j=1,my
            du(1,j,k)=wk1(j)*fmap(j)
            du(2,j,k)=wk2(j)*fmap(j)
         enddo

      enddo                     !!!  resuelve n vectores complex*8


      end

c----------------------------------------------------------------------
      subroutine deryr(uwk,fwk,n)
      use matrices,only: dt12,prem1

      implicit none
      include "ctes3D"

      integer n,i,j,k
      real*8 fwk(my),uwk(my)

      real*8  fmap,y2
      real*4  Re,alp,bet,a0,y,hy
      common /fis/ Re,alp,bet,a0,y(my),hy(my),fmap(my),y2(my)
      save   /fis/



c-------PREPARO EL TERMINO INDEPENDIENTE [d2]*{u} !!!d2 casi-pentadiagonal!!!
      fwk(1) =dt12(3,1)*uwk(1) + dt12(4,1)*uwk(2) + dt12(5,1)*uwk(3)
      fwk(2) =dt12(2,2)*uwk(1)+dt12(3,2)*uwk(2)+
     .     dt12(4,2)*uwk(3)+dt12(5,2)*uwk(4)
      do j=3,my-2
         fwk(j)=dt12(1,j)*uwk(1+j-3)
         do i=2,5
            fwk(j) = fwk(j) + dt12(i,j)*uwk(i+j-3)
         enddo
      enddo
      fwk(my-1)=dt12(1,my-1)*uwk(my-3)+dt12(2,my-1)*uwk(my-2)+
     .     dt12(3,my-1)*uwk(my-1)+dt12(4,my-1)*uwk(my)
      fwk(my)=  dt12(1,my)*uwk(my-2)  +dt12(2,my)*uwk(my-1)+
     .     dt12(3,my)*uwk(my)

      call banbks5(prem1,my,fwk)


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!poner estos bucles bien !!!!!!!!!!!!!!!!!!
      do j=1,my
         fwk(j)=fwk(j)*fmap(j)
      enddo


      end

c**************************************************
c     *    SUBRUTINAS PARA TRABAJAR CON 5-DIAGONALES
c     *    SUBROUTINES TO WORK WITH 5-DIAGONALS
c     *    (Numerical Recipes in FORTRAN + Sergio )
c     *
c     *               TRANSPUESTAS !!! (mas rapidas!!)
c     *
c     *    bandec5, banbks5
c**************************************************

      SUBROUTINE banbks5(a,n,b)
      INTEGER n
      REAL*8 a(5,n),b(n)
      INTEGER i,k

      do k=1,n-2
         b(k+1) = b(k+1)-a(4,k)*b(k)
         b(k+2) = b(k+2)-a(5,k)*b(k)
      enddo
      b(n) = b(n)- a(4,n-1)*b(n-1)

!     back substitution

      b(n) = b(n)*a(1,n)
      b(n-1) = (b(n-1)-a(2,n-1)*b(n))*a(1,n-1)
      do i=n-2,1,-1
         b(i) = (b(i)-a(2,i)*b(1+i)-a(3,i)*b(2+i))*a(1,i)
      enddo

      return
      END

      SUBROUTINE bandec5(a,n)
      INTEGER n
      REAL*8 a(5,n)
      INTEGER j,k

      do j=1,3
         a(j,1)=a(j+2,1)
      enddo
      do j=1,4
         a(j,2)=a(j+1,2)
      enddo


      do k=1,n-2
         a(1,k)   = 1d0/a(1,k)

         a(4,k)   = a(1,k+1)*a(1,k)

         a(1,k+1) = a(2,k+1)-a(4,k)*a(2,k)
         a(2,k+1) = a(3,k+1)-a(4,k)*a(3,k)
         a(3,k+1) = a(4,k+1)
!
         a(5,k)   = a(1,k+2)*a(1,k)

         a(1,k+2) = a(2,k+2)-a(5,k)*a(2,k)
         a(2,k+2) = a(3,k+2)-a(5,k)*a(3,k)
         a(3,k+2) = a(4,k+2)
         a(4,k+2) = a(5,k+2)
      enddo


      a(1,n-1) = 1d0/a(1,n-1)

      a(4,n-1)=a(1,n)*a(1,n-1)

      a(1,n) = a(2,n)-a(4,n-1)*a(2,n-1)
      a(2,n) = a(3,n)-a(4,n-1)*a(3,n-1)
      a(3,n) = a(4,n)

      a(1,n)=1d0/a(1,n)


      END


c
C
c**************************************************
c*    SUBRUTINA PARA RESOVER SIST LINEALES
c*    SUBROUTINE TO RESOVE LINEAR SYSTEMS
c*    (Numerical Recipes in FORTRAN)
c**************************************************
      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*8 big,dum,pivinv
      do 11 j=1,n
         ipiv(j)=0
 11   continue
      do 22 i=1,n
         big=0d0
         do 13 j=1,n
            if(ipiv(j).ne.1)then
               do 12 k=1,n
                  if (ipiv(k).eq.0) then
                     if (abs(a(j,k)).ge.big)then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  else if (ipiv(k).gt.1) then
                     pause 'singular matrix in gaussj'
                  endif
 12            continue
            endif
 13      continue
         ipiv(icol)=ipiv(icol)+1
         if (irow.ne.icol) then
            do 14 l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
 14         continue
            do 15 l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
 15         continue
         endif
         indxr(i)=irow
         indxc(i)=icol
         if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
         pivinv=1./a(icol,icol)
         a(icol,icol)=1d0
         do 16 l=1,n
            a(icol,l)=a(icol,l)*pivinv
 16      continue
         do 17 l=1,m
            b(icol,l)=b(icol,l)*pivinv
 17      continue
         do 21 ll=1,n
            if(ll.ne.icol)then
               dum=a(ll,icol)
               a(ll,icol)=0d0
               do 18 l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
 18            continue
               do 19 l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
 19            continue
            endif
 21      continue
 22   continue
      do 24 l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            do 23 k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
 23         continue
         endif
 24   continue

      return
      END


      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do 11 j=ii,i-1
               sum=sum-a(i,j)*b(j)
 11         continue
         else if (sum.ne.0d0) then
            ii=i
         endif
         b(i)=sum
 12   continue
      do 14 i=n,1,-1
         sum=b(i)
         do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
 13      continue
         b(i)=sum/a(i,i)
 14   continue
      return
      END



      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1d-20)
      INTEGER i,imax,j,k
      REAL*4 aamax,dum,sum,vv(NMAX)
      d=1d0
      do 12 i=1,n
         aamax=0d0
         do 11 j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11      continue
         if (aamax.eq.0d0) pause 'singular matrix in ludcmp'
         vv(i)=1d0/aamax
 12   continue
      do 19 j=1,n
         do 14 i=1,j-1
            sum=a(i,j)
            do 13 k=1,i-1
               sum=sum-a(i,k)*a(k,j)
 13         continue
            a(i,j)=sum
 14      continue
         aamax=0d0
         do 16 i=j,n
            sum=a(i,j)
            do 15 k=1,j-1
               sum=sum-a(i,k)*a(k,j)
 15         continue
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
 16      continue
         if (j.ne.imax)then
            do 17 k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
 17         continue
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.0d0)a(j,j)=TINY
         if(j.ne.n)then
            dum=1./a(j,j)
            do 18 i=j+1,n
               a(i,j)=a(i,j)*dum
 18         continue
         endif
 19   continue
      return
      END

