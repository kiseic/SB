c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine set_parameters()
*
      implicit NONE
*
      double precision R,eps
      common /const/ R,eps
*
      R=1d0
      eps=0.15d0
*
      R=R*(SIZE_PARAM)
*
      return
      end
*******************************************************************************
*
      subroutine test_cavity_symmetry(symmetry)
*
      implicit NONE
*
      integer CAVITY_SYM_CLASS,symmetry
*
      CAVITY_SYM_CLASS=2
*
C     CAVITY_SYM_CLASS=0: The cavity has no mirror symmetry.
C     CAVITY_SYM_CLASS=1: The cavity is symmetric with respect to the x-axis.
C     CAVITY_SYM_CLASS=2: The cavity is symmetric with respect to both x- and y-axis.
C     CAVITY_SYM_CLASS=4: The cavity is symmetric with respect to x=y and the y-axis.
*
      if(symmetry.ne.CAVITY_SYM_CLASS) then
         write(*,*) "Error in the parameter SYM."
         stop
      endif
*
      return
      end
*******************************************************************************
*
      subroutine def_bndry(nbe,ds,kappa,xl,yl,nx,ny)
*
      implicit NONE
*
      integer n,nbe,symmetry
      double precision pi,R,eps
      double precision ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe)
      double precision nx(1:nbe),ny(1:nbe)
      double precision cos1,sin1,cos2,sin2,u
      double precision phi,dphi,r0,r1,r2,X,Y,X1,Y1,X2,Y2,A,B
      double precision dx_ds,dy_ds,dum
      common /const/ R,eps
*
      pi=4d0*datan(1d0)
*
      call set_parameters()
*     
      dphi=0.5d0*pi/dble(nbe)
*
      do n=1,nbe
         phi = dphi*(dble(n)-0.5d0)
*
         cos1=dcos(phi)
         sin1=dsin(phi)
         cos2=dcos(2d0*phi)
         sin2=dsin(2d0*phi)
*
C        Cavity shape definitions:
         u=dsqrt( 1d0 + 2d0*eps*cos2 )
         r0=R*u ! r(phi)
         r1=-2d0*R*eps*sin2/u ! dr/dphi
         r2=-R*( 4d0*eps*cos2/u + 4d0*eps*eps*sin2*sin2/u/u/u ) ! d^2 r/dphi^2
*
C        Cavity boundary
         xl(n)=r0*cos1
         yl(n)=r0*sin1
*
         dum=dsqrt( r0*r0 + r1*r1 )
         ds(n)=dum*dphi ! Boundary element length
*
         X1=r1*cos1-r0*sin1
         Y1=r1*sin1+r0*cos1
         dx_ds=X1/dum
         dy_ds=Y1/dum
*
C        Normal vectors
         nx(n) =  dy_ds
         ny(n) = -dx_ds
*
         A=1d0/dum
         B=-r1*(r0+r2)/dum/dum/dum
*
         X2 = r2*cos1 - 2d0*r1*sin1 - r0*cos1
         Y2 = r2*sin1 + 2d0*r1*cos1 - r0*sin1
*
         X=A*( B*X1 + A*X2 )
         Y=A*( B*Y1 + A*Y2 )
*     
         kappa(n)=dsqrt(X*X+Y*Y) ! Curvature
*         
      enddo
*
      return
      end
*
*******************************************************************************
*
      subroutine indicator_func(x,y,flag)
*
      implicit NONE
*
      integer flag
      double precision R,eps,x,y,phi,rad,rad_c
      common /const/ R,eps
*
      call set_parameters()
*
      phi=datan(y/x)
      rad=dsqrt(x*x + y*y)
      rad_c=R*dsqrt( 1d0 + 2d0*eps*dcos(2d0*phi) )
*
      if(rad.LT.rad_c) then
         flag=1
      else
         flag=0
      endif
*
      return
      end
*
*******************************************************************************
