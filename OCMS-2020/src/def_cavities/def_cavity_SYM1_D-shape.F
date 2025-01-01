c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine set_parameters()
*
      implicit NONE
*
      double precision R,d
      common /const/ R,d
*
      R=1d0
      d=0.5d0
*
      R=R*(SIZE_PARAM)
      d=d*(SIZE_PARAM)
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
      CAVITY_SYM_CLASS=1
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
      integer n,nbe,N1,N2,m
      double precision R,d,ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe)
      double precision nx(1:nbe),ny(1:nbe)
      double precision pi,Lc,Ll,theta,arg
      double precision phi,phi0,dphi,y0,dy
      common /const/ R,d
*
      call set_parameters()
*     
      pi=4d0*datan(1d0)
*
      theta=dasin(d/R)
      arg=0.5d0*pi+theta
      Lc=arg*R   ! Length for the circular part
      Ll=dsqrt(R*R-d*d) ! Length for the linear part
*
      N2=aint(dble(nbe)*(Lc/(Lc+Ll)))
      N1=nbe-N2
*
C===============================
C     Linear part
C===============================
*     
      y0=dsqrt(R*R-d*d)
      dy=Ll/dble(N1)
*
      do n=1,N1
         xl(n)=d
         yl(n)=dy*(dble(n-1)+0.5d0)
         nx(n)=1d0
         ny(n)=0d0
         ds(n)=dy
         kappa(n)=0d0
      enddo
*
C==============================
C     Circular part
C==============================
      phi0=dacos(d/R)
      dphi=arg/dble(N2)
*
      do n=1,N2
         m=n+N1
         phi=phi0+dphi*(dble(n-1)+0.5d0)
         xl(m)=R*cos(phi)
         yl(m)=R*sin(phi)
         nx(m)=cos(phi)
         ny(m)=sin(phi)
         ds(m)=R*dphi
         kappa(m)=1d0/R
      enddo
*
      return
      end
*******************************************************************************

      subroutine indicator_func(x,y,flag)
*
      integer flag
      double precision R,d,x,y,rad
      common /const/ R,d
*
      call set_parameters()
*
      rad=dsqrt(x*x+y*y)
*
      if((rad.LT.R).AND.(x.LT.d)) then
         flag=1 ! Inside the D-shape
      else
         flag=0 ! Outside the D-shape
      endif
*
      return
      end
*
