c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine set_parameters()
*
      implicit NONE
*
      double precision pi,R,eps1,eps2,delta
      common /const/ R,eps1,eps2,delta
*
      pi=4d0*datan(1d0)
*
      R=1d0
      eps1=0.1d0
      eps2=0.075d0
      delta=pi*0.618d0
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
      CAVITY_SYM_CLASS=0
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
      integer n,nbe
      double precision pi,R,eps1,eps2,delta
      double precision ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe)
      double precision nx(1:nbe),ny(1:nbe)
      double precision rad,phi,dphi
      double precision cs1,sn1,cs2,sn2,dum
      double precision drdp,d2rdp2,dpds,d2pdpds
      double precision dxdp,dydp,d2xdp2,d2ydp2,X,Y
      common /const/ R,eps1,eps2,delta
*
      pi=4d0*datan(1d0)
*
      call set_parameters()
*     
      dphi=2d0*pi/dble(nbe)
*
      do n=1,nbe
         phi = dphi*(dble(n)-0.5d0)
*
         cs1=dcos(phi)
         cs2=dcos(2d0*phi+delta)
         sn1=dsin(phi)
         sn2=dsin(2d0*phi+delta)
*
C        Cavity Shape
*
         rad=R*(1d0 + eps1*cs1 + eps2*cs2)
         xl(n)=rad*cs1
         yl(n)=rad*sn1
*
C        Curvature
*
C        Terms dependent on the shape:
*
         drdp=-R*(eps1*sn1+2d0*eps2*sn2)
         d2rdp2=-R*(eps1*cs1+4d0*eps2*cs2)
*
C        Terms independent of the shape:
*
         dpds=1d0/dsqrt(rad*rad+drdp*drdp)
         d2pdpds=-drdp*(rad+d2rdp2)*dpds*dpds*dpds
*
         dxdp=drdp*cs1-rad*sn1
         dydp=drdp*sn1+rad*cs1
         d2xdp2=d2rdp2*cs1-2d0*drdp*sn1-rad*cs1
         d2ydp2=d2rdp2*sn1+2d0*drdp*cs1-rad*sn1
*
         X=dpds*(d2pdpds*dxdp+dpds*d2xdp2)
         Y=dpds*(d2pdpds*dydp+dpds*d2ydp2)
*     
         kappa(n)=dsqrt(X*X+Y*Y)
*
C        Normal Vector (nx,ny) & Boundary Element ds
*
         ds(n)=dphi/dpds
*
         nx(n)=dydp*dpds
         ny(n)=-dxdp*dpds
         
      enddo
*
      return
      end
*******************************************************************************
*
      subroutine indicator_func(x,y,flag)
*
      implicit NONE
*
      integer flag
      double precision R,eps1,eps2,delta,x,y,phi,rad,rad_c
      common /const/ R,eps1,eps2,delta
*
      call set_parameters()
*     
      phi=datan2(y,x)
      rad=dsqrt(x*x + y*y)
      rad_c=R*(1d0 + eps1*dcos(phi) + eps2*dcos(2d0*phi+delta))
*
      if(rad.LT.rad_c) then
         flag=1
      else
         flag=0
      endif
*******************************************************************************
      return
      end
