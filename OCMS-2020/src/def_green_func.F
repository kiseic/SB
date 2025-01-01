c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine gfn(nk,x1,y1,x2,y2,w)
*
      implicit NONE
*
      integer a,b,nu
      double precision x1,y1,x2,y2
      double precision dx,dy
      double complex imu,nk,z,v,w
*
      imu=dcmplx(0d0,1d0)
*
      dx=x1-x2
      dy=y1-y2
      z=nk*dsqrt(dx*dx+dy*dy)
      call hankel_1(0,z,v)
      w=-0.25d0*imu*v
*
      return
      end
*********************************************************
*
      subroutine gfn_redef(nk,xl,yl,x,y,ds,small,w)
*
      implicit NONE
*
      double precision pi,gamma,xl,yl,x,y,ds
      double precision dx,dy,dist,small
      double complex imu,w,nk
*
      imu=dcmplx(0d0,1d0)
      gamma=0.5772156649d0
      pi=4d0*datan(1d0)
*
      dx=xl-x
      dy=yl-y
      dist=dsqrt(dx*dx+dy*dy)
*
      if(dist.LT.small) then
         w=-1d0/2d0/pi*(1d0-cdlog(nk*ds/4d0)+imu*pi/2d0-gamma)
      else
         call gfn(nk,xl,yl,x,y,w)
      endif
*
      return
      end
*********************************************************
*
      subroutine dgfn(nk,x1,y1,x2,y2,nx,ny,w)
*
      implicit NONE
*
      double precision x1,x2,y1,y2,nx,ny
      double precision dx,dy,fact,dist
      double complex imu,nk,w,z,v
*
      imu=dcmplx(0d0,1d0)
*
      dx=x1-x2
      dy=y1-y2
      dist=dsqrt(dx*dx+dy*dy)
      fact=(nx*dx+ny*dy)/dist
      z=nk*dist
      call hankel_1(1,z,v)
      w=0.25d0*imu*nk*fact*v
*
      return
      end
*********************************************************
*
      subroutine dgfn_redef(nk,xl,yl,x,y,nx,ny,kappa,small,w)
*
      implicit NONE
*
      double precision xl,yl,x,y,nx,ny,kappa
      double precision pi,dx,dy,dist,small
      double complex w,nk
*
      pi=4d0*datan(1d0)
*
      dx=xl-x
      dy=yl-y
      dist=dsqrt(dx*dx+dy*dy)
*
      if(dist.LT.small) then
         w=kappa/4d0/pi
      else
         call dgfn(nk,xl,yl,x,y,nx,ny,w)
      endif
*
      return
      end
*********************************************************
