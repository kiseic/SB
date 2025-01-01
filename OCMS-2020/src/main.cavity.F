c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      program cavity_data
*
      implicit NONE
*
      integer nbe
*
      parameter (nbe=NBE)
*
      integer n,kmax,kx,ky,flag
      double precision ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe)
      double precision nx(1:nbe),ny(1:nbe)
      double precision xmin,xmax,ymin,ymax
      double precision max,x,y,dx,dy,xx,yy
      double precision dist,ddx,ddy
*
      call def_bndry(nbe,ds,kappa,xl,yl,nx,ny)
*
      open(100,FILE='data.cavity_boundary')
      open(200,FILE='data.cavity_domain')
*
      write(100,500) '# SYM=',SYM
*
 500  format(A,I1)
*
      do n=1,nbe-1
         write(100,*) n,xl(n),yl(n),nx(n),ny(n),ds(n),kappa(n)
*
C     Insert a newline if there is a discontinuity in the boundary data
         ddx=xl(n+1)-xl(n)
         ddy=yl(n+1)-yl(n)
         dist=dsqrt( ddx*ddx + ddy*ddy )
         if(dist.GT.(2d0*ds(n))) then
            write(100,*) ""
         endif
      enddo
*
      n=nbe
      write(100,*) n,xl(n),yl(n),nx(n),ny(n),ds(n),kappa(n)
*
      xmin=1d9
      xmax=-1d9
      ymin=1d9
      ymax=-1d9
*
      do n=1,nbe
         if(xmin.gt.xl(n)) then
            xmin=xl(n)
         endif
*
         if(xmax.lt.xl(n)) then
            xmax=xl(n)
         endif
*
         if(ymin.gt.yl(n)) then
            ymin=yl(n)
         endif
*
         if(ymax.lt.yl(n)) then
            ymax=yl(n)
         endif
      enddo
*
#if (SYM == 4 || SYM == 2)
      xmin=-xmax
      ymin=-ymax
#elif SYM == 1
      ymin=-ymax
#elif SYM ==0
C     Do nothing
#endif
*
      if(xmax.GT.(-xmin)) then
         xx=xmax
      else
         xx=-xmin
      endif
*
      if(ymax.GT.(-ymin)) then
         yy=ymax
      else
         yy=-ymin
      endif
*
      if(xx.GT.yy) then
         max=xx
      else
         max=yy
      endif
*
      max=max*1.2d0
*
      kmax=400
      dx=2d0*max/dble(kmax)
      dy=2d0*max/dble(kmax)
*
      do kx=-kmax,kmax
         do ky=-kmax,kmax 
            x=-max+dx*dble(kx)
            y=-max+dy*dble(ky)
*
            call indicator_func(x,y,flag)
*
            if(flag.EQ.1) then
               write(200,*) x,y
            endif
*         
         enddo
      enddo
*
      close(100)
      close(200)
*
C      stop
      end
**************************************************************************
