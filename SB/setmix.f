      subroutine setmix(ech,x,y,dscrad)

      integer m,mmax,mmin
      real*8  x,y,dscrad,rad,theta,pi
      real*8  amp,r0,rw,prad
      complex*16  ech,imu
      
      imu = dcmplx(0.00d0,1.00d0)
 
      pi = 4d0*datan(1d0)
      amp = 1.d-3 
      r0  = dscrad*0.7d0
      rw = 1d0
      mmax = 15
      mmin = 1

      rad = dsqrt(x*x+y*y)
      if (x.gt.0d0.and.y.ge.0d0) then
         theta=datan(dabs(y/x))
      elseif (x.lt.0d0.and.y.gt.0d0) then
         theta=pi-datan(dabs(y/x))
      elseif (x.lt.0d0.and.y.le.0d0) then
         theta=pi+datan(dabs(y/x))
      elseif (x.gt.0d0.and.y.lt.0d0) then
         theta=2d0*pi-datan(dabs(y/x))
      elseif (x.eq.0d0.and.y.gt.0d0) then
         theta=pi/2d0
      elseif (x.eq.0d0.and.y.lt.0d0) then
         theta=3d0*pi/2d0
      else
         theta=0d0
      endif

      ech=0.d0
      do m=mmin,mmax
      prad  = (rad - r0)/rw
      if ( dabs( prad ) .lt. 7d0 ) then
         prad  = prad*prad
         ech=ech+amp*dexp(-prad)*exp( imu*dble(m)*theta )
      else
         ez=0.00d0
      endif
      enddo
      
      return
      end
