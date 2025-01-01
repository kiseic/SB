c################################################################
c##                      xy座標と波数kx kyの設定                ##
c################################################################
      subroutine mkcord2dim(x,y,kx,ky,xmax,ymax,nx,ny)

      integer nx,ny,i
      real*8 xmin,xmax,ymin,ymax
      real*8 center,dlt,pi
      real*8 x(0:nx-1),y(0:ny-1)
      real*8 kx(0:nx-1),ky(0:ny-1)

      pi=4d0*datan(1d0)

      xmin=-xmax
      ymin=-ymax

      center=(xmax+xmin)/2d0
      dlt=(xmax-xmin)/dble(nx)
      do i=0,nx/2-1
         x(i)=center+dlt*dble(i)
         x(i+nx/2)=xmin+dlt*dble(i)
      enddo

      center=(ymax+ymin)/2d0
      dlt=(ymax-ymin)/dble(ny)
      do i=0,ny/2-1
         y(i)=center+dlt*dble(i)
         y(i+ny/2)=ymin+dlt*dble(i)
      enddo

      dlt=2d0*pi/(xmax-xmin)
      do i=0,nx/2
         kx(i)=dlt*dble(i)
      enddo
      do i=1,nx/2
         kx(nx-i)=-kx(i)
      enddo      

      dlt=2d0*pi/(ymax-ymin)
      do i=0,ny/2
         ky(i)=dlt*dble(i)
      enddo
      do i=1,ny/2
         ky(ny-i)=-ky(i)
      enddo

      return
      end

