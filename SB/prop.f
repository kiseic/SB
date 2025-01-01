      subroutine prop(echld,aa,bb,nx,ny)

      implicit none
      integer i,j,nx,ny,n1max
      complex*16 echld(0:nx-1,0:ny-1)
      complex*16 aa(0:nx-1,0:ny-1),bb(0:nx-1,0:ny-1)

      !---- FFT ---
      integer ip(0:nx/4+2)
      real*8  pw(0:nx/2),tt(0:8*ny-1)
      real*8  ae(0:2*nx-1,0:ny-1)
! aa = exp(imu*t*refs(i,j)/4d0-t*((back*conreg(i,j)+absorbs(i,j))/2d0))
! bb = exp(-(kx(i)*kx(i)+ky(j)*ky(j))*imu*t/2d0)

c----------------------------------------------------
c  Potential : Operation A
c----------------------------------------------------
      do j=0,ny-1
         do i=0,nx-1
            echld(i,j)  =  aa(i,j)*echld(i,j)
            ae(2*i,j)   =   dble( echld(i,j) )
            ae(2*i+1,j) =  dimag( echld(i,j) )
         enddo
      enddo

c----------------------------------------------------
c  Wavenumbers : Operation B
c----------------------------------------------------
      n1max = 2*nx
      ip(0) = 0
      
      call cdft2d(n1max,2*nx,ny,1,ae,tt,ip,pw)

      do j=0,ny-1
         do i=0,nx-1 
            echld(i,j)   = dcmplx( ae(2*i,j),ae(2*i+1,j) )
            echld(i,j)   =  bb(i,j)*echld(i,j)
            ae(2*i,j)    =  dble( echld(i,j) ) 
            ae(2*i+1,j)  = dimag( echld(i,j) ) 
         enddo
      enddo
!
!      ip(0) = 0
      
      call cdft2d(n1max,2*nx,ny,-1,ae,tt,ip,pw)
      
      do j=0,ny-1
         do i=0,nx-1
            ae(2*i,j)   = ae(2*i,j)/dble(nx*ny)
            ae(2*i+1,j) = ae(2*i+1,j)/dble(nx*ny)
            echld(i,j)  = dcmplx( ae(2*i,j),ae(2*i+1,j) )
c----------------------------------------------------
c  Potential : Operation A
c----------------------------------------------------
            echld(i,j)   = aa(i,j)*echld(i,j)
         enddo
      enddo

      return
      end
