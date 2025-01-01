      subroutine nlop(echld,p,ww,t,nx,ny,kpa,gpp,gpl,winf,delta0,w0)

      integer i,j,nx,ny
      real*8 kpa,gpl,gpp,winf,t,dww
      real*8 ww(0:nx-1,0:ny-1)
      real*8 w0(0:nx-1,0:ny-1)
      real*8 delta0
      complex*16 phase 

      complex*16 echld(0:nx-1,0:ny-1),p(0:nx-1,0:ny-1)
      complex*16 dp,ep
      phase = dcmplx( dcos(delta0*t), -dsin(delta0*t) )

      do j=0,ny-1
         do i=0,nx-1
            dp=-gpp*p(i,j)+kpa*ww(i,j)*echld(i,j)
            ep=echld(i,j)*dconjg(p(i,j))
            ep=ep+dconjg(echld(i,j))*p(i,j)
            dww=-gpl*(ww(i,j)-w0(i,j))-2d0*kpa*dble(ep)
            p(i,j)=( p(i,j)+dp*t ) * phase
            ww(i,j)=ww(i,j)+dww*t
         enddo
      enddo

      return
      end
