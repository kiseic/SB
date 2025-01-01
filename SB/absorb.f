c#####################################################
c    計算領域の端における吸収領域の設定
c#####################################################

      subroutine absorber(absorbs,x,y,r2,rho,cw,nx,ny)
      integer nx,ny,i,j
      real*8 absorbs(0:nx-1,0:ny-1),x(0:nx-1),y(0:ny-1)
      real*8 r2,rho,cw,rad,calcab
      
      do j=0,ny-1
         do i=0,nx-1
            rad=dsqrt(x(i)*x(i)+y(j)*y(j))
            if (rad.lt.r2) then
               calcab=rad-r2
               absorbs(i,j)=rho*exp(-(calcab*calcab)/(cw*cw))
            else
               absorbs(i,j)=rho
            endif
         enddo
      enddo

      return
      end
