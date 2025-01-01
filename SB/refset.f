c####################################################
c 屈折率分布の設定
c 規格化屈折率分布　refs(i,j)　を設定する。
c refs(i,j) = 1 : 共振器内　（共振器内では１にすること）
c           = neff = n_out^2/n_in^2　
c （外側では、外側屈折率n_outを（規格化前の）共振器屈折率n_inでわり、
c  ２乗したものを代入すること)
c  ただし、共振器境界はなめらかに変化するように、xdeで制御。
c     xdeの効果は、refs(i,j)を描画させて見てチェックすること
c
c####################################################

      subroutine refset(refs,x,y,delta,dscrad,xde,neff,conreg,nx,ny)
      integer i,j,nx,ny
      integer conreg(0:nx-1,0:ny-1)
      real*8 refs(0:nx-1,0:ny-1),x(0:nx-1),y(0:ny-1)
      real*8 delta,dscrad,xde,neff
      real*8 ax,ay,radius,radius2,rad,rad2
      do j=0,ny-1
         do i=0,nx-1
            ax=dabs(x(i))
            ay=dabs(y(j))
*     
            if(  (ay-ax/delta).lt.
     +           (-dscrad/dsqrt(2d0)*(1d0/delta-1d0))) then
               rad2=(ax-dscrad/dsqrt(2d0)*(1d0-delta) )**2+ay*ay
               radius2=(dscrad**2)/2d0*(1d0+delta*delta)
            else 
               rad2=ax**2+(ay+dscrad/dsqrt(2d0)*(1d0/delta-1d0))**2
               radius2=(dscrad**2)/2d0*(1d0+1d0/(delta*delta))
            endif
            if(rad2.lt.radius2) then 
               conreg(i,j)=1
            else
               conreg(i,j)=0
            endif
            rad=dsqrt(rad2)
            radius=dsqrt(radius2)
            if( rad. le. radius -xde ) then 
               refs(i,j) = 1d0 
            elseif( dabs(rad-radius).le.xde ) then 
               refs(i,j) = (1d0+neff)/2d0
     +              -((1d0-neff)/2d0) * dtanh( (rad-radius)*5d0 )
               !write(23,*) x(i),y(j)
            else
               refs(i,j)=neff
            endif
*     
         enddo
      enddo
      return
      end
