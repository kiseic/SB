c################################################################
c##                      時間tau(i)と周波数kt(i)の設定
c################################################################
c##  tau : time , kt : angular frequency                       ##
c##  tmin,tmax : time range                                    ##
c##  tau = (tmax-tmin)/hstep                                   ##
c################################################################

c---------------------------------------------------------------------------

      subroutine tcord(tau,kt,hstep,t)

      integer hstep,i
      real*8 t,tmax,tmin,center,dlt,pi
      real*8 tau(0:hstep-1),kt(0:hstep-1)

      pi=4d0*datan(1d0)

      tmax=dble(hstep)*t
      tmin=0d0

      center=(tmax+tmin)/2d0
      dlt=(tmax-tmin)/dble(hstep)
      do i=0,hstep/2-1
         tau(i)=center+dlt*dble(i)
         tau(i+hstep/2)=tmin+dlt*dble(i)
      enddo

      dlt=2d0*pi/(tmax-tmin)
      do i=0,hstep/2
         kt(i)=dlt*dble(i)
      enddo
      do i=1,hstep/2
         kt(hstep-i)=-kt(i)
      enddo

      return
      end
