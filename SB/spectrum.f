      subroutine power(hstep,kt,psref,nco)

      integer i,hstep,k,kpip(0:hstep/4),nco
      real*8 pw(0:hstep/2-1),kt(0:hstep-1)
      real*8 aaa(0:2*hstep-1)
      complex*16 psref(0:hstep-1)
      complex*16 powerspec(0:hstep-1)
       
      character*8 flname2

      kpip(0)=0
      
      do i=0,hstep-1
         aaa(2*i)=dble(psref(i))
         aaa(2*i+1)=dimag(psref(i))
      enddo

      call cdft(2*hstep,1,aaa,kpip,pw)

      do i=0,2*hstep-1
         aaa(i)=aaa(i)/dble(hstep)
      enddo

      do i=0,hstep-1
         powerspec(i)=dcmplx(aaa(2*i),aaa(2*i+1))
      enddo


      write(flname2,'(''fort.'',I3)') 100+nco
      open(10,file=flname2)
      do k=hstep-hstep/2/60,hstep-1
      write(10,*)kt(k),dble(abs(powerspec(k)*powerspec(k))*1d0)
      enddo

      do k=0,hstep/2/60-1
c         write(83,*)kt(k),dble(powerspec(k)*dconjg(powerspec(k)))
      write(10,*)kt(k),dble(abs(powerspec(k)*powerspec(k))*1d0)
      enddo
      close(10)

      return
      end



