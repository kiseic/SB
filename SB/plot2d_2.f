      subroutine plot2d(echld,nco,na,nx,ny,ww,echld_t,tw)
      integer i,j,nx,ny,nco,na
      complex*16 echld(0:nx-1,0:ny-1)
      real*8 echld_t(0:nx-1,0:ny-1)
      real*8 ww(0:nx-1,0:ny-1)
      real*8 tw
      character*8 flname1
      character*6 flname2
      character*8 flname3

      write(flname1,'(''wave.'',I3)') 100+nco
      open(10,file=flname1)
      write(10,*) "$ DATA=CONTOUR RETURN"
      write(10,*) "% nx=",100*na+1 
      write(10,*) "% ny=",100*na+1
      write(10,*) "% xmin=",-24.5*na 
      write(10,*) "% xmax=",24.5*na
      write(10,*) "% ymin=",-24.5*na 
      write(10,*) "% ymax=",24.5*na
      write(10,*) '% toplabel= ""'
      write(10,*) '% subtitle= ""'
      write(10,*) "% contfill"

      write(flname2,'(''ww.'',I3)') 100+nco
      open(13,file=flname2)
      write(13,*) "$ DATA=CONTOUR RETURN"
      write(13,*) "% nx=",100*na+1 
      write(13,*) "% ny=",100*na+1
      write(13,*) "% xmin=",-24.5*na 
      write(13,*) "% xmax=",24.5*na
      write(13,*) "% ymin=",-24.5*na 
      write(13,*) "% ymax=",24.5*na
      write(13,*) '% toplabel= ""'
      write(13,*) '% subtitle= ""'
      write(13,*) "% contfill"

      write(flname3,'(''smth.'',I3)') 100+nco
      open(15,file=flname3)
      write(15,*) "$ DATA=CONTOUR RETURN"
      write(15,*) "% nx=",100*na+1 
      write(15,*) "% ny=",100*na+1
      write(15,*) "% xmin=",-24.5*na 
      write(15,*) "% xmax=",24.5*na
      write(15,*) "% ymin=",-24.5*na 
      write(15,*) "% ymax=",24.5*na
      write(15,*) '% toplabel= ""'
      write(15,*) '% subtitle= ""'
      write(15,*) "% contfill"

      if (nco.eq.0) then
         tw=1d0
         do j=ny-100*na,ny-1,2
            do i=nx-100*na,nx-1,2
               write(10,777)dble(abs(echld(i,j)*echld(i,j)))
               write(13,777)dble(ww(i,j))
               write(15,777)dble(echld_t(i,j)/tw)
            enddo
            do i=0,100*na,2
               write(10,777)dble(abs(echld(i,j)*echld(i,j)))
               write(13,777)dble(ww(i,j))
               write(15,777)dble(echld_t(i,j)/tw)            
            enddo
         enddo
      
         do j=0,100*na,2
            do i=nx-100*na,nx-1,2
               write(10,777)dble(abs(echld(i,j)*echld(i,j)))
               write(13,777)dble(ww(i,j))
               write(15,777)dble(echld_t(i,j)/tw)
            enddo
            do i=0,100*na,2
               write(10,777)dble(abs(echld(i,j)*echld(i,j)))
               write(13,777)dble(ww(i,j))
               write(15,777)dble(echld_t(i,j)/tw)
            enddo
         enddo
      else
         do j=ny-100*na,ny-1,2
            do i=nx-100*na,nx-1,2
               write(10,777)dble(abs(echld(i,j)*echld(i,j)))
               write(13,777)dble(ww(i,j))
               write(15,777)dble(echld_t(i,j)/tw)
            enddo
            do i=0,100*na,2
               write(10,777)dble(abs(echld(i,j)*echld(i,j)))
               write(13,777)dble(ww(i,j))
               write(15,777)dble(echld_t(i,j)/tw)            
            enddo
         enddo
      
         do j=0,100*na,2
            do i=nx-100*na,nx-1,2
               write(10,777)dble(abs(echld(i,j)*echld(i,j)))
               write(13,777)dble(ww(i,j))
               write(15,777)dble(echld_t(i,j)/tw)
            enddo
            do i=0,100*na,2
               write(10,777)dble(abs(echld(i,j)*echld(i,j)))
               write(13,777)dble(ww(i,j))
               write(15,777)dble(echld_t(i,j)/tw)
            enddo
         enddo
      
      endif
      
 777  format(2e25.16e3)
      close(10)
      close(13)
      close(15)

      return
      end
