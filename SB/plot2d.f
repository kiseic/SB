      subroutine plot2d(echld,nco,na,nx,ny)
      integer i,j,nx,ny,nco,na
      complex*16 echld(0:nx-1,0:ny-1)
      character*8 flname1

      write(flname1,'(''wave.'',I3)') 100+nco
      open(10,file=flname1)
      write(10,*) "$ DATA=CONTOUR RETURN"
      write(10,*) "% nx=",100*na+1 
      write(10,*) "% ny=",100*na+1
      write(10,*) "% xmin=",-24.5*na 
      write(10,*) "% xmax=",24.5*na
      write(10,*) "% ymin=",-24.5*na 
      write(10,*) "% ymax=",24.5*na
      write(10,*) "% contfill"
       
      do j=ny-100*na,ny-1,2
         do i=nx-100*na,nx-1,2
            write(10,777)dble(abs(echld(i,j)*echld(i,j)))
         enddo
         do i=0,100*na,2
            write(10,777)dble(abs(echld(i,j)*echld(i,j)))
         enddo
      enddo
      
      do j=0,100*na,2
         do i=nx-100*na,nx-1,2
            write(10,777)dble(abs(echld(i,j)*echld(i,j)))
         enddo
         do i=0,100*na,2
            write(10,777)dble(abs(echld(i,j)*echld(i,j)))
         enddo
      enddo
 777  format(2e25.16e3)
      close(10)


      return
      end
