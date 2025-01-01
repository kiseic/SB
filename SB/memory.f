
!*******************************************************************
! memory 
!*******************************************************************
      subroutine memory(echld,p,ww,time,nco,nx,ny)
      integer i,j,nx,ny,nco
      real*8 time
      complex*16 echld(0:nx-1,0:ny-1)
      complex*16 p(0:nx-1,0:ny-1)
      real*8 ww(0:nx-1,0:ny-1)
      
      character(13) memory13
      character(14) memory14
      character(16) memory16
      character(17) memory17
      
      write(memory16,'(''memory-time.data'')')
      open(100,file=memory16)
      write(memory17,'(''memory-echld.data'')')
      open(101,file=memory17)
      write(memory13,'(''memory-p.data'')')
      open(102,file=memory13)
      write(memory14,'(''memory-ww.data'')')
      open(103,file=memory14)
      
      write(100,*) time,nco
      
      do j=0,ny
         do i=0,nx
            write(101,*) echld(i,j)
            write(102,*) p(i,j)
            write(103,*) ww(i,j)
         end do
      end do
      
      close(100)
      close(101)
      close(102)
      close(103)

      return
      end
      
!*******************************************************************
! load
!*******************************************************************

      subroutine load(echld,p,ww,time,nco,nx,ny)
      integer i,j,nx,ny,nco
      real*8 time
      complex*16 echld(0:nx-1,0:ny-1)
      complex*16 p(0:nx-1,0:ny-1)
      real*8 ww(0:nx-1,0:ny-1)
      
      open(200,file='memory-time.data')
      open(201,file='memory-echld.data')
      open(202,file='memory-p.data')
      open(203,file='memory-ww.data')
      
      read(200,*) time,nco
  
      do j=0,ny
         do i=0,nx
            read(201,*) echld(i,j)
            read(202,*) p(i,j)
            read(203,*) ww(i,j)
         end do
      end do

      close(200)
      close(201)
      close(202)
      close(203)
      
      return
      end
