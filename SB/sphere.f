      subroutine sphere(echld,x,y,xmax0,ymax0,nx,ny,nco)

      implicit none
      integer nx,ny,i,j,k,l,nco,nk
      parameter ( nk = 6 )
      integer hx(nk+3),hy(nk+3)

      real*8  xmin,xmax,ymin,ymax,stmax,srmax,srmin,srcen
      real*8  xcenter,ycenter,xdlt,ydlt,pi
      real*8  x(0:nx),y(0:ny),sx,sy,dlt
      real*8  xx(0:nx),yy(0:ny)
      real*8  pax,pay,xmax0,ymax0
      real*8  u(nk+3),v(nk+3)
      complex*16  zr(nk+3)
      complex*16  vxr(nk+3),vyr(nk+3)
      complex*16  val,tax,tay

      complex*16 echld(0:nx-1,0:ny-1)
      complex*16 ez(0:nx+3,0:ny+3)
      complex*16 sph(0:nx+3,0:ny+3)

      pi=4d0*datan(1d0)
      xmax = xmax0
      ymax = ymax0
      xmin = -xmax
      ymin = -ymax
 
      xcenter = ( xmax + xmin )/2.d0
      ycenter = ( ymax + ymin )/2.d0
      xdlt = ( xmax - xmin )/DBLE(nx)
      ydlt = ( ymax - ymin )/DBLE(ny)

      stmax = 2.d0*pi
      srmax = (xmax/DBLE(nx))*DBLE(nx-3)
      srmin = srmax/DBLE(nx-1)
      srcen = (srmax + srmin )/2.d0

      do i=0,nx-1
         xx(i) = x(i)
      enddo
      do j=0,ny-1
         yy(j) = y(j)
      enddo

      do j=0,ny-1
         do i=0,nx-1
            ez(i,j) = echld(i,j)
         enddo
      enddo

c      do i=0,nx/2-1
c         xx(i)=x(i+nx/2)
c         xx(i+nx/2)=x(i)
c      enddo
c      do j=0,ny/2-1
c         yy(j)=y(j+ny/2)
c         yy(j+ny/2)=y(j)
c      enddo
c
c      do j=0,ny/2-1
c         do i=0,nx/2-1
c            ez(i,j)=ezx(i+nx/2,j+ny/2)+ezy(i+nx/2,j+ny/2)
c            ez(i+nx/2,j)=ezx(i,j+ny/2)+ezy(i,j+ny/2)
c            ez(i,j+ny/2)=ezx(i+nx/2,j)+ezy(i+nx/2,j)
c            ez(i+nx/2,j+ny/2)=ezx(i,j)+ezy(i,j)
c         enddo
c      enddo

      do j = 0,ny-1
         do i = 0,nx-1

            dlt = (srmax - srmin)/DBLE(nx)
            if ( i .le. nx/2-1 ) then
               sx = srcen + dlt*DBLE(i)
            else 
               sx = srmin + dlt*DBLE(i-nx/2)
            endif
            dlt = stmax/DBLE(ny)
            if ( j .le. ny/2-1 ) then
               sy = pi + dlt*DBLE(j)
            else 
               sy = dlt*DBLE(j-ny/2)
            endif
            pax = sx*dcos(sy)
            pay = sx*dsin(sy)

            if (pax.ge.xcenter) then
               hx(1)=int((pax-xcenter)/xdlt)
               hx(2)=hx(1)+1
               hx(3)=hx(2)+1
               hx(5)=hx(3)+1

               if (hx(1).eq.0) then
                  hx(4)=nx-1
                  hx(6)=hx(4)-1
               else
                  hx(4)=hx(1)-1
                  if (hx(4).eq.0) then
                     hx(6)=nx-1
                  else
                     hx(6)=hx(4)-1
                  endif
               endif

               if (pay.ge.ycenter) then
                  hy(1)=int((pay-ycenter)/ydlt)
                  hy(2)=hy(1)+1
                  hy(3)=hy(2)+1
                  hy(5)=hy(3)+1

                  if (hy(1).eq.0) then
                     hy(4)=ny-1
                     hy(6)=hy(4)-1
                  else
                     hy(4)=hy(1)-1
                     if (hy(4).eq.0) then
                        hy(6)=ny-1
                     else
                        hy(6)=hy(4)-1
                     endif
                  endif

               else
                  hy(1)=ny/2+int(dabs(pay-ymin)/ydlt)
                  if (hy(1).eq.ny-1) then
                     hy(2)=0
                     hy(3)=hy(2)+1
                     hy(4)=hy(1)-1
                     hy(5)=hy(3)+1
                     hy(6)=hy(4)-1
                  else
                     hy(2)=hy(1)+1

                     if (hy(2).eq.ny-1) then
                        hy(3)=0
                        hy(5)=hy(3)+1
                     else
                        hy(3)=hy(2)+1
                        if (hy(3).eq.ny-1) then
                           hy(5)=0
                        else
                           hy(5)=hy(3)+1
                        endif
                     endif

                     hy(4)=hy(1)-1
                     hy(6)=hy(4)-1
                  endif
               endif
            else
               hx(1)=nx/2+int(dabs(pax-xmin)/xdlt)
               if (hx(1).eq.nx-1) then
                  hx(2)=0
                  hx(3)=hx(2)+1
                  hx(4)=hx(1)-1
                  hx(5)=hx(3)+1
                  hx(6)=hx(4)-1
               elseif (hx(1).eq.nx) then
                  hx(1)=0
                  hx(2)=hx(1)+1
                  hx(3)=hx(2)+1
                  hx(4)=nx-1
                  hx(5)=hx(3)+1
                  hx(6)=hx(4)-1
               else
                  hx(2)=hx(1)+1

                  if (hx(2).eq.nx-1) then
                     hx(3)=0
                     hx(5)=hx(3)+1
                  else
                     hx(3)=hx(2)+1
                     if (hx(3).eq.nx-1) then
                        hx(5)=0
                     else
                        hx(5)=hx(3)+1
                     endif
                  endif

                  hx(4)=hx(1)-1
                  hx(6)=hx(4)-1
               endif

               if (pay.ge.ycenter) then
                  hy(1)=int((pay-ycenter)/ydlt)
                  hy(2)=hy(1)+1
                  hy(3)=hy(2)+1
                  hy(5)=hy(3)+1

                  if (hy(1).eq.0) then
                     hy(4)=ny-1
                     hy(6)=hy(4)-1
                  else
                     hy(4)=hy(1)-1
                     if (hy(4).eq.0) then
                        hy(6)=ny-1
                     else
                        hy(6)=hy(4)-1
                     endif
                  endif

               else
                  hy(1)=ny/2+int(dabs(pay-ymin)/ydlt)
                  if (hy(1).eq.ny-1) then
                     hy(2)=0
                     hy(3)=hy(2)+1
                     hy(4)=hy(1)-1
                     hy(5)=hy(3)+1
                     hy(6)=hy(4)-1
                  else
                     hy(2)=hy(1)+1

                     if (hy(2).eq.ny-1) then
                        hy(3)=0
                        hy(5)=hy(3)+1
                     else
                        hy(3)=hy(2)+1
                        if (hy(3).eq.ny-1) then
                           hy(5)=0
                        else
                           hy(5)=hy(3)+1
                        endif
                     endif

                     hy(4)=hy(1)-1
                     hy(6)=hy(4)-1
                  endif
               endif
            endif


c###################################################################
c##                                                               ##
c##   4 times 4 = 16 points are decided!                          ##
c##                                                               ## 
c##   using those values, the value of point x is sought          ##
c##                                                               ##
c###################################################################

            do k=1,nk
               u(k)=xx(hx(k))
               v(k)=yy(hy(k))
            enddo

            do k=1,nk
               do l=1,nk
                  zr(l)=ez( hx(l),hy(k) )
               enddo
               call lagr(nk,pax,u,zr,val)
               vxr(k)=val
            enddo

            do l=1,nk
               do k=1,nk
                  zr(l)=ez( hx(l),hy(k) )
               enddo
               call lagr(nk,pay,v,zr,val)
               vyr(l)=val
            enddo
            call lagr(nk,pay,v,vxr,val)
            tax=val
            call lagr(nk,pax,u,vyr,val)
            tay=val
            sph(i,j)=(tax+tay)/2d0
         enddo
      enddo
      
      call mwrite(sph,srmax,srmin,srcen,nx,ny,nco)

      return
      end

c###########################################################
c##                                                       ##
c##   Lagrange interpolation method                       ##
c##                                                       ##
c###########################################################
      subroutine lagr(nf,zz,z,vz,pn)
      integer i,j,nf
      real*8  zz,z(nf+3)
      complex*16 vz(nf+3),pn,pp

      pn=0d0
      do i=1,nf
         pp=1d0
         do j=1,nf
            if (i.eq.j) then
            else
               pp=pp*(zz-z(j))/(z(i)-z(j))
            endif
         enddo
         pn=pn+vz(i)*pp
      enddo

      return
      end
c###########################################################
c##                                                       ##
c##   Plot m-spectrum                                     ##
c##                                                       ##
c###########################################################

      subroutine mwrite(sph,rmax,rmin,rcen,nx,ny,nco)

      integer nx,ny,ip(0:ny+2),nco
      real*8  a(0:2*ny-1),w(0:ny/2-1)
      real*8  m,r,rmax,rmin,rcen,dr,etotal
      complex*16 sph(0:nx+3,0:ny+3)
      complex*16 em(0:nx+3,0:ny+3)
      character*8 flname1

      write(flname1,'(''mome.'',I3)') 100+nco
      open(10,file=flname1)

      do i = 0,nx-1

         do j = 0,ny-1
            a(2*j)   = dble( sph(i,j) )
            a(2*j+1) = dimag( sph(i,j) )
         enddo
         ip(0) = 0
         call CDFT(2*ny,1,a,ip,w)
         do j = 0,2*ny-1
            a(j) = a(j)/DBLE(ny)
         enddo
         do j = 0,ny-1
            em(i,j) = DCMPLX(a(2*j),a(2*j+1))
         enddo
      enddo
      dr = (rmax - rmin)/DBLE(nx)
      do j = 0,ny-1
         etotal = 0.d0
         do i = 0,nx-1
            if ( i .le. nx/2-1 ) then
               r = rcen + dr*DBLE(i)
            else 
               r = rmin + dr*DBLE(i-nx/2)
            endif
            etotal = etotal + r*DBLE(em(i,j)*DCONJG(em(i,j)))*dr
         enddo
         
         if ( j .le. ny/2-1 ) then
            m = -DBLE(j)
         else
            m = DBLE(ny - j)
         endif

         write(10,*)m,etotal

      enddo
      close(10)
      return
      end

