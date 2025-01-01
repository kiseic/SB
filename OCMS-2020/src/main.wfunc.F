c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      program wave_function
*     
      implicit NONE
*
C Boundary zone width for the computation of the Hankel function
#define BZW_HANKEL (5d-3)
C Boundary zone width for the estimation of the maximum intensity
#define BZW_MAXINT (2d0*(BZW_HANKEL))
*
#define PLOT_NFP
#define PLOT_FFP
#define PLOT_WF
*
      integer nbe
      parameter (nbe=NBE)
*
      integer n,i,j,ix,iy,sign,nn,flag,imax,jx,jmax,nbe2
      integer a,b,ixmax,iymax,flag0
      integer info,lda,ldvl,ldvr,lwork
      double precision pi,nin,nout,rwork(nbe*4)
      double precision r,ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe),nx(1:nbe),ny(1:nbe)
      double precision kx,ky,dx,dy
      double precision err,minerr,x,y
      double precision lambda_in,ds_max,bndry
      double precision xmin,xmax,ymin,ymax,dist,dist_min
      double precision theta,dtheta,ex,ey,dum1,dum2,intensity,vmax
      double complex imu,k,nk_in,nk_out,alpha,z,wf,sa
      double complex P1(1:nbe,1:nbe),Q1(1:nbe,1:nbe)
      double complex P2(1:nbe,1:nbe),Q2(1:nbe,1:nbe)
      double complex T(1:nbe*2,1:nbe*2)
      double complex W(nbe*2),VL(nbe*2,nbe*2)
      double complex VR(nbe*2,nbe*2),work(nbe*4)
      double complex phi(1:nbe),dph_in(1:nbe),dph_out(1:nbe)
      character jobvl,jobvr
*
      pi=4d0*datan(1d0)
      imu=dcmplx(0d0,1d0)
*
C      OPEN(100,FILE='dat.W-dist')
*
#if defined(PLOT_WF)
      OPEN(200,FILE='dat.wfunc')
#endif
*
#if defined(PLOT_NFP)
      OPEN(300,FILE='dat.nearfield')
#endif
*
#if defined(PLOT_FFP)
      OPEN(400,FILE='dat.farfield')
C      OPEN(450,FILE='dat.ftheta')
#endif
*
#if (SYM == 4 || SYM == 2)
      read(5,*) nin,nout,a,b,kx,ky,xmin,xmax,ymin,ymax,ixmax,iymax
#elif SYM == 1
      read(5,*) nin,nout,b,kx,ky,xmin,xmax,ymin,ymax,ixmax,iymax
#elif SYM == 0
      read(5,*) nin,nout,kx,ky,xmin,xmax,ymin,ymax,ixmax,iymax
#endif
      call def_bndry(nbe,ds,kappa,xl,yl,nx,ny)
*
      ds_max=1d-100
      do n=1,nbe
         if(ds_max.LT.ds(n)) then
            ds_max=ds(n)
         endif
      enddo
*
      k = kx + imu*ky
*
#include "def_nk.F"
*
      lambda_in=2d0*pi/dble(nk_in)
*
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     B O U N D A R Y  I N T E G R A L  E Q U A T I O N
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
      sign=1
#if (SYM == 4 || SYM == 2)
      call def_matrix(nbe,nk_in,a,b,xl,yl,nx,ny,sign,ds,kappa,P1,Q1)
#elif SYM == 1
      call def_matrix(nbe,nk_in,b,xl,yl,nx,ny,sign,ds,kappa,P1,Q1)
#elif SYM ==0
      call def_matrix(nbe,nk_in,xl,yl,nx,ny,sign,ds,kappa,P1,Q1)
#endif
*     
      sign=-1
#if (SYM == 4 || SYM == 2)
      call def_matrix(nbe,nk_out,a,b,xl,yl,nx,ny,sign,ds,kappa,P2,Q2)
#elif SYM == 1
      call def_matrix(nbe,nk_out,b,xl,yl,nx,ny,sign,ds,kappa,P2,Q2)
#elif SYM == 0
      call def_matrix(nbe,nk_out,xl,yl,nx,ny,sign,ds,kappa,P2,Q2)
#endif
*     
      do i=1,nbe
         do j=1,nbe
#if defined(TM)
            T(i,j)=P1(i,j)
            T(i+nbe,j)=P2(i,j)
#elif defined(TE)
            T(i,j)=nin*nin*P1(i,j)
            T(i+nbe,j)=nout*nout*P2(i,j)
#else
            write(*,*) 'Error: polarization is not defined.'
            stop
#endif
            T(i,j+nbe)=Q1(i,j)
            T(i+nbe,j+nbe)=Q2(i,j)
         enddo
      enddo
*
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CALC: W.F. & NORMAL DERIVATIVE AT CAVITY BOUNDARY
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
      alpha=1d0+imu   ! alpha is an arbitrary (complex) number
*
      nbe2=nbe*2
*
      do n=1,nbe2
         T(n,n)=T(n,n)+alpha
      enddo
*    
      jobvl='N'
      jobvr='V'
      lda=nbe2
      ldvl=nbe2
      ldvr=nbe2
      lwork=2*nbe2
      call zgeev(jobvl,jobvr,nbe2,T,lda,W,VL,ldvl,VR,ldvr,
     &     work,lwork,rwork,info)
*
      minerr=1d10
      do n=1,nbe2
         dx=dble(W(n))-dble(alpha)
         dy=dimag(W(n))-dimag(alpha)
         err=dsqrt(dx*dx + dy*dy)
*
         if(err.LT.minerr) then
            nn=n
            minerr=err
         endif
      enddo
*
C      write(100,*) '# alpha = ',alpha
C      write(100,*) '# W(nn) = ',W(nn)
*
C      do n=1,nbe2
C         write(100,*) dble(W(n)),dimag(W(n))
C      enddo
C      close(100)
*
      do n=1,nbe
         phi(n)=VR(n+nbe,nn)   ! wavefunction at the boundary
#if defined(TM)
         dph_in(n)=VR(n,nn)     ! normal derivative (inside)
         dph_out(n)=VR(n,nn)    ! normal derivative (outside)
#elif defined(TE)
         dph_in(n)=nin*nin*VR(n,nn)    ! normal derivative (inside)
         dph_out(n)=nout*nout*VR(n,nn) ! normal derivative (outside)
#endif         
      enddo
*
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                  P L O T  F A R  F I E L D
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#if defined(PLOT_FFP)
*
      write(400,*) '#    [Angle (deg)]                [Intensity]'
*
      imax=aint(dble(nk_in)*dble(SIZE_PARAM))*20
*
      dtheta=2d0*pi/dble(imax)
*
      do i=0,imax-1
         theta=(dble(i)+0.5)*dtheta
         ex=dcos(theta)
         ey=dsin(theta)
*
         sa=dcmplx(0d0,0d0)
*
#if SYM == 4
         do n=1,nbe
            ! Region 1
            dum1 = ex*xl(n) + ey*yl(n)
            dum2 = ex*nx(n) + ey*ny(n)
            sa=sa+ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
*            
            ! Region 2
            dum1 = ex*yl(n) + ey*xl(n)
            dum2 = ex*ny(n) + ey*nx(n)
            sa=sa+dble(b)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
*     
            ! Region 3
            dum1 = ex*(-yl(n)) + ey*xl(n)
            dum2 = ex*(-ny(n)) + ey*nx(n)
            sa=sa+dble(a*b)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
*            
            ! Region 4
            dum1 = ex*(-xl(n)) + ey*yl(n)
            dum2 = ex*(-nx(n)) + ey*ny(n)
            sa=sa+dble(a)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
*            
            ! Region 5
            dum1 = ex*(-xl(n)) + ey*(-yl(n))
            dum2 = ex*(-nx(n)) + ey*(-ny(n))
            sa=sa+ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
*
            ! Region 6
            dum1 = ex*(-yl(n)) + ey*(-xl(n))
            dum2 = ex*(-ny(n)) + ey*(-nx(n))
            sa=sa+dble(b)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
*
            ! Region 7
            dum1 = ex*yl(n) + ey*(-xl(n))
            dum2 = ex*ny(n) + ey*(-nx(n))
            sa=sa+dble(a*b)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
*
            ! Region 8
            dum1 = ex*xl(n) + ey*(-yl(n))
            dum2 = ex*nx(n) + ey*(-ny(n))
            sa=sa+dble(a)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
*
         enddo
*
#elif SYM == 2
         do n=1,nbe
            dum1 = ex*xl(n) + ey*yl(n)
            dum2 = ex*nx(n) + ey*ny(n)
            sa=sa+ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))

            dum1 = ex*(-xl(n)) + ey*yl(n)
            dum2 = ex*(-nx(n)) + ey*ny(n)
            sa=sa+dble(a)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))

            dum1 = ex*(-xl(n)) + ey*(-yl(n))
            dum2 = ex*(-nx(n)) + ey*(-ny(n))
            sa=sa+dble(a*b)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))

            dum1 = ex*xl(n) + ey*(-yl(n))
            dum2 = ex*nx(n) + ey*(-ny(n))
            sa=sa+dble(b)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
         enddo
*
#elif SYM == 1
         do n=1,nbe
            dum1 = ex*xl(n) + ey*yl(n)
            dum2 = ex*nx(n) + ey*ny(n)
            sa=sa+ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))

            dum1 = ex*xl(n) + ey*(-yl(n))
            dum2 = ex*nx(n) + ey*(-ny(n))
            sa=sa+dble(b)*ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
         enddo
*
#elif SYM == 0
         do n=1,nbe
            dum1 = ex*xl(n) + ey*yl(n)
            dum2 = ex*nx(n) + ey*ny(n)
            sa=sa+ds(n)*exp(-imu*nk_out*dum1)*
     &           (imu*nk_out*dum2*phi(n)+dph_out(n))
         enddo
*
#endif
*
         sa=sa*(1d0+imu)/4d0/cdsqrt(pi*nk_out)
         intensity=cdabs(sa)*cdabs(sa)
         write(400,*) theta*180d0/pi,intensity
C         write(450,*) theta*180d0/pi,dble(sa),dimag(sa)
*
      enddo
*
      CLOSE(400)
C      CLOSE(450)
*
#endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                  P L O T  N E A R  F I E L D
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#if defined(PLOT_NFP)
*
      write(300,*) '#       [ s ]                     [ |phi|^2 ] 
     &        [ |dph|^2 ]'
*
      bndry=0.5d0*ds(1)
*
      do n=1,nbe
         if(n.EQ.1) then
            bndry=0.5d0*ds(1)
         else
            bndry=bndry+0.5d0*(ds(n-1)+ds(n))
         endif
*
         write(300,*) bndry,cdabs(phi(n))*cdabs(phi(n)),
     &        cdabs(dph_out(n))*cdabs(dph_out(n))
      enddo
*
      CLOSE(300)
*
#endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                  P L O T  W A V E F U N C T I O N
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#if defined(PLOT_WF)
*
      dx=(xmax-xmin)/dble(ixmax)
      dy=(ymax-ymin)/dble(iymax)
*
      write(200,*) '# content = "wavefunction"'
      write(200,*) "# NBE=",NBE
*
#if defined(TM)
      write(200,*) "# polarization = TM"
#elif defined(TE)
      write(200,*) "# polarization = TE"
#endif
*
      write(200,*) "# nin=",nin
      write(200,*) "# nout=",nout
#if (SYM == 4 || SYM == 2)
      write(200,*) "# a=",a
      write(200,*) "# b=",b
#elif SYM ==1
      write(200,*) "# b=",b
#endif
*
      write(200,*) "# kx=",kx
      write(200,*) "# ky=",ky
      write(200,*) "# lambda_in=",lambda_in
      write(200,*) "# ds_max=",ds_max
      write(200,*) "# 0.5*lambda_in/ds_max=",0.5d0*lambda_in/ds_max
#if defined(CUSTOMIZED)
      write(200,*) '# Dispersion relation = customized'
#else
      write(200,*) '# Dispersion relation = standard'
#endif
      write(200,*) "# BZW_HANKEL=",BZW_HANKEL
      write(200,*) "# BZW_MAXINT=",BZW_MAXINT
      write(200,*) "# nx=",ixmax
      write(200,*) "# ny=",iymax
      write(200,*) "# xmin=",xmin
      write(200,*) "# xmax=",xmax
      write(200,*) "# ymin=",ymin
      write(200,*) "# ymax=",ymax
*
      if(dimag(k).LT.0d0) then
         write(200,1000) ' # toplabel = "k =',
     &        dble(k),dimag(k),' i"'
      else
         write(200,1001) ' # toplabel = "k =',
     &        dble(k),' + ',dimag(k),'i"'
      endif
*
#if (SYM == 4 || SYM == 2)
      write(200,1002) ' # subtitle = "Parity: a=',a,", b=",b,'"'
#elif SYM == 1
      write(200,1003) ' # subtitle = "Parity: b=',b,'"'
#elif SYM == 0
      write(200,1004) ' # subtitle = ""'
#endif
*
 1000 format(A18,F14.10,F14.10,A3)
 1001 format(A18,F14.10,A3,F14.10,A3)
 1002 format(A25,I2,A4,I2,A1)
 1003 format(A25,I2,A1)
 1004 format(A25)
*
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      Estimate the Maximum Value (vmax)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
      vmax=0d0
      do iy=1,iymax
         do ix=1,ixmax
            x=xmin+dx*(dble(ix)-0.5)
            y=ymin+dy*(dble(iy)-0.5)
*
#if SYM == 4
            if((x.GE.0d0).AND.(y.GE.0d0)) then
               call far_from_bndry(x,y,flag0)
               call indicator_func(x,y,flag)
*
               if((flag.EQ.1).AND.(flag0.EQ.1)) then 
C              If (x,y) is inside the cavity and 
C              sufficiently far from the boundary
                  wf=dcmplx(0d0,0d0)
                  do n=1,nbe
                     call wfunc(a,b,nk_in,xl(n),yl(n),x,y,nx(n),ny(n),
     &                    ds(n),kappa(n),phi(n),dph_in(n),z)
                     wf=wf+z
                  enddo
                  intensity=cdabs(wf)*cdabs(wf)
                  if(intensity.GT.vmax) then
                     vmax=intensity
                  endif
               endif
            endif
#elif SYM == 2
            if((x.GE.0d0).AND.(y.GE.0d0)) then
               call far_from_bndry(x,y,flag0)
               call indicator_func(x,y,flag)
*
               if((flag.EQ.1).AND.(flag0.EQ.1)) then 
C              If (x,y) is inside the cavity and 
C              sufficiently far from the boundary
                  wf=dcmplx(0d0,0d0)
                  do n=1,nbe
                     call wfunc(a,b,nk_in,xl(n),yl(n),x,y,nx(n),ny(n),
     &                    ds(n),kappa(n),phi(n),dph_in(n),z)
                     wf=wf+z
                  enddo
                  intensity=cdabs(wf)*cdabs(wf)
                  if(intensity.GT.vmax) then
                     vmax=intensity
                  endif
               endif
            endif
#elif SYM == 1
            if(y.GE.0d0) then
               call far_from_bndry(x,y,flag0)
               call indicator_func(x,y,flag)
*
               if((flag.EQ.1).AND.(flag0.EQ.1)) then 
C              If (x,y) is inside the cavity and 
C              sufficiently far from the boundary
                  wf=dcmplx(0d0,0d0)
                  do n=1,nbe
                     call wfunc(b,nk_in,xl(n),yl(n),x,y,nx(n),ny(n),
     &                    ds(n),kappa(n),phi(n),dph_in(n),z)
                     wf=wf+z
                  enddo
                  intensity=cdabs(wf)*cdabs(wf)
                  if(intensity.GT.vmax) then
                     vmax=intensity
                  endif
               endif
            endif
#elif SYM == 0
            call far_from_bndry(x,y,flag0)
            call indicator_func(x,y,flag)
*
            if((flag.EQ.1).AND.(flag0.EQ.1)) then 
C           If (x,y) is inside the cavity and 
C           sufficiently far from the boundary
               wf=dcmplx(0d0,0d0)
               do n=1,nbe
                  call wfunc(nk_in,xl(n),yl(n),x,y,nx(n),ny(n),
     &                 ds(n),kappa(n),phi(n),dph_in(n),z)
                  wf=wf+z
               enddo
               intensity=cdabs(wf)*cdabs(wf)
               if(intensity.GT.vmax) then
                  vmax=intensity
               endif
            endif
#endif
         enddo
      enddo
*
      write(200,*) '# vmax=',vmax
C      stop
*
      write(200,*) "# ----------------------------------------"
      write(200,*) "#      [Re(phi)]           [Im(phi)]"
      write(200,*) "# ----------------------------------------"
*
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      Main Routine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do iy=1,iymax
         do ix=1,ixmax
            x=xmin+dx*(dble(ix)-0.5)
            y=ymin+dy*(dble(iy)-0.5)
*
            call indicator_func(x,y,flag)
*
            wf=dcmplx(0d0,0d0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      INSIDE THE CAVITY
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(flag.EQ.1) then
               do n=1,nbe
#if (SYM == 4 || SYM == 2)
                  call wfunc(a,b,nk_in,xl(n),yl(n),x,y,nx(n),ny(n),
     &                 ds(n),kappa(n),phi(n),dph_in(n),z)
#elif SYM == 1
                  call wfunc(b,nk_in,xl(n),yl(n),x,y,nx(n),ny(n),
     &                 ds(n),kappa(n),phi(n),dph_in(n),z)
#elif SYM == 0
                  call wfunc(nk_in,xl(n),yl(n),x,y,nx(n),ny(n),
     &                 ds(n),kappa(n),phi(n),dph_in(n),z)
#endif
                  wf=wf+z
               enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      OUTSIDE THE CAVITY
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            else if(flag.EQ.0) then
               do n=1,nbe
#if (SYM == 4 || SYM == 2)
                  call wfunc(a,b,nk_out,xl(n),yl(n),x,y,nx(n),ny(n),
     &                 ds(n),kappa(n),phi(n),dph_out(n),z)
#elif SYM == 1
                  call wfunc(b,nk_out,xl(n),yl(n),x,y,nx(n),ny(n),
     &                 ds(n),kappa(n),phi(n),dph_out(n),z)
#elif SYM == 0
                  call wfunc(nk_out,xl(n),yl(n),x,y,nx(n),ny(n),
     &                 ds(n),kappa(n),phi(n),dph_out(n),z)
#endif
                  wf=wf-z
               enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
            else
               write(*,*) 'Error: flag must be 0 or 1'
               stop
            endif
*
C            write(200,998) cdabs(wf)*cdabs(wf)
            write(200,998) dble(wf), dimag(wf)
*
         enddo
      enddo
 998  format(2E20.8) 
*
#endif
*
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c                           M A I N   E N D
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
C     stop ! commented out to avoid the message: The following floating-point exceptions are signalling: IEEE_DENORMAL
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
      subroutine far_from_bndry(x,y,flag0)
*
      integer nbe
      double precision bzw_maxint
      parameter (nbe=NBE)
      parameter (bzw_maxint=dble(BZW_MAXINT)*dble(SIZE_PARAM))
*
      integer n,flag0
      double precision x,y,dist_min,dist
      double precision ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe),nx(1:nbe),ny(1:nbe)
*
      call def_bndry(nbe,ds,kappa,xl,yl,nx,ny)
*
      dist_min=1d99
      do n=1,nbe
         dist = (xl(n)-x)*(xl(n)-x) + (yl(n)-y)*(yl(n)-y)
         dist = dsqrt(dist)
         if(dist_min.GT.dist) then
            dist_min=dist
         endif
      enddo
*
      if(dist_min.GE.bzw_maxint) then
         flag0=1
      else
         flag0=0
      endif
*
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
#if (SYM == 4 || SYM == 2)
*
      subroutine wfunc(a,b,nk,xl,yl,x,y,nx,ny,ds,kappa,phi,dph,w)
*
#elif SYM == 1
*
      subroutine wfunc(b,nk,xl,yl,x,y,nx,ny,ds,kappa,phi,dph,w)
*
#elif SYM == 0
*
      subroutine wfunc(nk,xl,yl,x,y,nx,ny,ds,kappa,phi,dph,w)
*
#else
      write(*,*) "Error in def_matrix.F: SYM value is not appropriate."
      stop
#endif     
*
      implicit NONE
*
      integer a,b
      double precision bzw_hankel
      double precision pi,ad,bd
      double precision xl,yl,x,y,nx,ny,kappa,ds
      double complex imu,nk,phi,dph
      double complex w,w1,w2,w3,w4,w5,w6,w7,w8
*
      imu=dcmplx(0d0,1d0)
      pi=4d0*datan(1d0)
      bzw_hankel=dble(BZW_HANKEL)*dble(SIZE_PARAM)
*
#if SYM == 4
C
C     When the cavity has mirror symmetries with respect to x=y and the y-axis.
C
      ad=dble(a) ! Parity index for the mirror symmetry w.r.t. the y-axis.
      bd=dble(b) ! Parity index for the mirror symmetry w.r.t. x=y.
*
      w=0d0
*
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTRIBUTION 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call dgfn_redef(nk,xl,yl,x,y,nx,ny,kappa,bzw_hankel,w1)
      call dgfn_redef(nk,yl,xl,x,y,ny,nx,kappa,bzw_hankel,w2)
      call dgfn_redef(nk,-yl,xl,x,y,-ny,nx,kappa,bzw_hankel,w3)
      call dgfn_redef(nk,-xl,yl,x,y,-nx,ny,kappa,bzw_hankel,w4)
      call dgfn_redef(nk,-xl,-yl,x,y,-nx,-ny,kappa,bzw_hankel,w5)
      call dgfn_redef(nk,-yl,-xl,x,y,-ny,-nx,kappa,bzw_hankel,w6)
      call dgfn_redef(nk,yl,-xl,x,y,ny,-nx,kappa,bzw_hankel,w7)
      call dgfn_redef(nk,xl,-yl,x,y,nx,-ny,kappa,bzw_hankel,w8)
*
      w=w+(w1 + bd*w2 + ad*bd*w3 + ad*w4
     &     + w5 + bd*w6 + ad*bd*w7 + ad*w8)*phi*ds
*     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTRIBUTION 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call gfn_redef(nk,xl,yl,x,y,ds,bzw_hankel,w1)
      call gfn_redef(nk,yl,xl,x,y,ds,bzw_hankel,w2)
      call gfn_redef(nk,-yl,xl,x,y,ds,bzw_hankel,w3)
      call gfn_redef(nk,-xl,yl,x,y,ds,bzw_hankel,w4)
      call gfn_redef(nk,-xl,-yl,x,y,ds,bzw_hankel,w5)
      call gfn_redef(nk,-yl,-xl,x,y,ds,bzw_hankel,w6)
      call gfn_redef(nk,yl,-xl,x,y,ds,bzw_hankel,w7)
      call gfn_redef(nk,xl,-yl,x,y,ds,bzw_hankel,w8)
*
      w=w-(w1 + bd*w2 + ad*bd*w3 + ad*w4
     &     + w5 + bd*w6 + ad*bd*w7 + ad*w8)*dph*ds
*
#elif SYM == 2
C
C     When the cavity has mirror symmetries with respect to both x- and y-axis.
C
      ad=dble(a) ! Parity index for the mirror symmetry w.r.t. the y-axis.
      bd=dble(b) ! Parity index for the mirror symmetry w.r.t. the x-axis.
*
      w=0d0
*
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTRIBUTION 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call dgfn_redef(nk,xl,yl,x,y,nx,ny,kappa,bzw_hankel,w1)
      call dgfn_redef(nk,-xl,yl,x,y,-nx,ny,kappa,bzw_hankel,w2)
      call dgfn_redef(nk,-xl,-yl,x,y,-nx,-ny,kappa,bzw_hankel,w3)
      call dgfn_redef(nk,xl,-yl,x,y,nx,-ny,kappa,bzw_hankel,w4)
*
      w=w+(w1 + ad*w2 + ad*bd*w3 + bd*w4)*phi*ds
*     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTRIBUTION 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call gfn_redef(nk,xl,yl,x,y,ds,bzw_hankel,w1)
      call gfn_redef(nk,-xl,yl,x,y,ds,bzw_hankel,w2)
      call gfn_redef(nk,-xl,-yl,x,y,ds,bzw_hankel,w3)
      call gfn_redef(nk,xl,-yl,x,y,ds,bzw_hankel,w4)
*
      w=w-(w1 + ad*w2 + ad*bd*w3 + bd*w4)*dph*ds
*
#elif SYM == 1
C
C     When the cavity has a mirror symmetry with respect to the x-axis.
C
      bd=dble(b) ! Parity index for the mirror symmetry w.r.t. the x-axis.
*
      w=0d0
*
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTRIBUTION 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call dgfn_redef(nk,xl,yl,x,y,nx,ny,kappa,bzw_hankel,w1)
      call dgfn_redef(nk,xl,-yl,x,y,nx,-ny,kappa,bzw_hankel,w4)
*
      w=w+(w1 + bd*w4)*phi*ds
*     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTRIBUTION 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call gfn_redef(nk,xl,yl,x,y,ds,bzw_hankel,w1)
      call gfn_redef(nk,xl,-yl,x,y,ds,bzw_hankel,w4)
*
      w=w-(w1 + bd*w4)*dph*ds
*
#elif SYM == 0
C
C     When the cavity has no mirror symmetries.
C
       w=0d0
*
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTRIBUTION 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call dgfn_redef(nk,xl,yl,x,y,nx,ny,kappa,bzw_hankel,w1)
*
      w=w+w1*phi*ds
*     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTRIBUTION 2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call gfn_redef(nk,xl,yl,x,y,ds,bzw_hankel,w1)
*
      w=w-w1*dph*ds
*
#endif
*
      return
      end
********************************************************************************
