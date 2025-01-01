*
      program husimi_distribution
*
      implicit NONE
*
c#define CALC_HUSIMI_DIST
c#define CALC_FLUX_DIST
#define CALC_INC_HUSIMI
c#define CALC_REFL_HUSIMI
c#define CALC_INC_HUSIMI_OUTSIDE
c#define CALC_REFL_HUSIMI_OUTSIDE
*
c#define PLOT_WHOLE_AREA
c#define PLOT_OPEN_AREA
#define PLOT_HALF_AREA      
*
#define NORMALIZED_TO_1
c#define NORMALIZED_BY_HUSIMI
*
      integer nbe,nbetot,ixmax,iymax
*
      parameter (nbe=NBE)
#if SYM == 4
      parameter (nbetot=NBE*8)
#elif SYM == 2
      parameter (nbetot=NBE*4)
#elif SYM == 1
      parameter (nbetot=NBE*2)
#elif SYM == 0
      parameter (nbetot=NBE)
#endif
*
       parameter (ixmax=IXMAX)
       parameter (iymax=IYMAX)
*
      integer n,nbe2,i,j,s,ix,iy,m,sign,nn,flag,a,b
      integer info,lda,ldvl,ldvr,lwork
      double precision nin,nout
      double precision rwork(nbe*4)
      double precision pi,sigma,k_re,ad,bd,fact,RC,pc
      double precision ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe),nx(1:nbe),ny(1:nbe)
      double precision Ls(1:nbetot),ds_a(1:nbetot),Hsum,C_husimi
      double precision kx,ky,dx,dy,Ltot,HD(1:ixmax,1:iymax)
      double precision dq,dp,dqdp,p,pmin,flux
      double precision err,minerr,xmax,xmin,ymax,ymin,x,y,Ld,q,q0,p0
      double complex imu,k,nk_in,nk_out,alpha,z,h,h1,h2,h1i,h2i
      double complex P1(1:nbe,1:nbe),Q1(1:nbe,1:nbe)
      double complex P2(1:nbe,1:nbe),Q2(1:nbe,1:nbe)
      double complex T(1:nbe*2,1:nbe*2)
      double complex W(nbe*2),VL(nbe*2,nbe*2)
      double complex VR(nbe*2,nbe*2),work(nbe*4)
      double complex phi(1:nbe),dph_in(1:nbe),dph_out(1:nbe)
      double complex phi_a(1:nbetot),dph_a(1:nbetot)
      double complex dum,dum2
      character jobvl,jobvr
*     
      pi=4d0*datan(1d0)
      imu=dcmplx(0d0,1d0)
*
      OPEN(100,FILE='dat.husimi')
*
      call def_bndry(nbe,ds,kappa,xl,yl,nx,ny)
*
#if (SYM == 4 || SYM == 2)
      read(*,*) nin,nout,a,b,kx,ky
#elif SYM == 1
      read(*,*) nin,nout,b,kx,ky
#elif SYM == 0
      read(*,*) nin,nout,kx,ky
#else
      write(*,*) "Error in the parameter SYM."
      stop
#endif
*      
      k = kx + imu*ky
*
#include "def_nk.F"
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
#elif SYM == 0
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
c     CALC: W.F. NORMAL DERIVATIVE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
      alpha=1d0+imu
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
c      do n=1,nbe2
c         write(*,*) dble(W(n)),imag(W(n))
c      enddo
c      stop
*
      minerr=1d10
      do n=1,nbe2
         dx=dble(W(n))-dble(alpha)
         dy=imag(W(n))-imag(alpha)
         err=dsqrt(dx*dx + dy*dy)
*
         if(err.LT.minerr) then
            nn=n
            minerr=err
         endif
      enddo
*
c      write(*,*) W(nn)
c      stop
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
c      do n=1,nbe
c         write(*,*) n,cdabs(phi(n)),cdabs(dph(n))
c      enddo
c      stop
*
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   H U S I M I  P L O T
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
      ad=dble(a)
      bd=dble(b)
*
      do n=1,nbe
*
#if SYM == 4
         ds_a(n)=ds(n)
         ds_a(n+nbe)=ds(nbe-n+1)
         ds_a(n+nbe*2)=ds(n)
         ds_a(n+nbe*3)=ds(nbe-n+1)
         ds_a(n+nbe*4)=ds(n)
         ds_a(n+nbe*5)=ds(nbe-n+1)
         ds_a(n+nbe*6)=ds(n)
         ds_a(n+nbe*7)=ds(nbe-n+1)
*
         phi_a(n)=phi(n)
         phi_a(n+nbe)=bd*phi(nbe-n+1)
         phi_a(n+nbe*2)=ad*bd*phi(n)
         phi_a(n+nbe*3)=ad*phi(nbe-n+1)
         phi_a(n+nbe*4)=phi(n)
         phi_a(n+nbe*5)=bd*phi(nbe-n+1)
         phi_a(n+nbe*6)=ad*bd*phi(n)
         phi_a(n+nbe*7)=ad*phi(nbe-n+1)
*
         dph_a(n)=dph_in(n)
         dph_a(n+nbe)=bd*dph_in(nbe-n+1)
         dph_a(n+nbe*2)=ad*bd*dph_in(n)
         dph_a(n+nbe*3)=ad*dph_in(nbe-n+1)
         dph_a(n+nbe*4)=dph_in(n)
         dph_a(n+nbe*5)=bd*dph_in(nbe-n+1)
         dph_a(n+nbe*6)=ad*bd*dph_in(n)
         dph_a(n+nbe*7)=ad*dph_in(nbe-n+1)
*
#elif SYM == 2
         ds_a(n)=ds(n)
         ds_a(n+nbe)=ds(nbe-n+1)
         ds_a(n+nbe*2)=ds(n)
         ds_a(n+nbe*3)=ds(nbe-n+1)
*
         phi_a(n)=phi(n)
         phi_a(n+nbe)=ad*phi(nbe-n+1)
         phi_a(n+nbe*2)=ad*bd*phi(n)
         phi_a(n+nbe*3)=bd*phi(nbe-n+1)
*
         dph_a(n)=dph_in(n)
         dph_a(n+nbe)=ad*dph_in(nbe-n+1)
         dph_a(n+nbe*2)=ad*bd*dph_in(n)
         dph_a(n+nbe*3)=bd*dph_in(nbe-n+1)
*
#elif SYM == 1
*
         ds_a(n)=ds(n)
         ds_a(n+nbe)=ds(nbe-n+1)
*
         phi_a(n)=phi(n)
         phi_a(n+nbe)=bd*phi(nbe-n+1)
*
         dph_a(n)=dph_in(n)
         dph_a(n+nbe)=bd*dph_in(nbe-n+1)
*
#elif SYM == 0
         ds_a(n)=ds(n)
         phi_a(n)=phi(n)
         dph_a(n)=dph_in(n)
#endif
*
      enddo
*
      Ld=0d0
      do n=1,nbe
         Ld=Ld+ds(n)
      enddo
*
#if SYM == 4
      Ltot=Ld*8d0
#elif SYM == 2
      Ltot=Ld*4d0
#elif SYM == 1
      Ltot=Ld*2d0
#elif SYM == 0
      Ltot=Ld
#endif
*
      Ls(1)=0.5d0*ds_a(1)
      do n=1,nbetot-1
         Ls(n+1)=Ls(n)+0.5d0*(ds_a(n)+ds_a(n+1))
      enddo
*
#if defined(CUSTOMIZED)
      k_re=dsqrt(2d0*kx+1d0)
#else
      k_re=dble(nin)*kx
#endif
*
#if defined(CALC_INC_HUSIMI_OUTSIDE)
      k_re=kx
#elif defined(CALC_REFL_HUSIMI_OUTSIDE)
      k_re=kx
#endif
*
      sigma=dsqrt(Ltot/2d0/k_re)
      pc=dble(nout)/dble(nin)
*
      xmax=Ltot
      xmin=0d0
*
****************************************************************************
#if defined(NORMALIZED_BY_HUSIMI)
      ymax=k_re
      ymin=-k_re
*
      dx=(xmax-xmin)/dble(ixmax)
      dy=(ymax-ymin)/dble(iymax)
*
      dq=1d0/dble(ixmax)
      dp=2d0/dble(iymax)
      dqdp=dq*dp
*
      C_husimi=0d0
      do iy=1,iymax
         do ix=1,ixmax
            x=xmin+dx*(dble(ix)-0.5)
            y=ymin+dy*(dble(iy)-0.5)
*
            h1=dcmplx(0d0,0d0)
            do n=1,nbetot
*
               q0=x
               p0=y
               q=Ls(n)
*
               call coherent_state(q0,p0,q,sigma,Ltot,10,z)
               h1=h1+ds_a(n)*dconjg(z)*phi_a(n)
            enddo
*
            HD(ix,iy)=cdabs(h1)*cdabs(h1)
*
            C_husimi = C_husimi + HD(ix,iy)*dqdp
*
         enddo
      enddo
***************************************************************************
#endif
*
#if defined(PLOT_WHOLE_AREA)
      ymax=k_re
      ymin=-k_re
#elif defined(PLOT_OPEN_AREA)
      ymax=k_re*pc
      ymin=-k_re*pc
#elif defined(PLOT_HALF_AREA)
c      ymax=k_re*pc
      ymax=k_re
      ymin=0d0
#else
      write(*,*) 'Error !!'
      stop
#endif
*
      dx=(xmax-xmin)/dble(ixmax)
      dy=(ymax-ymin)/dble(iymax)
*
      dq=1d0/dble(ixmax)
      dp=2d0/dble(iymax)
*
#if defined(PLOT_OPEN_AREA)
      dp=dp*pc
#endif
*
      dqdp=dq*dp
*
      write(100,*) '# content = "husimi"'
      write(100,*) "# NBE=",NBE
*
#if defined(TM)
      write(100,*) "# polarization = TM"
#elif defined(TE)
      write(100,*) "# polarization = TE"
#endif
*
      write(100,*) "# nin=",nin
      write(100,*) "# nout=",nout
#if (SYM == 4 || SYM == 2)
      write(100,*) "# a=",a
      write(100,*) "# b=",b
#elif SYM ==1
      write(100,*) "# b=",b
#endif
*
      write(100,*) "# kx=",kx
      write(100,*) "# ky=",ky
      write(100,*) "# nx=",ixmax
      write(100,*) "# ny=",iymax
      write(100,*) "# xmin=",0.0
      write(100,*) "# xmax=",1.0
#if defined(PLOT_WHOLE_AREA)
      write(100,*) "# ymin=",-1.0
      write(100,*) "# ymax=",1.0
#elif defined(PLOT_OPEN_AREA)
      write(100,*) "# ymin=",-pc
      write(100,*) "# ymax=",pc
#elif defined(PLOT_HALF_AREA)
      write(100,*) "# ymin=",0.0
      write(100,*) "# ymax=",1.0
#endif
*
#if defined(CUSTOMIZED)
      write(100,*) '# Dispersion relation = customized'
#else
      write(100,*) '# Dispersion relation = standard'
#endif      
*
      write(100,*) "# xyratio=1.0"
*
      if(dimag(k).LT.0d0) then
         write(100,1000) ' # toplabel = "k =',
     &        dble(k),dimag(k),' i"'
      else
         write(100,1001) ' # toplabel = "k =',
     &        dble(k),' + ',dimag(k),'i"'
      endif
*
#if (SYM == 4 || SYM == 2)
      write(100,1002) ' # subtitle = "Parity: a=',a,", b=",b,'"'
#elif SYM == 1
      write(100,1003) ' # subtitle = "Parity: b=',b,'"'
#elif SYM == 0
      write(100,1004) ' # subtitle = ""'
#endif
*
 1000 format(A18,F14.10,F14.10,A3)
 1001 format(A18,F14.10,A3,F14.10,A3)
 1002 format(A25,I2,A4,I2,A1)
 1003 format(A25,I2,A1)
 1004 format(A25)
*
cccccc FLUX CALC ccccccccccccccccccccccccccccccccccccccccccc
c      flux=0d0
c      do n=1,nbetot
c         flux=flux+dimag(k*phi_a(n)*dconjg(dph_a(n)))*ds_a(n)
c      enddo
c      write(*,*) flux
c      stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
*
      Hsum=0d0
      do iy=1,iymax
         do ix=1,ixmax
            x=xmin+dx*(dble(ix)-0.5)
            y=ymin+dy*(dble(iy)-0.5)
*
            h1=dcmplx(0d0,0d0)
            h2=dcmplx(0d0,0d0)
            h1i=dcmplx(0d0,0d0)
            h2i=dcmplx(0d0,0d0)
            do n=1,nbetot
*
               q0=x
               p0=y
               q=Ls(n)
*
               call coherent_state(q0,p0,q,sigma,Ltot,10,z)
               h1=h1+ds_a(n)*dconjg(z)*phi_a(n)
               h2=h2+ds_a(n)*dconjg(z)*dph_a(n)
               h1i=h1i+ds_a(n)*dconjg(z)*dconjg(phi_a(n))
               h2i=h2i+ds_a(n)*dconjg(z)*dconjg(dph_a(n))
            enddo
*
#if defined(CALC_HUSIMI_DIST)
            HD(ix,iy)=cdabs(h1)*cdabs(h1)
#elif defined(CALC_FLUX_DIST)
            HD(ix,iy)=dimag(dconjg(h1)*h2) 
#elif defined(CALC_INC_HUSIMI)
            fact=dsqrt(dble(nin)*dsqrt(1d0-y*y/k_re/k_re))
            dum = -fact*h1 + imu/fact/kx*h2
            dum2=dum
C            dum2= -fact*h1i + imu/fact/kx*h2i
            HD(ix,iy)=cdabs(dum)*cdabs(dum2)
#elif defined(CALC_REFL_HUSIMI)
            fact=dsqrt(dble(nin)*dsqrt(1d0-y*y/k_re/k_re))
            dum = -fact*h1 - imu/fact/kx*h2
            HD(ix,iy)=cdabs(dum)*cdabs(dum)
#else
            write(*,*) 'Error !!'
            stop
#endif
*
            Hsum = Hsum + HD(ix,iy)*dqdp
*
         enddo
      enddo
*
#if defined(PLOT_WHOLE_AREA)
      pmin=-1d0
#elif defined(PLOT_OPEN_AREA)
      pmin=-pc
#endif
*
      do iy=1,iymax
         do ix=1,ixmax
c            p=pmin+dp*(dble(iy)-0.5)
c            call reflectance(p,pc,RC)
#if defined(NORMALIZED_TO_1)
            HD(ix,iy)=HD(ix,iy)/Hsum
#elif defined(NORMALIZED_BY_HUSIMI)
            HD(ix,iy)=HD(ix,iy)/C_husimi
#else
            HD(ix,iy)=HD(ix,iy)
#endif
            write(100,998) HD(ix,iy)
         enddo
      enddo
 998  format(E20.10)
*
      CLOSE(100)
*     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c                           M A I N   E N D
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
C      stop ! commented out to avoid the message: The following floating-point exceptions are signalling: IEEE_DENORMAL
      end
*
*******************************************************************************
*
      subroutine reflectance(p,pc,RC)
*
      implicit NONE
*
      double precision p,pc,RC,dum
*
      if(dabs(p).GE.pc) then
         RC=1d0
      else
         dum=dsqrt((pc*pc-p*p)/(1d0-p*p))
         RC=(1d0-dum)/(1d0+dum)
         RC=RC*RC
      endif
*
      return
      end
*
*******************************************************************************
*
      subroutine coherent_state(q0,p0,q,sigma,Lx,lmax,z)
*
      implicit NONE
*
      integer l,lmax
      double precision pi,q0,p0,q,d,Lx
      double precision sigma,fact
      double complex imu,z
*
      imu=dcmplx(0d0,1d0)
      pi=4d0*datan(1d0)
*
      z=dcmplx(0d0,0d0)
      do l=-lmax,lmax
         d=q-dble(l)*Lx-q0
         z=z+cdexp(-d*d/2d0/sigma/sigma + imu*p0*d)
      enddo
*
      fact=dsqrt(dsqrt(1d0/pi/sigma/sigma))
      z=fact*z
*
*******************************************************************************
      return
      end
*
