c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      program determinant
*
      implicit NONE
*
      integer nbe
      double precision nin,nout
*     
      parameter (nbe=NBE)
*
      integer n,i,j,l,ix,iy,m,sign,a,b,ixmax,iymax,beta,symmetry
      double precision pi
      double precision ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe),nx(1:nbe),ny(1:nbe)
      double precision kx,ky,kx0,ky0,dx,dy,cx,cy,dwx,dwy
      double precision re_nk_max,lambda_min,ds_max,alpha
      double precision sum,lndet,det
      double complex imu,k,nk_in,nk_out
      double complex P1(1:nbe,1:nbe),Q1(1:nbe,1:nbe)
      double complex P2(1:nbe,1:nbe),Q2(1:nbe,1:nbe)
*
      pi=4d0*datan(1d0)
      imu=dcmplx(0d0,1d0)
*
C     Cavity Symmetry Check
#if (SYM > 4 || SYM == 3)
      write(*,*) "Error in the parameter SYM."
      stop
#endif
      symmetry=SYM
      call test_cavity_symmetry(symmetry)
*
#if (SYM == 4 || SYM == 2)
      read(5,*) nin,nout,a,b,cx,cy,dwx,dwy,dx,dy
#elif SYM == 1
      read(5,*) nin,nout,b,cx,cy,dwx,dwy,dx,dy
#elif SYM == 0
      read(5,*) nin,nout,cx,cy,dwx,dwy,dx,dy
#endif
*
      call def_bndry(nbe,ds,kappa,xl,yl,nx,ny)
*
      ds_max=0d0
      do i=1,nbe
         if(ds(i).GT.ds_max) then
            ds_max=ds(i)
         endif
      enddo
*
      ixmax=aint(2d0*dwx/dx)
      iymax=aint(2d0*dwy/dy)
      kx0=cx-dwx
      ky0=cy-dwy
*
      write(*,*) '# nbe = ',nbe
      write(*,*) '# nin = ',nin,' nout = ',nout
      write(*,*) '# a = ',a,' b = ',b
      write(*,*) '# cx = ',cx,' cy = ',cy
      write(*,*) '# dwx = ',dwx,' dwy = ',dwy
      write(*,*) '# dx = ',dx,' dy = ',dy
      write(*,*) '# ixmax = ',ixmax,' iymax = ',iymax
*
      k=dcmplx(cx+dwx,0d0)
*
#include "def_nk.F"
*
      re_nk_max=dble(nk_in)
*
      lambda_min=2d0*pi/re_nk_max
      write(*,*) '# ds_max = ',ds_max
      write(*,*) '# lambda_min = ',lambda_min
      write(*,*) '# 0.5*lambda_min/ds_max = ',0.5d0*lambda_min/ds_max
#if defined(CUSTOMIZED)
      write(*,*) '# Dispersion relation = customized'
#else
      write(*,*) '# Dispersion relation = standard'
#endif
*
*************************************************************************
*     Determine the value of 'beta'
*************************************************************************
      k=dcmplx(cx,cy)
*     
#include "def_nk.F"
      sign=1
#if (SYM == 4 || SYM == 2)
      call def_matrix(nbe,nk_in,a,b,xl,yl,nx,ny,sign,ds,kappa,P1,Q1)
#elif SYM == 1
      call def_matrix(nbe,nk_in,b,xl,yl,nx,ny,sign,ds,kappa,P1,Q1)
#elif SYM == 0
      call def_matrix(nbe,nk_in,xl,yl,nx,ny,sign,ds,kappa,P1,Q1)
#endif
*
#include "def_nk.F"
      sign=-1
#if (SYM == 4 || SYM == 2)
      call def_matrix(nbe,nk_out,a,b,xl,yl,nx,ny,sign,ds,kappa,P2,Q2)
#elif SYM == 1
      call def_matrix(nbe,nk_out,b,xl,yl,nx,ny,sign,ds,kappa,P2,Q2)
#elif SYM == 0
      call def_matrix(nbe,nk_out,xl,yl,nx,ny,sign,ds,kappa,P2,Q2)
#endif
*
#if defined(TE)
      do i=1,nbe
         do j=1,nbe
            P1(i,j)=nin*nin*P1(i,j)
            P2(i,j)=nout*nout*P2(i,j)
         enddo
      enddo
#endif
*
      call calc_detQ2(P1,P2,Q1,Q2,lndet)
      beta=aint(dlog(dabs(lndet))/dlog(10d0))
*
      write(*,*) '# beta = ',beta,'; lndet = ',lndet
*      
*************************************************************************
*     M A I N  L O O P
*************************************************************************
      do iy=1,iymax
         do ix=1,ixmax
            kx=kx0+dble(ix)*dx
            ky=ky0+dble(iy)*dy
*     
            k=dcmplx(kx,ky)
*
#include "def_nk.F"
            sign=1
#if (SYM == 4 || SYM == 2)
            call def_matrix(nbe,nk_in,a,b,xl,yl,nx,ny,sign,
     &           ds,kappa,P1,Q1)
#elif SYM == 1
            call def_matrix(nbe,nk_in,b,xl,yl,nx,ny,sign,
     &           ds,kappa,P1,Q1)
#elif SYM == 0
            call def_matrix(nbe,nk_in,xl,yl,nx,ny,sign,
     &           ds,kappa,P1,Q1)
#endif
*
#include "def_nk.F"
            sign=-1
#if (SYM == 4 || SYM == 2)
            call def_matrix(nbe,nk_out,a,b,xl,yl,nx,ny,sign,
     &           ds,kappa,P2,Q2)
#elif SYM == 1
            call def_matrix(nbe,nk_out,b,xl,yl,nx,ny,sign,
     &           ds,kappa,P2,Q2)
#elif SYM == 0
            call def_matrix(nbe,nk_out,xl,yl,nx,ny,sign,
     &           ds,kappa,P2,Q2)
#endif
*
#if defined(TE)
            do i=1,nbe
               do j=1,nbe
                  P1(i,j)=nin*nin*P1(i,j)
                  P2(i,j)=nout*nout*P2(i,j)
               enddo
            enddo
#endif
*
            call calc_detQ2(P1,P2,Q1,Q2,lndet)
            lndet=lndet*10d0**dble(-beta)
            det=exp(lndet)
*
*****************************************************************************
*     OUTPUT
*****************************************************************************
            write(*,200) dble(k),imag(k),det
         enddo
      enddo
*     
200   FORMAT(2F20.10,E20.8E3)
*******************************************************************************
      stop
      end
**************************************************************************
*
      subroutine calc_detQ2(P1,P2,Q1,Q2,lndet)
*
      implicit NONE
*
      integer i,j,l,info,ipiv(NBE)
      double precision lndet
      double complex work(NBE)
      double complex P1(1:NBE,1:NBE),Q1(1:NBE,1:NBE)
      double complex P2(1:NBE,1:NBE),Q2(1:NBE,1:NBE)
      double complex P1inv(1:NBE,1:NBE),Q2inv(1:NBE,1:NBE)
      double complex DUM1(1:NBE,1:NBE),DUM2(1:NBE,1:NBE)
      double complex MTRX(1:NBE,1:NBE)
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     M = Q2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,NBE
         do j=1,NBE
            MTRX(i,j)=Q2(i,j)
         enddo
      enddo
      call zgetrf(NBE,NBE,MTRX,NBE,ipiv,info)  
*      
      lndet=0d0
      do i=1,NBE
         lndet=lndet+dlog(cdabs(MTRX(i,i)))
      enddo
*
      return
      end
*
