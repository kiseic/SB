c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine set_parameters()
*
      implicit NONE
*
      double precision R,d,Xd
      common /const/ R,d,Xd
*
      R=1d0
      d=0.22d0
      Xd=0.391d0
*
      R=R*(SIZE_PARAM)
      d=d*(SIZE_PARAM)
      Xd=Xd*(SIZE_PARAM)
*
      return
      end
*******************************************************************************
*
      subroutine test_cavity_symmetry(symmetry)
*
      implicit NONE
*
      integer CAVITY_SYM_CLASS,symmetry
*
      CAVITY_SYM_CLASS=1
*
C     CAVITY_SYM_CLASS=0: The cavity has no mirror symmetry.
C     CAVITY_SYM_CLASS=1: The cavity is symmetric with respect to the x-axis.
C     CAVITY_SYM_CLASS=2: The cavity is symmetric with respect to both x- and y-axis.
C     CAVITY_SYM_CLASS=4: The cavity is symmetric with respect to x=y and the y-axis.
*
      if(symmetry.ne.CAVITY_SYM_CLASS) then
         write(*,*) "Error in the parameter SYM."
         stop
      endif
*
      return
      end
*******************************************************************************
*
      subroutine def_bndry(nbe,ds,kappa,xl,yl,nx,ny)
*
      implicit NONE
*
      integer n,nbe,N1,N2,m
      double precision R,d,Xd,ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe)
      double precision nx(1:nbe),ny(1:nbe)
      double precision pi,S1,S2,eta,x0,y0,phi,dphi
      common /const/ R,d,Xd
*     
      pi=4d0*datan(1d0)
*
      call set_parameters()
*
      S1=pi*R
      S2=pi*d
*
      N1=aint(dble(nbe)*(S1/(S1+S2)))
      N2=nbe-N1
*      
C==============================
C     Segment I (Outer Disk) 
C==============================
      dphi=pi/dble(N1)
*
      do n=1,N1
         phi=dphi*(dble(n-1)+0.5)
         xl(n)=R*cos(phi)
         yl(n)=R*sin(phi)
         nx(n)=cos(phi)
         ny(n)=sin(phi)
         ds(n)=R*dphi
         kappa(n)=1d0/R
      enddo
*
C===============================
C     Segment II (Inner Disk)
C===============================
*     
      x0=Xd
      y0=0d0
      dphi=pi/dble(N2)
*
      do n=1,N2
         phi=pi-dphi*(dble(n-1)+0.5)
         m=n+N1
         xl(m)=x0+d*dcos(phi)
         yl(m)=y0+d*dsin(phi)
         nx(m)=-cos(phi)
         ny(m)=-sin(phi)
         ds(m)=d*dphi
         kappa(m)=-1d0/d
      enddo
*
      return
      end
*
*******************************************************************************
*
      subroutine indicator_func(x,y,flag)
*
      integer flag
      double precision R,d,Xd,x,y,rad
      common /const/ R,d,Xd
*
      call set_parameters()
*
      rad=dsqrt((x-Xd)*(x-Xd)+y*y)
*
      if(rad.LT.d) then
         flag=0
      else
         rad=dsqrt(x*x+y*y)
*
         if(rad.LT.R) then
            flag=1
         else
            flag=0
         endif
      endif
*
      return
      end
*
