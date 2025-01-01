c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine set_parameters()
*
      implicit NONE
*
      double precision R,L,d
      common /const/ R,L,d
*
      R=0.5d0*(SIZE_PARAM)
      L=1d0*(SIZE_PARAM)
      d=0.2d0*(SIZE_PARAM)
*
      if(R.GT.L) then
         write(*,*) "Error: R must be less than L."
         stop
      endif
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
      CAVITY_SYM_CLASS=2
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
      integer n,m,nbe,N1,N2,N3,N4
      double precision pi,R,L,d,S,ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe)
      double precision nx(1:nbe),ny(1:nbe)
      double precision phi,dphi,cs,sn,dL,dx,x0
      common /const/ R,L,d
*
      pi=4d0*datan(1d0)
*     
      call set_parameters()
*
      S = pi*d/2d0 + 2d0*(L-d) + pi*R/2d0
      N1=aint(pi*d/4d0/S*dble(nbe))
      N2=aint(2d0*(L-d)/S*dble(nbe))
      N3=N1
      N4=nbe-N1-N2-N3
*
************************************************
*     Circular segment (outer,right) 
************************************************
      dphi=pi/4d0/dble(N1)
*
      do n=1,N1
         phi=dphi*dble(n-0.5d0)
         cs=dcos(phi)
         sn=dsin(phi)
*
         xl(n)=dsqrt(2d0)*(L-d)+d*cs
         yl(n)=d*sn
*
         ds(n)=d*dphi
         nx(n)=cs
         ny(n)=sn
         kappa(n)=1d0/d
      enddo
*
*********************************************
*     Linear segment
*********************************************
      dL=2d0*(L-d)/dble(N2)
*
      x0=dsqrt(2d0)*L-d/dsqrt(2d0)
      dx=dL/dsqrt(2d0)
*
      do n=1,N2
         m=N1+n
         xl(m)=x0-dx*(dble(n)-0.5d0)
         yl(m)=-xl(m)+dsqrt(2d0)*L
*
         ds(m)=dL
         nx(m)=1d0/dsqrt(2d0)
         ny(m)=1d0/dsqrt(2d0)
         kappa(m)=0d0
      enddo
*
************************************************
*     Circular segment (outer,top) 
************************************************
      dphi=pi/4d0/dble(N3)
*
      do n=1,N3
         m=N1+N2+n
         phi=pi/4d0+dphi*dble(n-0.5d0)
         cs=dcos(phi)
         sn=dsin(phi)
*
         xl(m)=d*cs
         yl(m)=dsqrt(2d0)*(L-d)+d*sn
*
         ds(m)=d*dphi
         nx(m)=cs
         ny(m)=sn
         kappa(m)=1d0/d
      enddo
*
************************************************
*     Circular segment (inner)
************************************************
      dphi=pi/2d0/dble(N4)
*
      do n=1,N4
         m=N1+N2+N3+n
*
         phi=dphi*dble(n-0.5d0)
         cs=dcos(phi)
         sn=dsin(phi)
*
         xl(m)=R*cs
         yl(m)=R*sn
*
         ds(m)=R*dphi
         nx(m)=-cs
         ny(m)=-sn
         kappa(m)=1d0/R
*
      enddo
*
      return
      end
*******************************************************************************
*
      subroutine indicator_func(x,y,flag)
*
      implicit NONE
*
      integer flag
      double precision R,L,d,x,y,xx,yy,buf
      double precision rad,rad2,dx,dy,xmax
      common /const/ R,L,d
*
      call set_parameters()
*
      if((x.GE.0).AND.(y.GE.0)) then
         xx=x
         yy=y
      else if((x.LT.0).AND.(y.GE.0)) then
         xx=-x
         yy=y
      else if((x.LT.0).AND.(y.LT.0)) then
         xx=-x
         yy=-y
      else
         xx=x
         yy=-y
      endif
*
      if(yy.GT.xx) then
         buf=xx
         xx=yy
         yy=buf
      endif
*
      rad = dsqrt( xx*xx + yy*yy )
*
      if(rad.LT.R) then
         flag=0
      else
*
         dx=xx-dsqrt(2d0)*(L-d)
         dy=yy
         rad2 = dsqrt( dx*dx + dy*dy )
*
         xmax=dsqrt(2d0)*L-d/dsqrt(2d0)
*     
         if((xx.LT.xmax).AND.(yy.LT.(-xx+dsqrt(2d0)*L))) then
            flag=1
         else if(rad2.LT.d) then
            flag=1
         else
            flag=0
         endif            
      endif
*
      return
      end
*
*******************************************************************************
