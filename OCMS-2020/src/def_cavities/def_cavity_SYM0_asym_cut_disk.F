c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine set_parameters()
*
      implicit NONE
*
      double precision pi,R,d1,d2,theta
      common /const/ R,d1,d2,theta
*
      pi=4d0*datan(1d0)
*
      R=1d0
      d1=0.5d0
      d2=0.9d0
      theta=(7d0-3d0*dsqrt(5d0))*pi
*
      R=R*(SIZE_PARAM)
      d1=d1*(SIZE_PARAM)
      d2=d2*(SIZE_PARAM)
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
      CAVITY_SYM_CLASS=0
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
      integer n,nn,nbe,nbe1,nbe2,nbe3,nbe4
      double precision pi,R,d1,d2,theta
      double precision ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe)
      double precision nx(1:nbe),ny(1:nbe)
      double precision phi0,phi1,phi2,phi3,xi,eta
      double precision x0,x1,x2,x3,y0,y1,y2,y3
      double precision S1,S2,S3,S4,Stot,arg,darg,dx,slope
      common /const/ R,d1,d2,theta
*
      pi=4d0*datan(1d0)
*
      call set_parameters()
*     
      xi=dasin(d1/R)
      eta=dacos(d2/R)
*
      phi0=pi-xi
      phi1=2d0*pi-(theta+eta)
      phi2=2d0*pi-(theta-eta)
      phi3=2d0*pi+xi
*
      x0=R*dcos(phi0)
      y0=R*dsin(phi0)
      x1=R*dcos(phi1)
      y1=R*dsin(phi1)
      x2=R*dcos(phi2)
      y2=R*dsin(phi2)
      x3=R*dcos(phi3)
      y3=R*dsin(phi3)
*
      S1=(phi1-phi0)*R
      S2=dsqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
      S3=(phi3-phi2)*R
      S4=dsqrt((x3-x0)*(x3-x0)+(y3-y0)*(y3-y0))
*
      Stot=S1+S2+S3+S4
*
      nbe1=nint(S1/Stot*dble(nbe))
      nbe2=nint(S2/Stot*dble(nbe))
      nbe3=nint(S3/Stot*dble(nbe))
      nbe4=nbe-(nbe1+nbe2+nbe3)
*
**********************************************
*     Circular part (left)                   *
**********************************************
      darg=(phi1-phi0)/dble(nbe1)
      do n=1,nbe1
         arg=phi0+darg*(dble(n)-0.5d0)
         xl(n)=R*dcos(arg)
         yl(n)=R*dsin(arg)
         ds(n)=R*darg
         kappa(n)=1d0/R
         nx(n)=dcos(arg)
         ny(n)=dsin(arg)
      enddo
*
**********************************************
*     Linear part (bottom)                   *
**********************************************
      dx=(x2-x1)/dble(nbe2)
      slope=(y2-y1)/(x2-x1)
      do n=1,nbe2
         nn=n+nbe1
         xl(nn)=x1+dx*(dble(n)-0.5d0)
         yl(nn)=slope*(xl(nn)-x1)+y1
         ds(nn)=dsqrt(1d0+slope*slope)*dx
         kappa(nn)=0d0
         nx(nn)=slope/dsqrt(1d0+slope*slope)
         ny(nn)=-1d0/dsqrt(1d0+slope*slope)
      enddo
*
**********************************************
*     Circular part (right)                  *
**********************************************
      darg=(phi3-phi2)/dble(nbe3)
      do n=1,nbe3
         nn=n+nbe1+nbe2
         arg=phi2+darg*(dble(n)-0.5d0)
         xl(nn)=R*dcos(arg)
         yl(nn)=R*dsin(arg)
         ds(nn)=R*darg
         kappa(nn)=1d0/R
         nx(nn)=dcos(arg)
         ny(nn)=dsin(arg)
      enddo
*
**********************************************
*     Linear part (top)                      *
**********************************************
      dx=(x0-x3)/dble(nbe4)
      slope=(y0-y3)/(x0-x3)
      do n=1,nbe4
         nn=n+nbe1+nbe2+nbe3
         xl(nn)=x3+dx*(dble(n)-0.5d0)
         yl(nn)=slope*(xl(nn)-x3)+y3
         ds(nn)=dsqrt(1d0+slope*slope)*(-dx)
         kappa(nn)=0d0
         nx(nn)=0d0
         ny(nn)=1d0
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
      double precision pi,x,y
      double precision rad,phi1,phi2,x1,x2,y1,y2,yy
      double precision R,d1,d2,theta
      common /const/ R,d1,d2,theta
*
      pi=4d0*datan(1d0)
*
      call set_parameters()
*
      phi1=-(theta+dacos(d2/R))
      phi2=-(theta-dacos(d2/R))
*
      x1=R*dcos(phi1)
      y1=R*dsin(phi1)
*
      x2=R*dcos(phi2)
      y2=R*dsin(phi2)
*
      flag=0
*
      rad=dsqrt(x*x+y*y)
      if(rad.LT.R) then
         yy=(y2-y1)/(x2-x1)*(x-x1)+y1
         if((y.LT.d1).AND.(yy.LT.y)) then
            flag=1
         endif
      endif
*
      return
      end
*
