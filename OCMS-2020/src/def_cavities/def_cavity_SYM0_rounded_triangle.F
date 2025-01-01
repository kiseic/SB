c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine set_parameters()
*
      implicit NONE
*
      double precision R,r1,r2,r3
      common /const/ R,r1,r2,r3
*
      R=1d0
      r1=0.1d0
      r2=0.2d0
      r3=0.3d0
*
      R=R*(SIZE_PARAM)
      r1=r1*(SIZE_PARAM)
      r2=r2*(SIZE_PARAM)
      r3=r3*(SIZE_PARAM)
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
      integer n,nn,nbe,nbe1,nbe2,nbe3,nbe4,nbe5,nbe6
      double precision pi,R,L,r1,r2,r3
      double precision theta,phi,chi,arg,darg,dx,slope
      double precision Xc1,Yc1,Xc2,Yc2,Xc3,Yc3
      double precision x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6
      double precision S1,S2,S3,S4,S5,S6,Stot
      double precision ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe)
      double precision nx(1:nbe),ny(1:nbe)
      common /const/ R,r1,r2,r3
*
      pi=4d0*datan(1d0)
*
      call set_parameters()
*     
      L=dsqrt(3d0)*R
      theta=dasin((r2-r1)/L)
      phi=dasin((r3-r2)/L)
      chi=dasin((r3-r1)/L)
*
      Xc1=-dsqrt(3d0)/2d0*R
      Yc1=-R/2d0
*
      Xc2=dsqrt(3d0)/2d0*R
      Yc2=-R/2d0
*
      Xc3=0d0
      Yc3=R
*
      x1=Xc1-r1*dcos(pi/6d0-chi)
      y1=Yc1+r1*dsin(pi/6d0-chi)
*
      x2=Xc1-r1*dsin(theta)
      y2=Yc1-r1*dcos(theta)
*
      x3=Xc2-r2*dsin(theta)
      y3=Yc2-r2*dcos(theta)
*
      x4=Xc2+r2*dcos(pi/6d0-phi)
      y4=Yc2+r2*dsin(pi/6d0-phi)
*
      x5=Xc3+r3*dcos(pi/6d0-phi)
      y5=Yc3+r3*dsin(pi/6d0-phi)
*
      x6=Xc3-r3*dcos(pi/6d0-chi)
      y6=Yc3+r3*dsin(pi/6d0-chi)
*      
      S1=(2d0*pi/3d0-theta-chi)*r1
      S2=dsqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))
      S3=(2d0*pi/3d0-phi+theta)*r2
      S4=dsqrt((x5-x4)*(x5-x4)+(y5-y4)*(y5-y4))
      S5=(2d0*pi/3d0+chi+phi)*r3
      S6=dsqrt((x1-x6)*(x1-x6)+(y1-y6)*(y1-y6))
*
      Stot=S1+S2+S3+S4+S5+S6
*
      nbe1=nint(nbe*S1/Stot)
      nbe2=nint(nbe*S2/Stot)
      nbe3=nint(nbe*S3/Stot)
      nbe4=nint(nbe*S4/Stot)
      nbe5=nint(nbe*S5/Stot)
      nbe6=nbe-(nbe1+nbe2+nbe3+nbe4+nbe5)
*
*****************************
*     Circular part S1      *
*****************************
      darg=(2d0*pi/3d0-theta-chi)/dble(nbe1)
      do n=1,nbe1
         arg=(5d0*pi/6d0+chi)+(dble(n)-0.5d0)*darg
         xl(n)=Xc1+r1*dcos(arg)
         yl(n)=Yc1+r1*dsin(arg)
         ds(n)=r1*darg
         kappa(n)=1d0/r1
         nx(n)=dcos(arg)
         ny(n)=dsin(arg)
      enddo
* 
*****************************
*     Linear part S2        *
*****************************
      dx=(x3-x2)/dble(nbe2)
      slope=(y3-y2)/(x3-x2)
      do n=1,nbe2
         nn=n+nbe1
         xl(nn)=x2+(dble(n)-0.5d0)*dx
         yl(nn)=slope*(xl(nn)-x2)+y2
         ds(nn)=dsqrt(1d0+slope*slope)*dx
         kappa(nn)=0d0
         nx(nn)=slope/dsqrt(1d0+slope*slope)
         ny(nn)=-1d0/dsqrt(1d0+slope*slope)
      enddo
*  
*****************************
*     Circular part S3      *
*****************************
      darg=(2d0*pi/3d0-phi+theta)/dble(nbe3)
      do n=1,nbe3
         nn=n+nbe2+nbe1
         arg=(3d0*pi/2d0-theta)+(dble(n)-0.5d0)*darg
         xl(nn)=Xc2+r2*dcos(arg)
         yl(nn)=Yc2+r2*dsin(arg)
         ds(nn)=r2*darg
         kappa(nn)=1d0/r2
         nx(nn)=dcos(arg)
         ny(nn)=dsin(arg)
      enddo
*
*****************************
*     Linear part S4        *
*****************************
      dx=(x5-x4)/dble(nbe4)
      slope=(y5-y4)/(x5-x4)
      do n=1,nbe4
         nn=n+nbe3+nbe2+nbe1
         xl(nn)=x4+(dble(n)-0.5d0)*dx
         yl(nn)=slope*(xl(nn)-x4)+y4
         ds(nn)=dsqrt(1d0+slope*slope)*(-dx)
         kappa(nn)=0d0
         nx(nn)=-slope/dsqrt(1d0+slope*slope)
         ny(nn)=1d0/dsqrt(1d0+slope*slope)
      enddo
*
*****************************
*     Circular part S5      *
*****************************
      darg=(2d0*pi/3d0+chi+phi)/dble(nbe5)
      do n=1,nbe5
         nn=n+nbe4+nbe3+nbe2+nbe1
         arg=(pi/6d0-phi)+(dble(n)-0.5d0)*darg
         xl(nn)=Xc3+r3*dcos(arg)
         yl(nn)=Yc3+r3*dsin(arg)
         ds(nn)=r3*darg
         kappa(nn)=1d0/r3
         nx(nn)=dcos(arg)
         ny(nn)=dsin(arg)
      enddo    
*
*****************************
*     Linear part S6        *
*****************************
      dx=(x1-x6)/dble(nbe6)
      slope=(y1-y6)/(x1-x6)
      do n=1,nbe6
         nn=n+nbe5+nbe4+nbe3+nbe2+nbe1
         xl(nn)=x6+(dble(n)-0.5d0)*dx
         yl(nn)=slope*(xl(nn)-x6)+y6
         ds(nn)=dsqrt(1d0+slope*slope)*(-dx)
         kappa(nn)=0d0
         nx(nn)=-slope/dsqrt(1d0+slope*slope)
         ny(nn)=1d0/dsqrt(1d0+slope*slope)
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
      double precision pi,R,r1,r2,r3,x,y,yy,rr1,rr2,rr3
      double precision L,theta,phi,chi,Xc1,Yc1,Xc2,Yc2,Xc3,Yc3
      double precision x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6
      common /const/ R,r1,r2,r3
*
      pi=4d0*datan(1d0)
*
      call set_parameters()
*
      L=dsqrt(3d0)*R
      theta=dasin((r2-r1)/L)
      phi=dasin((r3-r2)/L)
      chi=dasin((r3-r1)/L)
*
      Xc1=-dsqrt(3d0)/2d0*R
      Yc1=-R/2d0
*
      Xc2=dsqrt(3d0)/2d0*R
      Yc2=-R/2d0
*
      Xc3=0d0
      Yc3=R
*
      x1=Xc1-r1*dcos(pi/6d0-chi)
      y1=Yc1+r1*dsin(pi/6d0-chi)
*
      x2=Xc1-r1*dsin(theta)
      y2=Yc1-r1*dcos(theta)
*
      x3=Xc2-r2*dsin(theta)
      y3=Yc2-r2*dcos(theta)
*
      x4=Xc2+r2*dcos(pi/6d0-phi)
      y4=Yc2+r2*dsin(pi/6d0-phi)
*
      x5=Xc3+r3*dcos(pi/6d0-phi)
      y5=Yc3+r3*dsin(pi/6d0-phi)
*
      x6=Xc3-r3*dcos(pi/6d0-chi)
      y6=Yc3+r3*dsin(pi/6d0-chi)
*
C      write(*,*) x1,y1
C      write(*,*) x2,y2
C      write(*,*) x3,y3
C      write(*,*) x4,y4
C      write(*,*) x5,y5
C      write(*,*) x6,y6
C      stop
*
      flag=0
*
C     TEST L1
      yy=(y2-y1)/(x2-x1)*(x-x1)+y1
      if(yy.LE.y) then
C     TEST L2
         yy=(y3-y2)/(x3-x2)*(x-x2)+y2
         if(yy.LE.y) then
C     TEST L3
            yy=(y4-y3)/(x4-x3)*(x-x3)+y3
            if(yy.LE.y) then
C     TEST L4
               yy=(y5-y4)/(x5-x4)*(x-x4)+y4
               if(yy.GE.y) then
C     TEST L5
                  yy=(y6-y5)/(x6-x5)*(x-x5)+y5
                  if(yy.GE.y) then
C     TEST L6
                     yy=(y1-y6)/(x1-x6)*(x-x6)+y6
                     if(yy.GE.y) then
                        flag=1
                     endif
                  endif
               endif
            endif
         endif
      endif
*
      rr1=(x-Xc1)*(x-Xc1)+(y-Yc1)*(y-Yc1)
      rr1=dsqrt(rr1)
*
      rr2=(x-Xc2)*(x-Xc2)+(y-Yc2)*(y-Yc2)
      rr2=dsqrt(rr2)
*
      rr3=(x-Xc3)*(x-Xc3)+(y-Yc3)*(y-Yc3)
      rr3=dsqrt(rr3)
*
      if(rr1.LE.r1) then
         flag=1
      else if(rr2.LE.r2) then
         flag=1
      else if(rr3.LE.r3) then
         flag=1
      endif
*
      return
      end
