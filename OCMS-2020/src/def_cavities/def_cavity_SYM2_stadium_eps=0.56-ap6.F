c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine set_parameters()
*
      implicit NONE
*     
      double precision pi
      double precision eps,area0,size_param,R,L
      common /const/ R,L
*  
      pi=4d0*datan(1d0)
*
      eps=0.56d0 !eps=2L/2R=L/R
      area0=4d0+pi
      size_param=sqrt(area0/(4d0*eps+pi))
*
      R=(size_param) !Radius of circular part
      L=R*eps !Half the length of the straight part
      !write(*,*) "L/R=",L/R
      !write(*,*) "L=",L
      !write(*,*) "area=",4d0*L*R+pi*R**2
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
      integer n,nbe,Nr,Nl
      double precision pi,R,L,ds(1:nbe),kappa(1:nbe)
      double precision xl(1:nbe),yl(1:nbe)
      double precision nx(1:nbe),ny(1:nbe)
      double precision phi,dphi,cs1,sn1,dx,arc
      common /const/ R,L
*
      pi=4d0*datan(1d0)
*     
      call set_parameters()
*
      arc=pi/2d0*R
c      Nr=ceiling(L/(arc+L)*dble(nbe))
      Nr=aint(arc/(arc+L)*dble(nbe))
      Nl=nbe-Nr
*
*********************************************
*     Circular segment
*********************************************
      dphi=0.5d0*pi/dble(Nr)
*
      do n=1,Nr
         phi = dphi*(dble(n)-0.5d0)
*
         cs1=dcos(phi)
         sn1=dsin(phi)
*
         xl(n)=L+R*cs1
         yl(n)=R*sn1
*
         ds(n)=R*dphi
         nx(n)=cs1
         ny(n)=sn1
         kappa(n)=1d0/R
*
      enddo
************************************************
*     Line segment
************************************************
      dx=L/dble(Nl)
*
      do n=Nr+1,nbe
         xl(n)=L-(dble(n-Nr)-0.5d0)*dx
         yl(n)=R
*
         ds(n)=dx
         nx(n)=0d0
         ny(n)=1d0
         kappa(n)=0d0
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
      double precision R,L,x,y
      double precision rr,ll,dum1,dum2,dx,dy
      common /const/ R,L
*
      call set_parameters()
*
      rr=R
      ll=L
*
      dx=x-ll
      dy=y
      dum1=dsqrt(dx*dx + dy*dy)
*
      dx=x+ll
      dy=y
      dum2=dsqrt(dx*dx + dy*dy)
*
      if(dabs(x).LT.ll.AND.dabs(y).LT.rr) then
         flag=1
      else if(x.GT.ll.AND.dum1.LT.rr) then
         flag=1
      else if(x.LT.(-ll).AND.dum2.LT.rr) then
         flag=1
      else
         flag=0
      endif
*
      return
      end
*
*******************************************************************************
