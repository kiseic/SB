c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

#if (SYM == 4 || SYM ==2)
      subroutine def_matrix(nbe,nk,a,b,xl,yl,nx,ny,sign,ds,kappa,P,Q)
#elif SYM == 1
      subroutine def_matrix(nbe,nk,b,xl,yl,nx,ny,sign,ds,kappa,P,Q)
#elif SYM ==0
      subroutine def_matrix(nbe,nk,xl,yl,nx,ny,sign,ds,kappa,P,Q)
#else
      write(*,*) "Error in def_matrix.F: SYM value is not appropriate."
      stop
#endif
*
      implicit NONE
*
      integer nbe,i,j,sign,a,b
      double precision xl(1:nbe),yl(1:nbe),nx(1:nbe),ny(1:nbe)
      double precision ds(1:nbe),kappa(1:nbe)
      double precision pi,gamma,x1,x2,y1,y2,vx,vy,ad,bd
      double complex imu,nk,dum,w1,w2,w3,w4,w5,w6,w7,w8
      double complex P(1:nbe,1:nbe),Q(1:nbe,1:nbe)
*    
      pi=4d0*datan(1d0)
      imu=dcmplx(0d0,1d0)
      gamma=0.5772156649d0
*
#if SYM == 4
C
C     When the cavity has mirror symmetries with respect to x=y and the y-axis.
C
      ad=dble(a) ! Parity index for the mirror symmetry w.r.t. the y-axis.
      bd=dble(b) ! Parity index for the mirror symmetry w.r.t. x=y.
*
      do i=1,nbe
         do j=1,nbe
*     
            x2=xl(i)
            y2=yl(i)
            x1=xl(j)
            y1=yl(j)
            vx=nx(j)
            vy=ny(j)
*     
            if(i.NE.j) then
               call gfn(nk,x1,y1,x2,y2,w1)
               call gfn(nk,y1,x1,x2,y2,w2)
               call gfn(nk,-y1,x1,x2,y2,w3)
               call gfn(nk,-x1,y1,x2,y2,w4)
               call gfn(nk,-x1,-y1,x2,y2,w5)
               call gfn(nk,-y1,-x1,x2,y2,w6)
               call gfn(nk,y1,-x1,x2,y2,w7)
               call gfn(nk,x1,-y1,x2,y2,w8)
               P(i,j) = -2d0*( w1 + bd*w2 + ad*bd*w3 + ad*w4 
     &              + w5 + bd*w6 + ad*bd*w7 + ad*w8)*ds(j)
*
               call dgfn(nk,x1,y1,x2,y2,vx,vy,w1)
               call dgfn(nk,y1,x1,x2,y2,vy,vx,w2)
               call dgfn(nk,-y1,x1,x2,y2,-vy,vx,w3)
               call dgfn(nk,-x1,y1,x2,y2,-vx,vy,w4)
               call dgfn(nk,-x1,-y1,x2,y2,-vx,-vy,w5)
               call dgfn(nk,-y1,-x1,x2,y2,-vy,-vx,w6)
               call dgfn(nk,y1,-x1,x2,y2,vy,-vx,w7)
               call dgfn(nk,x1,-y1,x2,y2,vx,-vy,w8)
               Q(i,j) = 2d0*( w1 + bd*w2 + ad*bd*w3 + ad*w4 
     &              + w5 + bd*w6 + ad*bd*w7 + ad*w8)*ds(j)
               
            else
               dum=cdlog(nk*ds(j)/4d0)
               w1=ds(j)/pi*(1d0-dum+imu*pi/2d0-gamma)
               call gfn(nk,y1,x1,x2,y2,w2)
               call gfn(nk,-y1,x1,x2,y2,w3)
               call gfn(nk,-x1,y1,x2,y2,w4)
               call gfn(nk,-x1,-y1,x2,y2,w5)
               call gfn(nk,-y1,-x1,x2,y2,w6)
               call gfn(nk,y1,-x1,x2,y2,w7)
               call gfn(nk,x1,-y1,x2,y2,w8)
               P(i,j) = w1 - 2d0*( bd*w2 + ad*bd*w3 + ad*w4 
     &              + w5 + bd*w6 + ad*bd*w7 + ad*w8 )*ds(j)
*
               w1=kappa(j)/4d0/pi
               call dgfn(nk,y1,x1,x2,y2,vy,vx,w2)
               call dgfn(nk,-y1,x1,x2,y2,-vy,vx,w3)
               call dgfn(nk,-x1,y1,x2,y2,-vx,vy,w4)
               call dgfn(nk,-x1,-y1,x2,y2,-vx,-vy,w5)
               call dgfn(nk,-y1,-x1,x2,y2,-vy,-vx,w6)
               call dgfn(nk,y1,-x1,x2,y2,vy,-vx,w7)
               call dgfn(nk,x1,-y1,x2,y2,vx,-vy,w8)
               Q(i,j)=-dble(sign)+2d0*(w1 + bd*w2 + ad*bd*w3 + ad*w4 
     &              + w5 + bd*w6 + ad*bd*w7 + ad*w8)*ds(j)
            endif
*
         enddo
      enddo
*
#elif SYM == 2
C
C     When the cavity has mirror symmetries with respect to both x- and
C     y-axis.
C
      ad=dble(a) ! Parity index for the mirror symmetry w.r.t. the y-axis.
      bd=dble(b) ! Parity index for the mirror symmetry w.r.t. the x-axis.
*
      do i=1,nbe
         do j=1,nbe
*     
            x2=xl(i)
            y2=yl(i)
            x1=xl(j)
            y1=yl(j)
            vx=nx(j)
            vy=ny(j)
*     
            if(i.NE.j) then
               call gfn(nk,x1,y1,x2,y2,w1)
               call gfn(nk,-x1,y1,x2,y2,w2)
               call gfn(nk,-x1,-y1,x2,y2,w3)
               call gfn(nk,x1,-y1,x2,y2,w4)
               P(i,j) = -2d0*( w1 + ad*w2 + ad*bd*w3 + bd*w4 )*ds(j)
*
               call dgfn(nk,x1,y1,x2,y2,vx,vy,w1)
               call dgfn(nk,-x1,y1,x2,y2,-vx,vy,w2)
               call dgfn(nk,-x1,-y1,x2,y2,-vx,-vy,w3)
               call dgfn(nk,x1,-y1,x2,y2,vx,-vy,w4)
               Q(i,j) = 2d0*( w1 + ad*w2 + ad*bd*w3 + bd*w4 )*ds(j)
            else
               dum=cdlog(nk*ds(j)/4d0)
               w1=ds(j)/pi*(1d0-dum+imu*pi/2d0-gamma)
               call gfn(nk,-x1,y1,x2,y2,w2)
               call gfn(nk,-x1,-y1,x2,y2,w3)
               call gfn(nk,x1,-y1,x2,y2,w4)
               P(i,j) = w1 - 2d0*( ad*w2 + ad*bd*w3 + bd*w4 )*ds(j)
*
               w1=kappa(j)/4d0/pi
               call dgfn(nk,-x1,y1,x2,y2,-vx,vy,w2)
               call dgfn(nk,-x1,-y1,x2,y2,-vx,-vy,w3)
               call dgfn(nk,x1,-y1,x2,y2,vx,-vy,w4)
               Q(i,j)=-dble(sign)+2d0*(w1+ad*w2+ad*bd*w3+bd*w4)*ds(j)
            endif
*
         enddo
      enddo
*
#elif SYM == 1
C
C     When the cavity has a mirror symmetry with respect to the x-axis.
C
*
      bd=dble(b) ! Parity index for the mirror symmetry w.r.t. the x-axis.
*
      do i=1,nbe
         do j=1,nbe
*     
            x2=xl(i)
            y2=yl(i)
            x1=xl(j)
            y1=yl(j)
            vx=nx(j)
            vy=ny(j)
*     
            if(i.NE.j) then
               call gfn(nk,x1,y1,x2,y2,w1)
               call gfn(nk,x1,-y1,x2,y2,w4)
               P(i,j) = -2d0*( w1 + bd*w4 )*ds(j)
*
               call dgfn(nk,x1,y1,x2,y2,vx,vy,w1)
               call dgfn(nk,x1,-y1,x2,y2,vx,-vy,w4)
               Q(i,j) = 2d0*( w1 + bd*w4 )*ds(j)
            else
               dum=cdlog(nk*ds(j)/4d0)
               w1=ds(j)/pi*(1d0-dum+imu*pi/2d0-gamma)
               call gfn(nk,x1,-y1,x2,y2,w4)
               P(i,j) = w1 - 2d0*( bd*w4 )*ds(j)
*
               w1=kappa(j)/4d0/pi
               call dgfn(nk,x1,-y1,x2,y2,vx,-vy,w4)
               Q(i,j)=-dble(sign)+2d0*(w1+bd*w4)*ds(j)
            endif
*
         enddo
      enddo
*
#elif SYM == 0
C
C     When the cavity has no mirror symmetries.
C
*
      do i=1,nbe
         do j=1,nbe
*     
            x2=xl(i)
            y2=yl(i)
            x1=xl(j)
            y1=yl(j)
            vx=nx(j)
            vy=ny(j)
*     
            if(i.NE.j) then
               call gfn(nk,x1,y1,x2,y2,w1)
               P(i,j) = -2d0*w1*ds(j)
*
               call dgfn(nk,x1,y1,x2,y2,vx,vy,w1)
               Q(i,j) = 2d0*w1*ds(j)
            else
               dum=cdlog(nk*ds(j)/4d0)
               w1=ds(j)/pi*(1d0-dum+imu*pi/2d0-gamma)
               P(i,j) = w1
*
               w1=kappa(j)/4d0/pi
               Q(i,j)=-dble(sign)+2d0*w1*ds(j)
            endif
*
         enddo
      enddo   
*
#endif
      return
      end
