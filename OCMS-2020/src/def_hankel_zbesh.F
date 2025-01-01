c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

*
      subroutine hankel_1(nu,z,w)
*
      integer n
      parameter (n=2)
*
      integer nu,nz,ierr,m,kode
      double precision zr,zi,cyr(1:n),cyi(1:n),fnu
      double complex imu,z,w
*
      imu=dcmplx(0d0,1d0)
*
      kode=1
      m=1
      fnu=0
*
      zr=dble(z)
      zi=imag(z)
*      
      call ZBESH(zr,zi,fnu,kode,m,n,cyr,cyi,nz,ierr)
*
      if(nu.EQ.0) then
         w = cyr(1) + imu*cyi(1)
      else if(nu.EQ.1) then
         w = cyr(2) + imu*cyi(2)
      else
         write(*,*) 'Error in hankel_1().'
         stop
      endif
*
      return
      end
*
