c     Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
c     See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
c     Contact ocms@telecognix.com for further information.

C===============================================
C     Definition of the dispersion relation
C===============================================
#if defined(CUSTOMIZED)
C
C     Non-standard dispersion relation
C
C     (Ex) Schroedinger-Bloch model case
C          nk = 2*k + 1             (inside the cavity)
C          nk = 2*k + (nout/nin)^2  (outside the cavity)
C
      nk_in=cdsqrt(2d0*k+1d0)     ! Inside the cavity
      nk_out=cdsqrt(2d0*k+
     &     dble(nout)*dble(nout)
     &     /dble(nin)/dble(nin))  ! Outside the cavity
*
#else
C
C     Standard dispersion relation
C
C     nk = n * k
C
      nk_in=dble(nin)*k   ! Inside the cavity
      nk_out=dble(nout)*k ! Outside the cavity
*
#endif
