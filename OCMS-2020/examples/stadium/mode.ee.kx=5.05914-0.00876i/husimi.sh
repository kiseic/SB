#!/bin/bash -f

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.

#
#   husimi.sh for the cavity with symmetry class SYM=2
#
#   Shell script for running the executable file husimi.[TM/TE].nmax=...
#   for calculating the Husimi distribution for a given wave number.
#
#   See ocms/doc/users_guide.pdf for the details.
#
#===============================================================
#   P_A_R_A_M_E_T_E_R___S_E_T_T_I_N_G___S_T_A_R_T_S___H_E_R_E   
#===============================================================

# Executable file name (e.g., ../husimi.TM.nmax=50)
execfile=../bin/husimi.TM.nmax=50

# Parity indices (a,b = -1 or +1)
nin=3.3d0
nout=1d0

# Parity indices (a,b = -1(odd) or +1(even))
a=1
b=1

# Complex wave number of a resonant mode, k = kx + i ky
kx=5.05914d0
ky=-0.00876d0

#===============================================================
#   P_A_R_A_M_E_T_E_R___S_E_T_T_I_N_G___E_N_D_S___H_E_R_E
#===============================================================

# Run the execfile (DO NOT TOUCH)
echo $nin $nout $a $b $kx $ky | $execfile
