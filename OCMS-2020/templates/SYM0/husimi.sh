#!/bin/bash -f

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.

#
#   husimi.sh for the cavity with symmetry class SYM=0
#
#   Shell script for running the executable file husimi.[TM/TE].NBE=...
#   for calculating the Husimi distribution for a given wave number.
#
#   See ocms/doc/users_guide.pdf for the details.
#
#===============================================================
#   P_A_R_A_M_E_T_E_R___S_E_T_T_I_N_G___S_T_A_R_T_S___H_E_R_E   
#===============================================================

# Executable file name (e.g., ../husimi.TM.NBE=50)
execfile=

# Parity indices (a,b = -1 or +1)
nin=
nout=

# Complex wave number of a resonant mode, k = kx + i ky
kx=
ky=

#===============================================================
#   P_A_R_A_M_E_T_E_R___S_E_T_T_I_N_G___E_N_D_S___H_E_R_E
#===============================================================

# Run the execfile (DO NOT TOUCH)
echo $nin $nout $kx $ky | $execfile
