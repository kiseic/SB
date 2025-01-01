#!/bin/bash -f

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.

#
#   wfunc.sh for the cavity with symmetry class SYM=2
#
#   Shell script for running the executable file wfunc.[TM/TE].NBE=...
#   for calculating the wave function for a given wave number.
#
#   See ocms/doc/users_guide.pdf for the details.
#
#===============================================================
#   P_A_R_A_M_E_T_E_R___S_E_T_T_I_N_G___S_T_A_R_T_S___H_E_R_E   
#===============================================================

# Executable file name (e.g., ../wfunc.TM.NBE=50)
execfile=./../bin/wfunc.TM.NBE=50

# Refractive indices inside and outside the cavity
nin=2d0
nout=1d0

# Parity indices (a,b = -1(odd) or +1(even))
a=-1
b=1

# Complex wave number of a resonant mode, k = kx + i ky
kx=-0.024d0
ky=-0.005d0

# x and y ranges for outputting the data
xmin=-20d0
xmax=20d0
ymin=-10d0
ymax=10d0

# x and y coordinates grid points
ixmax=300
iymax=300

#===============================================================
#   P_A_R_A_M_E_T_E_R___S_E_T_T_I_N_G___E_N_D_S___H_E_R_E
#===============================================================

# Run the execfile (DO NOT TOUCH)
echo $nin $nout $a $b $kx $ky $xmin $xmax $ymin $ymax $ixmax $iymax | $execfile
