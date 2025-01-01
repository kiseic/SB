#!/bin/bash -f

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.

#
#   det.sh for the cavity with symmetry class SYM=2
#
#   Shell script for running the the executable file det.[TM/TE].NBE=...
#   for calculating the determinant value distribution for a given
#   resonance search domain in the complex wave number space.
#
#   See ocms/doc/users_guide.pdf for the details.
#
#===============================================================
#   P_A_R_A_M_E_T_E_R___S_E_T_T_I_N_G___S_T_A_R_T_S___H_E_R_E   
#===============================================================

# Executable file name (e.g., ./det.TM.NBE=50)
execfile=./../bin/det.TM.NBE=50
#det.TM.NBE=41￥paramsからでてくる￥という固有値のデータを計算するために必要なファイルをbinに入れる

# Refractive indices inside and outside the cavity
nin=2d0
nout=1d0

# Parity indices (a,b = -1(odd) or +1(even))
a=-1
b=1

# Resonance search domain:
#  (cx,cy) : center of the search domain
#  dwx,dwy : the half-widths of the search domain
#  dx,dy   : the grid spacings
cx=-0.02d0
cy=-0.005d0

dwx=0.01d0
dwy=0.01d0

 dx=0.001d0
 dy=0.001d0

#===============================================================
#   P_A_R_A_M_E_T_E_R___S_E_T_T_I_N_G___E_N_D_S___H_E_R_E
#===============================================================

# Output data file name (DO NOT TOUCH)
outputfile="dat.det.cx="$cx".cy="$cy."dx="$dx."dy="$dy

# Run the execfile (DO NOT TOUCH)
echo $nin $nout $a $b $cx $cy $dwx $dwy $dx $dy | $execfile > $outputfile
