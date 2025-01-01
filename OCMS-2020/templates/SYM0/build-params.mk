#
# build-params.mk
#

#=====================
#   Cavity symmetry   
#=====================
SYM = 0 # The cavity does not have any mirror symmetry

#==================
#   Cavity shape
#==================
CAVITY = def_cavity_SYM0_Sinai.F
#CAVITY = def_cavity_SYM0_asym_cut_disk.F
#CAVITY = def_cavity_SYM0_asym_limacon.F
#CAVITY = def_cavity_SYM0_rounded_triangle.F

#===========================
#   Cavity size parameter
#===========================
SIZE_PARAM = 1d0

#=============================
#   Polarization (TM or TE)
#=============================
PLZ = TM
#PLZ = TE

#=====================================
#   The number of boundary elements
#=====================================
NBE = 50

#======================================================
#   Grid points for plotting the Husimi distribution  
#   (parameters for main.husimi.F)
#======================================================
IXMAX = 200
IYMAX = 200
