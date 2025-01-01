#
# build-params.mk
#

#=====================
#   Cavity symmetry   
#=====================
SYM = 4 # The cavity is symmetric with respect to x=y and the y-axis

#==================
#   Cavity shape
#==================
CAVITY = def_cavity_SYM4_Sinai.F

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
