#
# build-params.mk
#

#=====================
#   Cavity symmetry   
#=====================
SYM = 1 # The cavity is symmetric only with respect to the x-axis

#==================
#   Cavity shape
#==================
CAVITY = def_cavity_SYM1_D-shape.F
#CAVITY = def_cavity_SYM1_annular.F
#CAVITY = def_cavity_SYM1_cardioid.F

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
