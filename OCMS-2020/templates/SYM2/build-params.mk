#
# build-params.mk
#

#=====================
#   Cavity symmetry   
#=====================
SYM = 2 # The cavity is symmetric with respect to both x- and y-axis

#==================
#   Cavity shape
#==================
CAVITY = def_cavity_SYM2_D2_deformed_circle.F
#CAVITY = def_cavity_SYM2_Sinai.F
#CAVITY = def_cavity_SYM2_ellipse.F
#CAVITY = def_cavity_SYM2_flattened_quadrupole.F
#CAVITY = def_cavity_SYM2_stadium.F

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
