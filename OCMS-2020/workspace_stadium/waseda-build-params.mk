#
# waseda-build-params.mk
#

#=======================
# Dispersion relation
#=======================
# When CUSTOMIZED is selected, the dispersion relation must be defined in
# defined in def_nk.F.
# STANDARD is the default and fallback mode.
#DISPERSION = STANDARD # Standard dispersion relation (i.e., omega=ck)
DISPERSION = CUSTOMIZED # Non-standard dispersion relation (defined in def_nk.F)
