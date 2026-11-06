#######################################################
# 2026-11-05
#
# Notes on lnL calculations for data, tree, and data+tree (SSE models)
#
# in R and Julia
#######################################################

# Location:

# ~/GitHub/PhyBEARS/ex/tutorials/_00_intro_to_lnLs_SSE_models_v1.R
#


# 1. Towards SSEs in R:
# ~/GitHub/PhyBEARS.jl/test/BiSSE_branchlikes_w_BD_v4_WORKING_n1.R

# See: 
# * Nee et al. 1994 for BD likelihood
# * Maddison et al. 2007 for BiSSE likelihood

# 2. This shows a match between BD likelihood and 
#    BiSSE lnL under matching parameters:
# ~/GitHub/PhyBEARS.jl/test/BiSSE_branchlikes_w_BD_v5_MATCHES_BD_YULE_n1.R

# Basically: to get a BiSSE model to match a BirthDeath lnL
# birthRate - has to be the same for all states
# deathRate - has to be 0 (or same for all states)
# q01 and q10 transition rates - 0
# tip data - all the same state (so that no character changes have to occur)




#######################################################
# Play with simulator in the castor R package
# (with some modifications by Nick)
#######################################################
# Originally designed for microbial macroevolution problems by
# crazy genius - Louca & co
# It lives here:
# ~/GitHub/PhyBEARS.jl/simulator/castor_1.7.2.000004.tar.gz


# Example scripts:
# ~/GitHub/PhyBEARS.jl/simulator/examples/tdsse_example_plus_converter_v1.R
# ~/GitHub/PhyBEARS.jl/simulator/examples/tdsse_example_plus_converter_v4.R

