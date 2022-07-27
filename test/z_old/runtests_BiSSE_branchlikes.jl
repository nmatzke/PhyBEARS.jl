using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DataFrames  # for DataFrame
using DifferentialEquations
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA


# List each PhyBEARS code file prefix here
using PhyBEARS.Example
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.TrUtils
using PhyBEARS.SSEs

# Sourced out of /GitHub/PhyBEARS.jl/test/BiSSE_branchlikes_w_BD_v4_WORKING_n1.R

# results from R:
"""
> R_result_branch_lnL = -3.128581
> R_result_total_lnL = -4.937608
> R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -6.579522

LnLs
# -4.937608 -3.128581

"""

include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.222222222, deathRate=0.1, d_val=0.0, e_val=0.0, a_val=0.0, j_val=0.0)
numstates = 2
n = 2
inputs = ModelLikes.setup_MuSSE_biogeo(numstates, tr; root_age_mult=1.5, in_params=in_params)

trdf = inputs.trdf
root_age = maximum(trdf[!, :node_age])
Es_interpolator = inputs.p_Ds_v5.sol_Es_v5;
prtQi(inputs)
prtCi(inputs)

# Checking Es
Julia_result_Es = Es_interpolator(1.0)
R_Es = R_result_EsDs[(2:(n+1))]
@test all(round.(Julia_result_Es; digits=4) .== round.(R_Es; digits=4))

# Now we do a downpass?
# Just kinda translating the R code over

numsteps = 100000
numstates = 2
