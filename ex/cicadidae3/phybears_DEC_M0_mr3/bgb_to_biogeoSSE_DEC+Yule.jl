using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames						# for DataFrame()
using DelimitedFiles				# for readdlm()
using NLopt									# seems to be the best gradient-free, box-constrained								
using RCall									# To call R code from Julia

# List each PhyBEARS code file prefix here
using PhyloBits.TrUtils			# for e.g. numstxt_to_df()
using PhyloBits.TreeTable
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.ModelLikes # e.g. setup_DEC_SSE2
using PhyBEARS.Uppass

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/ex/cicadidae3/phybears_DEC_M0_mr3/")
include("/GitHub/PhyBEARS.jl/ex/cicadidae3/phybears_DEC_M0_mr3/phybears_M0_mr3_DEC_v2.jl")
"""

setwd("/GitHub/PhyBEARS.jl/ex/cicadidae3/phybears_DEC_M0_mr3/")

# Input geography
lgdata_fn = "/GitHub/PhyBEARS.jl/ex/cicadidae3/phybears_DEC_M0_mr3/geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Input tree
trfn = "/GitHub/PhyBEARS.jl/ex/cicadidae3/phybears_DEC_M0_mr3/tree.newick"
tr = readTopology(trfn)
trdf = prt(tr)

#######################################################
# BioGeoBEARS results: DEC model
#######################################################
# BioGeoBEARS DEC on Cicadidae M0_unconstrained ancstates: global optim, 3 areas max. 
# d=0.001; e=7e−04; j=0; LnL=−310.93

# Basic tree info
numTips = sum(trdf.nodeType .== "tip")
numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")
ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
birthRate = yuleBirthRate = (numInternal-1) / ttl_tree_length

bmo = construct_BioGeoBEARS_model_object()
bmo.est[bmo.rownames .== "birthRate"] .= birthRate
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.0010
bmo.est[bmo.rownames .== "e"] .= 0.0007
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
numareas = 10
n = numstates_from_numareas(10,3,true)
#n = 176            # 10 areas, maxareas 3, 176 states

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
root_age_mult=1.5; max_range_size=3; include_null_range=true; max_range_size=NaN
max_range_size = 3 # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;


numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, sol_Es_v5=sol_Es_v5);



# Check the interpolator
p_Ds_v5.sol_Es_v5(1.0)
Es_interpolator(1.0)

# Do downpass - slow and fast calculation
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
# 8.8 seconds, lnLs match BioGeoBEARS
# (8.84, 16, -993.5629081577177, -9.03925286266909, -1002.6021610203868, -310.97299243972896)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
# (0.495, 16, -993.5629081577177, -9.03925286266909, -1002.6021610203868, -310.97299243972896)

#######################################################
# Showing the conversion from BioGeoBEARS lnL to a ClaSSE lnL (yule process, equal root state probs)
#######################################################
d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]

# Use the standard birth-death likelihood for just a phylogeny with a birth and death rate
ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
yuleBirthRate = (numInternal-1) / ttl_tree_length
yuleDeathRate = 0.0					# Yule process has 0 extinction
bd = bd_liks_trdf(trdf, yuleBirthRate, yuleDeathRate)
bd_lnL_noTopo = bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times

equal_root_prob2 = log(1/(numstates-include_null_range))
bgb_root_lnL = log(sum(d_root_orig)) + 1.0


# BioGeoBEARS lnL from Julia sum branch likelihoods
bgb_lnL  # from Julia result res
Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)
Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo
Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)

# Julia full lnL
Julia_total_lnLs1
Julia_sum_lq + rootstates_lnL
root_stateprobs = d_root_orig/sum(d_root_orig);
Julia_sum_lq + log(sum(root_stateprobs .* d_root_orig))

# Convert BioGeoBEARS bgb_lnl to total lnL
bgb_lnL - log(sum(d_root_orig)) - log(1/yuleBirthRate) + (bd.lnL - bd.lnl_topology) + log(sum(root_stateprobs .* d_root_orig))

