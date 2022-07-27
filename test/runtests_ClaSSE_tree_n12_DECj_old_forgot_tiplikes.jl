#######################################################
# This script shows:
# 1. The calculation of the BioGeoBEARS DEC and DEC+J
#    log-likelihoods, using numeric integration with
#    ClaSSE in Julia. (This requires taking the ClaSSE
#    branch log-likelihood, i.e. lq, and subtracting
#    the birthdeath-process lnL (BD lnL, but without
#    the logfactorial term), as well as the root state
#    frequencies.
#
# 2. Maximum likelihood, in a constrained ClaSSE model,
#    achieving the same ML parameters and lnL as 
#    BioGeoBEARS, for DEC and DEC+J.
#######################################################

using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames
using Optim                 # for e.g. L-BFGS-B Maximum Likelihood optimization,optimize

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DataFrames  # for DataFrame
using DifferentialEquations
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA


# List each PhyBEARS code file prefix here
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreeTable
using PhyBEARS.TreePass
using PhyBEARS.TrUtils
using PhyBEARS.SSEs

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_ClaSSE_tree_n12_DECj.jl")
"""
# 
# """
# # Run with:
# include("/GitHub/PhyBEARS.jl/test/runtests_ClaSSE_tree_n12_DECj.jl")
# """
# 
# @testset "Example" begin
# 	@test hello("runtests_ClaSSE_tree_n12_DECj.jl") == "runtests_ClaSSE_tree_n12_DECj.jl"
# #	@test domath(2.0) ≈ 7.0
# end
# 
# 
# #######################################################
# # Do a bunch of tests of the SSE calculation of 
# # Ds, Es, and likelihoods, on
# # branches, nodes, and trees,
# # under a variety of simple and more complex models
# #######################################################
# 
#@testset "runtests_ClaSSE_tree_n9_DECj.jl" begin
# 
#######################################################
# DEMONSTRATES MATCHING BETWEEN DIVERSITREE, BIOGEOBEARS, AND JULIA
# ON HAWAIIAN PSYCHOTRIA, 16-STATE DEC MODEL
#
# Run with:
# source("/GitHub/PhyBEARS.jl/Rsrc/compare_BGB_diversitree_DEC+J_v1.R")
# Truth:
DEC_lnL = -34.54313
DEC_R_result_branch_lnL = -67.6295
DEC_R_result_total_LnLs1 = -72.60212
DEC_R_result_total_LnLs1t = -71.48986
DEC_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -120.1545

# DEC+J
DECj_lnL = -20.94759
R_result_branch_lnL = -55.37332
R_result_total_LnLs1 = -58.83758
R_result_total_LnLs1t = -57.72533
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -96.34151

#######################################################


include("/GitHub/PhyBEARS.jl/src/TreePass.jl")
import .TreePass

# Repeat calculation in Julia
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes

# Read in the island areas using Parsers.jl (v10) rather than just manual input (v9)
include("/GitHub/PhyBEARS.jl/notes/Parsers.jl")
import .Parsers

# Island numbers (KOMH = 1234) in Rnodenums order:
#island_nums = [3, 3, 2, 2, 3, 3, 2, 1, 1, 3, 4, 2, 1, 1, 1, 1, 1, 1, 2]
lgdata_fn = "/GitHub/PhyBEARS.jl/Rsrc/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Psychotria tree
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
#in_params = (birthRate=0.3288164, deathRate=0.0, d_val=1e-12, e_val=1e-12, a_val=0.0, j_val=0.1142057)
bmo = construct_BioGeoBEARS_model_object()
bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0

numareas = 4
n = 16            # 4 areas, 16 states

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
#inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, in_params=in_params)
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo)
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs
inputs = Parsers.tipranges_to_tiplikes(inputs, geog_df);
res = inputs.res

prtQi(inputs)
inputs.setup.maxent01

# Update the tip likelihoods, with the geography data
inputs.res.likes_at_each_nodeIndex_branchTop
size(inputs.res.likes_at_each_nodeIndex_branchTop)
numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])

inputs.setup.observed_statenums
inputs.res.likes_at_each_nodeIndex_branchTop

inputs.setup.observed_statenums

inputs.res.likes_at_each_nodeIndex_branchTop


inputs.res.likes_at_each_nodeIndex_branchTop
res = inputs.res

res = inputs.res
trdf = inputs.trdf
p_Ds_v5 = inputs.p_Ds_v5
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

prtQi(inputs)
prtCi(inputs)
inputs.p_Ds_v5.params.mu_vals
p_Ds_v5.sol_Es_v5(1.0)

Es_interpolator(1.0)


# Parameters

# Do downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lqA, rootstates_lnLA, Julia_total_lnLs1A, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


bmo = construct_BioGeoBEARS_model_object()
bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 1.0e-12
bmo.est[bmo.rownames .== "e"] .= 1.0e-12
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.1142057

inputs = ModelLikes.setup_DEC_SSE2(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo)
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs
inputs = Parsers.tipranges_to_tiplikes(inputs, geog_df);
res = inputs.res

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

(total_calctime_in_sec, iteration_number, Julia_sum_lqA, rootstates_lnLA, Julia_total_lnLs1A, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)







in_params = (birthRate=0.3288164, deathRate=0.0, d_val=1e-12, e_val=1e-12, a_val=0.0, j_val=0.1142057)
numareas = 4
n = 16            # 4 areas, 16 states

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, in_params=in_params)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

(total_calctime_in_sec, iteration_number, Julia_sum_lqA, rootstates_lnLA, Julia_total_lnLs1A, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

prtCp(xy) .== prtCp(p_Ds_v5)


Rnames(res)
res.likes_at_each_nodeIndex_branchTop
res.normlikes_at_each_nodeIndex_branchTop
res.likes_at_each_nodeIndex_branchBot
res.normlikes_at_each_nodeIndex_branchBot

sum.(res.likes_at_each_nodeIndex_branchTop)
log.(sum.(res.likes_at_each_nodeIndex_branchTop))
sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))
Julia_sum_lq_old = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
nonroot_nodes = get_nonrootnodes(tr)
sum_likes_internal_branch_tops = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
Julia_sum_lq = Julia_sum_lq_old + sum_likes_internal_branch_tops


##############################################
# Consider the pure-birth log-likelihoods
# The Yule-process ML birthrate is just (# internal nodes - 1)/total_tree_length
# Get basic tree info
numTips = sum(trdf.nodeType .== "tip")
numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")


d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]

ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
yuleBirthRate = (numInternal-1) / ttl_tree_length
yuleDeathRate = 0.0					# Yule process has 0 extinction
bd = bd_liks_trdf(trdf, yuleBirthRate, yuleDeathRate)
bd_lnL_noTopo = bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times

# Convert to BioGeoBEARS lnL under Yule process assumption
# Check if the first state/geographic range is null
#if res.inputs.setup.states_list[1] == []
#	include_null_range = true
#end
include_null_range = true
numstates = length(res.normlikes_at_each_nodeIndex_branchTop[1])
equal_root_prob2 = log(1/(numstates-include_null_range)) 
bgb_root_lnL = log(sum(d_root_orig)) + 1.0

# res5t match
res5t = Julia_sum_lq + equal_root_prob2 + bgb_root_lnL - (1-log(1/yuleBirthRate))
# ...matches these in R: compare_BGB_diversitree_DEC_v1.R
# bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/(birthRate))
# bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 
# bgb1 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 + bgb_root_lnL
# (bgb1 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2
bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2
bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)
bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo
# -20.947589078592664
##############################################


@test abs(DECj_lnL - bgb_lnL) < 0.01




# Does the total of the branch log-likelihoods (lq) match?
@test round(R_result_branch_lnL; digits=4) == round(Julia_sum_lq; digits=4)

# Add the root probabilities

# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
root_stateprobs = d_root_orig/sum(d_root_orig)
rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL


# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE (these seem to be the defaults)
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
root_stateprobs = d_root_orig/sum(d_root_orig)
lambda = birthRate
e_root = Es_interpolator(root_age)


#d_root = d_root_orig ./ sum(root_stateprobs .* inputs.p_Ds_v5.params.Cijk_vals .* (1 .- e_root).^2)
# diversitree::rootfunc.classe
# lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
sum_of_lambdas = collect(repeat([0.0], n))
for i in 1:n
	sum_of_lambdas[i] = sum(inputs.p_Ds_v5.params.Cijk_vals[inputs.p_Ds_v5.p_TFs.Ci_eq_i[i]])
end
sum_of_lambdas
d_root = d_root_orig ./ sum(root_stateprobs .* sum_of_lambdas .* (1 .- e_root).^2)
rootstates_lnL = log(sum(root_stateprobs .* d_root))
# The above all works out to [0,1] for the Yule model with q01=q02=0.0

Julia_total_lnLs1t = Julia_sum_lq + rootstates_lnL

# Does the total lnL match R?
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
@test round(R_result_total_LnLs1; digits=4) == round(Julia_total_lnLs1; digits=4)

# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
@test round(R_result_total_LnLs1t; digits=4) == round(Julia_total_lnLs1t; digits=4)

# Does the total of branch likelihoods (lq) + node likelihoods match R?
# 
# R: R_result_sum_log_computed_likelihoods_at_each_node_x_lambda 
#    = sum(log(computed_likelihoods_at_each_node_x_lambda))
res.likes_at_each_nodeIndex_branchTop
log.(sum.(res.likes_at_each_nodeIndex_branchTop))

sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))

Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=4) == round(R_sum_lq_nodes; digits=4)








# The standard diversitree lnL calculation sums:
# 1. log-likelihoods at branch-bottoms (lq)
# 2. root-state likelihoods (e.g. rootstates_lnL)

# We can get this in R also:
# 1. sum(lq)
# 2. res1t = bisse_2areas(pars=bisse_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
#    (is sum(lq) + sum(rootstates_lnL)

# Weirdly, it seems like this R calculation extracts the total likelihood at the 
# branch bottoms (lq), but not the total likelihood at the nodes, after 
# the normlikes from above branches are combined at the node, and multiplied by
# the birthRate.
#
# This suggests another possible likelihood:
#
# 1. Julia: sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
#
# ...which matches...
#
# 2. sum(log(computed_likelihoods_at_each_node_x_lambda))
#
# (And which might be correct!  PhyBEARS will produce both.)
#




print("\nDifferences between Julia and R lnLs for\n/GitHub/PhyBEARS.jl/Rsrc/_compare_ClaSSE_calcs_v3_compare2julia.R\n calculation:\n")
print("R_result_branch_lnL (lq) - Julia_sum_lq: ")
print(R_result_branch_lnL - Julia_sum_lq)
print("\n")
print("R_result_total_LnLs1 (lq) - Julia_total_lnLs1: ")
print(R_result_total_LnLs1 - Julia_total_lnLs1)
print("\n")
print("R_result_total_LnLs1t (lq) - Julia_total_lnLs1t: ")
print(R_result_total_LnLs1t - Julia_total_lnLs1t)
print("\n")
print("R_result_total_lnL (lq) - Julia_sum_lq_nodes: ")
print(R_sum_lq_nodes - Julia_sum_lq_nodes)
print("\n")




##############################################
##############################################
# DEC ML on Psychotria
##############################################
##############################################


#in_params = (birthRate=0.3288164, deathRate=0.0, d_val=0.03505038, e_val=0.02832370, a_val=0.0, j_val=0.0)

bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0

numareas = 4
n = 16            # 4 areas, 16 states

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
#inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, in_params=in_params)
#(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo)
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs

# Update the tip likelihoods, with the geography data
inputs.res.likes_at_each_nodeIndex_branchTop
size(inputs.res.likes_at_each_nodeIndex_branchTop)
numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])

inputs.setup.observed_statenums
inputs.res.likes_at_each_nodeIndex_branchTop

# 2022-03-10: LEAVING OUT THIS PARSER CAUSED BIG PROBLEMS, LNLs COME OUT WRONG!
inputs = Parsers.tipranges_to_tiplikes(inputs, geog_df);
inputs.setup.observed_statenums

inputs.res.likes_at_each_nodeIndex_branchTop

# 2022-03-10: YOU ALSO HAVE TO MAKE SURE THE UPDATED RES IS STORED!!
inputs.res.likes_at_each_nodeIndex_branchTop
res = inputs.res

trdf = inputs.trdf
p_Ds_v5 = inputs.p_Ds_v5
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

prtQi(inputs)
prtCi(inputs)
inputs.p_Ds_v5.params.mu_vals
p_Ds_v5.sol_Es_v5(1.0)

Es_interpolator(1.0)


# Parameters

# Do downpass
res.likes_at_each_nodeIndex_branchTop
(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10)
res.likes_at_each_nodeIndex_branchTop

Rnames(res)
res.likes_at_each_nodeIndex_branchTop
res.normlikes_at_each_nodeIndex_branchTop
res.likes_at_each_nodeIndex_branchBot
res.normlikes_at_each_nodeIndex_branchBot

sum.(res.likes_at_each_nodeIndex_branchTop)
log.(sum.(res.likes_at_each_nodeIndex_branchTop))
sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))
Julia_sum_lq_old = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
nonroot_nodes = get_nonrootnodes(tr)
sum_likes_internal_branch_tops = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
Julia_sum_lq = Julia_sum_lq_old + sum_likes_internal_branch_tops


##############################################
# Consider the pure-birth log-likelihoods
# The Yule-process ML birthrate is just (# internal nodes - 1)/total_tree_length
# Get basic tree info
numTips = sum(trdf.nodeType .== "tip")
numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")


d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]

ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
yuleBirthRate = (numInternal-1) / ttl_tree_length
yuleDeathRate = 0.0					# Yule process has 0 extinction
bd = bd_liks_trdf(trdf, yuleBirthRate, yuleDeathRate)
bd_lnL_noTopo = bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times

# Convert to BioGeoBEARS lnL under Yule process assumption
# Check if the first state/geographic range is null
#if res.inputs.setup.states_list[1] == []
#	include_null_range = true
#end
include_null_range = true
numstates = length(res.normlikes_at_each_nodeIndex_branchTop[1])
equal_root_prob2 = log(1/(numstates-include_null_range)) 
bgb_root_lnL = log(sum(d_root_orig)) + 1.0

# res5t match
res5t = Julia_sum_lq + equal_root_prob2 + bgb_root_lnL - (1-log(1/yuleBirthRate))
# ...matches these in R: compare_BGB_diversitree_DEC_v1.R
# bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/(birthRate))
# bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 
# bgb1 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 + bgb_root_lnL
# (bgb1 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2
bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2
bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)
bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo
# # Julia ClaSSE:		-34.54312638071361
# R BioGeoBEARS:		-34.54313
@test abs(DEC_lnL - bgb_lnL) < 0.01
##############################################
##############################################
##############################################












#######################################################
# OK, let's do an ML inference
#######################################################
# 1. Figure out how ML inference works
# 
# https://julianlsolvers.github.io/Optim.jl/stable/#examples/generated/maxlikenlm/
#
# 
# 2. Write function to update parameters

Qdf_orig1 = prtQp(p_Ds_v5)
Cdf_orig1 = prtCp(p_Ds_v5)
Qdf_orig2 = prtQi(inputs)
Cdf_orig2 = prtCi(inputs)

bmo = construct_BioGeoBEARS_model_object()
d_rownum = collect(1:nrow(bmo))[bmo.rownames .== "d"]
e_rownum = collect(1:nrow(bmo))[bmo.rownames .== "e"]
bmo.est[d_rownum[1]] = 0.03505038
bmo.est[e_rownum[1]] = 0.0283237
bmo.est[d_rownum[1]] = 0.01
bmo.est[e_rownum[1]] = 0.001

bmo


#######################################################
# 2022-03-06: NOW PUT UPDATER INTO A FUNCTION!!
#######################################################

Qdf_orig1
Qdf_new1

# Set inputs to starting values
Rnames(p_Ds_v5)

# Yule birthRate
birthRate = 0.3288164



pars = [0.03505038, 0.02832370]
parnames = ["d", "e"]

# Assuming parameters are in the order of the bmo list of "free" parameters,
# 1. update the bmo
# 2. then update the Qij and Cijk arrays in p_Ds_v5 (the rate parameters arrays)
inputs.bmo.rownames[inputs.bmo.type.=="free"]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars

pars = [0.03505038, 0.02832370]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
prtQp(p_Ds_v5)
prtCp(p_Ds_v5);

pars = [0.01, 0.001]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
prtQp(p_Ds_v5)
prtCp(p_Ds_v5);

pars = [0.1, 0.2]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
prtQp(p_Ds_v5)
prtCp(p_Ds_v5);


function func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="lnL")
	mind = 0.0
	maxd = 5.0
	mine = 0.0
	maxe = 5.0
	minj = 0.0
	maxj = 2.99999
	nan_lnL = -100000
	
	# Get the Q, C
	res = inputs.res
	trdf = inputs.trdf

	inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;

	d_val = inputs.bmo.est[inputs.bmo.rownames.=="d"][1]
	inbounds = true
	if d_val < mind
		Julia_total_lnLs1A = nan_lnL
		inbounds = false
	end
	if d_val > maxd
		Julia_total_lnLs1A = nan_lnL
		inbounds = false
	end

	
	e_val = inputs.bmo.est[inputs.bmo.rownames.=="e"][1]
	if e_val < mine
		Julia_total_lnLs1A = nan_lnL
		inbounds = false
	end
	if e_val > maxe
		Julia_total_lnLs1A = nan_lnL
		inbounds = false
	end


	# Update
	p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
	
	# Check what updated params look like
	prtQp(p_Ds_v5)
	prtCp(p_Ds_v5)
	
	#res2 = deepcopy(res)
	sort!(trdf,"nodeIndex") # MAKE SURE this is sorted properly -- 
	# OTHERWISE I get a CRASH on 
	# iteration 1, node 19
	if inbounds == true
		(total_calctime_in_sec, iteration_number, Julia_sum_lqA, rootstates_lnLA, Julia_total_lnLs1A, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
	else
		rootstates_lnLA = nan_lnL
		Julia_sum_lqA = nan_lnL
		bgb_lnL = nan_lnL
	end
	
	txt = paste0(["pars[1]=", round(pars[1],digits=5), ", pars[2]=", round(pars[2],digits=5), ",	Julia_sum_lqA=", round(Julia_sum_lqA; digits=3), ", rootstates_lnLA=", round(rootstates_lnLA; digits=3), ",	Julia_total_lnLs1A=", Julia_total_lnLs1A, ", bgb_lnL=", round(bgb_lnL, digits=3)])
	print(txt) 
	print("\n")
	
	if returnval == "lnL"
		return(-Julia_total_lnLs1A)
	end
	if returnval == "bgb_lnL"
		return(-bgb_lnL)
	end
	if returnval == "inputs"
		return(inputs)
	end
	# Shouldn't get here
	return(NaN)
end # END function func_to_optimize(pars, parnames)




pars = [0.03505038, 0.02832370]
parnames = ["d", "e"]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
lnL = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func(pars)

pars = [0.03505038, 0.02832370]
parnames = ["d", "e"]
lnL = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func(pars)


pars = [0.03505038, 0.02832370]
parnames = ["d", "e"]
lower = [0.0, 0.0]
upper = [5.0, 5.0]
#func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL")
# http://julianlsolvers.github.io/Optim.jl/v0.9.3/user/minimization/
MLres = optimize(func, pars, LBFGS())

#MLres = NelderMead(func, pars)

p_Ds_v5.params.Qij_vals .= 0.1
iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


pars = MLres.minimizer
parnames = ["d", "e"]
lnL = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func(pars)



end # END @testset "runtests_BiSSE_tree_n3" begin

# Update Qmat
Qmat = (Qarray_ivals=p_Ds_v5.p_indices.Qarray_ivals, Qarray_jvals=p_Ds_v5.p_indices.Qarray_jvals, Qij_vals=p_Ds_v5.params.Qij_vals, Qarray_event_types=p_Ds_v5.p_indices.Qarray_event_types)
Qmat2 = update_Qij_vals(Qmat, areas_list, states_list, dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))), elist=repeat([1.0], length(areas_list)), amat=dmat )

