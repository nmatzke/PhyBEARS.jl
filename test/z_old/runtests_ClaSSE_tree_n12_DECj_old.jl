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
#using Optim                 # for e.g. LBFGS Maximum Likelihood optimization,optimize
# Optim really sucks, try LBFGSB: https://github.com/JuliaNLSolvers/Optim.jl/issues/953
# LBFGSB also sucks: https://github.com/Gnimuc/LBFGSB.jl
using NLopt									# seems to be the best gradient-free, box-constrained								

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DataFrames  # for DataFrame
using DifferentialEquations
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA


# List each PhyBEARS code file prefix here
using PhyBEARS.BGExample
using PhyBEARS.TrUtils
using PhyBEARS.StateSpace
using PhyBEARS.TreeTable
using PhyBEARS.TreePass
using PhyBEARS.Parsers
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
DEC_lnL = -34.54313;
DEC_R_result_branch_lnL = -67.6295;
DEC_R_result_total_LnLs1 = -72.60212;
DEC_R_result_total_LnLs1t = -71.48986;
DEC_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -120.1545;

# DEC+J
DECj_lnL = -20.94759;
DECj_R_result_branch_lnL = -55.37332;
DECj_R_result_total_LnLs1 = -58.83758;
DECj_R_result_total_LnLs1t = -57.72533;
DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -96.34151;

#######################################################


include("/GitHub/PhyBEARS.jl/src/TreePass.jl")
import .TreePass

# Repeat calculation in Julia
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes

# Island numbers (KOMH = 1234) in Rnodenums order:
#island_nums = [3, 3, 2, 2, 3, 3, 2, 1, 1, 3, 4, 2, 1, 1, 1, 1, 1, 1, 2]
lgdata_fn = "/GitHub/PhyBEARS.jl/Rsrc/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Psychotria tree
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")

# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.est[bmo.rownames .== "birthRate"] .= 0.32881638319078066
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
birthRate = 0.32881638319078066
numareas = 4
n = 16            # 4 areas, 16 states

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
#inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, in_params=in_params)
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;

numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);

# Check the interpolator
p_Ds_v5.sol_Es_v5(1.0)
Es_interpolator(1.0)

# Do downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DEC_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=2) == round(R_sum_lq_nodes; digits=2)


@test round(DEC_lnL, digits=2) == round(bgb_lnL, digits=2)
@test round(DEC_R_result_branch_lnL, digits=2) == round(Julia_sum_lq, digits=2)
@test round(DEC_R_result_total_LnLs1, digits=2) == round(Julia_total_lnLs1, digits=2)
@test round(DEC_R_result_total_LnLs1t, digits=2) == round(Julia_total_lnLs1t, digits=2)


print("\nDifferences between Julia and R lnLs for\n/GitHub/PhyBEARS.jl/Rsrc/_compare_ClaSSE_calcs_v3_compare2julia.R\n calculation:\n")
print("DEC_R_result_branch_lnL (lq) - Julia_sum_lq: ")
print(DEC_R_result_branch_lnL - Julia_sum_lq)
print("\n")
print("DEC_R_result_total_LnLs1 (lq) - Julia_total_lnLs1: ")
print(DEC_R_result_total_LnLs1 - Julia_total_lnLs1)
print("\n")
print("DEC_R_result_total_LnLs1t (lq) - Julia_total_lnLs1t: ")
print(DEC_R_result_total_LnLs1t - Julia_total_lnLs1t)
print("\n")
print("DEC_R_result_total_lnL (lq) - Julia_sum_lq_nodes: ")
print(R_sum_lq_nodes - Julia_sum_lq_nodes)
print("\n")

DEC_R_result_total_LnLs1 = -72.60212;
DEC_R_result_total_LnLs1t = -71.48986;






# DEC+J model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.32881638319078066
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 1.0e-12
bmo.est[bmo.rownames .== "e"] .= 1.0e-12
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.1142057

# Need to re-run the setup in order to create the j rows of Cijk_vals
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

@test round(DECj_lnL, digits=2) == round(bgb_lnL, digits=2)
@test round(DECj_R_result_branch_lnL, digits=2) == round(Julia_sum_lq, digits=2)
@test round(DECj_R_result_total_LnLs1, digits=2) == round(Julia_total_lnLs1, digits=2)
@test round(DECj_R_result_total_LnLs1t, digits=2) == round(Julia_total_lnLs1+log(1/birthRate), digits=2)



# Demonstrate the connections between SSE and BGB likelihoods
sum.(res.likes_at_each_nodeIndex_branchTop)
log.(sum.(res.likes_at_each_nodeIndex_branchTop))
sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))
Julia_sum_lq_old = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
nonroot_nodes = get_nonrootnodes(tr)
sum_likes_internal_branch_tops = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
Julia_sum_lq = Julia_sum_lq_old + sum_likes_internal_branch_tops

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=2) == round(R_sum_lq_nodes; digits=2)



##############################################
# Consider the pure-birth log-likelihoods
# The Yule-process ML birthrate is just (# internal nodes - 1)/total_tree_length
# Get basic tree info
numTips = sum(trdf.nodeType .== "tip")
numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")


d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]

ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
birthRate = yuleBirthRate = (numInternal-1) / ttl_tree_length
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
# bgb2 + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/(birthRate))
# bgb2 + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 
# bgb1 + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 + bgb_root_lnL
# (bgb1 + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

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
@test round(DECj_R_result_branch_lnL; digits=2) == round(Julia_sum_lq; digits=2)

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
@test round(DECj_R_result_total_LnLs1; digits=2) == round(Julia_total_lnLs1; digits=2)

# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
@test round(DECj_R_result_total_LnLs1t; digits=2) == round(Julia_total_lnLs1t; digits=2)

# Does the total of branch likelihoods (lq) + node likelihoods match R?
# 
# R: DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda 
#    = sum(log(computed_likelihoods_at_each_node_x_lambda))
res.likes_at_each_nodeIndex_branchTop
log.(sum.(res.likes_at_each_nodeIndex_branchTop))

sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=2) == round(R_sum_lq_nodes; digits=2)








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
print("DECj_R_result_branch_lnL (lq) - Julia_sum_lq: ")
print(DECj_R_result_branch_lnL - Julia_sum_lq)
print("\n")
print("DECj_R_result_total_LnLs1 (lq) - Julia_total_lnLs1: ")
print(DECj_R_result_total_LnLs1 - Julia_total_lnLs1)
print("\n")
print("DECj_R_result_total_LnLs1t (lq) - Julia_total_lnLs1t: ")
print(DECj_R_result_total_LnLs1t - Julia_total_lnLs1t)
print("\n")
print("DECj_R_result_total_lnL (lq) - Julia_sum_lq_nodes: ")
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
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;

# Update the tip likelihoods, with the geography data
inputs.res.likes_at_each_nodeIndex_branchTop
size(inputs.res.likes_at_each_nodeIndex_branchTop)
numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])

inputs.setup.observed_statenums
inputs.res.likes_at_each_nodeIndex_branchTop

inputs = Parsers.tipranges_to_tiplikes(inputs, geog_df);
inputs.setup.observed_statenums

inputs.res.likes_at_each_nodeIndex_branchTop


inputs.res.likes_at_each_nodeIndex_branchTop
res = inputs.res

trdf = inputs.trdf
p_Ds_v5 = inputs.p_Ds_v5
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);

prtQi(inputs)
prtCi(inputs)
inputs.p_Ds_v5.params.mu_vals
p_Ds_v5.sol_Es_v5(1.0)

Es_interpolator(1.0)


# Parameters

# Do downpass
res.likes_at_each_nodeIndex_branchTop
(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10);
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
# bgb2 + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/(birthRate))
# bgb2 + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 
# bgb1 + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 + bgb_root_lnL
# (bgb1 + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

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

#######################################################
# 2022-03-06: NOW PUT UPDATER INTO A FUNCTION!!
#######################################################

function func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="lnL", printlevel=1)

	# Get the results object (with old p_Ds_v5)
	# and treetable (trdf)
	res = inputs.res
	trdf = inputs.trdf
	
	# Set the result if parameters out of bounds
	nan_lnL = -Inf
	
	# Update the parameters
	inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
	
	# Get parameter names and lower/upper bounds
	parnames = inputs.bmo.rownames[inputs.bmo.type.=="free"];
	lower = inputs.bmo.min[bmo.type .== "free"]
	upper = inputs.bmo.max[bmo.type .== "free"]
	#bmo = inputs.bmo
	
	# Update
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est before bmo_updater_v1(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est[inputs.bmo.rownames .== "j"] .= 1.6
	#inputs.bmo.est[:] = bmo_updater_v1(inputs.bmo);
	bmo_updater_v1(inputs.bmo);
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est after bmo_updater_v1(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est .= bmo.est

	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights before p_Ds_v5_updater_v1:")
		print(round.(pdf_before.wt[1:5], digits=4))
	end
#	output = p_Ds_v5_updater_v1(p_Ds_v5, inputs);  # WORKS 2022-03-10
#	p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
#	p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
#	p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
	if printlevel >= 2
		pdf_after = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights after p_Ds_v5_updater_v1:")
		print(round.(pdf_before.wt[1:5], digits=4))
	end

	
	# Check what updated params look like
	#prtQp(p_Ds_v5)
	#prtCp(p_Ds_v5)

	# Double-check that lower/upper limits haven't been breached
	inbounds = true
	for i in 1:length(parnames)
		if pars[i] < lower[i]
			inbounds = false
		end
		if pars[i] > upper[i]
			inbounds = false
		end
	end # END for i in 1:length(parnames)
	
	#sort!(trdf,"nodeIndex") # MAKE SURE this is sorted properly -- 
	# OTHERWISE I get a CRASH on 
	# iteration 1, node 19
	if inbounds == true
		# Solve the Es
		prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
		# This solution is an interpolator
		sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
		p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
		
		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
	else
		Julia_sum_lq = nan_lnL
		rootstates_lnL = nan_lnL
		Julia_total_lnLs1 = nan_lnL
		bgb_lnL = nan_lnL
	end
	
	# Make a string listing the free parameters
	txt = ""
	for i in 1:length(parnames)
		txt = paste0([txt, parnames[i], "=", round(pars[i],digits=5), ",	"])
	end # END for i in 1:length(parnames)
	
	if printlevel >= 1
		txt = paste0([txt, "Julia_sum_lq=", round(Julia_sum_lq; digits=4), ", rootstates_lnL=", round(rootstates_lnL; digits=4), ",	Julia_total_lnLs1=", round(Julia_total_lnLs1, digits=4), ", bgb_lnL=", round(bgb_lnL, digits=4)])
		print("\n")
		print(txt) 
	end
	
	if returnval == "lnL"
		return(-Julia_total_lnLs1)
	end
	if returnval == "bgb_lnL"
		return(-bgb_lnL)
	end
	if returnval == "inputs"
		# Rates vectors
		inputs.p_Ds_v5.params.mu_vals[:] .= p_Ds_v5.params.mu_vals
		inputs.p_Ds_v5.params.Qij_vals[:] .= p_Ds_v5.params.Qij_vals
		inputs.p_Ds_v5.params.Cijk_weights[:] .= p_Ds_v5.params.Cijk_weights
		inputs.p_Ds_v5.params.Cijk_vals[:] .= p_Ds_v5.params.Cijk_vals
		inputs.p_Ds_v5.params.row_weightvals[:] .= p_Ds_v5.params.row_weightvals
		
		# results (res)
		inputs.res.regime[:] .= res.regime
		inputs.res.node_state[:] .= res.node_state
		inputs.res.node_Lparent_state[:] .= res.node_Lparent_state
		inputs.res.node_Rparent_state[:] .= res.node_Rparent_state
		#inputs.res.root_nodeIndex[:] .= res.root_nodeIndex
		#inputs.res.numNodes[:] .= res.numNodes
		#inputs.res.uppass_edgematrix[:] .= res.uppass_edgematrix
		inputs.res.thread_for_each_nodeOp[:] .= res.thread_for_each_nodeOp
		inputs.res.thread_for_each_branchOp[:] .= res.thread_for_each_branchOp
		inputs.res.calc_spawn_start[:] .= res.calc_spawn_start
		inputs.res.calc_start_time[:] .= res.calc_start_time
		inputs.res.calc_end_time[:] .= res.calc_end_time
		inputs.res.calc_duration[:] .= res.calc_duration
		inputs.res.calctime_iterations[:] .= res.calctime_iterations
		inputs.res.sumLikes_at_node_at_branchTop[:] .= res.sumLikes_at_node_at_branchTop
		inputs.res.lnL_at_node_at_branchTop[:] .= res.lnL_at_node_at_branchTop
		inputs.res.lq_at_branchBot[:] .= res.lq_at_branchBot
		inputs.res.like_at_branchBot[:] .= res.like_at_branchBot
		inputs.res.Es_at_each_nodeIndex_branchTop[:] .= res.Es_at_each_nodeIndex_branchTop
		inputs.res.Es_at_each_nodeIndex_branchBot[:] .= res.Es_at_each_nodeIndex_branchBot
		inputs.res.fakeX0s_at_each_nodeIndex_branchTop[:] .= res.fakeX0s_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchTop[:] .= res.likes_at_each_nodeIndex_branchTop
		inputs.res.normlikes_at_each_nodeIndex_branchTop[:] .= res.normlikes_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchBot[:] .= res.likes_at_each_nodeIndex_branchBot
		inputs.res.normlikes_at_each_nodeIndex_branchBot[:] .= res.normlikes_at_each_nodeIndex_branchBot
		
		return(inputs)
	end
	# Shouldn't get here
	return(NaN)
end # END function func_to_optimize(pars, parnames)


func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")

# The function used in NLopt.Opt requires a dummy "grad" input,
# as gradient-using version s of NLopt require this. Here, it's 
# just blank ([]).
# It also has to be modifiable, done with "!".
# (Kind of silly but anyways...
# https://discourse.julialang.org/t/using-nlopt-for-maximum-likelihood-estimation/31964 )
# ===================================
# "NLopt always expects two arguments objective function and 
#  if not given it silently fails because exception at the NLopt 
#  library calls are not caught but causes forced return. I also 
#  had the same problem. Make 
#  loglnorm function loglnorm(par, dummygrad) ll = 0.0 for i = 1:length(y)"
 
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)




#######################################################
# ML inference on DEC+J
#######################################################
# DEC+J model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 1e-12
bmo.est[bmo.rownames .== "e"] .= 1e-12
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.11
bmo.est
bmo.est[:] = bmo_updater_v1(bmo)
bmo.est

#
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;


inputs.bmo.type[bmo.rownames .== "j"] .= "free"
parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works
inputs.bmo

prtCp(p_Ds_v5)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);  # WORKS 2022-03-10
prtCp(p_Ds_v5)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
prtCp(p_Ds_v5)
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
prtCp(p_Ds_v5)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
pars = [0.9, 0.9, 0.9]
func(pars)



#######################################################
# Best optimizer so far - 2022-03-15
#######################################################
using NLopt
pars = [0.9, 0.9, 0.9]
func(pars)
opt = NLopt.Opt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2
opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
opt.lower_bounds
opt.upper_bounds
opt.ftol_abs = 0.00001  # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)

#######################################################






#LN_BOBYQA=34
opt = Opt(:LN_BOBYQA, 3)

opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
opt.min_objective = func
opt.ftol_abs = 0.001
opt.xtol_abs = 0.0001
opt.maxeval = 2000
opt.maxtime = 10000
opt.force_stop = 1
(minf,minx,ret) =  NLopt.optimize(opt, pars)




# Works on 1 run
#MLres = optimize(func, lower, upper, pars, Fminbox(GradientDescent()))

# Never ends, very annoying!
#MLres = optimize(func, lower, upper, pars, Fminbox(GradientDescent(), Optim.Options(x_tol=0.001)))
# https://github.com/Gnimuc/LBFGSB.jl

# copied from wrapper.jl, "typ_bnd"
function type_of_bounds(lb,ub)
	lb > ub && error("Inconsistent bounds")

	lbinf = lb==-Inf
	ubinf = ub==+Inf

	lbinf && ubinf && return 0
	!lbinf && ubinf && return 1
	!lbinf && !ubinf && return 2
	lbinf && !ubinf && return 3
end

function opt_bounds(n,m,lb,ub)
    optimizer=L_BFGS_B(n,m)
    bounds=zeros(3,n)
    bounds[2,:].=lb
    bounds[3,:].=ub
    bounds[1,:].= type_of_bounds.(lb,ub)
    return optimizer,bounds
end

optimizer,bounds = opt_bounds(length(pars),5,lower,upper)
bounds
tmp = lbfgsb(func, pars; m=10, lb=lower, ub=upper, factr=1e7, pgtol=1e-5, iprint=-1, maxfun=1000, maxiter=3)

vector_of_tuples = [(1.0e-12, 4.999999999999), (1.0e-12, 4.999999999999), (1.0e-5, 2.99999)]
for i in 1:length(lower)
	tmp = (Float64(lower[i]), Float64(upper[i]))
	push!(vector_of_tuples, tmp)
end
vector_of_tuples




using JuMP
using NLopt
using SCS
#model = Model(NLopt.Optimizer)
#set_optimizer_attribute(model, "algorithm", :LD_MMA)

dummy_gradient = []
opt=Opt(:LD_LBFGS, 3)
max_objective!(opt, func2(pars,dummy_gradient))
(maxf,maxx,ret) = optimize(opt, [0.0,0.0,0.0])


#LN_BOBYQA=34
opt = Opt(:LN_BOBYQA, 3)

opt.lower_bounds = lb::Union{AbstractVector,Real}
opt.upper_bounds = ub::Union{AbstractVector,Real}
opt.min_objective = func
opt.ftol_abs = 0.001
opt.xtol_abs = 0.0001
opt.maxeval = 2000
opt.maxtime = 10000
opt.nlopt_set_force_stop = 0
(minf,minx,ret) =  NLopt.optimize(opt, pars)




opt = Opt(:NLOPT_LD_SLSQP, 3)
opt.lower_bounds = lb::Union{AbstractVector,Real}
opt.upper_bounds = ub::Union{AbstractVector,Real}
opt.min_objective = func
opt.ftol_abs = 0.001
opt.xtol_abs = 0.0001
opt.maxeval = 2000
opt.maxtime = 10000
tmp2 = Optim.optimize(func, pars, NelderMead(), Optim.Options(x_tol=0.001))


opt = Opt(:NLOPT_LN_SBPLX, 3)


tmp = Optim.optimize(func, pars, NelderMead(), Optim.Options(x_tol=0.001))

using NLopt
pars = [0.9, 0.9, 0.9]
func(pars)
opt = NLoptOpt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2
opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
opt.lower_bounds
opt.upper_bounds
opt.ftol_abs = 0.00001  # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)


model = JuMP.Model(NLopt.Optimizer)
set_optimizer(model, NLopt.Optimizer)
@variable(model, pars[1:length(pars)])
#@constraint(model, [i=1:length(pars)], x0[i] >= lower[i])
#@constraint(model, [i=1:length(pars)], x0[i] <= upper[i])
@objective(model, Min, func2(x0,dummy_gradient))
set_optimizer_attribute(model, "lower_bounds", lb::Union{AbstractVector,Real})
set_optimizer_attribute(model, "upper_bounds", ub::Union{AbstractVector,Real})
set_optimizer_attribute(model, "algorithm", :LN_BOBYQA)
print(model)
JuMP.optimize!(model)
print(model)
solution_summary(model)

set_optimizer_attribute(model, "algorithm", :LN_BOBYQA)
set_optimizer_attribute(model, "min_objective", func)
set_optimizer_attribute(model, "lower_bounds", lb::Union{AbstractVector,Real})
set_optimizer_attribute(model, "upper_bounds", ub::Union{AbstractVector,Real})
set_optimizer_attribute(model, "variable", pars)
JuMP.optimize!(model)
solution_summary(model)



get_optimizer_attribute(model, "variable")

# Doesn't work -- crashes
tmp = bboptimize(func; SearchRange=vector_of_tuples, NumDimensions=3)

bounds = opt_bounds(3,5,lb,ub)
tmp = lbfgsb(func, pars, bounds; factr=1e7, pgtol=1e-5, iprint=-1, maxfun=1000, maxiter=3)

# http://julianlsolvers.github.io/Optim.jl/v0.9.3/user/minimization/
#MLres = optimize(func, pars, LBFGS())
#Terrible! MLres = optimize(func, pars, LBFGS(), Optim.Options(iterations=3, x_tol=0.001))
#MLres = optimize(func, pars, Optim.Options(iterations=3, x_tol=0.001))

# Worked but required 2 runs
# MLres = optimize(func, pars)
# pars = MLres.minimizer
# func(pars)
# MLres = optimize(func, pars)

# Get the inputs & res:
pars = MLres.minimizer
inputs.bmo.est[inputs.bmo.type .== "free"] .= MLres.minimizer
bmo_updater_v1!(inputs.bmo)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


inputs = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="inputs");
Rnames(inputs)


d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]

Julia_sum_lq_old = sum(inputs.res.lq_at_branchBot[1:(length(inputs.res.lq_at_branchBot)-1)])
nonroot_nodes = get_nonrootnodes_trdf(inputs.trdf)
sum_likes_internal_branch_tops = sum(log.(sum.(inputs.res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
Julia_sum_lq = Julia_sum_lq_old + sum_likes_internal_branch_tops

yuleBirthRate = inputs.bmo.est[inputs.bmo.rownames .== "birthRate"][1]
yuleDeathRate = inputs.bmo.est[inputs.bmo.rownames .== "deathRate"][1]
bd_ape = bd_liks_trdf(inputs.trdf, yuleBirthRate, yuleDeathRate)

include_null_range = inputs.setup.states_list[1] == []
numstates = length(inputs.res.normlikes_at_each_nodeIndex_branchTop[1])
equal_root_prob2 = log(1/(numstates-include_null_range)) 
bgb_root_lnL = log(sum(d_root_orig)) + 1.0

bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - (bd_ape.lnL - bd_ape.lnl_topology)

# res5t match
res5t = Julia_sum_lq + equal_root_prob2 + bgb_root_lnL - (1-log(1/yuleBirthRate))
# ...matches these in R: compare_BGB_diversitree_DEC_v1.R
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig)) + equal_root_prob2 + log(1/(birthRate))
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 
# bgb1 + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 + bgb_root_lnL
# (bgb1 + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2

bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2

bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)

bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo







#######################################################
# ML inference on DEC
#######################################################
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "fixed"
bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
bmo.est[:] = bmo_updater_v1(bmo) # works


inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;


parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works
inputs.bmo

prtCp(p_Ds_v5)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);  # WORKS 2022-03-10
prtCp(p_Ds_v5)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
prtCp(p_Ds_v5)
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
prtCp(p_Ds_v5)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
pars = [0.9, 0.9]

func(pars)
# http://julianlsolvers.github.io/Optim.jl/v0.9.3/user/minimization/
#MLres = optimize(func, pars, LBFGS())
#MLres = optimize(func, pars, LBFGS(), Optim.Options(iterations=3, x_tol=0.001))
#MLres = optimize(func, pars, Optim.Options(iterations=3, x_tol=0.001))
MLres = optimize(func, pars)

# Re-run to double-check
#pars = MLres.minimizer
#func(pars)
#MLres = optimize(func, pars)

# Get the inputs & res:
pars = MLres.minimizer
inputs.bmo.est[inputs.bmo.type .== "free"] .= MLres.minimizer
bmo_updater_v1!(inputs.bmo)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


inputs = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="inputs");
Rnames(inputs)


d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]

Julia_sum_lq_old = sum(inputs.res.lq_at_branchBot[1:(length(inputs.res.lq_at_branchBot)-1)])
nonroot_nodes = get_nonrootnodes_trdf(inputs.trdf)
sum_likes_internal_branch_tops = sum(log.(sum.(inputs.res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
Julia_sum_lq = Julia_sum_lq_old + sum_likes_internal_branch_tops

yuleBirthRate = inputs.bmo.est[inputs.bmo.rownames .== "birthRate"][1]
yuleDeathRate = inputs.bmo.est[inputs.bmo.rownames .== "deathRate"][1]
bd_ape = bd_liks_trdf(inputs.trdf, yuleBirthRate, yuleDeathRate)

include_null_range = inputs.setup.states_list[1] == []
numstates = length(inputs.res.normlikes_at_each_nodeIndex_branchTop[1])
equal_root_prob2 = log(1/(numstates-include_null_range)) 
bgb_root_lnL = log(sum(d_root_orig)) + 1.0

bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - (bd_ape.lnL - bd_ape.lnl_topology)

# res5t match
res5t = Julia_sum_lq + equal_root_prob2 + bgb_root_lnL - (1-log(1/yuleBirthRate))
# ...matches these in R: compare_BGB_diversitree_DEC_v1.R
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig)) + equal_root_prob2 + log(1/(birthRate))
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 
# bgb1 + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 + bgb_root_lnL
# (bgb1 + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2

bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2

bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)

bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo







#######################################################
# ML inference on DEC + birth-death
#######################################################
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "fixed"
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= 0.1
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
bmo.est[:] = bmo_updater_v1(bmo) # works


inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo)
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs


parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works
inputs.bmo

prtCp(p_Ds_v5)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);  # WORKS 2022-03-10
prtCp(p_Ds_v5)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
prtCp(p_Ds_v5)
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
prtCp(p_Ds_v5)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL")
pars = [0.9, 0.9, 0.3, 0.2]

func(pars)
# http://julianlsolvers.github.io/Optim.jl/v0.9.3/user/minimization/
#MLres = optimize(func, pars, LBFGS())
#MLres = optimize(func, pars, LBFGS(), Optim.Options(iterations=3, x_tol=0.001))
#MLres = optimize(func, pars, Optim.Options(iterations=3, x_tol=0.001))
MLres = optimize(func, pars)

# Re-run to double-check
#pars = MLres.minimizer
#func(pars)
#MLres = optimize(func, pars)

# Get the inputs & res:
pars = MLres.minimizer
inputs.bmo.est[inputs.bmo.type .== "free"] .= MLres.minimizer
bmo_updater_v1!(inputs.bmo)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)






#######################################################
# ML inference on DEC+J + birth-death
#######################################################
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= 0.1
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.1
bmo.est[:] = bmo_updater_v1(bmo) # works


inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo)
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs


parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works
inputs.bmo

prtCp(p_Ds_v5)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);  # WORKS 2022-03-10
prtCp(p_Ds_v5)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
prtCp(p_Ds_v5)
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
prtCp(p_Ds_v5)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL")
pars = [0.9, 0.9, 0.1, 0.3, 0.2]

func(pars)
# http://julianlsolvers.github.io/Optim.jl/v0.9.3/user/minimization/
#MLres = optimize(func, pars, LBFGS())
#MLres = optimize(func, pars, LBFGS(), Optim.Options(iterations=3, x_tol=0.001))
#MLres = optimize(func, pars, Optim.Options(iterations=3, x_tol=0.001))
MLres = optimize(func, pars)

# Re-run to double-check
pars = MLres.minimizer
func(pars)
MLres = optimize(func, pars)

# Get the inputs & res:
pars = MLres.minimizer
inputs.bmo.est[inputs.bmo.type .== "free"] .= MLres.minimizer
bmo_updater_v1!(inputs.bmo)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)





end # END @testset "runtests_BiSSE_tree_n3" begin
