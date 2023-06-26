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
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames						# for DataFrame()
using DelimitedFiles				# for readdlm()
using NLopt									# seems to be the best gradient-free, box-constrained								

# List each PhyBEARS code file prefix here
using PhyloBits
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
cd("/GitHub/PhyBEARS.jl/ex/cinthy/")
include("/GitHub/PhyBEARS.jl/ex/cinthy/Epic_ClaSSE_tree_n13_DECj_v2.jl")
#
Check against:
cd /GitHub/PhyBEARS.jl/ex/cinthy/
Rscript /GitHub/PhyBEARS.jl/ex/cinthy/epic_M0_ML_BSM_v1.R
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
#@testset "runtests_ClaSSE_tree_n13_DECj_WORKS.jl" begin
# 
#######################################################
# DEMONSTRATES MATCHING BETWEEN DIVERSITREE, BIOGEOBEARS, AND JULIA
# ON HAWAIIAN Epacridoideae, 16-STATE DEC MODEL
#
# Run with:
# source("/GitHub/PhyBEARS.jl/Rsrc/compare_BGB_diversitree_DEC+J_v1.R")
# Truth:
#DECj_lnL = -34.54313;
#DECj_R_result_branch_lnL = -67.6295;
#DECj_R_result_total_LnLs1 = -72.60212;
#DECj_R_result_total_LnLs1t = -71.48986;
#DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -120.1545;

# DEC+J
DECj_lnL = -71.95101;
DECj_R_result_branch_lnL = -383.0655;
DECj_R_result_total_LnLs1 = -387.4385;  # res1
DECj_R_result_total_LnLs1t = -385.2824; # res1t
#DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -96.34151;
#######################################################


#include("/GitHub/PhyBEARS.jl/src/TreePass.jl")
#import .TreePass

# Repeat calculation in Julia
#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes

# Island numbers (KOMH = 1234) in Rnodenums order:
#island_nums = [3, 3, 2, 2, 3, 3, 2, 1, 1, 3, 4, 2, 1, 1, 1, 1, 1, 1, 2]
lgdata_fn = "geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Epacridoideae tree
trfn = "tree.newick"
tr = PhyloBits.PNreadwrite.readTopology(trfn)

ML_yule_birthRate(tr)
ML_yule_birthRate_wRoot(tr)


# DEC model on Hawaiian Epacridoideae
bmo = construct_BioGeoBEARS_model_object()
birthRate = 0.11578
bmo.est[bmo.rownames .== "birthRate"] .= birthRate
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.0007341845
bmo.est[bmo.rownames .== "e"] .= 1.000000e-12 
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.05320481
bmo.type[bmo.rownames .== "j"] .= "free"
numareas = 3
n = 8            # 3 areas, 8 states

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
root_age_mult=1.5; max_range_size=NaN; include_null_range=false; max_range_size=NaN
max_range_size = NaN # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;

numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, sol_Es_v5=sol_Es_v5);

# Check the interpolator
p_Ds_v5.sol_Es_v5(1.0)
Es_interpolator(1.0)

# Do downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


# Do downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
#R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
#@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test abs(DECj_lnL - bgb_lnL) < 0.01
@test abs(DECj_R_result_branch_lnL - Julia_sum_lq) < 0.01
@test abs(DECj_R_result_total_LnLs1 - Julia_total_lnLs1) < 0.01
@test abs(DECj_R_result_total_LnLs1t - Julia_total_lnLs1t) < 0.01

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
#R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
#@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test abs(DECj_lnL - bgb_lnL) < 0.01
@test abs(DECj_R_result_branch_lnL - Julia_sum_lq) < 0.01
@test abs(DECj_R_result_total_LnLs1 - Julia_total_lnLs1) < 0.01
@test abs(DECj_R_result_total_LnLs1t - Julia_total_lnLs1t) < 0.01

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
#R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
#@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test abs(DECj_lnL - bgb_lnL) < 0.01
@test abs(DECj_R_result_branch_lnL - Julia_sum_lq) < 0.01
@test abs(DECj_R_result_total_LnLs1 - Julia_total_lnLs1) < 0.01
@test abs(DECj_R_result_total_LnLs1t - Julia_total_lnLs1t) < 0.01




# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
#R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
#@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)


print("\nDifferences between Julia and R lnLs for\n/GitHub/PhyBEARS.jl/ex/cinthy/epic_M0_ML_BSM_v1.R\n calculation:\n")
print("DECj_R_result_branch_lnL (lq) - Julia_sum_lq: ")
print(DECj_R_result_branch_lnL - Julia_sum_lq)
print("\n")
print("DECj_R_result_total_LnLs1 (lq) - Julia_total_lnLs1: ")
print(DECj_R_result_total_LnLs1 - Julia_total_lnLs1)
print("\n")
print("DECj_R_result_total_LnLs1t (lq) - Julia_total_lnLs1t: ")
print(DECj_R_result_total_LnLs1t - Julia_total_lnLs1t)
print("\n")
#print("DECj_R_result_total_lnL (lq) - Julia_sum_lq_nodes: ")
#print(R_sum_lq_nodes - Julia_sum_lq_nodes)
#print("\n")




#######################################################
# ML inference on DEC+J
#######################################################
# DEC+J model on Hawaiian Epacridoideae
bmo = construct_BioGeoBEARS_model_object()
birthRate = 0.11578
bmo.est[bmo.rownames .== "birthRate"] .= birthRate
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 7.260439e-04
bmo.est[bmo.rownames .== "e"] .= 1.000000e-12 
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 5.775812e-02
bmo.type[bmo.rownames .== "j"] .= "free"
numareas = 3
n = 8            # 3 areas, 8 states
bmo.est[:] = bmo_updater_v1(bmo, inputs.setup.bmo_rows);
bmo.est

max_range_size = 3 # replaces any background max_range_size=1
root_age_mult=1.5; max_range_size=NaN; include_null_range=false; max_range_size=NaN
max_range_size = NaN # replaces any background max_range_size=1
global inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;
res_ARCHIVE=deepcopy(res);

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



pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL", printlevel=1)
#pars = [0.9, 0.9, 0.9]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


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
#opt.ftol_abs = 0.001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################
# For saving to dataframe
optim_result = build_optim_result(opt, optf, optx, ret)


# Get the inputs & res:
pars = optx
#pars = [0.9747407112459348, 0.8, 0.11]
#pars = [100.0, 1.8, 0.11]
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
bmo_updater_v1!(inputs.bmo)
p_Ds_v5_updater_v1!(p_Es_v5, inputs);
# SEE runtests_ClaSSE_tree_n13_DECj_WORKS.jl
# save_everystep_EQ_false_CAN_MATTER_EVEN_ON_THE_Ds
#inputs.solver_options.save_everystep=false # CAN PRODUCE A -20.9 vs. -20.6 difference!
inputs.solver_options.save_everystep=true	# WORKS!! Can make a difference EVEN ON THE Ds!!

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
#R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
#@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
#R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
#@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
#R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
#@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)



# Calculate lnLs
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
# ...matches these in R: compare_BGB_diversitree_DECj_v1.R
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig)) + equal_root_prob2 + log(1/(birthRate))
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 
# bgb1 + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 + bgb_root_lnL
# (bgb1 + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
bgb2 = res5t - (bd_ape.lnL - bd_ape.lnl_topology) - equal_root_prob2
bgb2 = res5t - (bd_ape.lnl_numBirths + bd_ape.lnl_Births_above_root + bd_ape.lnl_numtips_wOneMinusDeathRate + bd_ape.lnl_branching_times) - equal_root_prob2

bgb2 = res5t - (bd_ape.lnl_numBirths + bd_ape.lnl_Births_above_root + bd_ape.lnl_branching_times) - equal_root_prob2

bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd_ape.lnL - bd_ape.lnl_topology)

bd_lnL_noTopo = bd_ape.lnL - bd_ape.lnl_topology
bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo



# Ancestral states estimation
uppass_ancstates_v7!(inputs.res, trdf, p_Ds_v7, solver_options);
resDECj = deepcopy(inputs.res);
geogfn = lgdata_fn
lnLs_tuple = (total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, total_loglikelihood=bgb_lnL)
optim_result = build_optim_result(opt, optf, optx, ret)
juliaRes_to_Rdata(resDECj, trdf, inputs, lnLs_tuple, optim_results, geogfn, trfn; outfns=NaN)
"""
# Then run, in R:
sourceall("~/HD/GitHub/PhyBEARS.jl/Rsrc/")
res = PhyBEARS_res_to_BGB_res(outfns=NaN)
"""


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
bmo.est[:] = bmo_updater_v1(bmo, inputs.setup.bmo_rows) # works


root_age_mult=1.5; max_range_size=NaN; include_null_range=false; max_range_size=NaN
max_range_size = NaN # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;



parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo); # works
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);  # WORKS 2022-03-10

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, sol_Es_v5=sol_Es_v5);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL")
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)

pars = [0.9, 0.9, 0.1, 0.3, 0.2]
func(pars)
func2(pars, [])

#######################################################
# Best optimizer so far - 2022-03-15
#######################################################
using NLopt
opt = NLopt.Opt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2
opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
opt.lower_bounds
opt.upper_bounds
#opt.ftol_abs = 0.00001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################

# Get the inputs & res:
pars = optx

# 2023-02-08_archived pars
optx = [1.0e-12, 1.0e-12, 0.10999899813382318, 0.34816245870282203, 0.0]
pars = optx

inputs.bmo.est[inputs.bmo.type .== "free"] .= optx
bmo_updater_v1!(inputs.bmo);
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

func(pars)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);


# inputs.solver_options.save_everystep = false
#	inputs.solver_options.saveat = nodetimes(trdf)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


# Ancestral states estimation
uppass_ancstates_v7!(inputs.res, trdf, p_Ds_v7, solver_options);
resDECj_BD = deepcopy(inputs.res);

Rnames(resDECj_BD)
vfft(resDECj_BD.likes_at_each_nodeIndex_branchTop)
vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchTop)
vfft(resDECj_BD.likes_at_each_nodeIndex_branchBot)
vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchBot)
vfft(resDECj_BD.anc_estimates_at_each_nodeIndex_branchBot)
vfft(resDECj_BD.anc_estimates_at_each_nodeIndex_branchTop)

rowSums(vfft(resDECj_BD.likes_at_each_nodeIndex_branchTop))
rowSums(vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchTop))
rowSums(vfft(resDECj_BD.likes_at_each_nodeIndex_branchBot))
rowSums(vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchBot))
rowSums(vfft(resDECj_BD.anc_estimates_at_each_nodeIndex_branchBot))
rowSums(vfft(resDECj_BD.anc_estimates_at_each_nodeIndex_branchTop))


# Write out to tables, for constructing a BioGeoBEARS R object
lnLs_tuple = (total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, total_loglikelihood=bgb_lnL)









# Birth EQUALS death
#######################################################
# ML inference on DEC+J + birth=death
#######################################################
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "birthRate"
bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= bmo.est[bmo.rownames .== "birthRate"]
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.1
bmo.est[:] = bmo_updater_v1(bmo, inputs.setup.bmo_rows) # works


root_age_mult=1.5; max_range_size=NaN; include_null_range=false; max_range_size=NaN
max_range_size = NaN # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;



parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo); # works
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);  # WORKS 2022-03-10

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, sol_Es_v5=sol_Es_v5);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL")
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)

pars = [0.9, 0.9, 0.1, 0.3]
func(pars)
func2(pars, [])

#######################################################
# Best optimizer so far - 2022-03-15
#######################################################
using NLopt
opt = NLopt.Opt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2
opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
opt.lower_bounds
opt.upper_bounds
#opt.ftol_abs = 0.00001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################

# Get the inputs & res:
pars = optx

# 2023-02-08_archived pars
#optx = [1.0e-12, 1.0e-12, 0.10999899813382318, 0.34816245870282203, 0.0]
pars = optx

inputs.bmo.est[inputs.bmo.type .== "free"] .= optx
bmo_updater_v1!(inputs.bmo);
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

func(pars)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);


# inputs.solver_options.save_everystep = false
#	inputs.solver_options.saveat = nodetimes(trdf)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

resDECj_BD = deepcopy(inputs.res);
Rnames(resDECj_BD)
vfft(resDECj_BD.likes_at_each_nodeIndex_branchTop)
vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchTop)
vfft(resDECj_BD.likes_at_each_nodeIndex_branchBot)
vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchBot)

rowSums(vfft(resDECj_BD.likes_at_each_nodeIndex_branchTop))
rowSums(vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchTop))
rowSums(vfft(resDECj_BD.likes_at_each_nodeIndex_branchBot))
rowSums(vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchBot))










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
bmo.est[:] = bmo_updater_v1(bmo, inputs.setup.bmo_rows) # works

#
root_age_mult=1.5; max_range_size=NaN; include_null_range=true; max_range_size=NaN
max_range_size = NaN # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;



inputs.bmo.type[bmo.rownames .== "j"] .= "fixed"
parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works
inputs.bmo

prtCp(p_Es_v5)
p_Ds_v5_updater_v1!(p_Es_v5, inputs);  # WORKS 2022-03-10
prtCp(p_Es_v5)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
prtCp(p_Es_v5)
p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);
prtCp(p_Es_v5)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)




(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)





pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL", printlevel=1)
pars = [0.9, 0.9]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


#######################################################
# Best optimizer so far - 2022-03-15
#######################################################
using NLopt
pars = [0.9, 0.9]
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
opt.ftol_abs = 0.00001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################




# Get the inputs & res:
pars = deepcopy(optx)
#pars = [0.03505038, 0.02832370]
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
#bmo_updater_v1!(inputs.bmo) # works
inputs.bmo.est .= bmo_updater_v1(inputs.bmo, inputs.setup.bmo_rows) # works
inputs.bmo
res = inputs.res;

p_Ds_v5_updater_v1!(p_Es_v5, inputs);  # WORKS 2022-03-10

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)




(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

# If you add BioGeoBEARS node likelihoods to Julia branch likelihoods...
Julia_total_lnLs1t = Julia_total_lnLs1 + log(1/birthRate)
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=1) == round(R_sum_lq_nodes; digits=1)

@test round(DECj_lnL, digits=1) == round(bgb_lnL, digits=1)
@test round(DECj_R_result_branch_lnL, digits=1) == round(Julia_sum_lq, digits=1)
@test round(DECj_R_result_total_LnLs1, digits=1) == round(Julia_total_lnLs1, digits=1)
@test round(DECj_R_result_total_LnLs1t, digits=1) == round(Julia_total_lnLs1t, digits=1)


# Calculate lnLs
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
# ...matches these in R: compare_BGB_diversitree_DECj_v1.R
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig)) + equal_root_prob2 + log(1/(birthRate))
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 
# bgb1 + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 + bgb_root_lnL
# (bgb1 + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2

bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2

bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)

bd_lnL_noTopo = bd_ape.lnL - bd_ape.lnl_topology
bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo







#######################################################
# ML inference on DEC + birth-death
#######################################################
Julia_sum_lq_ORIG = -68.50874661600467
Julia_rootstates_lnL_ORIG = -4.972290561776454
Julia_total_lnLs1_ORIG = -73.48103717778113
Julia_bgb_lnL_ORIG = -35.29504501641623



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
bmo.est[:] = bmo_updater_v1(bmo, inputs.setup.bmo_rows) # works


root_age_mult=1.5; max_range_size=NaN; include_null_range=false; max_range_size=NaN
max_range_size = NaN # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;



inputs.bmo.type[bmo.rownames .== "j"] .= "fixed"
parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works
inputs.bmo

p_Ds_v5_updater_v1!(p_Ds_v5, inputs);  # WORKS 2022-03-10

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, sol_Es_v5=sol_Es_v5);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

@test round(Julia_bgb_lnL_ORIG, digits=3) == round(bgb_lnL, digits=3)
@test round(Julia_sum_lq_ORIG, digits=3) == round(Julia_sum_lq, digits=3)
@test round(Julia_total_lnLs1_ORIG, digits=3) == round(Julia_total_lnLs1, digits=3)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

@test round(Julia_bgb_lnL_ORIG, digits=3) == round(bgb_lnL, digits=3)
@test round(Julia_sum_lq_ORIG, digits=3) == round(Julia_sum_lq, digits=3)
@test round(Julia_total_lnLs1_ORIG, digits=3) == round(Julia_total_lnLs1, digits=3)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

@test round(Julia_bgb_lnL_ORIG, digits=3) == round(bgb_lnL, digits=3)
@test round(Julia_sum_lq_ORIG, digits=3) == round(Julia_sum_lq, digits=3)
@test round(Julia_total_lnLs1_ORIG, digits=3) == round(Julia_total_lnLs1, digits=3)



pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
pars = [0.9, 0.9, 0.3, 0.2]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


#######################################################
# Best optimizer so far - 2022-03-15
#######################################################
using NLopt
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
opt.ftol_abs = 0.00001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################





# Get the inputs & res:
pars = optx
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
inputs_updater_v1!(inputs);
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

printlevel=1
returnval="bgb_lnL"
func(pars);

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);


# inputs.solver_options.save_everystep = false
#	inputs.solver_options.saveat = nodetimes(trdf)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

print("\n\nruntests_ClaSSE_tree_n13_DECj_WORKS.jl: Different optimizations on DEC+BD at home vs work...\n\n")

# Worked at home: 2023-02-11
Julia_sum_lq_ORIG = -68.63952019916809
Julia_rootstates_lnL_ORIG = -4.971193862177472
Julia_total_lnLs1_ORIG = -73.61071406134556
Julia_bgb_lnL_ORIG = -35.43075554544988

# Different optimizations, home and work
# Worked at work: 2023-02-10
Julia_sum_lq_ORIG_atWORK = -69.5224885849644
Julia_rootstates_lnL_ORIG_atWORK = -4.968808024086585
Julia_total_lnLs1_ORIG_atWORK = -74.49129660905099
Julia_bgb_lnL_ORIG_atWORK = -36.33474261227199

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


cutoff = 0.1
# @test Julia_bgb_lnL_ORIG - bgb_lnL) < cutoff
bgb_diffs = Julia_bgb_lnL_ORIG - bgb_lnL
bgb_diffs_atWORK = Julia_bgb_lnL_ORIG_atWORK - bgb_lnL

print("\nbgb_diffs: ")
print(bgb_diffs)
print("\nbgb_diffs_atWORK: ")
print(bgb_diffs_atWORK)

TF1 = abs(bgb_diffs) < cutoff
TF2 = abs(bgb_diffs_atWORK) < cutoff
@test TF1 || TF2

#@test Julia_sum_lq_ORIG - Julia_sum_lq) < cutoff

Julia_sum_lq_diffs = Julia_sum_lq_ORIG - Julia_sum_lq
Julia_sum_lq_diffs_atWORK = Julia_sum_lq_ORIG_atWORK - Julia_sum_lq

print("\nJulia_sum_lq_diffs: ")
print(Julia_sum_lq_diffs)
print("\nJulia_sum_lq_diffs_atWORK: ")
print(Julia_sum_lq_diffs_atWORK)

TF1 = abs(Julia_sum_lq_diffs) < cutoff
TF2 = abs(Julia_sum_lq_diffs_atWORK) < cutoff
@test TF1 || TF2


#@test Julia_total_lnLs1_ORIG - Julia_total_lnLs1) < cutoff
TF1 = abs(Julia_total_lnLs1_ORIG - Julia_total_lnLs1) < cutoff
TF2 = abs(Julia_total_lnLs1_ORIG_atWORK - Julia_total_lnLs1) < cutoff
@test TF1 || TF2

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

#@test Julia_bgb_lnL_ORIG - bgb_lnL) < cutoff
TF1 = abs(Julia_bgb_lnL_ORIG - bgb_lnL) < cutoff
TF2 = abs(Julia_bgb_lnL_ORIG_atWORK - bgb_lnL) < cutoff
@test TF1 || TF2

#@test Julia_sum_lq_ORIG - Julia_sum_lq) < cutoff
TF1 = abs(Julia_sum_lq_ORIG - Julia_sum_lq) < cutoff
TF2 = abs(Julia_sum_lq_ORIG_atWORK - Julia_sum_lq) < cutoff
@test TF1 || TF2


#@test Julia_total_lnLs1_ORIG - Julia_total_lnLs1) < cutoff
TF1 = abs(Julia_total_lnLs1_ORIG - Julia_total_lnLs1) < cutoff
TF2 = abs(Julia_total_lnLs1_ORIG_atWORK - Julia_total_lnLs1) < cutoff
@test TF1 || TF2

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

#@test Julia_bgb_lnL_ORIG - bgb_lnL) < cutoff
TF1 = abs(Julia_bgb_lnL_ORIG - bgb_lnL) < cutoff
TF2 = abs(Julia_bgb_lnL_ORIG_atWORK - bgb_lnL) < cutoff
@test TF1 || TF2

#@test Julia_sum_lq_ORIG - Julia_sum_lq) < cutoff
TF1 = abs(Julia_sum_lq_ORIG - Julia_sum_lq) < cutoff
TF2 = abs(Julia_sum_lq_ORIG_atWORK - Julia_sum_lq) < cutoff
@test TF1 || TF2

#@test Julia_total_lnLs1_ORIG - Julia_total_lnLs1) < cutoff
TF1 = abs(Julia_total_lnLs1_ORIG - Julia_total_lnLs1) < cutoff
TF2 = abs(Julia_total_lnLs1_ORIG_atWORK - Julia_total_lnLs1) < cutoff
@test TF1 || TF2

print("\n\n...END of runtests_ClaSSE_tree_n13_DECj_WORKS.jl: Different optimizations on DEC+BD at home vs work...\n\n")






# Calculate lnLs
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
# ...matches these in R: compare_BGB_diversitree_DECj_v1.R
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig)) + equal_root_prob2 + log(1/(birthRate))
bgb_lnL + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 
# bgb1 + bd_ape.lnL - bd_ape.lnl_topology + equal_root_prob2 + bgb_root_lnL
# (bgb1 + bd_ape.lnL - bd_ape.lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2

bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2

bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd_ape.lnL - bd_ape.lnl_topology)

bd_lnL_noTopo = bd_ape.lnL - bd_ape.lnl_topology
bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo






end # END @testset "runtests_BiSSE_tree_n3" begin
