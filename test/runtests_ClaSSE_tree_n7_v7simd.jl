using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Base.Threads  # <-- this is good for @spawn, Distributed.@spawn is BAD, it produces Futures that have to be scheduled etc]
using Combinatorics					# for e.g. combinations()
using DataFrames

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
using PhyBEARS.ModelLikes

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_ClaSSE_tree_n7_v7simd.jl")
"""
# 
# """
# # Run with:
# include("/GitHub/PhyBEARS.jl/test/runtests_BiSSE_tree_n3.jl")
# """
# 
# @testset "Example" begin
# 	@test hello("runtests_BiSSE_tree_n3.jl") == "Hello, runtests_BiSSE_tree_n3.jl"
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
@testset "runtests_ClaSSE_tree_n7_v7simd.jl" begin
# 
#######################################################
# Calculation of Es and Ds on a single branch
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, pure birth, no Q transitions, branchlength=1)
#
# Run with:
# source("/GitHub/PhyBEARS.jl/Rsrc/_compare_ClaSSE_calcs_n7_compare2julia.R")
# Truth:
R_result_branch_lnL = -7.337676
R_result_total_LnLs1 = -11.105443
R_result_total_LnLs1t = -5.932152
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -12.25723
#######################################################


include("/GitHub/PhyBEARS.jl/src/TreePass.jl")
import .TreePass

# Repeat calculation in Julia
#include("/GitHub/PhyBEARS.jl/src/ModelLikes.jl")
#import .ModelLikes
tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.2, deathRate=1.0, d_val=0.5, e_val=0.4, a_val=0.0, j_val=1.5)
numareas = 2
n = 3

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.2
bmo.est[bmo.rownames .== "deathRate"] .= 1.0
bmo.est[bmo.rownames .== "d"] .= 0.5
bmo.est[bmo.rownames .== "e"] .= 0.4
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 1.5
bmo.est[:] = bmo_updater_v1(bmo);

# Set up the model; NaN means areas & ranges are auto-generated
geog_df = NaN
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, ; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
prtCi(inputs)
# Change parameter inputs manually
# inputs.p_Ds_v5.params.Qij_vals[1] = 2*inputs.p_Ds_v5.params.Qij_vals[2]
# inputs.p_Ds_v5.params.Cijk_vals[1] = 2*inputs.p_Ds_v5.params.Cijk_vals[2]
# inputs.p_Ds_v5.params.mu_vals[2] = 2*inputs.p_Ds_v5.params.mu_vals[1]

res = inputs.res
res.likes_at_each_nodeIndex_branchTop
res.normlikes_at_each_nodeIndex_branchTop
# All tips equal, to match runtests_ClaSSE_tree_n7.jl
res.likes_at_each_nodeIndex_branchTop[1] = [0.0, 1.0, 0.0]			# state 2 for tip #1
res.normlikes_at_each_nodeIndex_branchTop[1] = [0.0, 1.0, 0.0]	# state 2 for tip #1
res.likes_at_each_nodeIndex_branchTop[2] = [1.0, 0.0, 0.0]			# state 2 for tip #2
res.normlikes_at_each_nodeIndex_branchTop[2] = [1.0, 0.0, 0.0]	# state 2 for tip #2
res.likes_at_each_nodeIndex_branchTop[4] = [0.0, 1.0, 0.0]			# state 2 for tip #3
res.normlikes_at_each_nodeIndex_branchTop[4] = [0.0, 1.0, 0.0]	# state 2 for tip #3
res.likes_at_each_nodeIndex_branchTop
res.normlikes_at_each_nodeIndex_branchTop
trdf = inputs.trdf
p_Es_v7 = inputs.p_Ds_v5
root_age = maximum(trdf[!, :node_age])
Es_tspan = inputs.Es_tspan
solver_options = inputs.solver_options


# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v7 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7)
# This solution is an interpolator
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v7;
p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms= p_Es_v7.terms, sol_Es_v5=sol_Es_v7);

prtQi(inputs)
prtCi(inputs)
inputs.p_Ds_v5.params.mu_vals
p_Ds_v7.params.mu_vals
p_Ds_v7.sol_Es_v5(1.0)

Es_interpolator(1.0)


# Parameters

# Do downpass
res = inputs.res;
(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=construct_SolverOpt(), max_iterations=10^10)

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
lambda = in_params.birthRate
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



#######################################################
# Test the GFlow calculation
#######################################################



end # END @testset "runtests_BiSSE_tree_n3" begin




