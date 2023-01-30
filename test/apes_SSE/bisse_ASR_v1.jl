

# List each PhyBEARS code file prefix here
using PhyBEARS.Example
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.TrUtils
using PhyBEARS.SSEs

using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using PhyloBits.TrUtils
using DataFrames
using CSV

using PhyBEARS
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.ModelLikes


# 
# """
# using darwins finches
# 
# 
# #######################################################
# # Do a bunch of tests of the SSE calculation of 
# # Ds, Es, and likelihoods, on
# # branches, nodes, and trees,
# # under a variety of simple and more complex models
# #######################################################
# 
# @testset "runtests_BiSSE_tree_n3.jl" begin
# 
#######################################################
# Calculation of Es and Ds on a single branch
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, pure birth, no Q transitions, branchlength=1)
#
# R code based on
# source("/GitHub/PhyBEARS.jl/wallis/TreeBig.R")

# Truth:
R_result_branch_lnL = -17.28156
R_result_total_LnLs1 = -18.44422
R_result_total_LnLs1t = -10.82559
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = NaN
#######################################################


#include("/GitHub/PhyBEARS.jl/src/TreePass.jl")
#import .TreePass

# Repeat calculation in Julia
#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes

tr = readTopology("((sp4:0.6248637277,sp5:0.6248637277):6.489662918,(sp6:0.1274213816,sp7:0.1274213816):6.987105264);")
in_params = (birthRate=0.222222222, deathRate=0.111111111, d_val=0.0, e_val=0.0, a_val=0.1, j_val=0.0)
# pars <- c(0.222222222, 0.222222222, 0.111111111, 0.05, 0.1, 0.15)

# phy$tip.state
# sp4 sp5 sp6 sp7 
#   1   1   0   0 

numstates = 2
n = 2

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
inputs = ModelLikes.setup_MuSSE_biogeo(numstates, tr; root_age_mult=1.5, in_params=in_params);
(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs;

# Change parameter inputs manually
inputs.p_Ds_v5.params.Cijk_vals[1] = 0.222222222
inputs.p_Ds_v5.params.Cijk_vals[2] = 0.222222222
inputs.p_Ds_v5.params.mu_vals[1] = 0.111111111
inputs.p_Ds_v5.params.mu_vals[2] = 0.05
inputs.p_Ds_v5.params.Qij_vals[1] = 0.1
inputs.p_Ds_v5.params.Qij_vals[2] = 0.15

prtQp(p_Ds_v5)
prtCp(p_Ds_v5)



inputs.res.likes_at_each_nodeIndex_branchTop
inputs.res.normlikes_at_each_nodeIndex_branchTop
res.likes_at_each_nodeIndex_branchTop[1] = [0.0, 1.0]
res.likes_at_each_nodeIndex_branchTop[2] = [0.0, 1.0]
res.likes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0]
res.likes_at_each_nodeIndex_branchTop[5] = [1.0, 0.0]
res.normlikes_at_each_nodeIndex_branchTop[1] = [0.0, 1.0]
res.normlikes_at_each_nodeIndex_branchTop[2] = [0.0, 1.0]
res.normlikes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0]
res.normlikes_at_each_nodeIndex_branchTop[5] = [1.0, 0.0]
trdf = inputs.trdf
p_Ds_v5 = inputs.p_Ds_v5;
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is a linear interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);


# Parameters - check
prtQi(inputs)
prtQp(p_Ds_v5)
prtCi(inputs)
prtCp(p_Ds_v5)

inputs.p_Ds_v5.params.mu_vals
p_Ds_v5.sol_Es_v5(1.0)

# Do downpass
(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10);



Rnames(res)
res.likes_at_each_nodeIndex_branchTop
res.normlikes_at_each_nodeIndex_branchTop
res.likes_at_each_nodeIndex_branchBot
res.normlikes_at_each_nodeIndex_branchBot

sum.(res.likes_at_each_nodeIndex_branchTop)
log.(sum.(res.likes_at_each_nodeIndex_branchTop))
sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))
Julia_sum_lq = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
sum(res.lq_at_branchBot)

# Does the total of the branch log-likelihoods (lq) match?
# @test round(R_result_branch_lnL; digits=4) == round(Julia_sum_lq; digits=4)

ttl_from_Julia = sum(res.lq_at_branchBot) +  sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))
@test round(R_result_total_LnLs1; digits=4) == round(ttl_from_Julia; digits=4)


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
d_root = d_root_orig ./ sum(root_stateprobs .* inputs.p_Ds_v5.params.Cijk_vals .* (1 .- e_root).^2)
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
# @test round(Julia_sum_lq_nodes; digits=4) == round(R_sum_lq_nodes; digits=4)


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




print("\nDifferences between Julia and R lnLs for\n/GitHub/PhyBEARS.jl/test/BiSSE_branchlikes_w_MLE_v6_WORKING.R\n calculation:\n")
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

# end # END @testset "runtests_BiSSE_tree_n3" begin




