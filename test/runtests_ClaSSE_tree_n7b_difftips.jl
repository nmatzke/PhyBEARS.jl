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

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_ClaSSE_tree_n7a_difftips.jl")
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
@testset "runtests_ClaSSE_tree_n7a_difftips.jl" begin
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
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.2, deathRate=1.0, d_val=0.5, e_val=0.4, a_val=0.0, j_val=1.5)
numareas = 2
n = 3

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, in_params=in_params)
(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs

# Change parameter inputs manually
# inputs.p_Ds_v5.params.Qij_vals[1] = 2*inputs.p_Ds_v5.params.Qij_vals[2]
# inputs.p_Ds_v5.params.Cijk_vals[1] = 2*inputs.p_Ds_v5.params.Cijk_vals[2]
# inputs.p_Ds_v5.params.mu_vals[2] = 2*inputs.p_Ds_v5.params.mu_vals[1]
inputs.p_Ds_v5.params.mu_vals .= [0.1, 0.5, 1.0]

res.likes_at_each_nodeIndex_branchTop[2] = [1.0, 0.0, 0.0]			# state 1 for tip #2
res.normlikes_at_each_nodeIndex_branchTop[2] = [1.0, 0.0, 0.0]	# state 1 for tip #2
inputs.res.likes_at_each_nodeIndex_branchTop
inputs.res.normlikes_at_each_nodeIndex_branchTop

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
(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10)

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
p_inputs = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);

# Version 5 ClaSSE standard calculation
res_nonFlow_v5 = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf, p_Ds_v5=p_inputs, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true)
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv5, rootstates_lnL_nFv5, Julia_total_lnLs1_nFv5) = res_nonFlow_v5
nf5_normlikes_top = deepcopy(res.normlikes_at_each_nodeIndex_branchTop);
nf6_normlikes_bot = deepcopy(res.normlikes_at_each_nodeIndex_branchBot);


# Version 6 ClaSSE standard calculation
res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_inputs, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnL_nFv6) = res_nonFlow_v6

nf6_normlikes_top = deepcopy(res.normlikes_at_each_nodeIndex_branchTop);
nf6_normlikes_bot = deepcopy(res.normlikes_at_each_nodeIndex_branchBot);


# Version 6 ClaSSE Gflow calculation
include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
import .Flow

# p_inputs -- add the sol_Es_v5
p_inputs = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);

# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
# The "I" requires "include NumericAlgebra"
G0 = Matrix{Float64}(I, n, n) ;

# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));

pG = (n=n, p_Ds_v5=p_inputs, A=A);
tspan = (0.0, 2.5);
prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG);

saveats = seq(tspan[1], tspan[2], tspan[2]/10000);
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, saveat=saveats, abstol = 1e-9, reltol = 1e-9);
#Gflow_to_01_GMRES  = solve(prob_Gs_v5, Tsit5(), save_everystep=true, saveat=saveats, abstol = 1e-9, reltol = 1e-9);

# What kind of solver?
#Gflow_to_01_GMRES  = solve(prob_Gs_v5, save_everystep=true, saveat=saveats, abstol = 1e-9, reltol = 1e-9, alg_hints=[:stiff]);
 
# Downpass with the Gflow system
Gflow = Gflow_to_01_GMRES;

res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);

gv6_normlikes_top = deepcopy(res.normlikes_at_each_nodeIndex_branchTop);
gv6_normlikes_bot = deepcopy(res.normlikes_at_each_nodeIndex_branchBot);

#(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = tmp_results1
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6) = res_Gflow_v6


res_nonFlow_v5
res_nonFlow_v6
res_Gflow_v6
res_nonFlow_v6[3:6] .- res_Gflow_v6[3:6]

nf5_normlikes_top
nf6_normlikes_top
gv6_normlikes_top
nf6_normlikes_bot
nf6_normlikes_bot 
gv6_normlikes_bot








# Higher resolution converges on a slightly different answer...
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-6, reltol = 1e-6);
Gflow = Gflow_to_01_GMRES;
res_Gflow_v6a = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true)
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-3, reltol = 1e-3);
Gflow = Gflow_to_01_GMRES;
res_Gflow_v6a = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true)
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-12, reltol = 1e-12);
Gflow = Gflow_to_01_GMRES;
res_Gflow_v6a = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true)



end # END @testset "runtests_BiSSE_tree_n3" begin




