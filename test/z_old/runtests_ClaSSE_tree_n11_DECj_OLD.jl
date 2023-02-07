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
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreeTable
using PhyBEARS.TreePass
using PhyBEARS.TrUtils
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.ModelLikes

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_ClaSSE_tree_n11_DECj.jl")
"""
# 
# """
# # Run with:
# include("/GitHub/PhyBEARS.jl/test/runtests_ClaSSE_tree_n11_DECj.jl")
# """# 
# @testset "Example" begin
# 	@test hello("runtests_ClaSSE_tree_n11_DECj.jl") == "runtests_ClaSSE_tree_n11_DECj.jl"
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
@testset "runtests_ClaSSE_tree_n11_DECj.jl" begin
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
DECj_R_result_branch_lnL = -55.37332
DECj_R_result_total_LnLs1 = -58.83758
DECj_R_result_total_LnLs1t = -57.72533
DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -96.34151

#######################################################


#include("/GitHub/PhyBEARS.jl/src/TreePass.jl")
#import .TreePass

# Repeat calculation in Julia
#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes

# Read in the island areas using Parsers.jl (v10) rather than just manual input (v9)
#include("/GitHub/PhyBEARS.jl/notes/Parsers.jl")
#import .Parsers

# Island numbers (KOMH = 1234) in Rnodenums order:
#island_nums = [3, 3, 2, 2, 3, 3, 2, 1, 1, 3, 4, 2, 1, 1, 1, 1, 1, 1, 2]
lgdata_fn = "/GitHub/PhyBEARS.jl/data/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Psychotria tree
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
in_params = (birthRate=0.3288164, deathRate=0.0, d_val=1e-12, e_val=1e-12, a_val=0.0, j_val=0.1142057)
numareas = 4
n = 16            # 4 areas, 16 states

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, in_params=in_params)
(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs

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

end # END @testset "runtests_BiSSE_tree_n3" begin




