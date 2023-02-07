using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
#using PhyloBits.PNtypes					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames						# for DataFrame()
using DelimitedFiles				# for readdlm()

# List each PhyBEARS code file prefix here
using PhyloBits.TrUtils			# for e.g. numstxt_to_df()
using PhyloBits.TreeTable
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
#using PhyBEARS.TreeTable
using PhyBEARS.TreePass
#using PhyBEARS.TrUtils
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.ModelLikes # e.g. setup_DEC_SSE2
using PhyBEARS.Uppass

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_ClaSSE_tree_n8_DEC.jl")
"""
# 
# """
# # Run with:
# include("/GitHub/PhyBEARS.jl/test/runtests_ClaSSE_tree_n8_DEC.jl")
# """
# 
# @testset "Example" begin
# 	@test hello("runtests_ClaSSE_tree_n8_DEC.jl") == "runtests_ClaSSE_tree_n8_DEC.jl"
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
@testset "runtests_ClaSSE_tree_n8_DEC.jl -- Psychotria DEC" begin
# 
#######################################################
# DEMONSTRATES MATCHING BETWEEN DIVERSITREE, BIOGEOBEARS, AND JULIA
# ON HAWAIIAN PSYCHOTRIA, 16-STATE DEC MODEL
#
# Run with:
# source("/GitHub/PhyBEARS.jl/R_examples/compare_BGB_diversitree_DEC_v1.R")
# Truth:
DEC_lnL = -34.54313
DEC_R_result_branch_lnL = -67.6295
DEC_R_result_total_LnLs1 = -72.60212
DEC_R_result_total_LnLs1t = -71.48986
DEC_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -120.1545
#######################################################


# Island numbers (KOMH = 1234) in Rnodenums order:
island_nums = [3, 3, 2, 2, 3, 3, 2, 1, 1, 3, 4, 2, 1, 1, 1, 1, 1, 1, 2]

"""
# Load geog data file
lgdata_fn = "/GitHub/PhyBEARS.jl/data/Psychotria/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
julian_dput(geog_df)
"""

geog_df = DataFrame(AbstractVector[Any["P_mariniana_Kokee2", "P_mariniana_Oahu", "P_mariniana_MauiNui", "P_hawaiiensis_Makaopuhi", "P_wawraeDL7428", "P_kaduana_PuuKukuiAS", "P_mauiensis_PepeAS", "P_hawaiiensis_WaikamoiL1", "P_mauiensis_Eke", "P_fauriei2", "P_hathewayi_1", "P_kaduana_HawaiiLoa", "P_greenwelliae07", "P_greenwelliae907", "P_grandiflora_Kal2", "P_hobdyi_Kuia", "P_hexandra_K1", "P_hexandra_M", "P_hexandra_Oahu"], Any["1", "0", "0", "0", "1", "0", "0", "0", "0", "0", "0", "0", "1", "1", "1", "1", "1", "1", "0"], Any["0", "1", "0", "0", "0", "0", "0", "0", "0", "1", "1", "1", "0", "0", "0", "0", "0", "0", "1"], Any["0", "0", "1", "0", "0", "1", "1", "1", "1", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"], Any["0", "0", "0", "1", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"]], DataFrames.Index(Dict(:M => 4, :H => 5, :tipnames => 1, :K => 2, :O => 3), [:tipnames, :K, :O, :M, :H]))

# Psychotria tree
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")


in_params = (birthRate=0.3288164, deathRate=0.0, d_val=0.03505038, e_val=0.02832370, a_val=0.0, j_val=0.0)
bmo = construct_bmo();
bmo.est[bmo.rownames.=="birthRate"] .= in_params.birthRate
bmo.est[bmo.rownames.=="deathRate"] .= in_params.deathRate
bmo.est[bmo.rownames.=="d"] .= in_params.d_val
bmo.est[bmo.rownames.=="e"] .= in_params.e_val
bmo.est[bmo.rownames.=="a"] .= in_params.a_val
bmo.est[bmo.rownames.=="j"] .= in_params.j_val
bmo.type[bmo.rownames.=="j"] .= "fixed"
numareas = 4
n = 16            # 4 areas, 16 states

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
max_range_size = NaN # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

prtQi(inputs)
prtCi(inputs)

# Redo trdf 
trdf = prt(tr)
R_order = sort(trdf, :Rnodenums).nodeIndex

# Look at tip likelihoods:
inputs.res.likes_at_each_nodeIndex_branchTop[R_order]
sum.(inputs.res.likes_at_each_nodeIndex_branchTop[R_order])
size(inputs.res.likes_at_each_nodeIndex_branchTop)
numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])


root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
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
vfft(res.likes_at_each_nodeIndex_branchTop)
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)
vfft(res.likes_at_each_nodeIndex_branchTop)

@test round(DEC_lnL; digits=2) == round(bgb_lnL; digits=2)


Rnames(res)
vfft(res.likes_at_each_nodeIndex_branchTop)
vfft(res.normlikes_at_each_nodeIndex_branchTop)
vfft(res.likes_at_each_nodeIndex_branchBot)
vfft(res.normlikes_at_each_nodeIndex_branchBot)

sum.(res.likes_at_each_nodeIndex_branchTop)
log.(sum.(res.likes_at_each_nodeIndex_branchTop))
sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))
Julia_sum_lq_old = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
nonroot_nodes = get_nonrootnodes(tr)
sum_likes_internal_branch_tops = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
Julia_sum_lq = Julia_sum_lq_old + sum_likes_internal_branch_tops

# Does the total of the branch log-likelihoods (lq) match?
@test round(DEC_R_result_branch_lnL; digits=2) == round(Julia_sum_lq; digits=2)

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
@test round(DEC_R_result_total_LnLs1; digits=2) == round(Julia_total_lnLs1; digits=2)

# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
@test round(DEC_R_result_total_LnLs1t; digits=2) == round(Julia_total_lnLs1t; digits=2)

# Does the total of branch likelihoods (lq) + node likelihoods match R?
# 
# R: DEC_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda 
#    = sum(log(computed_likelihoods_at_each_node_x_lambda))
res.likes_at_each_nodeIndex_branchTop
log.(sum.(res.likes_at_each_nodeIndex_branchTop))

sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))

Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = DEC_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
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




#######################################################
# Check ancestral state estimates
#######################################################

@testset "Psychotria ancstates DEC v5" begin

fn = "/GitHub/PhyBEARS.jl/data/Psychotria_DEC_ancstates_nodes.txt"
R_ancstates_nodes = numstxt_to_df(fn)

# Has NAs; these are auto-convert to NaN by numstxt_to_df
fn = "/GitHub/PhyBEARS.jl/data/Psychotria_DEC_ancstates_corners.txt"
R_ancstates_corners = numstxt_to_df(fn)

trdf = prt(tr)
R_order = sort(trdf, :Rnodenums).nodeIndex

uppass_ancstates_v5!(res, trdf, p_Ds_v5, solver_options; use_Cijk_rates_t=false, min_branchlength=1.0e-6);
rn(res)
Julia_ancstates_nodes_v5 = vvdf(deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order]))
Julia_ancstates_corners_v5 = vvdf(deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order]))

@test all(get_max_df_diffs_byCol(R_ancstates_nodes, Julia_ancstates_nodes_v5) .< 0.0002)
@test all(get_max_df_diffs_byCol(R_ancstates_corners, Julia_ancstates_corners_v5) .< 0.0002)

end # END @testset "Psychotria ancstates DEC" begin



@testset "Psychotria ancstates DEC v7" begin

fn = "/GitHub/PhyBEARS.jl/data/Psychotria_DEC_ancstates_nodes.txt"
R_ancstates_nodes = numstxt_to_df(fn)

# Has NAs; these are auto-convert to NaN by numstxt_to_df
fn = "/GitHub/PhyBEARS.jl/data/Psychotria_DEC_ancstates_corners.txt"
R_ancstates_corners = numstxt_to_df(fn)

trdf = prt(tr)
R_order = sort(trdf, :Rnodenums).nodeIndex

uppass_ancstates_v7!(res, trdf, p_Ds_v5, solver_options; use_Cijk_rates_t=false, min_branchlength=1.0e-6);
rn(res)
Julia_ancstates_nodes_v5 = vvdf(deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order]))
Julia_ancstates_corners_v5 = vvdf(deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order]))

@test all(get_max_df_diffs_byCol(R_ancstates_nodes, Julia_ancstates_nodes_v5) .< 0.0002)
@test all(get_max_df_diffs_byCol(R_ancstates_corners, Julia_ancstates_corners_v5) .< 0.0002)

end # END @testset "Psychotria ancstates DEC" begin



end # END @testset "runtests_BiSSE_tree_n3" begin




