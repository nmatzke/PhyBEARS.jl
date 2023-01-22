


using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using PhyloBits.TrUtils
using DataFrames
using CSV

# List each PhyBEARS code file prefix here
using PhyBEARS
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.ModelLikes
using PhyBEARS.TreePass
using PhyBEARS.SSEs


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
# source("/GitHub/PhyBEARS.jl/test/bisse_ASR_v1.R")

# Truth:
# > LnLs1
# [1] -9.574440 -6.670978
# > LnLs1t
# [1] -7.464283 -6.670978


R_result_branch_lnL = -6.670978
R_result_total_LnLs1 = -9.574440
R_result_total_LnLs1t = -7.464283
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -12.06325
#######################################################


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
rn(p_Ds_v5.p_TFs)


# Change parameter inputs manually
inputs.p_Ds_v5.params.Cijk_vals[1] = 0.222222222
inputs.p_Ds_v5.params.Cijk_vals[2] = 0.222222222
inputs.p_Ds_v5.params.mu_vals[1] = 0.111111111
inputs.p_Ds_v5.params.mu_vals[2] = 0.05
inputs.p_Ds_v5.params.Qij_vals[1] = 0.1
inputs.p_Ds_v5.params.Qij_vals[2] = 0.15

# Update the subs, after updating the individual values manually
inputs.p_Ds_v5.p_TFs.Qij_vals_sub_i
inputs.p_Ds_v5.p_TFs.Qji_vals_sub_j
update_Qij_vals_subs!(p_Ds_v5)
inputs.p_Ds_v5.p_TFs.Qij_vals_sub_i
inputs.p_Ds_v5.p_TFs.Qji_vals_sub_j

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


# Parameters
prtQi(inputs)
prtQp(p_Ds_v5)
prtCi(inputs)
prtCp(p_Ds_v5)

inputs.p_Ds_v5.params.mu_vals
p_Ds_v5.sol_Es_v5(1.0)

# Do downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=construct_SolverOpt(), return_lnLs=true, max_iterations=10^10)





Rnames(res)
res.likes_at_each_nodeIndex_branchTop
res.normlikes_at_each_nodeIndex_branchTop
res.likes_at_each_nodeIndex_branchBot
res.normlikes_at_each_nodeIndex_branchBot

# OLD:
# sum.(res.likes_at_each_nodeIndex_branchTop)
# log.(sum.(res.likes_at_each_nodeIndex_branchTop))
# sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))
# Julia_sum_lq = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
# sum(res.lq_at_branchBot)

# NEW:
Julia_sum_lq_old = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
nonroot_nodes = TreePass.get_nonrootnodes_trdf(trdf)
sum_likes_internal_branch_tops = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
Julia_sum_lq_new = Julia_sum_lq_old + sum_likes_internal_branch_tops

# Does the total of the branch log-likelihoods (lq) match?
@test round(R_result_branch_lnL; digits=4) == round(Julia_sum_lq; digits=4)

# OLD:
# ttl_from_Julia = sum(res.lq_at_branchBot) +  sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))
# NEW: Corresponding to R SSE script, "res1" result, which assumes diversitree options:
# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
root_stateprobs = d_root_orig/sum(d_root_orig)
rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL

@test round(R_result_total_LnLs1; digits=3) == round(Julia_total_lnLs1; digits=3)

# NEW: Corresponding to R SSE script, "res1t" result, which assumes diversitree options:
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

# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
@test round(R_result_total_LnLs1t; digits=4) == round(Julia_total_lnLs1t; digits=4)

# Does the total of branch likelihoods (lq) + node likelihoods match R?
# 
# R: R_result_sum_log_computed_likelihoods_at_each_node_x_lambda 
#    = sum(log(computed_likelihoods_at_each_node_x_lambda))
res.likes_at_each_nodeIndex_branchTop
log.(sum.(res.likes_at_each_nodeIndex_branchTop))

sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop)))


# This corresponds to:
# Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
# R_sum_lq_nodes = R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
#  sum(log(rowSums(EsDs[,3:4]))) + sum(attr(res1t,"intermediates")$lq)
# ...but is double-counting lnLs
Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
R_sum_lq_nodes = R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
@test round(Julia_sum_lq_nodes; digits=3) == round(R_sum_lq_nodes; digits=3)


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





#######################################################
# CHECK ancestral state inferences
# 
#######################################################
trtable = prt(tr)
tree_height = trtable.node_age[tr.root]
trtable = prt(tr)
ancnode = tr.root
decnode = trtable.leftNodeIndex[ancnode]
sisnode = trtable.rightNodeIndex[ancnode]
anctime = tree_height - trtable.node_age[ancnode]
dectime = tree_height - trtable.node_age[decnode]

# Solve the Ds, single branch
u0 = res.likes_at_each_nodeIndex_branchTop[tr.root]
u0 ./ sum(u0)
u0 = res.normlikes_at_each_nodeIndex_branchTop[tr.root]
u0 = [0.5, 0.5]

# Calculate bisse ancestral states:
ancnode = tr.root
lnode = trdf.leftNodeIndex[ancnode]
rnode = trdf.rightNodeIndex[ancnode]

left_likes = res.normlikes_at_each_nodeIndex_branchBot[lnode]
right_likes = res.normlikes_at_each_nodeIndex_branchBot[rnode]





#######################################################
# Estimate ancestral states
#######################################################

solver_options.solver = CVODE_BDF{:Newton, :GMRES, Nothing, Nothing}(0, 0, 0, false, 10, 5, 7, 3, 10, nothing, nothing, 0)
#solver_options.solver = Tsit5()
solver_options.solver = Vern9()
solver_options.abstol = 1e-6
solver_options.reltol = 1e-6
solver_options.save_everystep = false


include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
tspan = (anctime, dectime)

u0 = right_likes
prob_Ds_v5 = DifferentialEquations.ODEProblem(calcDs_4states2C, u0, tspan, p_Ds_v5);
sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v5(anctime)
sol_Ds_v5(anctime) ./ sum(sol_Ds_v5(anctime))
sol_Ds_v5(dectime)
lbranch_top = sol_Ds_v5(dectime) ./ sum(sol_Ds_v5(dectime))

uppass_likes = lbranch_top .* res.normlikes_at_each_nodeIndex_branchTop[lnode]
uppass_likes ./ sum(uppass_likes)

# Diversitree: asr.marginal
# 0.005395545 0.9946044545

# calcDs_4states2, calcDs_4states2A. 2E, 2F, 2G <- all the same on BiSSE
# 0.004611895520264392
# 0.9953881044797356

# calcDs_4states2B # <- closest!
# 0.005216130223249679
# 0.9947838697767504

# calcDs_4states2C # <-- same
# 0.005216130223249679
# 0.9947838697767504


include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")

ctable1 = prtCp(p_Ds_v5)
make_ctable_single_events(ctable1)


u0 = left_likes
prob_Ds_v5 = DifferentialEquations.ODEProblem(calcDs_4states2G, u0, tspan, p_Ds_v5);
sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v5(anctime)
sol_Ds_v5(anctime) ./ sum(sol_Ds_v5(anctime))
sol_Ds_v5(dectime)
rbranch_top = sol_Ds_v5(dectime) ./ sum(sol_Ds_v5(dectime))

uppass_likes = rbranch_top .* res.normlikes_at_each_nodeIndex_branchTop[rnode]
asr_at_node7 = uppass_likes ./ sum(uppass_likes)

# calcDs_4states2B, calcDs_4states2C   # <- closest! (and same)
#  0.9996842305583383
#  0.00031576944166156017



# Diversitree: asr.marginal
# 0.999620338 0.0003796623

diversitree_bisse_Rnode7_01 = [0.999620338, 0.0003796623]


asr_at_node7 .== diversitree_bisse_Rnode7_01

@test round(asr_at_node7[1]; digits=3) .== round(diversitree_bisse_Rnode7_01[1]; digits=3)
@test round(asr_at_node7[2]; digits=3) .== round(diversitree_bisse_Rnode7_01[2]; digits=3)

asr_at_node7[1] - diversitree_bisse_Rnode7_01[1]
asr_at_node7[2] - diversitree_bisse_Rnode7_01[2]


include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
u0 = left_likes
prob_Ds_v5 = DifferentialEquations.ODEProblem(calcDs_4states2F, u0, tspan, p_Ds_v5);
sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v5(anctime)
sol_Ds_v5(anctime) ./ sum(sol_Ds_v5(anctime))
sol_Ds_v5(dectime)
rbranch_top = sol_Ds_v5(dectime) ./ sum(sol_Ds_v5(dectime))

uppass_likes = rbranch_top .* res.normlikes_at_each_nodeIndex_branchTop[rnode]
uppass_likes ./ sum(uppass_likes)
# 0.9996842305583383
# 0.00031576944166156017

# 0.9996842305583383
# 0.00031576944166156017


#######################################################
# These all match, but this doesn't prove the 
# cladogenesis part works right when mu > 0 and 
# non-sympatry clado > 0
#######################################################

# v5 algorithm
include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
R_order = sort(trdf, :Rnodenums).nodeIndex
p_Ds_v7 = p_Ds_v5;
tmpnode = 6
u = res.uppass_probs_at_each_nodeIndex_branchBot[tmpnode]
p = p_Ds_v7;

p.params.Qij_vals
prtQp(p)
p.p_indices.Qarray_ivals
p.p_indices.Qarray_jvals

p.p_TFs.Qij_vals_sub_i
p.p_TFs.Qji_vals_sub_j


tree_age = trdf.node_age[tr.root]
t_end = tree_age - trdf.node_age[tmpnode]
t_start = tree_age - trdf.node_age[trdf.rightNodeIndex[tr.root]]
t = t_start
i = 2
j = 2
it = 1
du = repeat([0.0], 2)

# Individual calculations
res1 = calcDs_4states2D_print(du,u,p,t)
res2 = parameterized_ClaSSE_Ds_v7_simd_sums_2D_FWD_print(du,u,p,t)

@test all(res1 .== res2)

# Check that these make sense!
p_Ds_v5.p_TFs.Qij_vals_sub_i
p_Ds_v5.p_TFs.Qji_vals_sub_j
update_Qij_vals_subs!(p_Ds_v5)
p_Ds_v5.p_TFs.Qij_vals_sub_i
p_Ds_v5.p_TFs.Qji_vals_sub_j


uppass_ancstates_v5!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
rn(res)
res.uppass_probs_at_each_nodeIndex_branchBot[R_order,:]
res.uppass_probs_at_each_nodeIndex_branchTop[R_order,:]
res.anc_estimates_at_each_nodeIndex_branchBot[R_order,:]
res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:]

# Julia:
#  [0.43035780124527134, 0.5696421987547287]
#  [0.005216130223249679, 0.9947838697767503]
#  [0.9996911604384737, 0.0003088395615262636]

# R diversitrree
# t(st2)
#             [,1]         [,2]
# [1,] 0.430357148 0.5696428522  # <- matches res1 or res1t
# [2,] 0.005145312 0.9948546878  # <- closest of the below
# [3,] 0.999681941 0.0003180585  # <- closest of the below

R_bisse_anc_estimates = [0.430357148 0.5696428522; # <- matches res1 or res1t
0.005145312 0.9948546878;  # <- closest of the below
0.999681941 0.0003180585]

Julia_bisse_anc_estimates = res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:][5:7,:]
Julia_bisse_anc_estimates = mapreduce(permutedims, vcat, Julia_bisse_anc_estimates)
round.(R_bisse_anc_estimates; digits=3) .== round.(Julia_bisse_anc_estimates; digits=3)
@test all(round.(R_bisse_anc_estimates; digits=3) .== round.(Julia_bisse_anc_estimates; digits=3))

Julia_bisse_anc_estimates .- R_bisse_anc_estimates



# v7 algorithm
include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
R_order = sort(trdf, :Rnodenums).nodeIndex
p_Ds_v7 = p_Ds_v5;
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
rn(res)
res.uppass_probs_at_each_nodeIndex_branchBot[R_order,:]
res.uppass_probs_at_each_nodeIndex_branchTop[R_order,:]
res.anc_estimates_at_each_nodeIndex_branchBot[R_order,:]
res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:]

# Julia:
#  [0.43035780124527134, 0.5696421987547287]
#  [0.005216130223249679, 0.9947838697767503]
#  [0.9996911604384737, 0.0003088395615262636]

# R diversitrree
# t(st2)
#             [,1]         [,2]
# [1,] 0.430357148 0.5696428522  # <- matches res1 or res1t
# [2,] 0.005145312 0.9948546878  # <- closest of the below
# [3,] 0.999681941 0.0003180585  # <- closest of the below

R_bisse_anc_estimates = [0.430357148 0.5696428522; # <- matches res1 or res1t
0.005145312 0.9948546878;  # <- closest of the below
0.999681941 0.0003180585]

Julia_bisse_anc_estimates = res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:][5:7,:]
Julia_bisse_anc_estimates = mapreduce(permutedims, vcat, Julia_bisse_anc_estimates)
round.(R_bisse_anc_estimates; digits=3) .== round.(Julia_bisse_anc_estimates; digits=3)
@test all(round.(R_bisse_anc_estimates; digits=3) .== round.(Julia_bisse_anc_estimates; digits=3))

Julia_bisse_anc_estimates .- R_bisse_anc_estimates


@benchmark uppass_ancstates_v5!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)



@benchmark uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)





# Root state probabilities
res.normlikes_at_each_nodeIndex_branchTop[tr.root]
# 0.43035780124527134
# 0.5696421987547287
normlikes_root_states1 = 0.43035780124527134
normlikes_root_states2 = 0.5696421987547287

# attr(res1t,"intermediates")$root.p
# [1] 0.4303571 0.5696429
res1t_root_states1 = 0.4303571
res1t_root_states2 = 0.5696429


sqrt.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]) ./ sum(sqrt.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]))
# 2-element Vector{Float64}:
#  0.4650083587227312
#  0.5349916412772688
normlikes_sqrt_root_states1 = 0.4650083587227312
normlikes_sqrt_root_states2 = 0.5349916412772688

 
# Uppass through whole tree
include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
current_nodeIndex = tr.root
p_Ds_v7 = p_Ds_v5;
nodeOp_Cmat_uppass_v7!(res, current_nodeIndex, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

res.anc_estimates_at_each_nodeIndex_branchTop
res.uppass_probs_at_each_nodeIndex_branchTop

# end # END @testset "runtests_BiSSE_tree_n3" begin




