


using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using PhyloBits.TrUtils
using PhyloBits.TreeTable # for get_tree_height
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
using PhyBEARS.TimeDep
using PhyBEARS.Uppass


"""
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/bisse_ASR+Jv12_v1.jl")
"""


# 
# """
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

wd = "/GitHub/PhyBEARS.jl/test/"
cd(wd)


# These are archiving the Julia results;
# R's diversitree does not do BiSSE+J
R_rootstates_lnL = -2.9003085430341864
R_bgb_lnL = -1.6056875884138417

R_result_branch_lnL = -7.437234631321879
R_result_total_LnLs1 = -10.337543174356066
R_result_total_LnLs1t = -10.337543174356066
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -13.589787226916526
#######################################################

numareas = 2
tr = readTopology("((sp4:0.6248637277,sp5:0.6248637277):6.489662918,(sp6:0.1274213816,sp7:0.1274213816):6.987105264);")
oldest_possible_age = 2*get_tree_height(tr)
geog_df = DataFrame(tipnames=["sp4","sp5","sp6","sp7"],A=[1,1,0,0],B=[0,0,1,1]);

in_params = (birthRate=0.222222222, deathRate=0.111111111, d_val=0.0, e_val=0.0, a_val=0.1, j_val=0.25)
# pars <- c(0.222222222, 0.222222222, 0.111111111, 0.05, 0.1, 0.15)
bmo = construct_BioGeoBEARS_model_object();
bmo_rows = get_bmo_rows(bmo);
bmo.est[bmo.rownames.=="d"] .= 0.0;
bmo.est[bmo.rownames.=="e"] .= 0.0;
bmo.est[bmo.rownames.=="a"] .= in_params.a_val;
bmo.est[bmo.rownames.=="j"] .= in_params.j_val;
bmo.est[bmo.rownames.=="birthRate"] .= in_params.birthRate;
bmo.est[bmo.rownames.=="deathRate"] .= in_params.deathRate;

bmo.type[bmo.rownames.=="d"] .= "fixed"
bmo.type[bmo.rownames.=="e"] .= "fixed"
bmo.type[bmo.rownames.=="a"] .= "free"
bmo.type[bmo.rownames.=="j"] .= "free"
bmo.type[bmo.rownames.=="birthRate"] .= "free"
bmo.type[bmo.rownames.=="deathRate"] .= "free"



bmo.est .= bmo_updater_v2(bmo, bmo_rows);
bmo

# phy$tip.state
# sp4 sp5 sp6 sp7 
#   1   1   0   0 

numstates = 2
n = 2

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
#inputs = ModelLikes.setup_MuSSE_biogeo(numstates, tr; root_age_mult=1.5, in_params=in_params);
root_age_mult=1.5; max_range_size=1; include_null_range=false; manual_states_list=NaN; area_names=LETTERS(1:numareas)
inputs = ModelLikes.setup_DEC_SSE2(numstates, tr, geog_df; root_age_mult=1.5, max_range_size=1, include_null_range=false, bmo=bmo, manual_states_list=NaN, area_names=LETTERS(1:numareas));
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;
prtQi(inputs)
prtCi(inputs)

# Change parameter inputs manually
inputs.p_Ds_v5.params.mu_vals[1] = 0.111111111
inputs.p_Ds_v5.params.mu_vals[2] = 0.05
inputs.p_Ds_v5.params.Qij_vals[1] = 0.1
#inputs.p_Ds_v5.params.Qij_vals[2] = 0.15
inputs.p_Ds_v5.params.Qij_vals[2] = 0.1

# Update the subs, after updating the individual values manually
inputs.p_Ds_v5.p_TFs.Qij_vals_sub_i
inputs.p_Ds_v5.p_TFs.Qji_vals_sub_j
update_Qij_vals_subs!(p_Ds_v5)			# update
inputs.p_Ds_v5.p_TFs.Qij_vals_sub_i
inputs.p_Ds_v5.p_TFs.Qji_vals_sub_j

# Tip likelihoods
inputs.res.likes_at_each_nodeIndex_branchTop
inputs.res.normlikes_at_each_nodeIndex_branchTop
# res.likes_at_each_nodeIndex_branchTop[1] = [0.0, 1.0]
# res.likes_at_each_nodeIndex_branchTop[2] = [0.0, 1.0]
# res.likes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0]
# res.likes_at_each_nodeIndex_branchTop[5] = [1.0, 0.0]
# res.normlikes_at_each_nodeIndex_branchTop[1] = [0.0, 1.0]
# res.normlikes_at_each_nodeIndex_branchTop[2] = [0.0, 1.0]
# res.normlikes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0]
# res.normlikes_at_each_nodeIndex_branchTop[5] = [1.0, 0.0]
trdf = inputs.trdf
p_Ds_v5 = inputs.p_Ds_v5;
root_age = maximum(trdf[!, :node_age])

prtQi(inputs)
prtCi(inputs)


# Put in some fake files to modify this BiSSE+J run
files.times_fn = "bisse_ASR+Jv12_times.txt"
files.distances_fn = "bisse_ASR+Jv12_distances.txt"
files.area_of_areas_fn = "bisse_ASR+Jv12_area_of_areas.txt"

# Construct interpolators, times-at-which-to-interpolate QC
p = p_Ds_v5;
interpolators = files_to_interpolators(files, setup.numareas, setup.states_list, setup.v_rows, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals, trdf; oldest_possible_age=oldest_possible_age);

ts = [0.0, 1.0, 3.0]
interpolators.distances_interpolator(ts)
interpolators.area_of_areas_interpolator(ts)




p_Es_v12 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo, interpolators=interpolators);

# Add Q, C interpolators
p_Es_v12 = p = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v12, p_Es_v12.interpolators.times_for_SSE_interpolators);


# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is a linear interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);

p_Ds_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);



# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);


# Solve the Es
prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Es_v5(ts)
sol_Es_v7(ts)
sol_Es_v12(ts)

@test all(sol_Es_v7(ts).u[2] .== sol_Es_v5(ts).u[2])
@test all( (sol_Es_v7(ts).u[2] .- sol_Es_v12(ts).u[2]) .< 1e-8)
#@test all(sol_Es_v7(ts).u[3] .== sol_Es_v5(ts).u[3])
@test all( (sol_Es_v7(ts).u[3] .- sol_Es_v12(ts).u[3]) .< 1e-8)

p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


nodenum = 1
u = res.likes_at_each_nodeIndex_branchTop[nodenum]
du = repeat([0.0], length(u))


# We have to manually update one of the "a" values, to match the BiSSE example,
# as the default updater only allows 1 "a"
p_Es_v12.params.Qij_vals_t .= p_Es_v12.params.Qij_vals;

p = p_Es_v12;

t = 0.0
v12res = parameterized_ClaSSE_Es_v12_simd_sums_print(du,u,p,t)

du = repeat([0.0], length(u))
p = p_Ds_v5;
v5res = parameterized_ClaSSE_Es_v5_print(du,u,p,t)

v12res .- v5res
v12res .== v5res

#######################################################
# THESE MATCH, *IF*:

# a is equal across states, or
# #p.params.Qij_vals_t .= p.interpolators.Q_vals_interpolator(t)
# is commented out of 
#
# parameterized_ClaSSE_Es_v12_simd_sums_print
#######################################################



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
sum_Crates = sum_Cijk_rates_by_i(inputs.p_Ds_v5.params.Cijk_vals, inputs.p_Ds_v5.p_indices.Carray_ivals, inputs.p_Ds_v5.n)
d_root = d_root_orig ./ sum(root_stateprobs .* sum_Crates .* (1 .- e_root).^2)
rootstates_lnL = log(sum(root_stateprobs .* d_root))

# The above all works out to [0,1] for the Yule model with q01=q02=0.0
Julia_total_lnLs1t = Julia_sum_lq + rootstates_lnL

# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
# THE BELOW WORKS FOR BISSE, BUT NOT BISSE + J
# @test round(R_result_total_LnLs1t; digits=4) == round(Julia_total_lnLs1t; digits=4)

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
#solver_options.solver = Vern9()
solver_options.abstol = 1e-6
solver_options.reltol = 1e-6
solver_options.save_everystep = true  # This has to be true to calculate Es, as this is a continuous interpolator


#include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
tspan = (anctime, dectime)

u0 = right_likes
prob_Ds_v5 = DifferentialEquations.ODEProblem(calcDs_4states2D, u0, tspan, p_Ds_v7);
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

# calcDs_4states2D # <-- same
# 0.005216130223249679
# 0.9947838697767504


#include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")

ctable1 = prtCp(p_Ds_v5)
make_ctable_single_events(ctable1)


u0 = left_likes
prob_Ds_v5 = DifferentialEquations.ODEProblem(calcDs_4states2D, u0, tspan, p_Ds_v5);
sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v5(anctime)
sol_Ds_v5(anctime) ./ sum(sol_Ds_v5(anctime))
sol_Ds_v5(dectime)
rbranch_top = sol_Ds_v5(dectime) ./ sum(sol_Ds_v5(dectime))

uppass_likes = rbranch_top .* res.normlikes_at_each_nodeIndex_branchTop[rnode]
asr_at_node7 = uppass_likes ./ sum(uppass_likes)

# calcDs_4states2B, calcDs_4states2D   # <- closest! (and same)
#  0.9996842305583383
#  0.00031576944166156017



# Diversitree: asr.marginal
# 0.999620338 0.0003796623

# Diversitree biSSE, with unequal "a"
# diversitree_bisse_Rnode7_01 = [0.999620338, 0.0003796623]
julia_bisseJ_Rnode7_01 = [0.006446111346864297, 0.9935538886531358]

asr_at_node7 .== julia_bisseJ_Rnode7_01

#@test round(asr_at_node7[1]; digits=3) .== round(diversitree_bisse_Rnode7_01[1]; digits=3)
#@test round(asr_at_node7[2]; digits=3) .== round(diversitree_bisse_Rnode7_01[2]; digits=3)

asr_at_node7[1] - julia_bisseJ_Rnode7_01[1]
asr_at_node7[2] - julia_bisseJ_Rnode7_01[2]


#include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
u0 = left_likes
prob_Ds_v5 = DifferentialEquations.ODEProblem(calcDs_4states2D, u0, tspan, p_Ds_v5);
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
#include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
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
v5_anc_branchBot = deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order,:])
v5_anc_branchTop = deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:])

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

#R_bisse_anc_estimates = [0.430357148 0.5696428522; # <- matches res1 or res1t
#0.005145312 0.9948546878;  # <- closest of the below
#0.999681941 0.0003180585]

# Trivial
R_bisse_anc_estimates = res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:][5:7,]
R_bisse_anc_estimates = mapreduce(permutedims, vcat, R_bisse_anc_estimates)

Julia_bisse_anc_estimates = res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:][5:7,:]
Julia_bisse_anc_estimates = mapreduce(permutedims, vcat, Julia_bisse_anc_estimates)
round.(R_bisse_anc_estimates; digits=3) .== round.(Julia_bisse_anc_estimates; digits=3)

@test all(round.(R_bisse_anc_estimates; digits=3) .== round.(Julia_bisse_anc_estimates; digits=3))

Julia_bisse_anc_estimates .- R_bisse_anc_estimates



# v7 algorithm
#include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
R_order = sort(trdf, :Rnodenums).nodeIndex
p_Ds_v7 = p_Ds_v5;
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
rn(res)
res.uppass_probs_at_each_nodeIndex_branchBot[R_order,:]
res.uppass_probs_at_each_nodeIndex_branchTop[R_order,:]
v7_anc_branchBot = deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order,:])
v7_anc_branchTop = deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:])

#@test all( (v7_anc_branchBot[5] .- v5_anc_branchBot[5]) .< 1e-6) 
@test all( (v7_anc_branchBot[6] .- v5_anc_branchBot[6]) .< 1e-6) 
@test all( (v7_anc_branchBot[7] .- v5_anc_branchBot[7]) .< 1e-6) 

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

Julia_bisseJ_anc_estimates = res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:][5:7,]
Julia_bisseJ_anc_estimates = mapreduce(permutedims, vcat, Julia_bisseJ_anc_estimates)
# OLD:
# Julia_bisseJ_anc_estimates = [0.4761851826108468 0.5238148173891531; 0.9587285236937566 0.04127147630624328; 0.005800480504768555 0.9941995194952314]
Julia_bisseJ_anc_estimates = [0.476185 0.523815; 0.955845 0.0441553; 0.00605841 0.993942]



Julia_bisse_anc_estimates = res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:][5:7,:]
Julia_bisse_anc_estimates = mapreduce(permutedims, vcat, Julia_bisse_anc_estimates)
round.(Julia_bisseJ_anc_estimates; digits=3) .== round.(Julia_bisse_anc_estimates; digits=3)
@test all(round.(Julia_bisseJ_anc_estimates; digits=3) .== round.(Julia_bisse_anc_estimates; digits=3))

Julia_bisse_anc_estimates .- Julia_bisseJ_anc_estimates









# v12 algorithm
#include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
R_order = sort(trdf, :Rnodenums).nodeIndex
uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)
rn(res)
res.uppass_probs_at_each_nodeIndex_branchBot[R_order,:]
res.uppass_probs_at_each_nodeIndex_branchTop[R_order,:]
v12_anc_branchBot = deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order,:])
v12_anc_branchTop = deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:])

#@test all( (v7_anc_branchBot[5] .- v12_anc_branchBot[5]) .< 1e-6) 
@test all( (v7_anc_branchBot[6] .- v12_anc_branchBot[6]) .< 1e-6) 
@test all( (v7_anc_branchBot[7] .- v12_anc_branchBot[7]) .< 1e-6) 

#@test all( (v5_anc_branchBot[5] .- v12_anc_branchBot[5]) .< 1e-6) 
@test all( (v5_anc_branchBot[6] .- v12_anc_branchBot[6]) .< 1e-6) 
@test all( (v5_anc_branchBot[7] .- v12_anc_branchBot[7]) .< 1e-6) 


#@benchmark uppass_ancstates_v5!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

#@benchmark uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

#@benchmark uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)

# julia> @benchmark uppass_ancstates_v5!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
# BenchmarkTools.Trial: 1335 samples with 1 evaluation.
#  Range (min … max):  2.585 ms … 63.020 ms  ┊ GC (min … max): 0.00% … 74.05%
#  Time  (median):     2.758 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   3.740 ms ±  4.008 ms  ┊ GC (mean ± σ):  9.66% ±  8.99%
# 
#   █▂         ▁▁                                               
#   ███▇▇▇▆▆▇▇▆██▇▆▁▅▄▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▆ ▇
#   2.59 ms      Histogram: log(frequency) by time     25.2 ms <
# 
#  Memory estimate: 2.53 MiB, allocs estimate: 40329.
# 
# julia> @benchmark uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
# BenchmarkTools.Trial: 2357 samples with 1 evaluation.
#  Range (min … max):  1.455 ms … 72.653 ms  ┊ GC (min … max): 0.00% … 50.59%
#  Time  (median):     1.550 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   2.116 ms ±  2.962 ms  ┊ GC (mean ± σ):  5.30% ±  3.95%
# 
#   ██▆▄▂▁                                   ▁▂▂▁▁▁            ▁
#   █████████▆▇▇▅▆▆▇▇▅▅▅▆▆▆▆▆▅▇▃▅▆▆▆▆▄▆▅▃▆▆▆▇██████▆▇▃▅▄▄▄▃▁▃▅ █
#   1.46 ms      Histogram: log(frequency) by time     5.03 ms <
# 
#  Memory estimate: 739.73 KiB, allocs estimate: 9054.
# 
# julia> @benchmark uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)
# BenchmarkTools.Trial: 1529 samples with 1 evaluation.
#  Range (min … max):  2.364 ms … 45.462 ms  ┊ GC (min … max): 0.00% … 72.03%
#  Time  (median):     2.500 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   3.266 ms ±  3.362 ms  ┊ GC (mean ± σ):  8.26% ±  7.91%
# 
#   █▂         ▂                                                
#   █████▆▇▇▇▅▇█▇▆▆▅▃▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▅ ▇
#   2.36 ms      Histogram: log(frequency) by time       24 ms <
# 
#  Memory estimate: 2.16 MiB, allocs estimate: 22494.

#######################################################
# These match the v5 and v7 calculations, because u=0, x=0
#######################################################

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)
res.anc_estimates_at_each_nodeIndex_branchBot[R_order,]
res.anc_estimates_at_each_nodeIndex_branchTop[R_order,]

# (0.008, 2, -7.437234631321882, -2.900308543034186, -10.337543174356068, -1.6056875884138435)
# 
# julia> res.anc_estimates_at_each_nodeIndex_branchTop
# 7-element Vector{Vector{Float64}}:
#  [1.0, 0.0]
#  [1.0, 0.0]
#  [0.9587283714306406, 0.04127162856935951]
#  [0.0, 1.0]
#  [0.0, 1.0]
#  [0.005800402112135891, 0.9941995978878642]
#  [0.4761851826108467, 0.5238148173891533]
# 
# julia> res.anc_estimates_at_each_nodeIndex_branchTop[R_order,]
# 7-element Vector{Vector{Float64}}:
#  [1.0, 0.0]
#  [1.0, 0.0]
#  [0.0, 1.0]
#  [0.0, 1.0]
#  [0.4761851826108467, 0.5238148173891533]
#  [0.9587283714306406, 0.04127162856935951]
#  [0.005800402112135891, 0.9941995978878642]



#######################################################
# Change the parameters to be distance-dependent and area-dependent
#######################################################
# For comparison (won't be the same):

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
v7_ancstates_bots_v0 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order,])
v7_ancstates_tops_v0 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order,])

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

p_Es_v12_archive = deepcopy(p_Es_v12);
p_Ds_v12_archive = deepcopy(p_Ds_v12);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12_archive, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)
v12_ancstates_bots_v0 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order,])
v12_ancstates_tops_v0 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order,])

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12_archive, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


res_v0 = deepcopy(res)

prtQp(p_Ds_v5)
prtQp(p_Ds_v7)
prtQp(p_Ds_v12)

prtCp(p_Ds_v5)
prtCp(p_Ds_v7)
prtCp(p_Ds_v12)

# julia> prtCp(p_Ds_v5)
# 4×10 DataFrame
#  Row │ event   i      j      k      pair   wt        prob      rate       val        rates_t   
#      │ String  Int64  Int64  Int64  Int64  Float64   Float64   Float64    Float64    Float64   
# ─────┼─────────────────────────────────────────────────────────────────────────────────────────
#    1 │ y           1      1      1      1  0.916667  0.647059  0.143791   0.143791   0.143791
#    2 │ j           1      1      2      2  0.5       0.352941  0.0784314  0.0784314  0.0784314
#    3 │ j           2      2      1      2  0.5       0.352941  0.0784314  0.0784314  0.0784314
#    4 │ y           2      2      2      1  0.916667  0.647059  0.143791   0.143791   0.143791
# 
# julia> prtCp(p_Ds_v7)
# 4×10 DataFrame
#  Row │ event   i      j      k      pair   wt        prob      rate       val        rates_t   
#      │ String  Int64  Int64  Int64  Int64  Float64   Float64   Float64    Float64    Float64   
# ─────┼─────────────────────────────────────────────────────────────────────────────────────────
#    1 │ y           1      1      1      1  0.916667  0.647059  0.143791   0.143791   0.143791
#    2 │ j           1      1      2      2  0.5       0.352941  0.0784314  0.0784314  0.0784314
#    3 │ j           2      2      1      2  0.5       0.352941  0.0784314  0.0784314  0.0784314
#    4 │ y           2      2      2      1  0.916667  0.647059  0.143791   0.143791   0.143791
# 
# julia> prtCp(p_Ds_v12)
# 4×10 DataFrame
#  Row │ event   i      j      k      pair   wt        prob      rate       val        rates_t   
#      │ String  Int64  Int64  Int64  Int64  Float64   Float64   Float64    Float64    Float64   
# ─────┼─────────────────────────────────────────────────────────────────────────────────────────
#    1 │ y           1      1      1      1  0.916667  0.647059  0.143791   0.143791   0.143791
#    2 │ j           1      1      2      2  0.5       0.352941  0.0784314  0.0784314  0.0784314
#    3 │ j           2      2      1      2  0.5       0.352941  0.0784314  0.0784314  0.0784314
#    4 │ y           2      2      2      1  0.916667  0.647059  0.143791   0.143791   0.143791


# Update v12
# Different from v12, because different u

p_Es_v12.bmo.est[bmo.rownames.=="u"] .= 1.0;
p_Es_v12.bmo.est[bmo.rownames.=="x"] .= 0.0;

p_Es_v12.bmo.est .= bmo_updater_v2(p_Es_v12.bmo, p_Es_v12.setup.bmo_rows);
p_Es_v12.bmo.est

p_Es_v12.params.mu_vals
p_Es_v12.params.mu_vals_t
p_Ds_v5_updater_v1!(p_Es_v12, p_Es_v12; check_if_free_params_in_mat=true, printlevel=0);
p_Es_v12.params.mu_vals
p_Es_v12.params.mu_vals_t

# Add the mu[2]=0.05 back in, manually
p_Es_v12.params.mu_vals[2] = 0.05

# Add Q, C interpolators
p_Es_v12 = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v12, p_Es_v12.interpolators.times_for_SSE_interpolators);
# Solve the Es
prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);


(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)
v12_ancstates_bots_v2 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order,])
v12_ancstates_tops_v2 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order,])





# Different from v12, because different x

p_Es_v12.bmo.est[bmo.rownames.=="u"] .= 0.0;
p_Es_v12.bmo.est[bmo.rownames.=="x"] .= -1.0;

p_Es_v12.bmo.est .= bmo_updater_v2(p_Es_v12.bmo, p_Es_v12.setup.bmo_rows);
p_Es_v12.bmo.est

p_Ds_v5_updater_v1!(p_Es_v12, p_Es_v12; check_if_free_params_in_mat=true, printlevel=0);

# Make sure the base deathRate stays different
p_Es_v12.params.mu_vals
# Add the mu[2]=0.05 back in, manually
p_Es_v12.params.mu_vals[2] = 0.05

# Add Q, C interpolators
p_Es_v12 = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v12, p_Es_v12.interpolators.times_for_SSE_interpolators);
# Solve the Es
prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);


(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)
v12_ancstates_bots_v3 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order,])
v12_ancstates_tops_v3 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order,])






# Different from v12, because different u & x
p_Es_v12.bmo.est[bmo.rownames.=="u"] .= 1.0;
p_Es_v12.bmo.est[bmo.rownames.=="x"] .= -1.0;

p_Es_v12.bmo.est .= bmo_updater_v2(p_Es_v12.bmo, p_Es_v12.setup.bmo_rows);
p_Es_v12.bmo.est

p_Ds_v5_updater_v1!(p_Es_v12, p_Es_v12; check_if_free_params_in_mat=true, printlevel=0);
# Add the mu[2]=0.05 back in, manually
p_Es_v12.params.mu_vals[2] = 0.05

# Add Q, C interpolators
p_Es_v12 = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v12, p_Es_v12.interpolators.times_for_SSE_interpolators);
# Solve the Es
prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);


(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)
v12_ancstates_bots_v4 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order,])
v12_ancstates_tops_v4 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order,])







# Back to matching

p_Es_v12.bmo.est[bmo.rownames.=="u"] .= 0.0;
p_Es_v12.bmo.est[bmo.rownames.=="x"] .= 0.0;

p_Es_v12.bmo.est .= bmo_updater_v2(p_Es_v12.bmo, p_Es_v12.setup.bmo_rows);
p_Es_v12.bmo.est

p_Ds_v5_updater_v1!(p_Es_v12, p_Es_v12; check_if_free_params_in_mat=true, printlevel=0);
# Add the mu[2]=0.05 back in, manually
p_Es_v12.params.mu_vals[2] = 0.05

# Add Q, C interpolators
p_Es_v12 = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v12, p_Es_v12.interpolators.times_for_SSE_interpolators);
# Solve the Es
prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

rn(p_Ds_v12.interpolators)
p_Ds_v12.interpolators.Q_vals_interpolator[seq(0.0, 5.0, 0.5)]
p_Ds_v12.interpolators.C_rates_interpolator[seq(0.0, 5.0, 0.5)]
p_Ds_v12.interpolators.mu_vals_interpolator[seq(0.0, 5.0, 0.5)]
p_Ds_v12.params.mu_vals_t
p_Ds_v7.params.mu_vals

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)
v12_ancstates_bots_v5 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchBot[R_order,])
v12_ancstates_tops_v5 = deepcopy(res.anc_estimates_at_each_nodeIndex_branchTop[R_order,])

v7_anc_branchBot = [0.9594100108999198 0.040589989100080166; 0.9594100108999198 0.040589989100080166; 0.006534383820866188 0.9934656161791339; 0.006534383820866188 0.9934656161791339; NaN NaN; 0.5498905772628557 0.45010942273714427; 0.4024797879588382 0.5975202120411617;]

v7_anc_branchTop = [1.0 0.0; 1.0 0.0; 0.0 1.0; 0.0 1.0; 0.47618518261084697 0.523814817389153; 0.9558446837332132 0.04415531626678687; 0.006058412033613736 0.9939415879663863;]

v12_ancstates_bots_v5a = transpose(reshape(flat2(v12_ancstates_bots_v5), (2,7)))
v12_ancstates_tops_v5a = transpose(reshape(flat2(v12_ancstates_tops_v5), (2,7)))

v7_anc_branchBot[isnan.(v7_anc_branchBot)] .= 0.0
v12_ancstates_bots_v5a[isnan.(v12_ancstates_bots_v5a)] .= 0.0

v7_anc_branchTop
v12_ancstates_tops_v5a

(v7_anc_branchBot .- v12_ancstates_bots_v5a) .< 1.0e-3
(v7_anc_branchTop .- v12_ancstates_tops_v5a) .< 1.0e-3

@test all((v7_anc_branchBot .- v12_ancstates_bots_v5a) .< 1.0e-3)
@test all((v7_anc_branchTop .- v12_ancstates_tops_v5a) .< 1.0e-3)


# end # END @testset "runtests_BiSSE_tree_n3" begin




