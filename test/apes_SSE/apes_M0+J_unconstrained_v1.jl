
#######################################################
# Example inference on a example tree & geography dataset
#
# We will compare inference under 
# * SSE likelihood calculator v7 (constant rates through time, very fast)
# * SSE likelihood calculator v12 (changing rates through time, fast)
#   - Here, to start, we will only change the distances through time
#######################################################

using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using PhyloBits.TrUtils	# for vvdf
using PhyBEARS
using PhyBEARS.Uppass
using DataFrames
using CSV

# Change the working directory as needed
wd = expanduser("~/GitHub/PhyBEARS.jl/test/apes_SSE/")
cd(wd)


#######################################################
# Answers from BioGeoBEARS
#######################################################
R_bgb_lnL = -1.170587


# BioGeoBEARS ancestral states under DEC+J
tmp_bgb_ancstates = [];

bgb_ancstates_df = DataFrame(reshape(tmp_bgb_ancstates, (37, 16)), :auto)

# Psychotria tree from Ree & Smith 2008
trfn = "apes_tree.newick"
tr = readTopology(trfn)
trdf = prt(tr)
oldest_possible_age = 100.0

lgdata_fn = "geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = true
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0000001, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
bmo.type[bmo.rownames .== "j"] .= "free";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 1e-12;
bmo.est[bmo.rownames .== "e"] .= 1e-12;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 2.99999;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.est[bmo.rownames .== "x"] .= 0.0;
#bmo.est[bmo.rownames .== "xv"] .= 0.1;
bmo.max[bmo.rownames .== "xv"] .= 10.0;

bmo.est[:] .= bmo_updater_v1_OLD(bmo);

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

p_Es_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo);

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

include(expanduser("~/GitHub/PhyBEARS.jl/test/apes_SSE/uppass_clado_v7A.jl"))

# Truth
uppass_from_root_to_node6_branchBot = [0.0, 0.875, 0.0, 0.12500000]

# Reproduce
uppass_probs_just_below_node = [0.0, 0.3333333333, 0.3333333333, 0.3333333333]
uppass_probs_just_below_node = [0.25, 0.25, 0.25, 0.25]
Ldownpass_likes = res.normlikes_at_each_nodeIndex_branchBot[3]
Rdownpass_likes = res.normlikes_at_each_nodeIndex_branchBot[4]
Ldownpass_likes_ones = [1.0, 1.0, 1.0, 1.0]
Rdownpass_likes_ones = [1.0, 1.0, 1.0, 1.0]

relprob_each_split_scenario = nodeOp_Cmat_get_condprobs_v7A(uppass_probs_just_below_node, Ldownpass_likes_ones, Rdownpass_likes, p_Ds_v7; use_Cijk_rates_t=false)

ctable1 = prtCp(p)
ctable = make_ctable_single_events(ctable1)

uppass_lprobs = repeat([0.0], n)
uppass_rprobs = repeat([0.0], n)
for statei in 1:n
	uppass_lprobs[statei] = sum(relprob_each_split_scenario[ctable.j .== statei])
	#uppass_rprobs[statei] = sum(relprob_each_split_scenario[ctable.k .== statei]) # discarded, non-target corner
end
uppass_lprobs


relprob_each_split_scenario = nodeOp_Cmat_get_condprobs_v7A(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes_ones, p_Ds_v7; use_Cijk_rates_t=false)

ctable1 = prtCp(p)
ctable = make_ctable_single_events(ctable1)

uppass_lprobs = repeat([0.0], n)
uppass_rprobs = repeat([0.0], n)
for statei in 1:n
	#uppass_lprobs[statei] = sum(relprob_each_split_scenario[ctable.j .== statei])
	uppass_rprobs[statei] = sum(relprob_each_split_scenario[ctable.k .== statei]) # discarded, non-target corner
end
uppass_rprobs

uppass_from_root_to_node6_branchBot



# All ancestral states:
R_order = sort(trdf, :Rnodenums).nodeIndex
uppass_ancstates_v7A!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

# Branch bottoms ("corners")
round.(vvdf(res.uppass_probs_at_each_nodeIndex_branchBot[R_order]), digits=4)
round.(vvdf(res.uppass_probs_at_each_nodeIndex_branchTop[R_order]), digits=4)


round.(vvdf(res.anc_estimates_at_each_nodeIndex_branchBot[R_order]), digits=3)
# Branch tops ("corners")
round.(vvdf(res.anc_estimates_at_each_nodeIndex_branchTop[R_order]), digits=3)

R_bgb_lnL
bgb_lnL
@test ( (R_bgb_lnL - bgb_lnL) < 0.0001)





#######################################################
# Maximum likelihood inference
#######################################################
#bmo.type[bmo.rownames .== "xv"] .= "free"
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "fixed"
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
#bmo.type[bmo.rownames .== "deathRate"] .= "birthRate"
#bmo.type[bmo.rownames .== "x"] .= "free"
#bmo.type[bmo.rownames .== "x"] .= "free"

bmo[bmo.type .== "free",:]

pars = bmo.est[bmo.type .== "free"]
parnames = bmo.rownames[bmo.type .== "free"]
#func = x -> func_to_optimize_v7(x, parnames, inputs, p_Ds_v7; returnval="bgb_lnL", printlevel=1)
func = x -> func_to_optimize_v7(x, parnames, inputs, p_Ds_v7; returnval="lnL", printlevel=1)
#pars = [0.04, 0.01, 0.001, 0.34, 0.0]
#pars = [0.04, 0.01, 0.001, 0.34]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


using NLopt
opt = NLopt.Opt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]
opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
#opt.ftol_abs = 0.0001 # tolerance on log-likelihood
#opt.ftol_rel = 0.01 # tolerance on log-likelihood
#opt.xtol_abs = 0.00001 # tolerance on parameters
#opt.xtol_rel = 0.001 # tolerance on parameters
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################



# Get the inputs & res:
pars = optx
pars
#pars = [0.9747407112459348, 0.8, 0.11]
#pars = [100.0, 1.8, 0.11]
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
bmo_updater_v2(inputs.bmo, inputs.setup.bmo_rows);
inputs.bmo.est[:] = bmo_updater_v2(inputs.bmo, inputs.setup.bmo_rows);

p_Ds_v5_updater_v1!(p_Ds_v7, inputs);

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

# Root ancestral states:
Rnames(res)
round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]; digits=3)
#  0.0  0.0  0.0  0.0  0.0  0.002  0.003  0.947  0.0  0.003  0.0  0.001  0.024  0.02  0.0  0.0

# All ancestral states:
R_order = sort(trdf, :Rnodenums).nodeIndex
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
rn(res)

# Show ancestral state probability estimates

# Branch bottoms ("corners")
round.(vvdf(res.anc_estimates_at_each_nodeIndex_branchBot[R_order]), digits=3)
# Branch tops ("corners")
round.(vvdf(res.anc_estimates_at_each_nodeIndex_branchTop[R_order]), digits=3)

bgb_ancstates_df

bgb_ancstates_df .- round.(vvdf(res.anc_estimates_at_each_nodeIndex_branchTop[R_order]), digits=4)


round.(vvdf(res.likes_at_each_nodeIndex_branchBot[R_order]), digits=4)