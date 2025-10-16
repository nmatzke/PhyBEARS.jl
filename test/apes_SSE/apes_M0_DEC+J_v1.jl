
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
using DifferentialEquations
using Test						# for @test, @testset
using PhyloBits
using PhyloBits.TrUtils	# for vvdf
using PhyloBits.PNreadwrite # for readTopology
using PhyloBits.TreeTable # for ML_yule_birthRate
using PhyBEARS
using PhyBEARS.Uppass
using PhyBEARS.Parsers
using PhyBEARS.StateSpace		# for numstates_from_numareas
using PhyBEARS.Optimizers		# for bmo_updater_v1_SLOW
using DataFrames
using CSV


# Change the working directory as needed
wd = expanduser("~/GitHub/PhyBEARS.jl/test/apes_SSE/")
cd(wd)

"""

include(expanduser("~/GitHub/PhyBEARS.jl/test/apes_SSE/apes_M0_DEC+J_v1.jl"))
"""

#######################################################
# DEMONSTRATES MATCHING BETWEEN BIOGEOBEARS AND JULIA
# ON SIMPLE great ape phylogeny, 4-STATE DEC MODEL
#
# Run with:
# source(expanduser("~/GitHub/PhyBEARS.jl/Rsrc/compare_BGB_diversitree_DEC+J_v1.R"))
# Truth:
R_bgb_lnL = -1.170587






# BioGeoBEARS ancestral states under DEC+J
bgb_ancstates_AT_branchBots = [0, 0, 0, 0, NaN, 0, 0, 5.483716e-31, 1, 2.000002e-24, 4.500005e-24, NaN, 0.940298, 0.9850746, 1, 7.463047e-39, 1, 1, NaN, 4.486683e-07, 5.890815e-13, 1.492539e-14, 1.492539e-14, 1.791046e-13, 8.059706e-13, NaN, 0.05970153, 0.01492539];
bgb_ancstates_AT_nodes = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.4029848, 0.940298, 0.9850746, 1, 0, 1, 1, 0.4029853, 4.486684e-07, 6.338576e-13, 0, 0, 0, 0, 0.1940299, 0.05970153, 0.01492539];
bgb_ancstates_AT_branchBots_df = DataFrame(reshape(bgb_ancstates_AT_branchBots, (7, 4)), :auto)
bgb_ancstates_AT_nodes_df = DataFrame(reshape(bgb_ancstates_AT_nodes, (7, 4)), :auto)

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

bmo.est .= bmo_updater_v1_SLOW(bmo);

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

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

p_Ds_v5_updater_v1!(p_Ds_v7, inputs);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

# Save results
results_fn = expanduser("~/GitHub/PhyBEARS.jl/test/test_results.txt") # lnLs and times

txtvec = ["test: ", "apes_M0_DEC+J_v1.jl", "apes", "areas:2", "states:4", "DEC+J+Yule", "1 like", total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, R_bgb_lnL]

append_txtvec(results_fn, txtvec)



@testset "Apes DEC+J lnL vs. BioGeoBEARS lnL" begin
	@test abs(R_bgb_lnL - bgb_lnL) < 1e-4
end


# All ancestral states:
R_order = sort(trdf, :Rnodenums).nodeIndex
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

df1 = df1bot = bgb_ancstates_AT_branchBots_df
df2 = df2bot = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_order])
write_trdf_ancstates(results_fn, ["botstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")

df1 = df1top = bgb_ancstates_AT_nodes_df
df2 = df2top = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_order])
write_trdf_ancstates(results_fn, ["topstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")

compare_dfs(df1bot, df2bot; tol=1e-4)
get_max_df_diffs_byCol(df1bot, df2bot)

compare_dfs(df1top, df2top; tol=1e-4)
get_max_df_diffs_byCol(df1top, df2top)

@testset "Apes DEC+J ancstates vs Julia ancstates (with set parameters)" begin
	@test all(flat2(compare_dfs(df1bot, df2bot; tol=1e-4) .== 1.0))
	@test all(flat2(compare_dfs(df1top, df2top; tol=1e-4) .== 1.0))
end


#######################################################
# Maximum likelihood inference
#######################################################
inputs.bmo.type[inputs.bmo.rownames .== "birthRate"] .= "free"
inputs.bmo.type[inputs.bmo.rownames .== "deathRate"] .= "fixed"
inputs.bmo.est[inputs.bmo.rownames .== "deathRate"] .= 0.0
inputs.bmo[inputs.bmo.type .== "free",:]

pars = inputs.bmo.est[inputs.bmo.type .== "free"]
parnames = inputs.bmo.rownames[inputs.bmo.type .== "free"]
#func = x -> func_to_optimize_v7(x, parnames, inputs, p_Ds_v7; returnval="bgb_lnL", printlevel=1)
func = x -> func_to_optimize_v7(x, parnames, inputs, p_Ds_v7; returnval="bgb_lnL", printlevel=1)
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

# (enforce the birthRate)
#inputs.bmo.est[inputs.bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr)

bmo_updater_v2(inputs.bmo, inputs.setup.bmo_rows);
inputs.bmo.est[:] = bmo_updater_v2(inputs.bmo, inputs.setup.bmo_rows);

p_Ds_v5_updater_v1!(p_Ds_v7, inputs);

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

# Save results
txtvec = ["test: ", "apes_M0_DEC+J_v1.jl", "apes", "areas:2", "states:4", "DEC+J+Yule", "ML", total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, R_bgb_lnL]

append_txtvec(results_fn, txtvec)


# Root ancestral states:
Rnames(res)
round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]; digits=3)

# All ancestral states:
# Show ancestral state probability estimates
R_order = sort(trdf, :Rnodenums).nodeIndex
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

df1 = df1bot = bgb_ancstates_AT_branchBots_df
df2 = df2bot = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_order])
write_trdf_ancstates(results_fn, ["botstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")

df1 = df1top = bgb_ancstates_AT_nodes_df
df2 = df2top = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_order])
write_trdf_ancstates(results_fn, ["topstates"], trdf[R_order,:], df1, df2; mode="a", delim="\t")

compare_dfs(df1bot, df2bot; tol=1e-4)
get_max_df_diffs_byCol(df1bot, df2bot)

compare_dfs(df1top, df2top; tol=1e-4)
get_max_df_diffs_byCol(df1top, df2top)

@testset "Apes DEC+J ancstates vs. Julia ancstates (with ML parameters)" begin
	@test all(flat2(compare_dfs(df1bot, df2bot; tol=1e-4) .== 1.0))
	@test all(flat2(compare_dfs(df1top, df2top; tol=1e-4) .== 1.0))
end

#######################################################
# NOTE: 
# 
# ML_yule_birthRate(tr), following APE's Yule, gets a birthRate of 0.222222
#
# However, ML optimization here, optimizing on bgb_lnl, gives a birthRate of 0.333333
#
# This is a matter of counting, or not, an extra speciation rate (the root), i.e.
# conditioning on survival.
#
#######################################################






