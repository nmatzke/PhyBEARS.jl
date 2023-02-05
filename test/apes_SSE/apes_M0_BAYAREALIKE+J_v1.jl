
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
wd = "/GitHub/PhyBEARS.jl/test/apes_SSE/"
cd(wd)

"""

include("/GitHub/PhyBEARS.jl/test/apes_SSE/apes_M0_BAYAREALIKE+J_v1.jl")
"""

#######################################################
# DEMONSTRATES MATCHING BETWEEN BIOGEOBEARS AND JULIA
# ON SIMPLE great ape phylogeny, 4-STATE BAYAREALIKE MODEL
#
# Run with:
# source("/GitHub/PhyBEARS.jl/Rsrc/compare_BGB_diversitree_BAYAREALIKE+J_v1.R")
# Truth:
R_bgb_lnL = -1.386306






# BioGeoBEARS ancestral states under BAYAREALIKE+J
bgb_ancstates_AT_branchBots = [0, 0, 0, 0, NaN, 0, 0, 5.00005e-20, 1, 2.00002e-14, 4.500045e-14, NaN, 0.999995, 1, 1, 2.750131e-29, 1, 1, NaN, 5.000075e-06, 5.002075e-11, 2.000013e-21, 2.000013e-21, 3.200016e-20, 6.000033e-14, NaN, 6.000033e-14, 3.200016e-20];
bgb_ancstates_AT_nodes = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.4999975, 0.999995, 1, 1, 0, 1, 1, 0.5000025, 5.000075e-06, 5.002575e-11, 0, 0, 0, 0, 6.000033e-14, 3.200017e-20, 2.000013e-21];
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

# BAYAREALIKE-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0000001, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
bmo.type[bmo.rownames .== "j"] .= "free";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 1.0e-07;
bmo.est[bmo.rownames .== "e"] .= 1.0e-07;
bmo.est[bmo.rownames .== "mx01y"] .= 0.9999;
bmo.est[bmo.rownames .== "mx01v"] .= 0.9999;
bmo.est[bmo.rownames .== "j"] .= 0.99999;
bmo.est[bmo.rownames .== "s"] .= 0.0;
bmo.est[bmo.rownames .== "v"] .= 0.0;
bmo.type[bmo.rownames .== "ysv"] .= "1-j";
bmo.type[bmo.rownames .== "ys"] .= "ysv*1";
bmo.type[bmo.rownames .== "y"] .= "ysv*1";
bmo.type[bmo.rownames .== "s"] .= "ysv*0";
bmo.type[bmo.rownames .== "v"] .= "ysv*0";
bmo.max[bmo.rownames .== "j"] .= 0.99999;
bmo.max[bmo.rownames .== "ysv"] .= 1.0;
bmo.max[bmo.rownames .== "ys"] .= 1.0;
bmo.max[bmo.rownames .== "y"] .= 1.0;
bmo.est[bmo.rownames .== "mx01v"] .= 0.0001;
bmo.est[bmo.rownames .== "mx01y"] .= 0.9999;
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


@testset "Apes BAYAREALIKE lnL" begin
	@test abs(R_bgb_lnL - bgb_lnL) < 1e-4
end


# All ancestral states:
R_order = sort(trdf, :Rnodenums).nodeIndex
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

df1 = df1bot = bgb_ancstates_AT_branchBots_df
df2 = df2bot = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_order])
df1 = df1top = bgb_ancstates_AT_nodes_df
df2 = df2top = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_order])

compare_dfs(df1bot, df2bot; tol=1e-4)
get_max_df_diffs_byCol(df1bot, df2bot)

compare_dfs(df1top, df2top; tol=1e-4)
get_max_df_diffs_byCol(df1top, df2top)

@testset "Apes BAYAREALIKE ancstates" begin
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

# Root ancestral states:
Rnames(res)
round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]; digits=3)

# All ancestral states:
# Show ancestral state probability estimates
R_order = sort(trdf, :Rnodenums).nodeIndex
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

df1 = df1bot = bgb_ancstates_AT_branchBots_df
df2 = df2bot = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_order])
df1 = df1top = bgb_ancstates_AT_nodes_df
df2 = df2top = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_order])

compare_dfs(df1bot, df2bot; tol=1e-4)
get_max_df_diffs_byCol(df1bot, df2bot)

compare_dfs(df1top, df2top; tol=1e-4)
get_max_df_diffs_byCol(df1top, df2top)

@testset "Apes BAYAREALIKE+J ancstates vs. Julia ML ancstates" begin
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






