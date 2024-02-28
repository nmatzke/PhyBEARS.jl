
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
wd = "/GitHub/PhyBEARS.jl/test/apes_SSE/"
cd(wd)

"""
include("/GitHub/PhyBEARS.jl/test/apes_SSE/fossils_apes_M0_DEC_v1.jl")
"""

#######################################################
# DEMONSTRATES MATCHING BETWEEN BIOGEOBEARS AND JULIA
# ON SIMPLE great ape phylogeny, 4-STATE DEC MODEL
#
# Run with:
# source("/GitHub/PhyBEARS.jl/ex/compare_BGB_diversitree_DEC+J_v1.R")
# Truth:
R_bgb_lnL = -4.481012

# BioGeoBEARS ancestral states under DEC+J
bgb_ancstates_AT_branchBots = [0, 0, 0, 0, NaN, 0, 0, 9.55885872371469e-14, 0.999999999997088, 1.02736516865682e-13, 2.3942137600093e-13, NaN, 0.0212357703981079, 0.0324086154040224, 0.999999999998852, 1.85939277741373e-12, 0.999999999999754, 0.999999999999244, NaN, 0.757828224601766, 0.630413600097194, 1.05227375864171e-12, 1.05227375864171e-12, 1.43663791560109e-13, 5.17042461743951e-13, NaN, 0.220936005000126, 0.337177784498784];
bgb_ancstates_AT_nodes = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 6.42693382782259e-14, 2.27415872607374e-14, 9.55885855108166e-14, 1, 0, 1, 1, 0.757828232249181, 0.630413602214152, 1.85939274383416e-12, 0, 0, 0, 0, 0.242171767750754, 0.369586397785825, 0.999999999998045];
bgb_ancstates_AT_branchBots_df = DataFrame(reshape(bgb_ancstates_AT_branchBots, (7, 4)), :auto)
bgb_ancstates_AT_nodes_df = DataFrame(reshape(bgb_ancstates_AT_nodes, (7, 4)), :auto)

# Ape tree, with or without fossils
trfn = "apes_tree.newick"
#trstr = "(((chimp:1.0,(human:0.5):0.5):1.0,gorilla:2.0):1.0,orang:3.0);"
#tr = readTopology(trstr)
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
bmo.type[bmo.rownames .== "j"] .= "fixed";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 0.1010557;
bmo.est[bmo.rownames .== "e"] .= 1e-12;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.est[bmo.rownames .== "x"] .= 0.0;

bmo.est[:] .= bmo_updater_v1_SLOW(bmo);

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

prtCp(p_Ds_v7)
p_Ds_v5_updater_v1!(p_Ds_v7, inputs);
prtCp(p_Ds_v7)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


@testset "Apes DEC lnL, regular tree" begin
	@test abs(R_bgb_lnL - bgb_lnL) < 1e-5
end

resNF = deepcopy(res);
Julia_sum_lqNF = deepcopy(Julia_sum_lq)
rootstates_lnL_NF = deepcopy(rootstates_lnL)
Julia_total_lnLs1_NF = deepcopy(Julia_total_lnLs1)
bgb_lnL_NF = deepcopy(bgb_lnL)
trdfNF = deepcopy(trdf)
p_Ds_v7_NF = deepcopy(p_Ds_v7);

prtCp(p_Ds_v7_NF)



# Ape tree, with direct ancestor
#trfn = "apes_tree.newick"
trstr = "(((chimp:1.0,(human:0.5):0.5):1.0,gorilla:2.0):1.0,orang:3.0);"
tr = readTopology(trstr)
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
bmo.type[bmo.rownames .== "j"] .= "fixed";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 0.1010557;
bmo.est[bmo.rownames .== "e"] .= 1e-12;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.est[bmo.rownames .== "x"] .= 0.0;

bmo.est[:] .= bmo_updater_v1_SLOW(bmo);

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

prtCp(p_Ds_v7)
p_Ds_v5_updater_v1!(p_Ds_v7, inputs);
prtCp(p_Ds_v7)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

# Save results
results_fn = "/GitHub/PhyBEARS.jl/test/test_results.txt" # lnLs and times

txtvec = ["test: ", "fossils_apes_M0_DEC_v1.jl", "apes_wFossilAncNode", "areas:2", "states:4", "DEC+Yule", "1 like", total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, R_bgb_lnL]

append_txtvec(results_fn, txtvec)



prtQp(p_Ds_v7)
prtQp(p_Ds_v7_NF)
prtCp(p_Ds_v7)
prtCp(p_Ds_v7_NF)

@testset "Apes DEC lnL, tree with a direct ancestor node" begin
	@test abs(R_bgb_lnL - bgb_lnL) < 1e-5

	@test abs(Julia_sum_lq - Julia_sum_lqNF) < 1e-5
	@test abs(Julia_total_lnLs1 - Julia_total_lnLs1_NF) < 1e-5
	@test abs(rootstates_lnL - rootstates_lnL_NF) < 1e-5
	@test abs(bgb_lnL - bgb_lnL_NF) < 1e-5

end



Julia_sum_lq - Julia_sum_lqNF
Julia_total_lnLs1 - Julia_total_lnLs1_NF
rootstates_lnL - rootstates_lnL_NF
bgb_lnL - bgb_lnL_NF

vfft(res.normlikes_at_each_nodeIndex_branchTop)
vfft(resNF.normlikes_at_each_nodeIndex_branchTop)

vfft(res.likes_at_each_nodeIndex_branchTop)
vfft(resNF.likes_at_each_nodeIndex_branchTop)

ind = [1,2,4,5,6,7,8]
vfft(res.likes_at_each_nodeIndex_branchBot[ind,:]) ./ vfft(resNF.likes_at_each_nodeIndex_branchBot)


# This is the difference?
resNF.likes_at_each_nodeIndex_branchTop[3][2]
log(resNF.likes_at_each_nodeIndex_branchTop[3][2])

trdf
trdfNF


# All ancestral states:
R_order = sort(trdf, :Rnodenums).nodeIndex
index_branchBots = [1,8,3,4,5,6,7]
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
trdf2 = trdf[index_branchBots,:]
sort!(trdf2, :Rnodenums)
trdf2

df1 = df1bot = bgb_ancstates_AT_branchBots_df
df2 = df2bot = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_order][index_branchBots])
write_trdf_ancstates(results_fn, ["botstates"], trdf2, df1, df2; mode="a", delim="\t")

df1 = df1top = bgb_ancstates_AT_nodes_df
df2 = df2top = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_order][1:7])
write_trdf_ancstates(results_fn, ["topstates"], trdf2, df1, df2; mode="a", delim="\t")

compare_dfs(df1bot, df2bot; tol=1e-4)
get_max_df_diffs_byCol(df1bot, df2bot)

compare_dfs(df1top, df2top; tol=1e-4)
get_max_df_diffs_byCol(df1top, df2top)

@testset "Apes DEC ancstates, with a direct ancestor" begin
	@test all(flat2(compare_dfs(df1bot, df2bot; tol=1e-4) .== 1.0))
	@test all(flat2(compare_dfs(df1top, df2top; tol=1e-4) .== 1.0))
end


uppass_ancstates_v7!(resNF, trdfNF, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
R_order = sort(trdfNF, :Rnodenums).nodeIndex
df1 = df1bot = bgb_ancstates_AT_branchBots_df
df2 = df2bot = vfft(resNF.anc_estimates_at_each_nodeIndex_branchBot[R_order])
df1 = df1top = bgb_ancstates_AT_nodes_df
df2 = df2top = vfft(resNF.anc_estimates_at_each_nodeIndex_branchTop[R_order])

@testset "Apes DEC ancstates, with a direct ancestor" begin
	@test all(flat2(compare_dfs(df1bot, df2bot; tol=1e-4) .== 1.0))
	@test all(flat2(compare_dfs(df1top, df2top; tol=1e-4) .== 1.0))
end








#######################################################
# Ape tree, with hooknode
#######################################################

#trfn = "apes_tree.newick"
trstr = "(((chimp:1.0,(human:0.5,fossil:1.0e-6):0.5):1.0,gorilla:2.0):1.0,orang:3.0);"
tr = readTopology(trstr)
prt(tr)

oldest_possible_age = 100.0

lgdata_fn = "geog_wFossil.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = true
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0000001, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
bmo.type[bmo.rownames .== "j"] .= "fixed";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 0.1010557;
bmo.est[bmo.rownames .== "e"] .= 1e-12;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.est[bmo.rownames .== "x"] .= 0.0;

bmo.est[:] .= bmo_updater_v1_SLOW(bmo);

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

# the Fossil hooknode is node #3;
# replacing with all 1s, as this is about comparison to previous trees
hooknode_num = 3
res.likes_at_each_nodeIndex_branchTop[hooknode_num] .= 1.0
res.normlikes_at_each_nodeIndex_branchTop[hooknode_num] .= 1.0


p_Es_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo);

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

prtCp(p_Ds_v7)
p_Ds_v5_updater_v1!(p_Ds_v7, inputs);
prtCp(p_Ds_v7)

# 2023-02-27
#p_Ds_v7.params.psi_vals .= 0.00000001
#modify_tiplikes_sampling_fossils_v7!(inputs, p_Ds_v7, geog_df)

solver_options=inputs.solver_options; max_iterations=10^5; return_lnLs=true
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

prtQp(p_Ds_v7)
prtQp(p_Ds_v7_NF)
prtCp(p_Ds_v7)
prtCp(p_Ds_v7_NF)

# You added a hook tip with 1s at any state state, this multiplies
# the total marginal likelihood by 4
Julia_sum_lq
# -7.309498852764715

Julia_sum_lqNF
# -8.695784978077587

# Multiply the raw likelihood by 4
log(exp(Julia_sum_lqNF) * 4)
# -7.309490616957696

(Julia_sum_lq+log(1/4)) - Julia_sum_lqNF
(Julia_total_lnLs1+log(1/4)) - Julia_total_lnLs1_NF
rootstates_lnL - rootstates_lnL_NF
(bgb_lnL+log(1/4)) - bgb_lnL_NF


txtvec = ["test: ", "fossils_apes_M0_DEC_v1.jl", "apes_wFossilHookNode", "areas:2", "states:4", "DEC+Yule", "1 like", total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, bgb_lnL_NF-log(1/4)]

append_txtvec(results_fn, txtvec)





ind = [1,2,5,6,7,8,9]
vfft(res.likes_at_each_nodeIndex_branchBot[ind])
vfft(resNF.likes_at_each_nodeIndex_branchBot)

vfft(res.normlikes_at_each_nodeIndex_branchTop[ind])
vfft(resNF.normlikes_at_each_nodeIndex_branchTop)



vfft(res.likes_at_each_nodeIndex_branchTop[ind])
vfft(resNF.likes_at_each_nodeIndex_branchTop)

vfft(res.normlikes_at_each_nodeIndex_branchTop[ind])
vfft(resNF.normlikes_at_each_nodeIndex_branchTop)

res_wHookTip = deepcopy(res);

p_Ds_v7.params.psi_vals # All psi values are 0.0

@testset "Apes DEC lnL, after adding a fossil hooktip with all 1s, and adding a log(1/4) correction to the lnL (with psi=0.0; modify_tiplikes_sampling_fossils_v7!() is NOT RUN.)" begin
	@test abs(R_bgb_lnL - (bgb_lnL+log(1/4))) < 1e-5

	@test abs((Julia_sum_lq+log(1/4)) - Julia_sum_lqNF) < 1e-5
	@test abs((Julia_total_lnLs1+log(1/4)) - Julia_total_lnLs1_NF) < 1e-5
	@test abs(rootstates_lnL - rootstates_lnL_NF) < 1e-5
	@test abs((bgb_lnL+log(1/4)) - bgb_lnL_NF) < 1e-5

end

# Ancestral states
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
ind = [1,2,5,6,7,8,9] # cut the hooknode/tip from the Julia-ordered table
R_order = sort(trdf, :Rnodenums).nodeIndex
Rind = [1,2,4,5,6,7,8] # cut the hooknode/tip from the R-ordered table

df1 = df1bot = bgb_ancstates_AT_branchBots_df
df2 = df2bot = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_order][Rind])
write_trdf_ancstates(results_fn, ["botstates"], trdf[R_order,:][Rind,:], df1, df2; mode="a", delim="\t")

df1 = df1top = bgb_ancstates_AT_nodes_df
df2 = df2top = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_order][Rind])
write_trdf_ancstates(results_fn, ["topstates"], trdf[R_order,:][Rind,:], df1, df2; mode="a", delim="\t")

@testset "Apes DEC ancstates, after adding a fossil hooktip with all 1s" begin
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








#######################################################
# By adding:
# a sampling event to the tip likelihood
# a psi > 0 to the rates, changing the likelihoods on the branches
#
# ...we should get a predictable change in the lnL
#######################################################

#trfn = "apes_tree.newick"
trstr = "(((chimp:1.0,(human:0.5,fossil:1.0e-6):0.5):1.0,gorilla:2.0):1.0,orang:3.0);"
tr = readTopology(trstr)
prt(tr)

oldest_possible_age = 100.0

lgdata_fn = "geog_wFossil.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = true
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0000001, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
bmo.type[bmo.rownames .== "j"] .= "fixed";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 0.1010557;
bmo.est[bmo.rownames .== "e"] .= 1e-12;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.est[bmo.rownames .== "x"] .= 0.0;

bmo.est[:] .= bmo_updater_v1_SLOW(bmo);

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

# the Fossil hooknode is node #3;
# replacing with all 1s, as this is about comparison to previous trees
hooknode_num = 3
res.likes_at_each_nodeIndex_branchTop[hooknode_num] .= 1.0
res.normlikes_at_each_nodeIndex_branchTop[hooknode_num] .= 1.0


# Modify the tip likelihoods
p_Ds_v5.params.psi_vals
inputs.bmo.est[inputs.setup.bmo_rows.psiRate] = 0.00000000001;
inputs.bmo.est .= bmo_updater_v1(inputs.bmo, inputs.setup.bmo_rows)
p_Ds_v5_updater_v1!(inputs.p_Ds_v5, inputs)
p_Ds_v5.params.psi_vals

orig_likes = deepcopy(inputs.res.likes_at_each_nodeIndex_branchTop)
modify_tiplikes_sampling_fossils_v7!(inputs, p_Ds_v5, geog_df)
inputs.res.likes_at_each_nodeIndex_branchTop
inputs.res.normlikes_at_each_nodeIndex_branchTop



# 2023-02-24_HERE



p_Es_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo);

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

prtCp(p_Ds_v7)
p_Ds_v5_updater_v1!(p_Ds_v7, inputs);
prtCp(p_Ds_v7)

solver_options=inputs.solver_options; max_iterations=10^5; return_lnLs=true
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


# Node 3 has 1.0 vs. 0.25 here:
vfft(res_wHookTip.normlikes_at_each_nodeIndex_branchTop)
vfft(res.normlikes_at_each_nodeIndex_branchTop)

# Node 3 has 1.0 vs. 4.0 here:
res_wHookTip.sumLikes_at_node_at_branchTop
res.sumLikes_at_node_at_branchTop


prtQp(p_Ds_v7)
prtQp(p_Ds_v7_NF)
prtCp(p_Ds_v7)
prtCp(p_Ds_v7_NF)

p_Ds_v7_NF.params.psi_vals
p_Ds_v7.params.psi_vals

nonroot_nodes = 1:nrow(trdf)
keepTF = nonroot_nodes .!= tr.root
nonroot_nodes = nonroot_nodes[keepTF]
sum_edge_lengths = sum(trdf.brlen[nonroot_nodes])

# You added a hook tip with 1s at any state state, this originally multiplied
# the total marginal likelihood by 4, but now with psi it's different
# (When psi = 1.0, we add log(1.0) + -psi*9.0 to the lnL)
Julia_sum_lq
# -17.695774724445098

Julia_sum_lqNF
# -8.695784978077587

Julia_sum_lqNF - 9
# -17.695784978077587

(Julia_sum_lq+log(1)) - Julia_sum_lqNF
(Julia_total_lnLs1+log(1)) - Julia_total_lnLs1_NF
rootstates_lnL - rootstates_lnL_NF
(bgb_lnL+log(1)) - bgb_lnL_NF
R_bgb_lnL - bgb_lnL


# Where did the x4 multiplier go?!?
res_wHookTip.lnL_at_node_at_branchTop
res.lnL_at_node_at_branchTop

res_wHookTip.like_at_branchBot
res.like_at_branchBot

res_wHookTip.likes_at_each_nodeIndex_branchTop
res.likes_at_each_nodeIndex_branchTop

# All branch bottoms different!
# Node 3 has 1.0 vs. 0.25 here
vfft(res_wHookTip.normlikes_at_each_nodeIndex_branchTop)
vfft(res.normlikes_at_each_nodeIndex_branchTop)

log.(vfft(res_wHookTip.likes_at_each_nodeIndex_branchBot)) .- log.(vfft(res.likes_at_each_nodeIndex_branchBot))


res_likes_at_each_nodeIndex_branchBot_mod = vfft(deepcopy(res.likes_at_each_nodeIndex_branchBot))
res_likes_at_each_nodeIndex_branchBot_mod[[1,5,7],:] = res_likes_at_each_nodeIndex_branchBot_mod[[1,5,7],:] .* exp(1.0);
res_likes_at_each_nodeIndex_branchBot_mod[[2,4],:] = res_likes_at_each_nodeIndex_branchBot_mod[[2,4],:] .* exp(1/2);
res_likes_at_each_nodeIndex_branchBot_mod[[6],:] = res_likes_at_each_nodeIndex_branchBot_mod[[6],:] .* exp(2.0);
res_likes_at_each_nodeIndex_branchBot_mod[[8],:] = res_likes_at_each_nodeIndex_branchBot_mod[[8],:] .* exp(3.0);
1.38629 

vfft(res_wHookTip.likes_at_each_nodeIndex_branchBot)
vfft(res.likes_at_each_nodeIndex_branchBot)
res_likes_at_each_nodeIndex_branchBot_mod

# Node 3 has 1.0 vs. 0.25 here:
vfft(res_wHookTip.normlikes_at_each_nodeIndex_branchTop)
vfft(res.normlikes_at_each_nodeIndex_branchTop)

# Node 3 has 1.0 vs. 4.0 here:
res_wHookTip.sumLikes_at_node_at_branchTop
res.sumLikes_at_node_at_branchTop

#######################################################
# The change in lnL is due to:
#
# 1. One fossil-sampling event, contributing log(psi) 
#    to the log-likelihood
# 
# 2. 9 branch-length units of non-sampling, contributing
#    -psi*9 to the lnL
# 
#######################################################

# psi=1.0
# difference: -3.202615406138527
psi = 1.0
log(psi) + -psi * sum_edge_lengths
# -9.000001000000001


# psi=0.1
# difference: -3.202615406138527
psi = 0.1
log(psi) + -psi * sum_edge_lengths
# -3.2025850929940455

# psi=0.2
# difference: -3.409511488102442
psi = 0.2
log(psi) + -psi * sum_edge_lengths


ind = [1,2,5,6,7,8,9]
vfft(res.likes_at_each_nodeIndex_branchBot[ind])
vfft(resNF.likes_at_each_nodeIndex_branchBot)

vfft(res.normlikes_at_each_nodeIndex_branchTop[ind])
vfft(resNF.normlikes_at_each_nodeIndex_branchTop)



vfft(res.likes_at_each_nodeIndex_branchTop[ind])
vfft(resNF.likes_at_each_nodeIndex_branchTop)

vfft(res.normlikes_at_each_nodeIndex_branchTop[ind])
vfft(resNF.normlikes_at_each_nodeIndex_branchTop)














psi = 1.0

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

# the Fossil hooknode is node #3;
# replacing with all 1s, as this is about comparison to previous trees
hooknode_num = 3
res.likes_at_each_nodeIndex_branchTop[hooknode_num] .= 1.0

p_Ds_v5.params.psi_vals
inputs.bmo.est[inputs.setup.bmo_rows.psiRate] = psi;
inputs.bmo.est .= bmo_updater_v1(inputs.bmo, inputs.setup.bmo_rows)
p_Ds_v5_updater_v1!(inputs.p_Ds_v5, inputs)
p_Ds_v5.params.psi_vals

orig_likes = deepcopy(inputs.res.likes_at_each_nodeIndex_branchTop)
modify_tiplikes_sampling_fossils_v7!(inputs, p_Ds_v5, geog_df)
inputs.res.likes_at_each_nodeIndex_branchTop

p_Es_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo);

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

psi_modifier_to_lnL = log(psi) + -(psi * sum_edge_lengths)



txtvec = ["test: ", "fossils_apes_M0_DEC_v1.jl", "apes_sampledFossilBranch", "areas:2", "states:4", "DEC+Yule", "1 like", total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, bgb_lnL_NF-log(1/4)+psi_modifier_to_lnL]

append_txtvec(results_fn, txtvec)




@testset "Apes DEC lnL, after adding a fossil hooktip with all 1s with a log(1/4) correction, & subtracting a psi*9 lnL" begin
	@test abs((bgb_lnL + log(1/4) - psi_modifier_to_lnL) - R_bgb_lnL) < 1e-4
	@test abs((Julia_sum_lq + log(1/4) - psi_modifier_to_lnL) - Julia_sum_lqNF) < 1e-4
	@test abs((Julia_total_lnLs1 + log(1/4) - psi_modifier_to_lnL) - Julia_total_lnLs1_NF) < 1e-4
	@test abs(rootstates_lnL - rootstates_lnL_NF) < 1e-4
	@test abs((bgb_lnL + log(1/4) - psi_modifier_to_lnL) - bgb_lnL_NF) < 1e-4
end



# Ancestral states
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
ind = [1,2,5,6,7,8,9] # cut the hooknode/tip from the Julia-ordered table
R_order = sort(trdf, :Rnodenums).nodeIndex
Rind = [1,2,4,5,6,7,8] # cut the hooknode/tip from the R-ordered table

df1 = df1bot = bgb_ancstates_AT_branchBots_df
df2 = df2bot = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_order][Rind])
write_trdf_ancstates(results_fn, ["botstates"], trdf[R_order,:][Rind,:], df1, df2; mode="a", delim="\t")

df1 = df1top = bgb_ancstates_AT_nodes_df
df2 = df2top = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_order][Rind])
write_trdf_ancstates(results_fn, ["topstates"], trdf[R_order,:][Rind,:], df1, df2; mode="a", delim="\t")

@testset "Apes DEC ancstates, after adding a fossil hooktip with all 1s" begin
	@test all(flat2(compare_dfs(df1bot, df2bot; tol=1e-4) .== 1.0))
	@test all(flat2(compare_dfs(df1top, df2top; tol=1e-4) .== 1.0))
end











#######################################################
# Different psi
#######################################################

psi = exp(1.0) # 2.718281828459045

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

# the Fossil hooknode is node #3;
# replacing with all 1s, as this is about comparison to previous trees
hooknode_num = 3
res.likes_at_each_nodeIndex_branchTop[hooknode_num] .= 1.0

p_Ds_v5.params.psi_vals
inputs.bmo.est[inputs.setup.bmo_rows.psiRate] = psi;
inputs.bmo.est .= bmo_updater_v1(inputs.bmo, inputs.setup.bmo_rows)
p_Ds_v5_updater_v1!(inputs.p_Ds_v5, inputs)
p_Ds_v5.params.psi_vals

orig_likes = deepcopy(inputs.res.likes_at_each_nodeIndex_branchTop)
modify_tiplikes_sampling_fossils_v7!(inputs, p_Ds_v5, geog_df)
inputs.res.likes_at_each_nodeIndex_branchTop

p_Es_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo);

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

psi_modifier_to_lnL = log(psi) + -(psi * sum_edge_lengths)

txtvec = ["test: ", "fossils_apes_M0_DEC_v1.jl", "apes_sampledFossilBranch_DiffPsi", "areas:2", "states:4", "DEC+Yule", "1 like", total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, bgb_lnL_NF-log(1/4)+psi_modifier_to_lnL]

append_txtvec(results_fn, txtvec)




@testset "Apes DEC lnL, after adding a fossil hooktip with all 1s with a log(1/4) correction, & subtracting a psi*9 lnL" begin
	@test abs((bgb_lnL + log(1/4) - psi_modifier_to_lnL) - R_bgb_lnL) < 1e-2
	@test abs((Julia_sum_lq + log(1/4) - psi_modifier_to_lnL) - Julia_sum_lqNF) < 1e-2
	@test abs((Julia_total_lnLs1 + log(1/4) - psi_modifier_to_lnL) - Julia_total_lnLs1_NF) < 1e-2
	@test abs(rootstates_lnL - rootstates_lnL_NF) < 1e-4
	@test abs((bgb_lnL + log(1/4) - psi_modifier_to_lnL) - bgb_lnL_NF) < 1e-2
end


# Ancestral states
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
ind = [1,2,5,6,7,8,9] # cut the hooknode/tip from the Julia-ordered table
R_order = sort(trdf, :Rnodenums).nodeIndex
Rind = [1,2,4,5,6,7,8] # cut the hooknode/tip from the R-ordered table

df1 = df1bot = bgb_ancstates_AT_branchBots_df
df2 = df2bot = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_order][Rind])
write_trdf_ancstates(results_fn, ["botstates"], trdf[R_order,:][Rind,:], df1, df2; mode="a", delim="\t")

df1 = df1top = bgb_ancstates_AT_nodes_df
df2 = df2top = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_order][Rind])
write_trdf_ancstates(results_fn, ["topstates"], trdf[R_order,:][Rind,:], df1, df2; mode="a", delim="\t")

@testset "Apes DEC ancstates, after adding a fossil hooktip with all 1s" begin
	@test all(flat2(compare_dfs(df1bot, df2bot; tol=1e-4) .== 1.0))
	@test all(flat2(compare_dfs(df1top, df2top; tol=1e-4) .== 1.0))
end







