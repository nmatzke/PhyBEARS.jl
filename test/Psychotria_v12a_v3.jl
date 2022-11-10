
#######################################################
# Example inference on a simulated tree & geography dataset
# (from Wallis's ss8_sim_001)
# 2022-11-09
#
# We will compare inference under 
# * SSE likelihood calculator v7 (constant rates through time, very fast)
# * SSE likelihood calculator v12 (changing rates through time, very fast)
#   - Here, to start, we will only change the distances through time
#######################################################

using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using PhyBEARS
using DataFrames
using CSV

# Change the working directory as needed
wd = "/GitHub/PhyBEARS.jl/data/"
cd(wd)

# This simulation has 50 living species
trfn = "Psychotria_tree.newick"
tr = readTopology(trfn)
trdf = prt(tr)
oldest_possible_age = 100.0

lgdata_fn = "Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = true
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
bmo.type[bmo.rownames .== "j"] .= "free";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 0.034;
bmo.est[bmo.rownames .== "e"] .= 0.028;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.est[bmo.rownames .== "x"] .= 0.0;
bmo.est .= bmo_updater_v1(bmo);

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;
numstates = setup.numstates

p = p_Ds_v5
solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-12;
solver_options.reltol = 1e-12;

#######################################################
# Read in and parse distances and area-of-areas
#######################################################
#files.times_fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/v12a_times.txt"
#files.distances_fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/v12a_distances.txt"
#files.area_of_areas_fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/v12a_area_of_areas.txt"

interpolators = files_to_interpolators(files, setup.numareas, setup.states_list, setup.v_rows, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals; oldest_possible_age=100.0);


p_Es_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo, interpolators=interpolators);


p_Es_v10 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo, interpolators=interpolators);

p = p_Es_v10

# Add Q, C interpolators
temptimes = reduce(vcat, [interpolators.interpolator_times, sort(unique(trdf.node_age))])
times = sort(unique(temptimes))
p = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v10, times);

# Interpolators
p.interpolators.Q_vals_interpolator(0.0)[1:3]
p.interpolators.Q_vals_interpolator(1.0)[1:3]
p.interpolators.Q_vals_interpolator(2.0)[1:3]
p.interpolators.Q_vals_interpolator(3.0)[1:3]

p.interpolators.C_rates_interpolator(0.0)[1:3]
p.interpolators.C_rates_interpolator(1.0)[1:3]
p.interpolators.C_rates_interpolator(2.0)[1:3]
p.interpolators.C_rates_interpolator(3.0)[1:3]

p_Es_v12 = p
p_Es_v10 = p


# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")


prob_Es_v10 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v10_simd_sums, p_Es_v10.uE, Es_tspan, p_Es_v10);
# This solution is an interpolator
sol_Es_v10 = solve(prob_Es_v10, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v10;

p_Ds_v7 = (n=p_Es_v10.n, params=p_Es_v10.params, p_indices=p_Es_v10.p_indices, p_TFs=p_Es_v10.p_TFs, uE=p_Es_v10.uE, terms=p_Es_v10.terms, setup=p_Es_v10.setup,states_as_areas_lists=p_Es_v10.states_as_areas_lists, bmo=p_Es_v10.bmo, sol_Es_v5=sol_Es_v10);

# Check the interpolator
p_Ds_v7.sol_Es_v5(1.0)
Es_interpolator(1.0)


prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
# This solution is an interpolator
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v12;



p_Ds_v10 = (n=p_Es_v10.n, params=p_Es_v10.params, p_indices=p_Es_v10.p_indices, p_TFs=p_Es_v10.p_TFs, uE=p_Es_v10.uE, terms=p_Es_v10.terms, setup=p_Es_v10.setup, states_as_areas_lists=p_Es_v10.states_as_areas_lists, use_distances=p_Es_v10.use_distances, bmo=p_Es_v10.bmo, interpolators=p_Es_v10.interpolators, sol_Es_v10=sol_Es_v10);

p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);


# Calculate the Ds & total lnL via downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)




(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)



(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)



p_Ds_v12.interpolators.Q_vals_interpolator(0.0)[1:3]
p_Ds_v12.interpolators.Q_vals_interpolator(1.0)[1:3]
p_Ds_v12.interpolators.Q_vals_interpolator(2.0)[1:3]
p_Ds_v12.interpolators.Q_vals_interpolator(3.0)[1:3]

p_Ds_v12.interpolators.C_rates_interpolator(0.0)[1:3]
p_Ds_v12.interpolators.C_rates_interpolator(1.0)[1:3]
p_Ds_v12.interpolators.C_rates_interpolator(2.0)[1:3]
p_Ds_v12.interpolators.C_rates_interpolator(3.0)[1:3]

