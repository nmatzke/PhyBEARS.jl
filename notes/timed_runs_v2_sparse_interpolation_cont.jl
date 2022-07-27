#######################################################
# Goal: a script that runs a bunch of variant likelihoods
# at different error tolerances
#######################################################

using PhyloNetworks
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DoubleFloats
using LSODA						# for 
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyBEARS.Parsers
using PhyBEARS.Gmaps
using Base.Threads				# for @spawn
using DataFrames
using CSV

"""
# Run with:
cd /GitHub/PhyBEARS.jl
JULIA_NUM_THREADS=23 julia
cd("/GitHub/PhyBEARS.jl/notes/")
include("/GitHub/PhyBEARS.jl/notes/timed_runs_v2_sparse_interpolation_cont.jl")
"""


"""
cd("/GitHub/PhyBEARS.jl/notes/")
include("/GitHub/PhyBEARS.jl/notes/timed_runs_v2_sparse_interpolation_cont.jl")

trfn = "/GitHub/PhyBEARS.jl/data/Klaus_Matzke_2020_PodoArau_197sp.newick"
tr = readTopology(trfn)

geogfn = "/GitHub/PhyBEARS.jl/data/Podocarpaceae_197_9areas_5Araucariaceae.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)

include_null_range=false
tol=1e-6
i=0

i = 0; tres = Dict()		# tres=time results, i=incrementer

i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
#i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
#i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
#i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15, Gseg_times=NaN, num_incs=100)
#i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-17, Gseg_times=NaN, num_incs=100)

i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9, Gseg_times=NaN, num_incs=100)
#i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15, Gseg_times=NaN, num_incs=100)
#i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-17, Gseg_times=NaN, num_incs=100)

tdf3 = indexed_Dict_to_DF(tres)
CSV.write("tdf3a_cont_sparse.txt", tdf3, delim="\t")

trfn = "/GitHub/PhyBEARS.jl/data/Klaus_Matzke_2020_PodoArau_197sp.newick"
tr = readTopology(trfn)

geogfn = "/GitHub/PhyBEARS.jl/data/Podocarpaceae_197_9areas_5Araucariaceae.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)

include_null_range=false
tol=1e-6
i=0

i = 0; tres = Dict()		# tres=time results, i=incrementer

i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
i+=1; tres[i]=timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
i+=1; tres[i]=timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
i+=1; tres[i]=timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
#i+=1; tres[i]=timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
#i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3)
i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9)
i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12)
#i+=1; tres[i]=timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15)

i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15, Gseg_times=NaN, num_incs=100)
#i+=1; tres[i]=timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-17, Gseg_times=NaN, num_incs=100)

i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-3, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-9, Gseg_times=NaN, num_incs=100)
#i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-12, Gseg_times=NaN, num_incs=100)
i+=1; tres[i]=tres[i-1]
i+=1; tres[i]=tres[i-1]
#i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-15, Gseg_times=NaN, num_incs=100)
#i+=1; tres[i]=timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-17, Gseg_times=NaN, num_incs=100)

tdf3 = indexed_Dict_to_DF(tres)
CSV.write("tdf3b_cont_sparse.txt", tdf3, delim="\t")




"""

#######################################################
# Gold standard run: Traditional SSE
# CVODE_BDF(linear_solver=:GMRES);
# User sets the error tolerance
#######################################################
function timed_run_goldstandard(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# Some reasonable defaults (includes death & jump processes)
	if isnan(bmo) == true
		bmo = construct_BioGeoBEARS_model_object()
		bmo.type[bmo.rownames .== "j"] .= "free"
		bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
		bmo.est[bmo.rownames .== "deathRate"] .= 0.02
		bmo.est[bmo.rownames .== "d"] .= 0.01
		bmo.est[bmo.rownames .== "e"] .= 0.01
		bmo.est[bmo.rownames .== "a"] .= 0.0
		bmo.est[bmo.rownames .== "j"] .= 0.05
		bmo.est[:] = bmo_updater_v1(bmo);
	end # END if isnan(bmo) == true
	

	# Set up the model
	tmpt = Dates.now()
	inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
	(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
	p_Ds_v5 = inputs.p_Ds_v5;
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

	push!(calctimes, (Dates.now() - tmpt).value/1000.0);

	saveats =  trdf.node_age[trdf.nodeType .!= "tip"];
	sort!(saveats);

	solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
	solver_options.saveat = saveats;
	solver_options.save_everystep = false;
	solver_options.abstol = tol;
	solver_options.reltol = tol;	
	solver_options.dense = true;	

	# Solve the Es
	# The Es solver must be continuous!
	tmpt = Dates.now()
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0);
	
	# Add the Es interpolator to p_Ds_v5
	p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	tmpt = Dates.now()
	res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6, iteration_number_nFv6, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6
	push!(calctimes, mean(res.calc_duration))
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_nonFlow_v6 = vcat(flat2(res_nonFlow_v6), calctimes)
	
	# Add the algorithm choice
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [string(sol_Es_v5.alg)])
	res_nonFlow_v6 = vcat(res_nonFlow_v6, ["CVODE_BDF"])
	ttl_t = res_nonFlow_v6[7] + res_nonFlow_v6[8] + res_nonFlow_v6[10]
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [ttl_t])
	
	# names of saved info
	col_names = ["downpass_t", "num_iters", "branch_lnL", "root_lnL", "ttl_LnL", "BGB_lnL", "setup_t", "Esolve_t", "Ds_1calc_t", "downpass_t2", "Es_alg", "Ds_alg", "ttl_t"]
	
	df = DataFrame(reshape(res_nonFlow_v6, 1, length(res_nonFlow_v6)), col_names)
	#rename_df(df, newnames)

	return(df)
end # END timed_run_goldstandard

#######################################################
# Gold standard run: Traditional SSE
# CVODE_BDF(linear_solver=:GMRES);
# User sets the error tolerance
#######################################################
function timed_run_trad_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# Some reasonable defaults (includes death & jump processes)
	if isnan(bmo) == true
		bmo = construct_BioGeoBEARS_model_object()
		bmo.type[bmo.rownames .== "j"] .= "free"
		bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
		bmo.est[bmo.rownames .== "deathRate"] .= 0.02
		bmo.est[bmo.rownames .== "d"] .= 0.01
		bmo.est[bmo.rownames .== "e"] .= 0.01
		bmo.est[bmo.rownames .== "a"] .= 0.0
		bmo.est[bmo.rownames .== "j"] .= 0.05
		bmo.est[:] = bmo_updater_v1(bmo);
	end # END if isnan(bmo) == true

	# Set up the model
	tmpt = Dates.now()
	inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
	(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
	p_Ds_v5 = inputs.p_Ds_v5;
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

	push!(calctimes, (Dates.now() - tmpt).value/1000.0);

	saveats =  trdf.node_age[trdf.nodeType .!= "tip"];
	sort!(saveats);

	solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
	solver_options.saveat = saveats;
	solver_options.save_everystep = true;
	solver_options.abstol = tol;
	solver_options.reltol = tol;	
	solver_options.dense = true;	

	# Solve the Es
	tmpt = Dates.now()
	# The Es solver must be continuous!
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, save_at=[], save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0);
	
	# Add the Es interpolator to p_Ds_v5
	p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	tmpt = Dates.now()
	res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6_solverFree!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6, iteration_number_nFv6, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6, most_common_solver) = res_nonFlow_v6
	push!(calctimes, mean(res.calc_duration))
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	#push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_nonFlow_v6 = vcat(flat2(res_nonFlow_v6[1:6]), calctimes)
	
	# Add the algorithm choice
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [string(sol_Es_v5.alg)])
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [string(most_common_solver)])
	ttl_t = res_nonFlow_v6[7] + res_nonFlow_v6[8] + res_nonFlow_v6[10]
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [ttl_t])

	# names of saved info
	col_names = ["downpass_t", "num_iters", "branch_lnL", "root_lnL", "ttl_LnL", "BGB_lnL", "setup_t", "Esolve_t", "Ds_1calc_t", "downpass_t2", "Es_alg", "Ds_alg", "ttl_t"]

	df = DataFrame(reshape(res_nonFlow_v6, 1, length(res_nonFlow_v6)), col_names)
	#rename_df(df, newnames)
	
	return(df)
end # END timed_run_trad_solverFree




#######################################################
# Gold standard run: Traditional SSE
# CVODE_BDF(linear_solver=:GMRES);
# Parallel version
# User sets the error tolerance
#######################################################
function timed_run_goldstandard_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
	# Error check
	if Base.Threads.nthreads() == 1
		txt = "STOP ERROR in timed_run_goldstandard_parallel(). Julia must be run with more than 1 thread for the parallelized function to work. Use e.g. 'JULIA_NUM_THREADS=23 julia' to start Julia, and Base.Threads.nthreads() to check that the number of threads is > 1."
		throw(txt)
	end
	
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# Some reasonable defaults (includes death & jump processes)
	if isnan(bmo) == true
		bmo = construct_BioGeoBEARS_model_object()
		bmo.type[bmo.rownames .== "j"] .= "free"
		bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
		bmo.est[bmo.rownames .== "deathRate"] .= 0.02
		bmo.est[bmo.rownames .== "d"] .= 0.01
		bmo.est[bmo.rownames .== "e"] .= 0.01
		bmo.est[bmo.rownames .== "a"] .= 0.0
		bmo.est[bmo.rownames .== "j"] .= 0.05
		bmo.est[:] = bmo_updater_v1(bmo);
	end # END if isnan(bmo) == true
	

	# Set up the model
	tmpt = Dates.now()
	inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
	(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
	p_Ds_v5 = inputs.p_Ds_v5;
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

	push!(calctimes, (Dates.now() - tmpt).value/1000.0);

	saveats =  trdf.node_age[trdf.nodeType .!= "tip"];
	sort!(saveats);

	solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
	solver_options.saveat = saveats;
	solver_options.save_everystep = false;
	solver_options.abstol = tol;
	solver_options.reltol = tol;	
	solver_options.dense = true;	

	# Solve the Es
	# The Es solver must be continuous!
	tmpt = Dates.now()
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0);
	
	# Add the Es interpolator to p_Ds_v5
	p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	tmpt = Dates.now()
	res_nonFlow_v6 = iterative_downpass_parallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6, iteration_number_nFv6, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6
	push!(calctimes, mean(res.calc_duration))
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_nonFlow_v6 = vcat(flat2(res_nonFlow_v6), calctimes)
	
	# Add the algorithm choice
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [string(sol_Es_v5.alg)])
	res_nonFlow_v6 = vcat(res_nonFlow_v6, ["CVODE_BDF"])
	ttl_t = res_nonFlow_v6[7] + res_nonFlow_v6[8] + res_nonFlow_v6[10]
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [ttl_t])
	
	# names of saved info
	col_names = ["downpass_t", "num_iters", "branch_lnL", "root_lnL", "ttl_LnL", "BGB_lnL", "setup_t", "Esolve_t", "Ds_1calc_t", "downpass_t2", "Es_alg", "Ds_alg", "ttl_t"]
	
	df = DataFrame(reshape(res_nonFlow_v6, 1, length(res_nonFlow_v6)), col_names)
	#rename_df(df, newnames)

	return(df)
end # END timed_run_goldstandard_parallel


#######################################################
# Traditional SSE, solver is decided by solve()
# Parallel version
# User sets the error tolerance
#######################################################
function timed_run_trad_solverFree_parallel(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
	# Error check
	if Base.Threads.nthreads() == 1
		txt = "STOP ERROR in timed_run_trad_solverFree_parallel(). Julia must be run with more than 1 thread for the parallelized function to work. Use e.g. 'JULIA_NUM_THREADS=23 julia' to start Julia, and Base.Threads.nthreads() to check that the number of threads is > 1."
		throw(txt)
	end
	
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# Some reasonable defaults (includes death & jump processes)
	if isnan(bmo) == true
		bmo = construct_BioGeoBEARS_model_object()
		bmo.type[bmo.rownames .== "j"] .= "free"
		bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
		bmo.est[bmo.rownames .== "deathRate"] .= 0.02
		bmo.est[bmo.rownames .== "d"] .= 0.01
		bmo.est[bmo.rownames .== "e"] .= 0.01
		bmo.est[bmo.rownames .== "a"] .= 0.0
		bmo.est[bmo.rownames .== "j"] .= 0.05
		bmo.est[:] = bmo_updater_v1(bmo);
	end # END if isnan(bmo) == true

	# Set up the model
	tmpt = Dates.now()
	inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
	(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
	p_Ds_v5 = inputs.p_Ds_v5;
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

	push!(calctimes, (Dates.now() - tmpt).value/1000.0);

	saveats =  trdf.node_age[trdf.nodeType .!= "tip"];
	sort!(saveats);

	solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
	solver_options.saveat = saveats;
	solver_options.save_everystep = true;
	solver_options.abstol = tol;
	solver_options.reltol = tol;	
	solver_options.dense = true;	

	# Solve the Es
	tmpt = Dates.now()
	# The Es solver must be continuous!
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0);
	
	# Add the Es interpolator to p_Ds_v5
	p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	tmpt = Dates.now()
	res_nonFlow_v6 = iterative_downpass_parallel_ClaSSE_v6_solverFree!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6, iteration_number_nFv6, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6, most_common_solver) = res_nonFlow_v6
	push!(calctimes, mean(res.calc_duration))
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	#push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_nonFlow_v6 = vcat(flat2(res_nonFlow_v6[1:6]), calctimes)
	
	# Add the algorithm choice
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [string(sol_Es_v5.alg)])
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [string(most_common_solver)])
	ttl_t = res_nonFlow_v6[7] + res_nonFlow_v6[8] + res_nonFlow_v6[10]
	res_nonFlow_v6 = vcat(res_nonFlow_v6, [ttl_t])

	# names of saved info
	col_names = ["downpass_t", "num_iters", "branch_lnL", "root_lnL", "ttl_LnL", "BGB_lnL", "setup_t", "Esolve_t", "Ds_1calc_t", "downpass_t2", "Es_alg", "Ds_alg", "ttl_t"]

	df = DataFrame(reshape(res_nonFlow_v6, 1, length(res_nonFlow_v6)), col_names)
	#rename_df(df, newnames)
	
	return(df)
end # END timed_run_trad_solverFree_parallel




#######################################################
# Gold standard run: Gflow SSE
# CVODE_BDF(linear_solver=:GMRES);
# User sets the error tolerance
#######################################################
function timed_run_Gflow_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# Some reasonable defaults (includes death & jump processes)
	if isnan(bmo) == true
		bmo = construct_BioGeoBEARS_model_object()
		bmo.type[bmo.rownames .== "j"] .= "free"
		bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
		bmo.est[bmo.rownames .== "deathRate"] .= 0.02
		bmo.est[bmo.rownames .== "d"] .= 0.01
		bmo.est[bmo.rownames .== "e"] .= 0.01
		bmo.est[bmo.rownames .== "a"] .= 0.0
		bmo.est[bmo.rownames .== "j"] .= 0.05
		bmo.est[:] = bmo_updater_v1(bmo);
	end # END if isnan(bmo) == true

	# Set up the model
	tmpt = Dates.now()
	inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
	(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
	p_Ds_v5 = inputs.p_Ds_v5;
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

	push!(calctimes, (Dates.now() - tmpt).value/1000.0);

	saveats =  trdf.node_age[trdf.nodeType .!= "tip"];
	sort!(saveats);

	solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
	solver_options.saveat = saveats;
	solver_options.save_everystep = false;
	solver_options.abstol = tol;
	solver_options.reltol = tol;	
	solver_options.dense = true;	

	# Solve the Es
	# The Es solver must be continuous!
	tmpt = Dates.now()
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0);
	
	# Add the Es interpolator to p_Ds_v5
	p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	tmpt = Dates.now()
	# Version 7/2 ClaSSE Gflow calculations
	G0 = Matrix{Float64}(I, n, n) ;
	# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
	tmpzero = repeat([0.0], n^2);
	A = reshape(tmpzero, (n,n));
	pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
	tspan = (0.0, 1.01 * maximum(trdf.node_age))
	#prob_Gs_v5_sub_i = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE_sub_i, G0, tspan, pG);
	prob_Gs_v5_sub_i = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE, G0, tspan, pG);
	#Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, solver_options.solver, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=true);
	Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=true);

	push!(calctimes, (Dates.now() - tmpt).value/1000.0)

	res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_Gflow_v6 = vcat(flat2(res_Gflow_v6), calctimes)
	# Add the algorithm choice
	res_Gflow_v6 = vcat(res_Gflow_v6, [string(sol_Es_v5.alg)])
	res_Gflow_v6 = vcat(res_Gflow_v6, [string(Gflow_to_01_GMRES.alg)])
	ttl_t = res_Gflow_v6[7] + res_Gflow_v6[8] + res_Gflow_v6[9] + res_Gflow_v6[10]
	res_Gflow_v6 = vcat(res_Gflow_v6, [ttl_t])
	
	col_names = ["downpass_t", "num_iters", "branch_lnL", "root_lnL", "ttl_LnL", "BGB_lnL", "setup_t", "Esolve_t", "Ds_1calc_t", "downpass_t2", "Es_alg", "Ds_alg", "ttl_t"]

	df = DataFrame(reshape(res_Gflow_v6, 1, length(res_Gflow_v6)), col_names)
	#rename_df(df, newnames)
	
	return(df)
end # END timed_run_Gflow_CVODE






#######################################################
# Auto-solver run: Gflow SSE
# CVODE_BDF(linear_solver=:GMRES);
# User sets the error tolerance
#######################################################
function timed_run_Gflow_solverFree(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# Some reasonable defaults (includes death & jump processes)
	if isnan(bmo) == true
		bmo = construct_BioGeoBEARS_model_object()
		bmo.type[bmo.rownames .== "j"] .= "free"
		bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
		bmo.est[bmo.rownames .== "deathRate"] .= 0.02
		bmo.est[bmo.rownames .== "d"] .= 0.01
		bmo.est[bmo.rownames .== "e"] .= 0.01
		bmo.est[bmo.rownames .== "a"] .= 0.0
		bmo.est[bmo.rownames .== "j"] .= 0.05
		bmo.est[:] = bmo_updater_v1(bmo);
	end # END if isnan(bmo) == true

	# Set up the model
	tmpt = Dates.now()
	inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
	(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
	p_Ds_v5 = inputs.p_Ds_v5;
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

	push!(calctimes, (Dates.now() - tmpt).value/1000.0);

	saveats =  trdf.node_age[trdf.nodeType .!= "tip"];
	sort!(saveats);

	solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
	solver_options.saveat = saveats;
	solver_options.save_everystep = false;
	solver_options.abstol = tol;
	solver_options.reltol = tol;	
	solver_options.dense = true;	

	# Solve the Es
	# The Es solver must be continuous!
	tmpt = Dates.now()
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0);
	
	# Add the Es interpolator to p_Ds_v5
	p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	tmpt = Dates.now()
	# Version 7/2 ClaSSE Gflow calculations
	G0 = Matrix{Float64}(I, n, n) ;
	# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
	tmpzero = repeat([0.0], n^2);
	A = reshape(tmpzero, (n,n));
	pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
	tspan = (0.0, 1.01 * maximum(trdf.node_age))
	prob_Gs_v5_sub_i = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE_sub_i, G0, tspan, pG);
	#Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, save_at=solver_options.saveat, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=solver_options.dense);
	Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=true);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)

	res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_Gflow_v6 = vcat(flat2(res_Gflow_v6), calctimes)
	# Add the algorithm choice
	res_Gflow_v6 = vcat(res_Gflow_v6, [string(sol_Es_v5.alg)])
	res_Gflow_v6 = vcat(res_Gflow_v6, [string(Gflow_to_01_GMRES.alg)])
	ttl_t = res_Gflow_v6[7] + res_Gflow_v6[8] + res_Gflow_v6[9] + res_Gflow_v6[10]
	res_Gflow_v6 = vcat(res_Gflow_v6, [ttl_t])
	
	col_names = ["downpass_t", "num_iters", "branch_lnL", "root_lnL", "ttl_LnL", "BGB_lnL", "setup_t", "Esolve_t", "Ds_1calc_t", "downpass_t2", "Es_alg", "Ds_alg", "ttl_t"]
	
	df = DataFrame(reshape(res_Gflow_v6, 1, length(res_Gflow_v6)), col_names)
	#rename_df(df, newnames)
	
	
	return(df)
end # END timed_run_Gflow_solverFree






#######################################################
# (Float64 version)
# Gold standard run: GflowArrays SSE
# CVODE_BDF(linear_solver=:GMRES);
# User sets the error tolerance
#######################################################
function timed_run_GflowArrays_CVODE(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6, Gseg_times=NaN, num_incs=NaN)
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# Some reasonable defaults (includes death & jump processes)
	if isnan(bmo) == true
		bmo = construct_BioGeoBEARS_model_object()
		bmo.type[bmo.rownames .== "j"] .= "free"
		bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
		bmo.est[bmo.rownames .== "deathRate"] .= 0.02
		bmo.est[bmo.rownames .== "d"] .= 0.01
		bmo.est[bmo.rownames .== "e"] .= 0.01
		bmo.est[bmo.rownames .== "a"] .= 0.0
		bmo.est[bmo.rownames .== "j"] .= 0.05
		bmo.est[:] = bmo_updater_v1(bmo);
	end # END if isnan(bmo) == true

	# Set up the model
	tmpt = Dates.now()
	inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
	(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
	p_Ds_v5 = inputs.p_Ds_v5;
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

	push!(calctimes, (Dates.now() - tmpt).value/1000.0);

	saveats =  trdf.node_age[trdf.nodeType .!= "tip"];
	sort!(saveats);

	solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
	solver_options.saveat = saveats;
	solver_options.save_everystep = false;
	solver_options.abstol = tol;
	solver_options.reltol = tol;	
	solver_options.dense = true;	

	# Solve the Es
	# The Es solver must be continuous!
	tmpt = Dates.now()
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0);
	
	# Add the Es interpolator to p_Ds_v5
	p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	tmpt = Dates.now()
	# Version 7/2 ClaSSE Gflow calculations
	G0 = Matrix{Float64}(I, n, n) ;
	# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
	tmpzero = repeat([0.0], n^2);
	A = reshape(tmpzero, (n,n));
	pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
	tspan = (0.0, 1.01 * maximum(trdf.node_age))
	prob_Gs_v5_sub_i = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE_sub_i, G0, tspan, pG);
	#Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, solver_options.solver, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=true);
	#Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=true);
	
	# Calculate the number of increments
	root_age = maximum(trdf.node_age)
	if isnan(Gseg_times) == true
		if isnan(num_incs) == true
			num_incs = 10.0
		end
		Gseg_times = seq(0.0, root_age*1.01, (root_age*1.01)/num_incs);
	end # END if isnan(Gseg_times) == true
	
	# Calculate array of Gflow matrices with float64 matrix multiplication
	(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap = Gmaps.construct_Gmap_interpolator(pG, Gseg_times; abstol=tol, reltol=tol, use_double=false);
	Gflow_via_Gmap = t -> Gmaps.interp_from_Gmap(t, Gmap)


	algorithms = collect(repeat([Any[]], length(Gseg_times)))
	for i in 1:length(Gseg_times)
		algorithms[i] = [string(Gflows_dict[i].alg)]
	end
	algorithms
	most_common_solver = get_a_most_common_value(algorithms)
	
	
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)

	res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_via_Gmap, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_Gflow_v6 = vcat(flat2(res_Gflow_v6), calctimes)
	# Add the algorithm choice
	res_Gflow_v6 = vcat(res_Gflow_v6, [string(sol_Es_v5.alg)])
	res_Gflow_v6 = vcat(res_Gflow_v6, [string(most_common_solver)])
	ttl_t = res_Gflow_v6[7] + res_Gflow_v6[8] + res_Gflow_v6[9] + res_Gflow_v6[10]
	res_Gflow_v6 = vcat(res_Gflow_v6, [ttl_t])
	
	col_names = ["downpass_t", "num_iters", "branch_lnL", "root_lnL", "ttl_LnL", "BGB_lnL", "setup_t", "Esolve_t", "Ds_1calc_t", "downpass_t2", "Es_alg", "Ds_alg", "ttl_t"]

	df = DataFrame(reshape(res_Gflow_v6, 1, length(res_Gflow_v6)), col_names)
	#rename_df(df, newnames)
	
	return(df)
end # END timed_run_Gflow_CVODE






#######################################################
# (Double64 version)
# Gold standard run: GflowArrays SSE
# CVODE_BDF(linear_solver=:GMRES);
# User sets the error tolerance
#######################################################
function timed_run_GflowArrays_CVODE_dbl(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6, Gseg_times=NaN, num_incs=NaN)
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# Some reasonable defaults (includes death & jump processes)
	if isnan(bmo) == true
		bmo = construct_BioGeoBEARS_model_object()
		bmo.type[bmo.rownames .== "j"] .= "free"
		bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
		bmo.est[bmo.rownames .== "deathRate"] .= 0.02
		bmo.est[bmo.rownames .== "d"] .= 0.01
		bmo.est[bmo.rownames .== "e"] .= 0.01
		bmo.est[bmo.rownames .== "a"] .= 0.0
		bmo.est[bmo.rownames .== "j"] .= 0.05
		bmo.est[:] = bmo_updater_v1(bmo);
	end # END if isnan(bmo) == true

	# Set up the model
	tmpt = Dates.now()
	inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
	(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
	p_Ds_v5 = inputs.p_Ds_v5;
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

	push!(calctimes, (Dates.now() - tmpt).value/1000.0);

	saveats =  trdf.node_age[trdf.nodeType .!= "tip"];
	sort!(saveats);

	solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
	solver_options.saveat = saveats;
	solver_options.save_everystep = false;
	solver_options.abstol = tol;
	solver_options.reltol = tol;	
	solver_options.dense = true;	

	# Solve the Es
	# The Es solver must be continuous!
	tmpt = Dates.now()
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0);
	
	# Add the Es interpolator to p_Ds_v5
	p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	tmpt = Dates.now()
	# Version 7/2 ClaSSE Gflow calculations
	G0 = Matrix{Float64}(I, n, n) ;
	# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
	tmpzero = repeat([0.0], n^2);
	A = reshape(tmpzero, (n,n));
	pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
	tspan = (0.0, 1.01 * maximum(trdf.node_age))
	prob_Gs_v5_sub_i = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE_sub_i, G0, tspan, pG);
	#Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, solver_options.solver, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=true);
	#Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, save_at=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=true);
	
	# Calculate the number of increments
	root_age = maximum(trdf.node_age)
	if isnan(Gseg_times) == true
		if isnan(num_incs) == true
			num_incs = 10.0
		end
		Gseg_times = seq(0.0, root_age*1.01, (root_age*1.01)/num_incs);
	end # END if isnan(Gseg_times) == true
	
	# Calculate array of Gflow matrices with float64 matrix multiplication
	(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap = Gmaps.construct_Gmap_interpolator(pG, Gseg_times; abstol=tol, reltol=tol, use_double=true);
	Gflow_via_Gmap = t -> Gmaps.interp_from_Gmap(t, Gmap)

	algorithms = collect(repeat([Any[]], length(Gseg_times)))
	for i in 1:length(Gseg_times)
		algorithms[i] = [string(Gflows_dict[i].alg)]
	end
	algorithms
	most_common_solver = get_a_most_common_value(algorithms)
	
	
	
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)

	res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_via_Gmap, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_Gflow_v6 = vcat(flat2(res_Gflow_v6), calctimes)
	# Add the algorithm choice
	res_Gflow_v6 = vcat(res_Gflow_v6, [string(sol_Es_v5.alg)])
	res_Gflow_v6 = vcat(res_Gflow_v6, [string(most_common_solver)])
	ttl_t = res_Gflow_v6[7] + res_Gflow_v6[8] + res_Gflow_v6[9] + res_Gflow_v6[10]
	res_Gflow_v6 = vcat(res_Gflow_v6, [ttl_t])
	
	col_names = ["downpass_t", "num_iters", "branch_lnL", "root_lnL", "ttl_LnL", "BGB_lnL", "setup_t", "Esolve_t", "Ds_1calc_t", "downpass_t2", "Es_alg", "Ds_alg", "ttl_t"]

	df = DataFrame(reshape(res_Gflow_v6, 1, length(res_Gflow_v6)), col_names)
	#rename_df(df, newnames)
	
	return(df)
end # END timed_run_Gflow_CVODE


