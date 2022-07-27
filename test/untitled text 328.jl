#######################################################
# Goal: a script that runs a bunch of variant likelihoods
# at different error tolerances
#######################################################

using PhyloNetworks
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using LSODA						# for 
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyBEARS.Parsers
using PhyBEARS.Gmaps


"""
# Run with:
cd /GitHub/PhyBEARS.jl
JULIA_NUM_THREADS=23 julia
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/speedtests_PodoArau_wExtinction+J.jl")
"""




trfn = "/GitHub/PhyBEARS.jl/data/Psychotria_tree.newick"
tr = readTopology(trfn)

geogfn = "/GitHub/PhyBEARS.jl/data/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)
timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-3)
timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-6)
timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-9)
timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-12)
timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-15)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-3)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-6)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-9)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-12)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-15)


trfn = "/GitHub/PhyBEARS.jl/data/Cyrtandra.newick"
tr = readTopology(trfn)

geogfn = "/GitHub/PhyBEARS.jl/data/Cyrtandra_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(geogfn)

include_null_range=false
tol=1e-6

timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-3)
timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-6)
timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-9)
timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-12)
timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-15)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-3)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-6)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-9)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-12)
timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-15)


#######################################################
# Gold standard run: traditional SSE
# CVODE_BDF(linear_solver=:GMRES);
# User sets the error tolerance
#######################################################
function timed_run_Gflow(trfn, geogfn; include_null_range=false, tol=1e-6)
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# DEC model on Hawaiian Psychotria
	bmo = construct_BioGeoBEARS_model_object()
	bmo.type[bmo.rownames .== "j"] .= "free"
	bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
	bmo.est[bmo.rownames .== "deathRate"] .= 0.02
	bmo.est[bmo.rownames .== "d"] .= 0.01
	bmo.est[bmo.rownames .== "e"] .= 0.01
	bmo.est[bmo.rownames .== "a"] .= 0.0
	bmo.est[bmo.rownames .== "j"] .= 0.05
	bmo.est[:] = bmo_updater_v1(bmo);

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
	solver_options.dense = false;	

	# Solve the Es
	tmpt = Dates.now()
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_at=solver_options.saveat, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=solver_options.dense);
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
	Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, save_at=solver_options.saveat, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=solver_options.dense);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)

	res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_Gflow_v6 = vcat(flat2(res_Gflow_v6), calctimes)
	
	return(res_Gflow_v6)
end

function timed_run_goldstandard(trfn, geogfn; include_null_range=false, tol=1e-6)
	tr = readTopology(trfn)
	geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
	
	numareas = Rncol(geog_df)-1
	max_range_size = numareas
	n = numstates_from_numareas(numareas, max_range_size, include_null_range)
	
	# Just push a series of times onto calctime
	calctimes = []
	
	# DEC model on Hawaiian Psychotria
	bmo = construct_BioGeoBEARS_model_object()
	bmo.type[bmo.rownames .== "j"] .= "free"
	bmo.est[bmo.rownames .== "birthRate"] .= 0.032881638319078066
	bmo.est[bmo.rownames .== "deathRate"] .= 0.02
	bmo.est[bmo.rownames .== "d"] .= 0.01
	bmo.est[bmo.rownames .== "e"] .= 0.01
	bmo.est[bmo.rownames .== "a"] .= 0.0
	bmo.est[bmo.rownames .== "j"] .= 0.05
	bmo.est[:] = bmo_updater_v1(bmo);

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
	solver_options.dense = false;	

	# Solve the Es
	tmpt = Dates.now()
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_at=solver_options.saveat, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=solver_options.dense);
	push!(calctimes, (Dates.now() - tmpt).value/1000.0);
	
	# Add the Es interpolator to p_Ds_v5
	p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	tmpt = Dates.now()
	res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6, iteration_number_nFv6, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	push!(calctimes, (Dates.now() - tmpt).value/1000.0)
	
	res_nonFlow_v6 = vcat(flat2(res_nonFlow_v6), calctimes)
	
	return(res_nonFlow_v6)
end




