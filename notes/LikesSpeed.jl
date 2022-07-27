module LikesSpeed
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/
using PhyBEARS.Parsers 		# for e.g. getranges_from_LagrangePHYLIP
using PhyBEARS.StateSpace 	# for e.g. numstates_from_numareas
using PhyBEARS.TrUtils 		# for e.g. flat2
using PhyBEARS.SSEs 
using PhyBEARS.ModelLikes
using PhyBEARS.TreePass
using PhyBEARS.Flow 
using DataFrames			# for e.g. DataFrame()
using PhyloNetworks		# for e.g. readTopology()
using Dates						# for e.g. DateTime, Dates.now()
using Distributed			# for e.g. @spawn
using Random					# for MersenneTwister()
using DifferentialEquations # for ODEProblem
using Sundials						# for CVODE_BDF()
using Dates						# for Dates.now()
using Statistics 			# for mean()
export timed_run_goldstandard




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
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, saveat=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
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



function timed_run_LapackBand(trfn, geogfn; bmo=NaN, include_null_range=false, tol=1e-6)
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
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, saveat=[], save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
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




end # END module LikesSpeed
