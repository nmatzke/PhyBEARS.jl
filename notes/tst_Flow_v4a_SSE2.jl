# 4a: Trying with a breakedge
# 2021-07-18_WORKS!!

#######################################################
# A test module to design FlowLikes.jl, not in scope Main
# 
# Development workflow recommended by:
#
# https://docs.julialang.org/en/v1/manual/workflow-tips/
#
# Setup:
"""
cd /GitHub/PhyBEARS.jl/notes/
julia
"""


"""
cd("/GitHub/PhyBEARS.jl/notes/")
include("tst_Flow_v4c_WORKS.jl")

"""

module Design_FlowLikes
	cd("/GitHub/PhyBEARS.jl/notes/")
	include("ModelLikes.jl")
	import .ModelLikes
	#using .Tmp

	cd("/GitHub/PhyBEARS.jl/notes/")
	include("Flow.jl")
	import .Flow

	using LinearAlgebra  # for factorize(); "I" in: Matrix{Float64}(I, 2, 2)
											 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
	using Profile     # for @profile
	using DataFrames  # for DataFrame
	using PhyloNetworks
	using PhyBEARS.TrUtils # for flat2() (similar to unlist)
	using PhyBEARS.StateSpace
	using PhyBEARS.TreePass
	using PhyBEARS.SSEs
	
	using DifferentialEquations
	using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
	#Pkg.add(PackageSpec(url="https://github.com/JuliaDiffEq/deSolveDiffEq.jl"))
	#using deSolveDiffEq 
	# https://docs.juliadiffeq.org/stable/solvers/ode_solve/index.html

	using Profile     # for @profile
	using DataFrames  # for DataFrame
	using PhyloNetworks
#	using RCall       # for df_to_Rdata, reval, g = globalEnv

	
	tr = readTopology("(((chimp:5):5,human:10):10,gorilla:20);");
	trdf = prt(tr);
	
	tipnames = ["chimp", "human", "gorilla"]
	A = ["1"; "0"; "0"]
	B = ["0"; "1"; "0"]
	geog_df = DataFrame(tipnames=tipnames, A=A, B=B);
	geog_df

	# DEC model on Hawaiian Psychotria
	bmo = construct_BioGeoBEARS_model_object()
	bmo.est[bmo.rownames .== "birthRate"] .= 0.32881638319078066
	bmo.est[bmo.rownames .== "deathRate"] .= 0.0
	bmo.est[bmo.rownames .== "d"] .= 0.03505038
	bmo.est[bmo.rownames .== "e"] .= 0.02832370
	bmo.est[bmo.rownames .== "a"] .= 0.0
	bmo.est[bmo.rownames .== "j"] .= 0.0
	birthRate = 0.32881638319078066
	numareas = 2
	n = 4            # 4 areas, 16 states

	# Set up a DEC-like model; will calculate Es over 150% of root depth
	inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true)
	(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs;
	
	#inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
	#(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;

	res = inputs.res
	trdf = inputs.trdf
	root_age = maximum(trdf[!, :node_age])
	solver_options = inputs.solver_options
	solver_options.save_everystep
	p_Ds_v5 = inputs.p_Ds_v5
	Rnames(p_Ds_v5)
	
	# Look at the model parameters (Q and C matrix)
	Rcbind(p_Ds_v5.p_indices.Qarray_ivals, p_Ds_v5.p_indices.Qarray_jvals, p_Ds_v5.params.Qij_vals)
	Rcbind(p_Ds_v5.p_indices.Carray_ivals, p_Ds_v5.p_indices.Carray_jvals, p_Ds_v5.p_indices.Carray_kvals, p_Ds_v5.params.Cijk_vals)
	
	
	# Change the d/e parameters
	d = 0.034
	e = 0.029
	p_Ds_v5.params.Qij_vals[1:2] .= d .+ 0.0
	p_Ds_v5.params.Qij_vals[3:4] .= e .+ 0.0
	Rcbind(p_Ds_v5.p_indices.Qarray_ivals, p_Ds_v5.p_indices.Qarray_jvals, p_Ds_v5.params.Qij_vals)

	# Change the mu (extinction) parameters
	p_Ds_v5.params.mu_vals .= 0.2


	# Check that the E interpolator flatted out at large times
	#	sol_Es = p_Ds_v5.sol_Es_v5
	#	sol_Es(seq(1.0, root_age, 1.0) )
	print("\nSolving the Es once, for the whole tree timespan...")
	
	# Solve the Es
	n = inputs.p_Ds_v5.n
	uE = collect(repeat([0.0], n))
	Es_tspan = (0.0, 1.4*root_age)
	p_Es_v5 = inputs.p_Ds_v5
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, uE, Es_tspan, p_Es_v5)

	print(Es_tspan)
	
	# This solution is a linear interpolator
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
	
	print("...solving Es has finished, created interpolator 'sol_Es_v5'.\n")


	# Run numerical integration of the Ds, for ground truth
	# Creates an interpolator, ground_truth_Ds_interpolator, to produce
	# Ds at any timepoint
	n = inputs.p_Ds_v5.n
	u0 = collect(repeat([0.0], n))
	u0[2] = 1.0
	tspan = (0.0, 1.2*root_age)
	
	
	# p_inputs -- add the sol_Es_v5
	p_inputs = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
	# Start with an identity matrix
	# The "I" requires "include NumericAlgebra"
	G0 = Matrix{Float64}(I, n, n) 

	# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))


	pG = (n=n, p_Ds_v5=p_inputs, A=A)
	tspan = (0.0, 20.0)
	prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG)

	Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9);
# 	Gflow_to_01_Tsit5  = solve(prob_Gs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
# 	Gflow_to_01_Lsoda  = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
# 	
# 	# Check that the different interpolators match
# 	Gflow_to_01_GMRES(0.001)
# 	Gflow_to_01_Tsit5(0.001)
# 	Gflow_to_01_Lsoda(0.001)
# 
# 	Gflow_to_01_GMRES(0.01)
# 	Gflow_to_01_Tsit5(0.01)
# 	Gflow_to_01_Lsoda(0.01)
# 
# 	Gflow_to_01_GMRES(0.1)
# 	Gflow_to_01_Tsit5(0.1)
# 	Gflow_to_01_Lsoda(0.1)
# 
# 	Gflow_to_01_GMRES(1.0)
# 	Gflow_to_01_Tsit5(1.0)
# 	Gflow_to_01_Lsoda(1.0)
# 
# 	Gflow_to_01_GMRES(10.0)
# 	Gflow_to_01_Tsit5(10.0)
# 	Gflow_to_01_Lsoda(10.0)
# 
# 
# 	Gflow_to_01_GMRES(0.0)
# 	Gflow_to_01_GMRES(0.5)
# 	Gflow_to_01_GMRES(1.0)
	
	# OK, we now have an equation that calculates the flow, G, down any timespan

	# Downpass in the standard, non-flow way
	# Insert the pre-calculated Es into the non-flow lnL calculator
	p_inputs = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
	
	res_nonFlow = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_inputs, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
	(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = res_nonFlow

	nonFlow_normlikes = res.normlikes_at_each_nodeIndex_branchTop;

	
	# Downpass with the Gflow system
	Gflow = Gflow_to_01_GMRES;

	res_wFlow = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);

	#(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = tmp_results1
	(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1) = res_wFlow
	
	wFlow_normlikes = res.normlikes_at_each_nodeIndex_branchTop;


	tmp_results2 = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf, p_Ds_v5=p_inputs, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true)
	(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1) = tmp_results2

	(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1) = tmp_results1
	
	tmp_results2
	res_nonFlow
	res_wFlow
	
	prtCp(p_Ds_v5)
	prtQp(p_Ds_v5)
	p_Ds_v5.params.mu_vals

	# Close-ish match with extinction = 0

	
	nonFlow_normlikes
	wFlow_normlikes
		
	
end

