module Optimizers
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

print("PhyBEARS: Loading Optimizers.jl dependencies...")
using Distributed				# for @everywhere
using DataFrames  			# for e.g. DataFrame()
using Dates							# for e.g. DateTime, Dates.now()
using Base.Threads			# for e.g. Base.Threads.@spawn  (NOTE: Distributed.@spawn is different, deprecated, and BAD - returns a Future)
using Random						# for MersenneTwister()
using DifferentialEquations # for ODEProblem
using LSODA							# for lsoda()
using Sundials					# for CVODE_BDF(linear_solver=:GMRES)

using PhyloBits.TrUtils		 	# for paste0
using PhyloBits.TreeTable 	# for nodetimes
using PhyBEARS.TimeDep			# for construct_QC_interpolators
using PhyBEARS.StateSpace 	# for bmo_updater_v1
using PhyBEARS.SSEs				# for parameterized_ClaSSE_Es_v7_simd_sums
using PhyBEARS.TreePass
print("...done.\n")


export func_to_optimize, func_to_optimize_nonparallel_v7, func2_EXAMPLE, func_EXAMPLE, func_to_optimize, func2_v5, func_v5, update_Qij_vals, update_Qij_vals2!, p_Ds_v5_updater_v1_OLD, bmo_updater_v1_SLOW, bmo_updater_v1, bmo_updater_v2, update_maxent01, update_Cijk_vals, update_Cijk_vals2_noUpdate, update_Qij_vals_subs!, update_Cijk_vals2!, p_Ds_v5_updater_v1!, inputs_updater_v1!, inputs_updater_v2!, bmo_updater_v1!, func_to_optimize_v7, func_to_optimize_v7c, func_to_optimize_v12







#######################################################
# OK, let's do an ML inference
#######################################################
# 1. Figure out how ML inference works
# 
# https://julianlsolvers.github.io/Optim.jl/stable/#examples/generated/maxlikenlm/
#
# 
# 2. Write function to update parameters

#######################################################
# 2022-03-06: NOW PUT UPDATER INTO A FUNCTION!!
#######################################################
function func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="lnL", printlevel=1)

	# Get the results object (with old p_Ds_v5)
	# and treetable (trdf)
	res = inputs.res
	trdf = inputs.trdf
	solver_options = inputs.solver_options
	Es_tspan = inputs.Es_tspan
	
	# Set the result if parameters out of bounds
	nan_lnL = -1000000000
	
	# Update the parameters
	inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
	
	# Get parameter names and lower/upper bounds
	parnames = inputs.bmo.rownames[inputs.bmo.type.=="free"];
	lower = inputs.bmo.min[inputs.bmo.type .== "free"]
	upper = inputs.bmo.max[inputs.bmo.type .== "free"]
	#bmo = inputs.bmo
	
	# Update
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est before inputs_updater_v1!(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est[inputs.bmo.rownames .== "j"] .= 1.6
	#inputs.bmo.est = bmo_updater_v1(inputs.bmo);
	inputs_updater_v1!(inputs);
	
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est after inputs_updater_v1!(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est .= bmo.est

	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights before p_Ds_v5_updater_v1!:")
		print(round.(pdf_before.wt[1:5], digits=4))

		print("\n")
		pdf_before = prtQp(p_Ds_v5)
		print("\nfunc_to_optimize: Qmat before p_Ds_v5_updater_v1!:")
		print(round.(pdf_before.val[1:5], digits=4))

	end
#	output = p_Ds_v5_updater_v1_OLD(p_Ds_v5, inputs);  # WORKS 2022-03-10
#	p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
#	p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
#	p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
	p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
	if printlevel >= 2
		pdf_after = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights after p_Ds_v5_updater_v1:")
		print(round.(pdf_after.wt[1:5], digits=4))

		pdf_after = prtQp(p_Ds_v5)
		print("\nfunc_to_optimize: Qmat after p_Ds_v5_updater_v1:")
		print(round.(pdf_after.val[1:5], digits=4))
	end

	
	# Check what updated params look like
	#prtQp(p_Ds_v5)
	#prtCp(p_Ds_v5)

	# Double-check that lower/upper limits haven't been breached
	inbounds = true
	for i in 1:length(parnames)
		if pars[i] < lower[i]
			inbounds = false
		end
		if pars[i] > upper[i]
			inbounds = false
		end
	end # END for i in 1:length(parnames)
	
	#sort!(trdf,"nodeIndex") # MAKE SURE this is sorted properly -- 
	# OTHERWISE I get a CRASH on 
	# iteration 1, node 19
	if inbounds == true
		# Solve the Es
		prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, inputs.Es_tspan, p_Ds_v5)
		# This solution is an interpolator
		sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
		p_Ds_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
		
		
# SEE runtests_ClaSSE_tree_n13_DECj_WORKS.jl
# save_everystep_EQ_false_CAN_MATTER_EVEN_ON_THE_Ds
#inputs.solver_options.save_everystep=false # CAN PRODUCE A -20.9 vs. -20.6 difference!
# inputs.solver_options.save_everystep=true	# WORKS!! Can make a difference EVEN ON THE Ds!!

# Slightly faster for solver to just calculate Ds at nodes
#		inputs.solver_options.save_everystep = false
#		inputs.solver_options.saveat = nodetimes(trdf)

		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
	else
		Julia_sum_lq = nan_lnL
		rootstates_lnL = nan_lnL
		Julia_total_lnLs1 = nan_lnL
		bgb_lnL = nan_lnL
	end
	
	# Make a string listing the free parameters
	txt = ""
	for i in 1:length(parnames)
		txt = paste0([txt, parnames[i], "=", round(pars[i],digits=5), ",	"])
	end # END for i in 1:length(parnames)
	
	if printlevel >= 1
		txt = paste0([txt, "Julia_sum_lq=", round(Julia_sum_lq; digits=4), ", rootstates_lnL=", round(rootstates_lnL; digits=4), ",	Julia_total_lnLs1=", round(Julia_total_lnLs1, digits=4), ", bgb_lnL=", round(bgb_lnL, digits=4)])
		print(txt) 
		print("\n")
	end
	
	if returnval == "lnL"
		return(-Julia_total_lnLs1)
	end
	if returnval == "bgb_lnL"
		return(-bgb_lnL)
	end
	if returnval == "inputs"
		# Rates vectors
		inputs.p_Ds_v5.params.mu_vals[:] .= p_Ds_v5.params.mu_vals
		inputs.p_Ds_v5.params.psi_vals[:] .= p_Ds_v5.params.psi_vals
		inputs.p_Ds_v5.params.Qij_vals[:] .= p_Ds_v5.params.Qij_vals
		inputs.p_Ds_v5.params.Cijk_weights[:] .= p_Ds_v5.params.Cijk_weights
		inputs.p_Ds_v5.params.Cijk_vals[:] .= p_Ds_v5.params.Cijk_vals
		inputs.p_Ds_v5.params.row_weightvals[:] .= p_Ds_v5.params.row_weightvals
		
		# results (res)
		inputs.res.regime[:] .= res.regime
		inputs.res.node_state[:] .= res.node_state
		inputs.res.node_Lparent_state[:] .= res.node_Lparent_state
		inputs.res.node_Rparent_state[:] .= res.node_Rparent_state
		#inputs.res.root_nodeIndex[:] .= res.root_nodeIndex
		#inputs.res.numNodes[:] .= res.numNodes
		#inputs.res.uppass_edgematrix[:] .= res.uppass_edgematrix
		inputs.res.thread_for_each_nodeOp[:] .= res.thread_for_each_nodeOp
		inputs.res.thread_for_each_branchOp[:] .= res.thread_for_each_branchOp
		inputs.res.calc_spawn_start[:] .= res.calc_spawn_start
		inputs.res.calc_start_time[:] .= res.calc_start_time
		inputs.res.calc_end_time[:] .= res.calc_end_time
		inputs.res.calc_duration[:] .= res.calc_duration
		inputs.res.calctime_iterations[:] .= res.calctime_iterations
		inputs.res.sumLikes_at_node_at_branchTop[:] .= res.sumLikes_at_node_at_branchTop
		inputs.res.lnL_at_node_at_branchTop[:] .= res.lnL_at_node_at_branchTop
		inputs.res.lq_at_branchBot[:] .= res.lq_at_branchBot
		inputs.res.like_at_branchBot[:] .= res.like_at_branchBot
		inputs.res.Es_at_each_nodeIndex_branchTop[:] .= res.Es_at_each_nodeIndex_branchTop
		inputs.res.Es_at_each_nodeIndex_branchBot[:] .= res.Es_at_each_nodeIndex_branchBot
		inputs.res.fakeX0s_at_each_nodeIndex_branchTop[:] .= res.fakeX0s_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchTop[:] .= res.likes_at_each_nodeIndex_branchTop
		inputs.res.normlikes_at_each_nodeIndex_branchTop[:] .= res.normlikes_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchBot[:] .= res.likes_at_each_nodeIndex_branchBot
		inputs.res.normlikes_at_each_nodeIndex_branchBot[:] .= res.normlikes_at_each_nodeIndex_branchBot
		
		return(inputs)
	end
	# Shouldn't get here
	return(NaN)
end # END function func_to_optimize(pars, parnames)


function func_to_optimize_nonparallel_v7(pars, parnames, inputs, p_Ds_v5; returnval="lnL", printlevel=1)

	# Get the results object (with old p_Ds_v5)
	# and treetable (trdf)
	res = inputs.res
	trdf = inputs.trdf
	solver_options = inputs.solver_options
	Es_tspan = inputs.Es_tspan
	
	# Set the result if parameters out of bounds
	nan_lnL = -1000000000
	
	# Update the parameters
	inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
	
	# Get parameter names and lower/upper bounds
	parnames = inputs.bmo.rownames[inputs.bmo.type.=="free"];
	lower = inputs.bmo.min[inputs.bmo.type .== "free"]
	upper = inputs.bmo.max[inputs.bmo.type .== "free"]
	#bmo = inputs.bmo
	
	# Update
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est before inputs_updater_v1!(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est[inputs.bmo.rownames .== "j"] .= 1.6
	#inputs.bmo.est = bmo_updater_v1(inputs.bmo);
	inputs_updater_v1!(inputs);
	
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est after inputs_updater_v1!(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est .= bmo.est

	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights before p_Ds_v5_updater_v1!:")
		print(round.(pdf_before.wt[1:5], digits=4))

		print("\n")
		pdf_before = prtQp(p_Ds_v5)
		print("\nfunc_to_optimize: Qmat before p_Ds_v5_updater_v1!:")
		print(round.(pdf_before.val[1:5], digits=4))

	end
	#	output = p_Ds_v5_updater_v1_OLD(p_Ds_v5, inputs);  # WORKS 2022-03-10
	#	p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
	#	p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
	#	p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
	p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
	if printlevel >= 2
		pdf_after = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights after p_Ds_v5_updater_v1:")
		print(round.(pdf_after.wt[1:5], digits=4))

		pdf_after = prtQp(p_Ds_v5)
		print("\nfunc_to_optimize: Qmat after p_Ds_v5_updater_v1:")
		print(round.(pdf_after.val[1:5], digits=4))
	end

	
	# Check what updated params look like
	#prtQp(p_Ds_v5)
	#prtCp(p_Ds_v5)

	# Double-check that lower/upper limits haven't been breached
	inbounds = true
	for i in 1:length(parnames)
		if pars[i] < lower[i]
			inbounds = false
		end
		if pars[i] > upper[i]
			inbounds = false
		end
	end # END for i in 1:length(parnames)
	
	#sort!(trdf,"nodeIndex") # MAKE SURE this is sorted properly -- 
	# OTHERWISE I get a CRASH on 
	# iteration 1, node 19
	if inbounds == true
		# Solve the Es
		prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, inputs.Es_tspan, p_Ds_v5)
		# This solution is an interpolator
		sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
		p_Ds_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
		
# SEE runtests_ClaSSE_tree_n13_DECj_WORKS.jl
# save_everystep_EQ_false_CAN_MATTER_EVEN_ON_THE_Ds
#inputs.solver_options.save_everystep=false # CAN PRODUCE A -20.9 vs. -20.6 difference!
# inputs.solver_options.save_everystep=true	# WORKS!! Can make a difference EVEN ON THE Ds!!

# Slightly faster for solver to just calculate Ds at nodes
#		inputs.solver_options.save_everystep = false
#		inputs.solver_options.saveat = nodetimes(trdf)

		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
	else
		Julia_sum_lq = nan_lnL
		rootstates_lnL = nan_lnL
		Julia_total_lnLs1 = nan_lnL
		bgb_lnL = nan_lnL
	end
	
	# Make a string listing the free parameters
	txt = ""
	for i in 1:length(parnames)
		txt = paste0([txt, parnames[i], "=", round(pars[i],digits=5), ",	"])
	end # END for i in 1:length(parnames)
	
	if printlevel >= 1
		txt = paste0([txt, "Julia_sum_lq=", round(Julia_sum_lq; digits=4), ", rootstates_lnL=", round(rootstates_lnL; digits=4), ",	Julia_total_lnLs1=", round(Julia_total_lnLs1, digits=4), ", bgb_lnL=", round(bgb_lnL, digits=4)])
		print("\n")
		print(txt) 
	end
	
	if returnval == "lnL"
		return(-Julia_total_lnLs1)
	end
	if returnval == "bgb_lnL"
		return(-bgb_lnL)
	end
	if returnval == "inputs"
		# Rates vectors
		inputs.p_Ds_v5.params.mu_vals[:] .= p_Ds_v5.params.mu_vals
		inputs.p_Ds_v5.params.psi_vals[:] .= p_Ds_v5.params.psi_vals
		inputs.p_Ds_v5.params.Qij_vals[:] .= p_Ds_v5.params.Qij_vals
		inputs.p_Ds_v5.params.Cijk_weights[:] .= p_Ds_v5.params.Cijk_weights
		inputs.p_Ds_v5.params.Cijk_vals[:] .= p_Ds_v5.params.Cijk_vals
		inputs.p_Ds_v5.params.row_weightvals[:] .= p_Ds_v5.params.row_weightvals
		
		# results (res)
		inputs.res.regime[:] .= res.regime
		inputs.res.node_state[:] .= res.node_state
		inputs.res.node_Lparent_state[:] .= res.node_Lparent_state
		inputs.res.node_Rparent_state[:] .= res.node_Rparent_state
		#inputs.res.root_nodeIndex[:] .= res.root_nodeIndex
		#inputs.res.numNodes[:] .= res.numNodes
		#inputs.res.uppass_edgematrix[:] .= res.uppass_edgematrix
		inputs.res.thread_for_each_nodeOp[:] .= res.thread_for_each_nodeOp
		inputs.res.thread_for_each_branchOp[:] .= res.thread_for_each_branchOp
		inputs.res.calc_spawn_start[:] .= res.calc_spawn_start
		inputs.res.calc_start_time[:] .= res.calc_start_time
		inputs.res.calc_end_time[:] .= res.calc_end_time
		inputs.res.calc_duration[:] .= res.calc_duration
		inputs.res.calctime_iterations[:] .= res.calctime_iterations
		inputs.res.sumLikes_at_node_at_branchTop[:] .= res.sumLikes_at_node_at_branchTop
		inputs.res.lnL_at_node_at_branchTop[:] .= res.lnL_at_node_at_branchTop
		inputs.res.lq_at_branchBot[:] .= res.lq_at_branchBot
		inputs.res.like_at_branchBot[:] .= res.like_at_branchBot
		inputs.res.Es_at_each_nodeIndex_branchTop[:] .= res.Es_at_each_nodeIndex_branchTop
		inputs.res.Es_at_each_nodeIndex_branchBot[:] .= res.Es_at_each_nodeIndex_branchBot
		inputs.res.fakeX0s_at_each_nodeIndex_branchTop[:] .= res.fakeX0s_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchTop[:] .= res.likes_at_each_nodeIndex_branchTop
		inputs.res.normlikes_at_each_nodeIndex_branchTop[:] .= res.normlikes_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchBot[:] .= res.likes_at_each_nodeIndex_branchBot
		inputs.res.normlikes_at_each_nodeIndex_branchBot[:] .= res.normlikes_at_each_nodeIndex_branchBot
		
		return(inputs)
	end
	# Shouldn't get here
	return(NaN)
end # END function func_to_optimize_nonparallel_v7(pars, parnames)

# CUT "_EXAMPLE", put in main script to use
func_EXAMPLE = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")

# The function used in NLopt.Opt requires a dummy "grad" input,
# as gradient-using version s of NLopt require this. Here, it's 
# just blank ([]).
# It also has to be modifiable, done with "!".
# (Kind of silly but anyways...
# https://discourse.julialang.org/t/using-nlopt-for-maximum-likelihood-estimation/31964 )
# ===================================
# "NLopt always expects two arguments objective function and 
#  if not given it silently fails because exception at the NLopt 
#  library calls are not caught but causes forced return. I also 
#  had the same problem. Make 
#  loglnorm function loglnorm(par, dummygrad) ll = 0.0 for i = 1:length(y)"
 
function func2_EXAMPLE(pars, dummy_gradient!)
	return func(pars)
end # END function func_to_optimize



function func_to_optimize_parallel_v7(pars, parnames, inputs, p_Ds_v5; returnval="lnL", printlevel=1)

	# Get the results object (with old p_Ds_v5)
	# and treetable (trdf)
	res = inputs.res
	trdf = inputs.trdf
	solver_options = inputs.solver_options
	Es_tspan = inputs.Es_tspan
	
	# Set the result if parameters out of bounds
	nan_lnL = -1000000000
	
	# Update the parameters
	inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
	
	# Get parameter names and lower/upper bounds
	parnames = inputs.bmo.rownames[inputs.bmo.type.=="free"];
	lower = inputs.bmo.min[inputs.bmo.type .== "free"]
	upper = inputs.bmo.max[inputs.bmo.type .== "free"]
	#bmo = inputs.bmo
	
	# Update
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est before inputs_updater_v1!): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est[inputs.bmo.rownames .== "j"] .= 1.6
	#inputs.bmo.est = bmo_updater_v1(inputs.bmo);
	inputs_updater_v1!(inputs);
	
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est after inputs_updater_v1!(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est .= bmo.est

	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights before p_Ds_v5_updater_v1!:")
		print(round.(pdf_before.wt[1:5], digits=4))

		print("\n")
		pdf_before = prtQp(p_Ds_v5)
		print("\nfunc_to_optimize: Qmat before p_Ds_v5_updater_v1!:")
		print(round.(pdf_before.val[1:5], digits=4))

	end
	#	output = p_Ds_v5_updater_v1_OLD(p_Ds_v5, inputs);  # WORKS 2022-03-10
	#	p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
	#	p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
	#	p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
	p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
	if printlevel >= 2
		pdf_after = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights after p_Ds_v5_updater_v1:")
		print(round.(pdf_after.wt[1:5], digits=4))

		pdf_after = prtQp(p_Ds_v5)
		print("\nfunc_to_optimize: Qmat after p_Ds_v5_updater_v1:")
		print(round.(pdf_after.val[1:5], digits=4))
	end

	
	# Check what updated params look like
	#prtQp(p_Ds_v5)
	#prtCp(p_Ds_v5)

	# Double-check that lower/upper limits haven't been breached
	inbounds = true
	for i in 1:length(parnames)
		if pars[i] < lower[i]
			inbounds = false
		end
		if pars[i] > upper[i]
			inbounds = false
		end
	end # END for i in 1:length(parnames)
	
	#sort!(trdf,"nodeIndex") # MAKE SURE this is sorted properly -- 
	# OTHERWISE I get a CRASH on 
	# iteration 1, node 19
	if inbounds == true
		# Solve the Es
		prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, inputs.Es_tspan, p_Ds_v5)
		# This solution is an interpolator
		sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
		p_Ds_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
		
		# Slightly faster for solver to just calculate Ds at nodes
# SEE runtests_ClaSSE_tree_n13_DECj_WORKS.jl
# save_everystep_EQ_false_CAN_MATTER_EVEN_ON_THE_Ds
#inputs.solver_options.save_everystep=false # CAN PRODUCE A -20.9 vs. -20.6 difference!
# inputs.solver_options.save_everystep=true	# WORKS!! Can make a difference EVEN ON THE Ds!!

# Slightly faster for solver to just calculate Ds at nodes
#		inputs.solver_options.save_everystep = false
#		inputs.solver_options.saveat = nodetimes(trdf)
		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_parallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
	else
		Julia_sum_lq = nan_lnL
		rootstates_lnL = nan_lnL
		Julia_total_lnLs1 = nan_lnL
		bgb_lnL = nan_lnL
	end
	
	# Make a string listing the free parameters
	txt = ""
	for i in 1:length(parnames)
		txt = paste0([txt, parnames[i], "=", round(pars[i],digits=5), ",	"])
	end # END for i in 1:length(parnames)
	
	if printlevel >= 1
		txt = paste0([txt, "Julia_sum_lq=", round(Julia_sum_lq; digits=4), ", rootstates_lnL=", round(rootstates_lnL; digits=4), ",	Julia_total_lnLs1=", round(Julia_total_lnLs1, digits=4), ", bgb_lnL=", round(bgb_lnL, digits=4)])
		print("\n")
		print(txt) 
	end
	
	if returnval == "lnL"
		return(-Julia_total_lnLs1)
	end
	if returnval == "bgb_lnL"
		return(-bgb_lnL)
	end
	if returnval == "inputs"
		# Rates vectors
		inputs.p_Ds_v5.params.mu_vals[:] .= p_Ds_v5.params.mu_vals
		inputs.p_Ds_v5.params.psi_vals[:] .= p_Ds_v5.params.psi_vals
		inputs.p_Ds_v5.params.Qij_vals[:] .= p_Ds_v5.params.Qij_vals
		inputs.p_Ds_v5.params.Cijk_weights[:] .= p_Ds_v5.params.Cijk_weights
		inputs.p_Ds_v5.params.Cijk_vals[:] .= p_Ds_v5.params.Cijk_vals
		inputs.p_Ds_v5.params.row_weightvals[:] .= p_Ds_v5.params.row_weightvals
		
		# results (res)
		inputs.res.regime[:] .= res.regime
		inputs.res.node_state[:] .= res.node_state
		inputs.res.node_Lparent_state[:] .= res.node_Lparent_state
		inputs.res.node_Rparent_state[:] .= res.node_Rparent_state
		#inputs.res.root_nodeIndex[:] .= res.root_nodeIndex
		#inputs.res.numNodes[:] .= res.numNodes
		#inputs.res.uppass_edgematrix[:] .= res.uppass_edgematrix
		inputs.res.thread_for_each_nodeOp[:] .= res.thread_for_each_nodeOp
		inputs.res.thread_for_each_branchOp[:] .= res.thread_for_each_branchOp
		inputs.res.calc_spawn_start[:] .= res.calc_spawn_start
		inputs.res.calc_start_time[:] .= res.calc_start_time
		inputs.res.calc_end_time[:] .= res.calc_end_time
		inputs.res.calc_duration[:] .= res.calc_duration
		inputs.res.calctime_iterations[:] .= res.calctime_iterations
		inputs.res.sumLikes_at_node_at_branchTop[:] .= res.sumLikes_at_node_at_branchTop
		inputs.res.lnL_at_node_at_branchTop[:] .= res.lnL_at_node_at_branchTop
		inputs.res.lq_at_branchBot[:] .= res.lq_at_branchBot
		inputs.res.like_at_branchBot[:] .= res.like_at_branchBot
		inputs.res.Es_at_each_nodeIndex_branchTop[:] .= res.Es_at_each_nodeIndex_branchTop
		inputs.res.Es_at_each_nodeIndex_branchBot[:] .= res.Es_at_each_nodeIndex_branchBot
		inputs.res.fakeX0s_at_each_nodeIndex_branchTop[:] .= res.fakeX0s_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchTop[:] .= res.likes_at_each_nodeIndex_branchTop
		inputs.res.normlikes_at_each_nodeIndex_branchTop[:] .= res.normlikes_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchBot[:] .= res.likes_at_each_nodeIndex_branchBot
		inputs.res.normlikes_at_each_nodeIndex_branchBot[:] .= res.normlikes_at_each_nodeIndex_branchBot
		
		return(inputs)
	end
	# Shouldn't get here
	return(NaN)
end # END function func_to_optimize(pars, parnames)

# CUT "_EXAMPLE", put in main script to use
func_EXAMPLE = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")

# The function used in NLopt.Opt requires a dummy "grad" input,
# as gradient-using version s of NLopt require this. Here, it's 
# just blank ([]).
# It also has to be modifiable, done with "!".
# (Kind of silly but anyways...
# https://discourse.julialang.org/t/using-nlopt-for-maximum-likelihood-estimation/31964 )
# ===================================
# "NLopt always expects two arguments objective function and 
#  if not given it silently fails because exception at the NLopt 
#  library calls are not caught but causes forced return. I also 
#  had the same problem. Make 
#  loglnorm function loglnorm(par, dummygrad) ll = 0.0 for i = 1:length(y)"
 
function func2_EXAMPLE(pars, dummy_gradient!)
	return func(pars)
end # END function func_to_optimize_parallel










function func_to_optimize_v7(pars, parnames, inputs, p_Ds_v5; returnval="lnL", printlevel=1)

	# Get the results object (with old p_Ds_v5)
	# and treetable (trdf)
	res = inputs.res
	trdf = inputs.trdf
	solver_options = inputs.solver_options
	
	# Set the result if parameters out of bounds
	nan_lnL = -Inf
	
	# Update the parameters
	inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
	
	# Get parameter names and lower/upper bounds
	parnames = inputs.bmo.rownames[inputs.bmo.type.=="free"];
	lower = inputs.bmo.min[inputs.bmo.type .== "free"]
	upper = inputs.bmo.max[inputs.bmo.type .== "free"]
	#bmo = inputs.bmo
	
	# Update
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est before bmo_updater_v2(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est[inputs.bmo.rownames .== "j"] .= 1.6
	#inputs.bmo.est = bmo_updater_v1(inputs.bmo);
	inputs.bmo.est = bmo_updater_v2(inputs.bmo, inputs.setup.bmo_rows);
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est after bmo_updater_v2(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est .= bmo.est

	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights before p_Ds_v5_updater_v2:")
		print(round.(pdf_before.wt[1:5], digits=4))
	end
#	output = p_Ds_v5_updater_v1_OLD(p_Ds_v5, inputs);  # WORKS 2022-03-10
#	p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
#	p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
#	p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
	p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
	if printlevel >= 2
		pdf_after = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights after p_Ds_v5_updater_v1:")
		print(round.(pdf_before.wt[1:5], digits=4))
	end

	
	# Check what updated params look like
	#prtQp(p_Ds_v5)
	#prtCp(p_Ds_v5)

	# Double-check that lower/upper limits haven't been breached
	inbounds = true
	for i in 1:length(parnames)
		if pars[i] < lower[i]
			inbounds = false
		end
		if pars[i] > upper[i]
			inbounds = false
		end
	end # END for i in 1:length(parnames)
	
	#sort!(trdf,"nodeIndex") # MAKE SURE this is sorted properly -- 
	# OTHERWISE I get a CRASH on 
	# iteration 1, node 19
	if inbounds == true
		# Solve the Es
		prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, inputs.Es_tspan, p_Ds_v5)
		# This solution is an interpolator
		sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
		p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
		
		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
	else
		Julia_sum_lq = nan_lnL
		rootstates_lnL = nan_lnL
		Julia_total_lnLs1 = nan_lnL
		bgb_lnL = nan_lnL
	end
	
	# Make a string listing the free parameters
	txt = ""
	for i in 1:length(parnames)
		txt = paste0([txt, parnames[i], "=", round(pars[i],digits=5), ",	"])
	end # END for i in 1:length(parnames)
	
	if printlevel >= 1
		txt = paste0([txt, "Julia_sum_lq=", round(Julia_sum_lq; digits=4), ", rootstates_lnL=", round(rootstates_lnL; digits=4), ",	Julia_total_lnLs1=", round(Julia_total_lnLs1, digits=4), ", bgb_lnL=", round(bgb_lnL, digits=4)])
		print("\n")
		print(txt) 
	end
	
	if returnval == "lnL"
		return(-Julia_total_lnLs1)
	end
	if returnval == "bgb_lnL"
		return(-bgb_lnL)
	end
	if returnval == "inputs"
		# Rates vectors
		inputs.p_Ds_v5.params.mu_vals[:] .= p_Ds_v5.params.mu_vals
		inputs.p_Ds_v5.params.psi_vals[:] .= p_Ds_v5.params.psi_vals
		inputs.p_Ds_v5.params.Qij_vals[:] .= p_Ds_v5.params.Qij_vals
		inputs.p_Ds_v5.params.Cijk_weights[:] .= p_Ds_v5.params.Cijk_weights
		inputs.p_Ds_v5.params.Cijk_vals[:] .= p_Ds_v5.params.Cijk_vals
		inputs.p_Ds_v5.params.row_weightvals[:] .= p_Ds_v5.params.row_weightvals
		
		# results (res)
		inputs.res.regime[:] .= res.regime
		inputs.res.node_state[:] .= res.node_state
		inputs.res.node_Lparent_state[:] .= res.node_Lparent_state
		inputs.res.node_Rparent_state[:] .= res.node_Rparent_state
		#inputs.res.root_nodeIndex[:] .= res.root_nodeIndex
		#inputs.res.numNodes[:] .= res.numNodes
		#inputs.res.uppass_edgematrix[:] .= res.uppass_edgematrix
		inputs.res.thread_for_each_nodeOp[:] .= res.thread_for_each_nodeOp
		inputs.res.thread_for_each_branchOp[:] .= res.thread_for_each_branchOp
		inputs.res.calc_spawn_start[:] .= res.calc_spawn_start
		inputs.res.calc_start_time[:] .= res.calc_start_time
		inputs.res.calc_end_time[:] .= res.calc_end_time
		inputs.res.calc_duration[:] .= res.calc_duration
		inputs.res.calctime_iterations[:] .= res.calctime_iterations
		inputs.res.sumLikes_at_node_at_branchTop[:] .= res.sumLikes_at_node_at_branchTop
		inputs.res.lnL_at_node_at_branchTop[:] .= res.lnL_at_node_at_branchTop
		inputs.res.lq_at_branchBot[:] .= res.lq_at_branchBot
		inputs.res.like_at_branchBot[:] .= res.like_at_branchBot
		inputs.res.Es_at_each_nodeIndex_branchTop[:] .= res.Es_at_each_nodeIndex_branchTop
		inputs.res.Es_at_each_nodeIndex_branchBot[:] .= res.Es_at_each_nodeIndex_branchBot
		inputs.res.fakeX0s_at_each_nodeIndex_branchTop[:] .= res.fakeX0s_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchTop[:] .= res.likes_at_each_nodeIndex_branchTop
		inputs.res.normlikes_at_each_nodeIndex_branchTop[:] .= res.normlikes_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchBot[:] .= res.likes_at_each_nodeIndex_branchBot
		inputs.res.normlikes_at_each_nodeIndex_branchBot[:] .= res.normlikes_at_each_nodeIndex_branchBot
		
		return(inputs)
	end
	# Shouldn't get here
	return(NaN)
end # END function func_to_optimize_v7(pars, parnames)





# Distributed calculation; requires multiple cores and processes,
# "@everywhere using" instead of "using..."
function func_to_optimize_v7c(pars, parnames, inputs, p_Ds_v5; returnval="lnL", printlevel=1)

	# Get the results object (with old p_Ds_v5)
	# and treetable (trdf)
	res = inputs.res
	trdf = inputs.trdf
	solver_options = inputs.solver_options
	Es_tspan = inputs.Es_tspan
	
	# Set the result if parameters out of bounds
	nan_lnL = -1000000000
	
	# Update the parameters
	inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
	
	# Get parameter names and lower/upper bounds
	parnames = inputs.bmo.rownames[inputs.bmo.type.=="free"];
	lower = inputs.bmo.min[inputs.bmo.type .== "free"]
	upper = inputs.bmo.max[inputs.bmo.type .== "free"]
	#bmo = inputs.bmo
	
	# Update
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est before bmo_updater_v2(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est[inputs.bmo.rownames .== "j"] .= 1.6
	#inputs.bmo.est = bmo_updater_v1(inputs.bmo);
	inputs_updater_v2!(inputs);
	
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est after bmo_updater_v2(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est .= bmo.est

	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights before p_Ds_v5_updater_v1:")
		print(round.(pdf_before.wt[1:5], digits=4))

		print("\n")
		pdf_before = prtQp(p_Ds_v5)
		print("\nfunc_to_optimize: Qmat before p_Ds_v5_updater_v1:")
		print(round.(pdf_before.val[1:5], digits=4))

	end
	#	output = p_Ds_v5_updater_v1_OLD(p_Ds_v5, inputs);  # WORKS 2022-03-10
	#	p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
	#	p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
	#	p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
	p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
	if printlevel >= 2
		pdf_after = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights after p_Ds_v5_updater_v1:")
		print(round.(pdf_after.wt[1:5], digits=4))

		pdf_after = prtQp(p_Ds_v5)
		print("\nfunc_to_optimize: Qmat after p_Ds_v5_updater_v1:")
		print(round.(pdf_after.val[1:5], digits=4))
	end

	
	# Check what updated params look like
	#prtQp(p_Ds_v5)
	#prtCp(p_Ds_v5)

	# Double-check that lower/upper limits haven't been breached
	inbounds = true
	for i in 1:length(parnames)
		if pars[i] < lower[i]
			inbounds = false
		end
		if pars[i] > upper[i]
			inbounds = false
		end
	end # END for i in 1:length(parnames)
	
	#sort!(trdf,"nodeIndex") # MAKE SURE this is sorted properly -- 
	# OTHERWISE I get a CRASH on 
	# iteration 1, node 19
	if inbounds == true
		# Solve the Es
		prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, inputs.Es_tspan, p_Ds_v5)
		# This solution is an interpolator
		sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
		Distributed.@everywhere p_Ds_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
		
# SEE runtests_ClaSSE_tree_n13_DECj_WORKS.jl
# save_everystep_EQ_false_CAN_MATTER_EVEN_ON_THE_Ds
#inputs.solver_options.save_everystep=false # CAN PRODUCE A -20.9 vs. -20.6 difference!
# inputs.solver_options.save_everystep=true	# WORKS!! Can make a difference EVEN ON THE Ds!!

# Slightly faster for solver to just calculate Ds at nodes
#		inputs.solver_options.save_everystep = false
#		inputs.solver_options.saveat = nodetimes(trdf)
		
		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_parallel_ClaSSE_v7c!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
	else
		Julia_sum_lq = nan_lnL
		rootstates_lnL = nan_lnL
		Julia_total_lnLs1 = nan_lnL
		bgb_lnL = nan_lnL
	end
	
	# Make a string listing the free parameters
	txt = ""
	for i in 1:length(parnames)
		txt = paste0([txt, parnames[i], "=", round(pars[i],digits=5), ",	"])
	end # END for i in 1:length(parnames)
	
	if printlevel >= 1
		txt = paste0([txt, "Julia_sum_lq=", round(Julia_sum_lq; digits=4), ", rootstates_lnL=", round(rootstates_lnL; digits=4), ",	Julia_total_lnLs1=", round(Julia_total_lnLs1, digits=4), ", bgb_lnL=", round(bgb_lnL, digits=4)])
		print("\n")
		print(txt) 
	end
	
	if returnval == "lnL"
		return(-Julia_total_lnLs1)
	end
	if returnval == "bgb_lnL"
		return(-bgb_lnL)
	end
	if returnval == "inputs"
		# Rates vectors
		inputs.p_Ds_v5.params.mu_vals[:] .= p_Ds_v5.params.mu_vals
		inputs.p_Ds_v5.params.psi_vals[:] .= p_Ds_v5.params.psi_vals
		inputs.p_Ds_v5.params.Qij_vals[:] .= p_Ds_v5.params.Qij_vals
		inputs.p_Ds_v5.params.Cijk_weights[:] .= p_Ds_v5.params.Cijk_weights
		inputs.p_Ds_v5.params.Cijk_vals[:] .= p_Ds_v5.params.Cijk_vals
		inputs.p_Ds_v5.params.row_weightvals[:] .= p_Ds_v5.params.row_weightvals
		
		# results (res)
		inputs.res.regime[:] .= res.regime
		inputs.res.node_state[:] .= res.node_state
		inputs.res.node_Lparent_state[:] .= res.node_Lparent_state
		inputs.res.node_Rparent_state[:] .= res.node_Rparent_state
		#inputs.res.root_nodeIndex[:] .= res.root_nodeIndex
		#inputs.res.numNodes[:] .= res.numNodes
		#inputs.res.uppass_edgematrix[:] .= res.uppass_edgematrix
		inputs.res.thread_for_each_nodeOp[:] .= res.thread_for_each_nodeOp
		inputs.res.thread_for_each_branchOp[:] .= res.thread_for_each_branchOp
		inputs.res.calc_spawn_start[:] .= res.calc_spawn_start
		inputs.res.calc_start_time[:] .= res.calc_start_time
		inputs.res.calc_end_time[:] .= res.calc_end_time
		inputs.res.calc_duration[:] .= res.calc_duration
		inputs.res.calctime_iterations[:] .= res.calctime_iterations
		inputs.res.sumLikes_at_node_at_branchTop[:] .= res.sumLikes_at_node_at_branchTop
		inputs.res.lnL_at_node_at_branchTop[:] .= res.lnL_at_node_at_branchTop
		inputs.res.lq_at_branchBot[:] .= res.lq_at_branchBot
		inputs.res.like_at_branchBot[:] .= res.like_at_branchBot
		inputs.res.Es_at_each_nodeIndex_branchTop[:] .= res.Es_at_each_nodeIndex_branchTop
		inputs.res.Es_at_each_nodeIndex_branchBot[:] .= res.Es_at_each_nodeIndex_branchBot
		inputs.res.fakeX0s_at_each_nodeIndex_branchTop[:] .= res.fakeX0s_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchTop[:] .= res.likes_at_each_nodeIndex_branchTop
		inputs.res.normlikes_at_each_nodeIndex_branchTop[:] .= res.normlikes_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchBot[:] .= res.likes_at_each_nodeIndex_branchBot
		inputs.res.normlikes_at_each_nodeIndex_branchBot[:] .= res.normlikes_at_each_nodeIndex_branchBot
		
		return(inputs)
	end
	# Shouldn't get here
	return(NaN)
end # END function func_to_optimize_v7c(pars, parnames)




func_v5 = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")

# The function used in NLopt.Opt requires a dummy "grad" input,
# as gradient-using version s of NLopt require this. Here, it's 
# just blank ([]).
# It also has to be modifiable, done with "!".
# (Kind of silly but anyways...
# https://discourse.julialang.org/t/using-nlopt-for-maximum-likelihood-estimation/31964 )
# ===================================
# "NLopt always expects two arguments objective function and 
#  if not given it silently fails because exception at the NLopt 
#  library calls are not caught but causes forced return. I also 
#  had the same problem. Make 
#  loglnorm function loglnorm(par, dummygrad) ll = 0.0 for i = 1:length(y)"
 
function func2_v5(pars, dummy_gradient!)
	return func_v5(pars)
end # END function func2(pars, dummy_gradient!)



function func_to_optimize_v12(pars, parnames, inputs, p_Ds_v12; returnval="lnL", printlevel=1)

	# Get the results object (with old p_Ds_v12)
	# and treetable (trdf)
	res = inputs.res
	trdf = inputs.trdf
	solver_options = inputs.solver_options
	Es_tspan = inputs.Es_tspan

	
	# Set the result if parameters out of bounds
	nan_lnL = -Inf
	
	# Update the parameters
	inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
	
	# Get parameter names and lower/upper bounds
	parnames = inputs.bmo.rownames[inputs.bmo.type.=="free"];
	lower = inputs.bmo.min[inputs.bmo.type .== "free"]
	upper = inputs.bmo.max[inputs.bmo.type .== "free"]
	#bmo = inputs.bmo
	
	# Update
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est before bmo_updater_v2(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est[inputs.bmo.rownames .== "j"] .= 1.6
	#inputs.bmo.est = bmo_updater_v1(inputs.bmo);
	inputs.bmo.est = bmo_updater_v2(inputs.bmo, inputs.setup.bmo_rows);
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est after bmo_updater_v2(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est .= bmo.est

	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v12)
		print("\nfunc_to_optimize: clado weights before p_Ds_v5_updater_v2:")
		print(round.(pdf_before.wt[1:5], digits=4))
	end
#	output = p_Ds_v5_updater_v1_OLD(p_Ds_v12, inputs);  # WORKS 2022-03-10
#	p_Ds_v12.params.Cijk_weights[:] .= output.Cijk_weights
#	p_Ds_v12.params.Cijk_vals[:] .= output.Cijk_vals
#	p_Ds_v12.params.row_weightvals[:] .= output.row_weightvals
	p_Ds_v5_updater_v1!(p_Ds_v12, inputs);
	if printlevel >= 2
		pdf_after = prtCp(p_Ds_v12)
		print("\nfunc_to_optimize: clado weights after p_Ds_v5_updater_v1:")
		print(round.(pdf_before.wt[1:5], digits=4))
	end

	
	# Check what updated params look like
	#prtQp(p_Ds_v12)
	#prtCp(p_Ds_v12)

	# Double-check that lower/upper limits haven't been breached
	inbounds = true
	for i in 1:length(parnames)
		if pars[i] < lower[i]
			inbounds = false
		end
		if pars[i] > upper[i]
			inbounds = false
		end
	end # END for i in 1:length(parnames)
	
	#sort!(trdf,"nodeIndex") # MAKE SURE this is sorted properly -- 
	# OTHERWISE I get a CRASH on 
	# iteration 1, node 19
	if inbounds == true
		p_Es_v12 = TimeDep.construct_QC_interpolators(p_Ds_v12, p_Ds_v12.interpolators.times_for_SSE_interpolators);
		
		# Solve the Es
		prob_Es_v12 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, inputs.Es_tspan, p_Es_v12)
		# This solution is an interpolator
		sol_Es_v12 = solve(prob_Es_v12, inputs.solver_options.solver, save_everystep=true, abstol=inputs.solver_options.abstol, reltol=inputs.solver_options.reltol);
		p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

		# Calculate the Ds
		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
	else
		Julia_sum_lq = nan_lnL
		rootstates_lnL = nan_lnL
		Julia_total_lnLs1 = nan_lnL
		bgb_lnL = nan_lnL
	end
	
	# Make a string listing the free parameters
	txt = ""
	for i in 1:length(parnames)
		txt = paste0([txt, parnames[i], "=", round(pars[i],digits=5), ",	"])
	end # END for i in 1:length(parnames)
	
	if printlevel >= 1
		txt = paste0([txt, "Julia_sum_lq=", round(Julia_sum_lq; digits=4), ", rootstates_lnL=", round(rootstates_lnL; digits=4), ",	Julia_total_lnLs1=", round(Julia_total_lnLs1, digits=4), ", bgb_lnL=", round(bgb_lnL, digits=4)])
		print("\n")
		print(txt) 
	end
	
	if returnval == "lnL"
		return(-Julia_total_lnLs1)
	end
	if returnval == "bgb_lnL"
		return(-bgb_lnL)
	end
	if returnval == "inputs"
		# Rates vectors
		inputs.p_Ds_v12.params.mu_vals[:] .= p_Ds_v12.params.mu_vals
		inputs.p_Ds_v12.params.psi_vals[:] .= p_Ds_v12.params.psi_vals
		inputs.p_Ds_v12.params.Qij_vals[:] .= p_Ds_v12.params.Qij_vals
		inputs.p_Ds_v12.params.Cijk_weights[:] .= p_Ds_v12.params.Cijk_weights
		inputs.p_Ds_v12.params.Cijk_vals[:] .= p_Ds_v12.params.Cijk_vals
		inputs.p_Ds_v12.params.row_weightvals[:] .= p_Ds_v12.params.row_weightvals
		
		# results (res)
		inputs.res.regime[:] .= res.regime
		inputs.res.node_state[:] .= res.node_state
		inputs.res.node_Lparent_state[:] .= res.node_Lparent_state
		inputs.res.node_Rparent_state[:] .= res.node_Rparent_state
		#inputs.res.root_nodeIndex[:] .= res.root_nodeIndex
		#inputs.res.numNodes[:] .= res.numNodes
		#inputs.res.uppass_edgematrix[:] .= res.uppass_edgematrix
		inputs.res.thread_for_each_nodeOp[:] .= res.thread_for_each_nodeOp
		inputs.res.thread_for_each_branchOp[:] .= res.thread_for_each_branchOp
		inputs.res.calc_spawn_start[:] .= res.calc_spawn_start
		inputs.res.calc_start_time[:] .= res.calc_start_time
		inputs.res.calc_end_time[:] .= res.calc_end_time
		inputs.res.calc_duration[:] .= res.calc_duration
		inputs.res.calctime_iterations[:] .= res.calctime_iterations
		inputs.res.sumLikes_at_node_at_branchTop[:] .= res.sumLikes_at_node_at_branchTop
		inputs.res.lnL_at_node_at_branchTop[:] .= res.lnL_at_node_at_branchTop
		inputs.res.lq_at_branchBot[:] .= res.lq_at_branchBot
		inputs.res.like_at_branchBot[:] .= res.like_at_branchBot
		inputs.res.Es_at_each_nodeIndex_branchTop[:] .= res.Es_at_each_nodeIndex_branchTop
		inputs.res.Es_at_each_nodeIndex_branchBot[:] .= res.Es_at_each_nodeIndex_branchBot
		inputs.res.fakeX0s_at_each_nodeIndex_branchTop[:] .= res.fakeX0s_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchTop[:] .= res.likes_at_each_nodeIndex_branchTop
		inputs.res.normlikes_at_each_nodeIndex_branchTop[:] .= res.normlikes_at_each_nodeIndex_branchTop
		inputs.res.likes_at_each_nodeIndex_branchBot[:] .= res.likes_at_each_nodeIndex_branchBot
		inputs.res.normlikes_at_each_nodeIndex_branchBot[:] .= res.normlikes_at_each_nodeIndex_branchBot
		
		return(inputs)
	end
	# Shouldn't get here
	return(NaN)
end # END function func_to_optimize_v12(pars, parnames)








"""
# Update Qij_vals
# Takes dmat, amat, elist as inputs
# dmat is d * various dispersal multipliers
# amat is a * various dispersal multipliers
# elist is "e" for each area

# Update Qij_vals
numareas = 3
areas_list = collect(1:numareas)
states_list = areas_list_to_states_list(areas_list, 3, true)
numstates = length(states_list)
amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
elist = repeat([0.123], numstates)
allowed_event_types=["d","e"]

Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
Qarray_event_types = Qmat.Qarray_event_types
Qmat1_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

# Update!
dmat = reshape(repeat([0.5], numareas^2), (numareas,numareas))
Qmat2 = update_Qij_vals(Qmat, areas_list, states_list, dmat, elist, amat )
Qmat2

Qarray_ivals = Qmat2.Qarray_ivals
Qarray_jvals = Qmat2.Qarray_jvals
Qij_vals = Qmat2.Qij_vals
Qarray_event_types = Qmat2.Qarray_event_types
Qmat2_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

Qmat1_df
Qmat2_df

"""
function update_Qij_vals(Qmat, areas_list, states_list, dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))), elist=repeat([1.0], length(areas_list)), amat=dmat; return_df=false)

	numstates = length(states_list)
	statenums = collect(1:numstates)
	range_sizes = length.(states_list)
	#areas_list = sort(unique(flat2(states_list)))
	numareas = length(areas_list)
	
	if (Rclass(Qmat) == "DataFrame")
		Qarray_ivals = Qmat[!,:i]
		Qarray_jvals = Qmat[!,:j]
		Qij_vals = Qmat[!,:val]
		Qarray_event_types = Qmat[!,:event]
	else
		Qarray_ivals = Qmat.Qarray_ivals
		Qarray_jvals = Qmat.Qarray_jvals
		Qij_vals = Qmat.Qij_vals
		Qarray_event_types = Qmat.Qarray_event_types
	end
	
	# Update the "d" events (anagenetic range expansion)
	TF = Qarray_event_types .== "d"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			end_areanums_not_found_in_start_areas = setdiff(ending_areanums, starting_areanums)
			
			# Add up the d events
			tmp_d_sum = 0.0
			for k in 1:ancsize
				# Because there is only 1 end_areanums_not_found_in_start_areas
				tmp_d_sum += dmat[starting_areanums[k], end_areanums_not_found_in_start_areas[1]][]
			end
			# Store the result
			rates[z] = tmp_d_sum
		end
		Qij_vals[TF] = rates
	end # End update of d event weights


	# Update the "a" events (anagenetic range expansion)
	TF = Qarray_event_types .== "a"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanum = ancstate # because it has size=1 by def.
			ending_areanum = decstate   # because it has size=1 by def.

			# Store the result
			rates[z] = amat[starting_areanum,ending_areanum]
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Update the "e" events (anagenetic range extinction/
	# local extirpation/range loss)
	TF = Qarray_event_types .== "e"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			start_areanums_not_found_in_end_areas = setdiff(starting_areanums, ending_areanums)

			# Store the result
			rates[z] = elist[start_areanums_not_found_in_end_areas[1]]
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Return results
	if return_df == true
		Qmat2 = DataFrames.DataFrame(Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals, Qarray_event_types=Qarray_event_types)
	else
		Qmat2 = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals, Qarray_event_types=Qarray_event_types)
	end
	return Qmat2
end # end function update_Qij_vals

"""
# update_Qij_vals_subs!(p_Ds_v5)
# 
# Updates: Qij_vals_sub_i, Qij_vals_sub_j. 
# 
# For example: after manually changing the Qij_vals
#
# Example:
p_Ds_v5 = (n = 2, params = (mu_vals = [0.111111111, 0.05], Qij_vals = [0.1, 0.15], Qij_vals_t = [0.1, 0.1], Cijk_weights = [1.0, 1.0], Cijk_probs = [1.0, 1.0], Cijk_vals = [0.222222222, 0.222222222], Cijk_rates = [0.222222222, 0.222222222], Cijk_rates_t = [0.222222222, 0.222222222], row_weightvals = [1.0, 1.0]), p_indices = (Qarray_ivals = [1, 2], Qarray_jvals = [2, 1], Qarray_event_types = ["a", "a"], Carray_ivals = [1, 2], Carray_jvals = [1, 2], Carray_kvals = [1, 2], Carray_pair = [1, 1], Carray_event_types = ["y", "y"]), p_TFs = (Qi_eq_i = Vector{Bool}[[1, 0], [0, 1]], Qi_eq_i_index = [[1], [2]], Ci_eq_i = Vector{Bool}[[1, 0], [0, 1]], Ci_eq_i_index = [[1], [2]], Qi_sub_i = [[1], [2]], Qj_sub_i = [[2], [1]], Qi_sub_j = [[2], [1]], Qij_vals_sub_i = [[0.1], [0.1]], Qij_vals_sub_i_t = [[0.1], [0.1]], Qji_vals_sub_j = [[0.15], [0.1]], Qji_vals_sub_j_t = [[0.1], [0.1]], Ci_sub_i = [[1], [2]], Cj_sub_i = [[1], [2]], Ck_sub_i = [[1], [2]], Ci_sub_j = [[1], [2]], Cj_sub_j = [[1], [2]], Ck_sub_j = [[1], [2]], Qij_singleNum_sub_i = Vector{Any}[[3], [2]], Cij_singleNum_sub_i = Vector{Any}[[1], [4]], Cik_singleNum_sub_i = Vector{Any}[[1], [4]], Cijk_not_y_sub_i = Vector{Bool}[[0], [0]], Cijk_pair_sub_i = [[1], [1]], Cijk_rates_sub_i = [[0.222222222], [0.222222222]], Cijk_rates_sub_i_t = [[0.222222222], [0.222222222]], Cjik_rates_sub_j = [[0.222222222], [0.222222222]], Cjik_rates_sub_j_t = [[0.222222222], [0.222222222]]), uE = [0.0, 0.0])

# Manual changes
p_Ds_v5.params.Cijk_vals[1] = 0.222222222
p_Ds_v5.params.Cijk_vals[2] = 0.222222222
p_Ds_v5.params.mu_vals[1] = 0.111111111
p_Ds_v5.params.mu_vals[2] = 0.05
p_Ds_v5.params.Qij_vals[1] = 0.1
p_Ds_v5.params.Qij_vals[2] = 0.15

# View the subs
inputs.p_Ds_v5.p_TFs.Qij_vals_sub_i
inputs.p_Ds_v5.p_TFs.Qji_vals_sub_j

# Update
update_Qij_vals_subs!(p_Ds_v5)

# View the updated subs
inputs.p_Ds_v5.p_TFs.Qij_vals_sub_i
inputs.p_Ds_v5.p_TFs.Qji_vals_sub_j
"""
function update_Qij_vals_subs!(p_Ds_v5)

	# Update Qij_vals_sub_i
	for i in 1:p_Ds_v5.n
		# 2023-01-22 Qi_sub_i replaces Qi_eq_i
		# 2023-01-30_THIS WAS A HORRIBLE MISTAKE, GO BACK

		# For this:
		# 6×5 DataFrame
		#  Row │ event   i      j      val       vals_t  
		#      │ String  Int64  Int64  Float64   Float64 
		# ─────┼─────────────────────────────────────────
		#    1 │ d           2      4  0.101056      0.0
		#    2 │ d           3      4  0.101056      0.0
		#    3 │ e           2      1  1.0e-12       0.0
		#    4 │ e           3      1  1.0e-12       0.0
		#    5 │ e           4      2  1.0e-12       0.0
		#    6 │ e           4      3  1.0e-12       0.0

		# Qi_sub_i, when i=2, gives "2 2"
		# Qi_eq_i, when i=2, gives TFs for rows 1 & 3
		# Qi_sub_j, when j=2, gives "4" (the i resulting when j==2)
		# But Qj_eq_j, when j=2, gives TFs for row 5

		#p_Ds_v5.p_TFs.Qij_vals_sub_i[i] .= p_Ds_v5.params.Qij_vals[p_Ds_v5.p_TFs.Qi_sub_i[i]]
		p_Ds_v5.p_TFs.Qij_vals_sub_i[i] .= p_Ds_v5.params.Qij_vals[p_Ds_v5.p_TFs.Qi_eq_i[i]]
		# 2023-01-22 Qi_sub_i replaces Qi_eq_i  (PS: Yes, do Qji_vals not Qij_vals)
		# 2023-01-30 NO! Use Qj_eq_j, see above
		p_Ds_v5.p_TFs.Qji_vals_sub_j[i] .= p_Ds_v5.params.Qij_vals[p_Ds_v5.p_TFs.Qj_eq_j[i]]
	end
	
end # END update_Qij_vals_subs!(p_Ds_v5)



"""
# Update Qij_vals2!
# Takes dmat, amat, elist as inputs
# dmat is d * various dispersal multipliers
# amat is a * various dispersal multipliers
# elist is "e" for each area

# Update Qij_vals2

# Assumes a previously set-up p_Ds_v5 for a 3-area system (see below for run-from-scratch)
numareas = 3
areas_list = collect(1:numareas)
states_list = areas_list_to_states_list(areas_list, 3, true)
numstates = length(states_list)
amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
elist = repeat([0.123], numstates)
allowed_event_types=["d","e"]

Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
Qarray_event_types = Qmat.Qarray_event_types
Qmat1_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

# Update!
dmat = reshape(repeat([0.5], numareas^2), (numareas,numareas))
Qmat2 = update_Qij_vals2!(Qmat, areas_list, states_list, dmat, elist, amat )
Qmat2

Qarray_ivals = Qmat2.Qarray_ivals
Qarray_jvals = Qmat2.Qarray_jvals
Qij_vals = Qmat2.Qij_vals
Qarray_event_types = Qmat2.Qarray_event_types
Qmat2_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

Qmat1_df
Qmat2_df


# Setup for a 2-area, MuSSE (ie BiSSE) model
p_Ds_v5 = (n = 2, params = (mu_vals = [0.111111111, 0.05], Qij_vals = [0.1, 0.15], Qij_vals_t = [0.1, 0.1], Cijk_weights = [1.0, 1.0], Cijk_probs = [1.0, 1.0], Cijk_vals = [0.222222222, 0.222222222], Cijk_rates = [0.222222222, 0.222222222], Cijk_rates_t = [0.222222222, 0.222222222], row_weightvals = [1.0, 1.0]), p_indices = (Qarray_ivals = [1, 2], Qarray_jvals = [2, 1], Qarray_event_types = ["a", "a"], Carray_ivals = [1, 2], Carray_jvals = [1, 2], Carray_kvals = [1, 2], Carray_pair = [1, 1], Carray_event_types = ["y", "y"]), p_TFs = (Qi_eq_i = Vector{Bool}[[1, 0], [0, 1]], Qi_eq_i_index = [[1], [2]], Ci_eq_i = Vector{Bool}[[1, 0], [0, 1]], Ci_eq_i_index = [[1], [2]], Qi_sub_i = [[1], [2]], Qj_sub_i = [[2], [1]], Qi_sub_j = [[2], [1]], Qij_vals_sub_i = [[0.1], [0.1]], Qij_vals_sub_i_t = [[0.1], [0.1]], Qji_vals_sub_j = [[0.15], [0.1]], Qji_vals_sub_j_t = [[0.1], [0.1]], Ci_sub_i = [[1], [2]], Cj_sub_i = [[1], [2]], Ck_sub_i = [[1], [2]], Ci_sub_j = [[1], [2]], Cj_sub_j = [[1], [2]], Ck_sub_j = [[1], [2]], Qij_singleNum_sub_i = Vector{Any}[[3], [2]], Cij_singleNum_sub_i = Vector{Any}[[1], [4]], Cik_singleNum_sub_i = Vector{Any}[[1], [4]], Cijk_not_y_sub_i = Vector{Bool}[[0], [0]], Cijk_pair_sub_i = [[1], [1]], Cijk_rates_sub_i = [[0.222222222], [0.222222222]], Cijk_rates_sub_i_t = [[0.222222222], [0.222222222]], Cjik_rates_sub_j = [[0.222222222], [0.222222222]], Cjik_rates_sub_j_t = [[0.222222222], [0.222222222]]), uE = [0.0, 0.0])

numareas = p_Ds_v5.n
areas_list = collect(1:numareas)
states_list = areas_list_to_states_list(areas_list, 1, false)
numstates = length(states_list)
amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
elist = repeat([0.123], numstates)
allowed_event_types=["a"]
return_df=false

Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=allowed_event_types)
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
Qarray_event_types = Qmat.Qarray_event_types
Qmat1_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

# Update!
dmat = reshape(repeat([0.5], numareas^2), (numareas,numareas))
Qmat2 = update_Qij_vals2!(p_Ds_v5, areas_list, states_list, dmat, elist, amat; return_df=return_df)
Qmat2

Qarray_ivals = Qmat2.p_indices.Qarray_ivals
Qarray_jvals = Qmat2.p_indices.Qarray_jvals
Qij_vals = Qmat2.params.Qij_vals
Qarray_event_types = Qmat2.p_indices.Qarray_event_types
Qmat2_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

Qmat1_df
Qmat2_df


# Update just the subs, update_Qij_vals_subs!

# Manual changes
p_Ds_v5.params.Cijk_vals[1] = 0.222222222
p_Ds_v5.params.Cijk_vals[2] = 0.222222222
p_Ds_v5.params.mu_vals[1] = 0.111111111
p_Ds_v5.params.mu_vals[2] = 0.05
p_Ds_v5.params.Qij_vals[1] = 0.1
p_Ds_v5.params.Qij_vals[2] = 0.15

# View the subs
inputs.p_Ds_v5.p_TFs.Qij_vals_sub_i
inputs.p_Ds_v5.p_TFs.Qji_vals_sub_j

# Update
update_Qij_vals_subs!(p_Ds_v5)

# View the updated subs
inputs.p_Ds_v5.p_TFs.Qij_vals_sub_i
inputs.p_Ds_v5.p_TFs.Qji_vals_sub_j

"""
function update_Qij_vals2!(p_Ds_v5, areas_list, states_list, dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))), elist=repeat([1.0], length(areas_list)), amat=dmat; return_df=false)

	numstates = length(states_list)
	statenums = collect(1:numstates)
	range_sizes = length.(states_list)
	#areas_list = sort(unique(flat2(states_list)))
	numareas = length(areas_list)
	
	Qarray_ivals = p_Ds_v5.p_indices.Qarray_ivals
	Qarray_jvals = p_Ds_v5.p_indices.Qarray_jvals
	Qij_vals = p_Ds_v5.params.Qij_vals
	Qarray_event_types = p_Ds_v5.p_indices.Qarray_event_types
	
	# Update the "d" events (anagenetic range expansion)
	TF = Qarray_event_types .== "d"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			end_areanums_not_found_in_start_areas = setdiff(ending_areanums, starting_areanums)
			
			# Add up the d events
			tmp_d_sum = 0.0
			for k in 1:ancsize
				# Because there is only 1 end_areanums_not_found_in_start_areas -- 
				# (all d dispersals add just 1 area)
				tmp_d_sum += dmat[starting_areanums[k], end_areanums_not_found_in_start_areas[1]][]
			end
			# Store the result
			rates[z] = tmp_d_sum
		end
		Qij_vals[TF] = rates
	end # End update of d event weights


	# Update the "a" events (anagenetic range expansion)
	TF = Qarray_event_types .== "a"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanum = ancstate # because it has size=1 by def.
			ending_areanum = decstate   # because it has size=1 by def.

			# Store the result
			rates[z] = amat[starting_areanum,ending_areanum][]  # 2023-01-22 needs the [] to de-Matrixify
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Update the "e" events (anagenetic range extinction/
	# local extirpation/range loss)
	TF = Qarray_event_types .== "e"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			# Because only 1 area is lost at a time under this model
			start_areanums_not_found_in_end_areas = setdiff(starting_areanums, ending_areanums)

			# Store the result
			rates[z] = elist[start_areanums_not_found_in_end_areas[1]]
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Update Qij_vals_sub_i
	for i in 1:p_Ds_v5.n
		# 2023-01-22 Qi_sub_i replaces Qi_eq_i
		# 2023-01-30_THIS WAS A HORRIBLE MISTAKE, GO BACK

		# For this:
		# 6×5 DataFrame
		#  Row │ event   i      j      val       vals_t  
		#      │ String  Int64  Int64  Float64   Float64 
		# ─────┼─────────────────────────────────────────
		#    1 │ d           2      4  0.101056      0.0
		#    2 │ d           3      4  0.101056      0.0
		#    3 │ e           2      1  1.0e-12       0.0
		#    4 │ e           3      1  1.0e-12       0.0
		#    5 │ e           4      2  1.0e-12       0.0
		#    6 │ e           4      3  1.0e-12       0.0

		# Qi_sub_i, when i=2, gives "2 2"
		# Qi_eq_i, when i=2, gives TFs for rows 1 & 3
		# Qi_sub_j, when j=2, gives "4" (the i resulting when j==2)
		# But Qj_eq_j, when j=2, gives TFs for row 5
				
		#p_Ds_v5.p_TFs.Qij_vals_sub_i[i] .= Qij_vals[p_Ds_v5.p_TFs.Qi_sub_i[i]]
		p_Ds_v5.p_TFs.Qij_vals_sub_i[i] .= Qij_vals[p_Ds_v5.p_TFs.Qi_eq_i[i]]
		# 2023-01-22 Qi_sub_i replaces Qi_eq_i  (PS: Yes, do Qji_vals not Qij_vals)
		# 2023-01-30 NO! Use Qj_eq_j, see above
		p_Ds_v5.p_TFs.Qji_vals_sub_j[i] .= Qij_vals[p_Ds_v5.p_TFs.Qj_eq_j[i]]
	end


	# Return results
	if return_df == true
		Qmat2 = DataFrames.DataFrame(Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals, Qarray_event_types=Qarray_event_types, Qij_vals_sub_i=Qij_vals_sub_i)
		return Qmat2
	else
		p_Ds_v5.params.Qij_vals[:] .= Qij_vals
		return p_Ds_v5
	end
	return NaN  # shouldn't get here
end # end function update_Qij_vals2!

"""
# Update the BioGeoBEARS_model_object after changing e.g. "j"
# Full updater (as of 2023-02-02)
#
# bmo_updater_v1_SLOW -- full updater, slow due to lots of == 
# bmo_updater_v1! -- very fast, but requires inputs.setup.bmo_rows, and doesn't update "u" etc.
# bmo_updater_v2 -- 2x slower, requires inputs.setup.bmo_rows, but *DOES* update "u" etc.
# 
"""
function bmo_updater_v1_SLOW(bmo)
	# Update cladogenesis parameters
	j_wt = bmo.est[bmo.rownames .== "j"][1]
	ysv = bmo.est[bmo.rownames .== "ysv"][1]
	ys = bmo.est[bmo.rownames .== "ys"][1]
	y = bmo.est[bmo.rownames .== "y"][1]
	s = bmo.est[bmo.rownames .== "s"][1]
	v = bmo.est[bmo.rownames .== "v"][1]
	
	# Update y, s, v (sympatry, subset sympatry, vicariance)
	ysv_func = bmo.type[bmo.rownames .== "ysv"][1] # YOU NEED THE [1] IN HERE
	if ysv_func == "3-j"
		ysv = 3.0-j_wt
		ys = ysv*2/3
		y = ysv*1/3
		s = ysv*1/3
		v = ysv*1/3
	end
	if ysv_func == "2-j" # assumes "s" is 0.0
		ysv = 2.0-j_wt
		ys = ysv*1/2
		y = ysv*1/2
		v = ysv*1/2
	end
	if ysv_func == "1-j"
		ysv = 1.0-j_wt
		ys = ysv*1.0
		y_val = ysv*1.0
	end
	
	bmo.est[bmo.rownames .== "ysv"] .= ysv
	bmo.est[bmo.rownames .== "ys"] .= ys
	bmo.est[bmo.rownames .== "y"] .= y
	bmo.est[bmo.rownames .== "s"] .= s
	bmo.est[bmo.rownames .== "v"] .= v
	
	
	# Update u_e and u_mu based on u?
	type_eq_u_TF = bmo.type .== "u"
	rownames_eq_u_e_TF = bmo.rownames .== "u_e"
	rownames_eq_u_mu_TF = bmo.rownames .== "u_mu"
	u_e_from_u_TF = type_eq_u_TF .+ rownames_eq_u_e_TF
	if sum(u_e_from_u_TF) > 0.0
		u = bmo.est[bmo.rownames .== "u"][1]
		bmo.est[bmo.rownames .== "u_e"] .= u
	end
	u_mu_from_u_TF = type_eq_u_TF .+ rownames_eq_u_mu_TF
	if sum(u_mu_from_u_TF) > 0.0
		u = bmo.est[bmo.rownames .== "u"][1]
		bmo.est[bmo.rownames .== "u_mu"] .= u
	end
	
	# Special relationships between parameters
	# Update mx01's based on mx01?

	# Update deathRate based on birthRate?
	if bmo.type[bmo.rownames .== "deathRate"][1] == "birthRate"
		bmo.est[bmo.rownames .== "deathRate"] .= bmo.est[bmo.rownames .== "birthRate"]
	end

	# Update e based on d?
	if bmo.type[bmo.rownames .== "e"] == "d"
		bmo.est[bmo.rownames .== "e"] = bmo.est[bmo.rownames .== "d"]
	end

	# Update deathRate based on e (single-area extirpation rate)
	# NOTE: This also needs a mu of 0.0 for ranges with more than 1 area
	# 
	# You can achieve this with an area matrix with areas of size 1.0,
	# and a "u" fixed to -100.0. Thus, the e for 1 area is:
	#
	# 1.0 * 1.0^-100.0 = 1.0
	#
	# And for 2 areas, it is
	# 1.0 * 2.0^-100.0 = 7.888609052210118e-31
	if bmo.type[bmo.rownames .== "deathRate"][1] == "e"
		bmo.est[bmo.rownames .== "deathRate"] .= bmo.est[bmo.rownames .== "e"]
	end
	
	# If e also equates the birthRate (forces it to be higher)
	if bmo.type[bmo.rownames .== "e"][1] == "birthRate"
		bmo.est[bmo.rownames .== "e"] .= bmo.est[bmo.rownames .== "birthRate"]
	end


	
	# Update xv based on x?
	# i.e., an x of -1 gives an xv of 1
	if bmo.type[bmo.rownames .== "xv"][1] == "-x"
		#u = bmo.est[bmo.rownames .== "e"]
		#bmo.est[bmo.rownames .== "u_e"] = u
		bmo.est[bmo.rownames .== "xv"] = -1.0 * bmo.est[bmo.rownames .== "x"]
	end
	
	
	return bmo.est
end # END bmo_updater_v1_SLOW


"""
# FAST UPDATER WITH PRE-SPECIFIED bmo_rows (from inputs.setup.bmo_rows)
# bmo_updater_v2 *DOES* update u, etc.

# bmo_updater_v1_SLOW -- full updater, slow due to lots of == 
# bmo_updater_v1! -- very fast, but requires inputs.setup.bmo_rows, and doesn't update "u" etc.
# bmo_updater_v2 -- 2x slower, requires inputs.setup.bmo_rows, but *DOES* update "u" etc.

# Example:
numareas = 2
tr = readTopology("((sp4:0.6248637277,sp5:0.6248637277):6.489662918,(sp6:0.1274213816,sp7:0.1274213816):6.987105264);")
geog_df = DataFrame(tipnames=["sp4","sp5","sp6","sp7"],A=[1,1,0,0],B=[0,0,1,1]);

# in_params = (birthRate=0.222222222, deathRate=0.111111111, d_val=0.0, e_val=0.0, a_val=0.1, j_val=0.0)
# pars <- c(0.222222222, 0.222222222, 0.111111111, 0.05, 0.1, 0.15)
bmo = construct_BioGeoBEARS_model_object();
bmo_rows = get_bmo_rows(bmo)

bmo.est[bmo.rownames.=="d"] .= 0.0;
bmo.est[bmo.rownames.=="e"] .= 0.0;
bmo.est[bmo.rownames.=="j"] .= 0.5;
bmo
bmo.est .= bmo_updater_v2(bmo, bmo_rows)
bmo
"""
function bmo_updater_v1(bmo, bmo_rows)
	# Update cladogenesis parameters
	j_wt = bmo.est[bmo_rows.j][1]
	ysv = bmo.est[bmo_rows.ysv][1]
	ys = bmo.est[bmo_rows.ys][1]
	y = bmo.est[bmo_rows.y][1]
	s = bmo.est[bmo_rows.s][1]
	v = bmo.est[bmo_rows.v][1]
	
	# Update
	ysv_func = bmo.type[bmo_rows.ysv]  # NO [1] is needed here
	if ysv_func == "3-j"
		ysv = 3.0-j_wt
		ys = ysv*2/3
		y = ysv*1/3
		s = ysv*1/3
		v = ysv*1/3
	end
	if ysv_func == "2-j" # assumes "s" is 0.0
		ysv = 2.0-j_wt
		ys = ysv*1/2
		y = ysv*1/2
		v = ysv*1/2
	end
	if ysv_func == "1-j"
		ysv = 1.0-j_wt
		ys = ysv*1.0
		y_val = ysv*1.0
	end
	
	bmo.est[bmo_rows.ysv] = ysv
	bmo.est[bmo_rows.ys] = ys
	bmo.est[bmo_rows.y] = y
	bmo.est[bmo_rows.s] = s
	bmo.est[bmo_rows.v] = v


	# Update deathRate based on birthRate?
	if bmo.type[bmo_rows.deathRate] == "birthRate"
		bmo.est[bmo_rows.deathRate] = bmo.est[bmo_rows.birthRate]
	end

	# Update e based on d?
	if bmo.type[bmo_rows.e] == "d"
		bmo.est[bmo_rows.e] = bmo.est[bmo_rows.d]
	end

	# Update deathRate based on e?
	if bmo.type[bmo_rows.deathRate] == "e"
		bmo.est[bmo_rows.deathRate] = bmo.est[bmo_rows.e]
	end

	# Update e based on birthRate?
	if bmo.type[bmo_rows.e] == "birthRate"
		bmo.est[bmo_rows.e] = bmo.est[bmo_rows.birthRate]
	end

	
	# Update mx01's based on mx01?
	
	return bmo.est
end # END bmo_updater_v1


"""
# FAST UPDATER WITH PRE-SPECIFIED bmo_rows (from inputs.setup.bmo_rows)
# bmo_updater_v2 *DOES* update u, etc.

# bmo_updater_v1_SLOW -- full updater, slow due to lots of == 
# bmo_updater_v1! -- very fast, but requires inputs.setup.bmo_rows, and doesn't update "u" etc.
# bmo_updater_v2 -- 2x slower, requires inputs.setup.bmo_rows, but *DOES* update "u" etc.

# Example:
numareas = 2
tr = readTopology("((sp4:0.6248637277,sp5:0.6248637277):6.489662918,(sp6:0.1274213816,sp7:0.1274213816):6.987105264);")
geog_df = DataFrame(tipnames=["sp4","sp5","sp6","sp7"],A=[1,1,0,0],B=[0,0,1,1]);

# in_params = (birthRate=0.222222222, deathRate=0.111111111, d_val=0.0, e_val=0.0, a_val=0.1, j_val=0.0)
# pars <- c(0.222222222, 0.222222222, 0.111111111, 0.05, 0.1, 0.15)
bmo = construct_BioGeoBEARS_model_object();
bmo_rows = get_bmo_rows(bmo)

bmo.est[bmo.rownames.=="d"] .= 0.0;
bmo.est[bmo.rownames.=="e"] .= 0.0;
bmo.est[bmo.rownames.=="j"] .= 0.5;
bmo
bmo.est .= bmo_updater_v2(bmo, bmo_rows)
bmo
"""
function bmo_updater_v2(bmo, bmo_rows)
	# Update cladogenesis parameters
	j_wt = bmo.est[bmo_rows.j][1]
	ysv = bmo.est[bmo_rows.ysv][1]
	ys = bmo.est[bmo_rows.ys][1]
	y = bmo.est[bmo_rows.y][1]
	s = bmo.est[bmo_rows.s][1]
	v = bmo.est[bmo_rows.v][1]
	
	# Update
	ysv_func = bmo.type[bmo_rows.ysv]  # NO [1] is needed here
	if ysv_func == "3-j"
		ysv = 3.0-j_wt
		ys = ysv*2/3
		y = ysv*1/3
		s = ysv*1/3
		v = ysv*1/3
	end
	if ysv_func == "2-j" # assumes "s" is 0.0
		ysv = 2.0-j_wt
		ys = ysv*1/2
		y = ysv*1/2
		v = ysv*1/2
	end
	if ysv_func == "1-j"
		ysv = 1.0-j_wt
		ys = ysv*1.0
		y_val = ysv*1.0
	end
	
	bmo.est[bmo_rows.ysv] = ysv
	bmo.est[bmo_rows.ys] = ys
	bmo.est[bmo_rows.y] = y
	bmo.est[bmo_rows.s] = s
	bmo.est[bmo_rows.v] = v

	# Does this make it slower?? YES, 2x slowdown!
	
	# Make the single-area extinction rate equal to 
	# mu + e
	#
	# (if mu = 0.0, then extinction only happens when e>0.0)
	# mu: intrinsic extinction rate for a range
	# e: range-loss rate
	#
	# 2023-03-12_doesn't work, because this calculation has to be 
	# done elsewhere: TimeDep.update_mus_time_t!
	#
	# 
	#if ((bmo.type[bmo_rows.deathRate] == "mu+e") || (bmo.type[bmo_rows.deathRate] == "e+mu"))
	#	#u = bmo.est[bmo_rows.u]
	#	#bmo.est[bmo_rows.u_e] = u
	#	bmo.est[bmo_rows.u_e] = bmo.est[bmo_rows.u]
	#end
	
		

	# Update u_e and u_mu based on u?
	if bmo.type[bmo_rows.u_e] == "u"
		#u = bmo.est[bmo_rows.u]
		#bmo.est[bmo_rows.u_e] = u
		bmo.est[bmo_rows.u_e] = bmo.est[bmo_rows.u]
	end
	if bmo.type[bmo_rows.u_mu] == "u"
		#u = bmo.est[bmo_rows.u]
		#bmo.est[bmo_rows.u_mu] = u
		bmo.est[bmo_rows.u_mu] = bmo.est[bmo_rows.u]
	end	
	
	# Update deathRate based on birthRate?
	if bmo.type[bmo_rows.deathRate] == "birthRate"
		bmo.est[bmo_rows.deathRate] = bmo.est[bmo_rows.birthRate]
	end

	# Update e based on d?
	if bmo.type[bmo_rows.e] == "d"
		bmo.est[bmo_rows.e] = bmo.est[bmo_rows.d]
	end


	# Update deathRate based on e?
	if bmo.type[bmo_rows.deathRate] == "e"
		bmo.est[bmo_rows.deathRate] = bmo.est[bmo_rows.e]
	end

	# Update e based on birthRate?
	if bmo.type[bmo_rows.e] == "birthRate"
		bmo.est[bmo_rows.e] = bmo.est[bmo_rows.birthRate]
	end
	
	# Update mx01's based on mx01?


	# Update xv based on x?
	# i.e., an x of -1 gives an xv of 1
	if bmo.type[bmo_rows.xv] == "-x"
		#u = bmo.est[bmo_rows.u]
		#bmo.est[bmo_rows.u_e] = u
		bmo.est[bmo_rows.xv] = -1.0 * bmo.est[bmo_rows.x]
	end

	
	return bmo.est
end # END bmo_updater_v2



"""
# p_Ds_v5 updater #1
# p_Ds_v5 contains the actual rates that go straight into the lnL calculations

pars = [0.03505038, 0.02832370]
parnames = ["d", "e"]

# Assuming parameters are in the order of the bmo list of "free" parameters,
# 1. update the bmo
# 2. then update the Qij and Cijk arrays in p_Ds_v5 (the rate parameters arrays)
inputs.bmo.rownames[inputs.bmo.type.=="free"]

pars = [0.03505038, 0.02832370]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
output = p_Ds_v5_updater_v1_OLD(p_Ds_v5, inputs);  # WORKS 2022-03-10
p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)

pars = [0.01, 0.001]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
output = p_Ds_v5_updater_v1_OLD(p_Ds_v5, inputs);  # WORKS 2022-03-10
p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)

pars = [0.1, 0.2]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
output = p_Ds_v5_updater_v1_OLD(p_Ds_v5, inputs);  # WORKS 2022-03-10
p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)


lnL = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func(pars)
"""
function p_Ds_v5_updater_v1_OLD(p_Ds_v5, inputs; check_if_free_params_in_mat=true)
	
	if check_if_free_params_in_mat == true
		free_param_names = inputs.bmo.rownames[inputs.bmo.type .== "free"]
		# Are all of the free param names %in% the p_Ds_v5
		#TF = in([unique(p_Ds_v5.p_indices.Qarray_event_types); unique(p_Ds_v5.p_indices.Carray_event_types)]).(free_param_names)
		TF1 = in(free_param_names).(["j"])[1]
		TF2 = in(unique(p_Ds_v5.p_indices.Carray_event_types)).(["j"])[1]
		if TF1 != TF2
			txt = paste(["ERROR in p_Ds_v5_updater_v1_OLD(): j's must be in both places. j_in_bmo:", TF1, "; j_in_params:", TF2])
			print("\n")
			print(txt)
			print("\n")
			throw(txt)
		end # END if TF1 != TF2
	end # END if check_if_free_params_in_mat == true
	
	# Extract the parameters
	d = inputs.bmo.est[inputs.bmo.rownames .== "d"][1]
	a = inputs.bmo.est[inputs.bmo.rownames .== "a"][1]
	e = inputs.bmo.est[inputs.bmo.rownames .== "e"][1]
	x = inputs.bmo.est[inputs.bmo.rownames .== "x"][1]
	w = inputs.bmo.est[inputs.bmo.rownames .== "w"][1]
	n = inputs.bmo.est[inputs.bmo.rownames .== "n"][1]
	x2 = 1.0
	x3 = 1.0 
	u = inputs.bmo.est[inputs.bmo.rownames .== "u"][1]
	
	# Update the dmat and elist
	inputs.setup.dmat_base .= d
	inputs.setup.amat_base .= a
	inputs.setup.elist_base .= e
	
	# Create the final dmat (dispersal rate multipliers) by multiplying the various input dmats
	inputs.setup.dmat .= inputs.setup.dmat_base .* inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	inputs.setup.amat .= inputs.setup.amat_base .* inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	inputs.setup.elist .= inputs.setup.elist_base .* inputs.setup.area_of_areas.^u
	
	# jmat does not include the "d" or "j" parameter
	inputs.setup.jmat .= inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	
	# Update the mus
	#p_Ds_v5.params.mu_vals[:] .= inputs.bmo.est[inputs.bmo.rownames .== "deathRate"][1]
	p_Ds_v5.params.mu_vals[:] .= inputs.bmo.est[inputs.setup.bmo_rows.deathRate]
	p_Ds_v5.params.psi_vals[:] .= inputs.bmo.est[inputs.setup.bmo_rows.psiRate]
	
	# Now update the p_Ds_v5 (the rates) for Q matrix
	p_Ds_v5 = update_Qij_vals2!(p_Ds_v5, inputs.setup.areas_list, inputs.setup.states_list, inputs.setup.dmat, inputs.setup.elist, inputs.setup.amat; return_df=false);
	
	#print("\n")
	#print(prtQp(p_Ds_v5))
	
	# Now update the p_Ds_v5 (the rates) for C matrix
	pdf_before = prtCp(p_Ds_v5)
	print("\np_Ds_v5_updater_v1, before:")
	print(round.(pdf_before.wt[1:5], digits=4))
	output = update_Cijk_vals2(p_Ds_v5, inputs.setup.areas_list, inputs.setup.states_list, inputs.bmo, inputs.setup.maxent01, inputs.setup.jmat);
	p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
	p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
	p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
	pdf_after = prtCp(p_Ds_v5)
	print("\np_Ds_v5_updater_v1, after:")
	print(round.(pdf_before.wt[1:5], digits=4))
		
	return output
end # END function p_Ds_v5_updater_v1_OLD()



"""
######################################
# update_maxent01(bmo)
######################################
#
# The "maxent01" object contains a series of tables specifying the relative weights
# of daughters of different range sizes.  This was constructed in BioGeoBEARS to 
# allow a single parameter to control the decisions about daughter range sizes.
# 
# The "maxent" part only refers to the fact that a Maximum Entropy function applied
# to discrete numerical data (e.g. numbers 1-6 on a die) can put flat or skewed 
# weights/probabilities on the various possible outcomes. 
#
# Basically, there is a mx01 parameter for each cladogenesis process: y,s,v, and j.
#
# If the mx01 parameter = 0.0, then the smaller daughter with a rangesize of 1
# has 100% probability. (Thus, DEC has mx01y, mx01s, mx01v = 0.0.)
#
# If the mx01 parameter = 0.5, then the smaller daughter with a rangesize of 1
# has 100% probability. (Thus, DIVALIKE has mx01v = 0.5.)
#
# If the mx01 parameter = 1.0 then the smaller daughter will have the 
# largest possible rangesize, given the ancestor range. (Thus, BAYAREALIKE has mx01y=0.999.)
# 
# update_maxent01(bmo) re-calculates the maxent01 list of tables, given the mx01 parameters
# in bmo (the BioGeoBEARS_model_object)
"""
function update_maxent01(bmo)
	# Generate the tables (weights for smaller daughter rangesize, for each ancestral rangesize)
	maxent01symp = relative_probabilities_of_subsets(total_numareas, bmo.est[bmo.rownames .== "mx01y"][1])
	maxent01sub = relative_probabilities_of_subsets(total_numareas, bmo.est[bmo.rownames .== "mx01s"][1])
	maxent01jump = relative_probabilities_of_subsets(total_numareas, bmo.est[bmo.rownames .== "mx01j"][1])
	maxent01vic = relative_probabilities_of_vicariants(total_numareas, bmo.est[bmo.rownames .== "mx01v"][1])
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
	return maxent01
end # END function update_maxent01(bmo)



#######################################################
# Update the Cijk_vals
#######################################################
function update_Cijk_vals(Carray, areas_list, states_list, maxent01, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); row_weightvals=null )

	if (Rclass(Carray) == "DataFrame")
		if (row_weightvals == null)
			txt = "STOP ERROR in update_Cijk_vals(): Input 'Carray' is a DataFrame, so row_weightvals cannot be null."
			println("\n")
			println(txt)
			println("\n")
			error(txt)
		end
		Carray_ivals = Carray[!,:i]
		Carray_jvals = Carray[!,:j]
		Carray_kvals = Carray[!,:k]
		Carray_pair = Carray[!,:pair]
		Cijk_weights = Carray[!,:wt]
		Cijk_probs = Carray[!,:prob] # abcd
		Cijk_rates = Carray[!,:rate]
		Cijk_vals = Carray[!,:val]
		row_weightvals = row_weightvals
	else
		Carray_event_types = Carray.Carray_event_types;
		Carray_ivals = Carray.Carray_ivals;
		Carray_jvals = Carray.Carray_jvals;
		Carray_kvals = Carray.Carray_kvals;
		Carray_pair = Carray. Carray_pair;
		Cijk_weights = Carray.Cijk_weights;
		Cijk_probs = Carray.Cijk_probs
		Cijk_rates = Carray.Cijk_rates
		Cijk_vals = Carray.Cijk_vals;
		row_weightvals = Carray.row_weightvals
	end # END if (Rclass(Carray) == "DataFrame")



	numstates = length(states_list)
	
	maxent01symp = maxent01.maxent01symp
	maxent01sub = maxent01.maxent01sub
	maxent01vic = maxent01.maxent01vic
	maxent01jump = maxent01.maxent01jump
	
	# Weights
	y_wt = Cparams.y
	s_wt = Cparams.s
	v_wt = Cparams.v
	j_wt = Cparams.j
	
	# Update the "y" events (narrow sympatry / range-copying)
	TF = Carray_event_types .== "y"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			weights[z] = y_wt * maxent01symp[ancsize, lsize] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of y event weights


	# Update the "s" events (subset sympatry)
	TF = Carray_event_types .== "s"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			weights[z] = s_wt * maxent01sub[ancsize, lsize] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	# Update the "v" events (vicariance)
	TF = Carray_event_types .== "v"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of v event weights
	
	# Update the "j" events (jump dispersal)
	TF = Carray_event_types .== "j"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			
			# j events, modified by distance / multipliers (via input "dmat") if needed
			try_jump_dispersal_based_on_dist = true
			normalize_by_number_of_dispersal_events = true
			jweight_for_cell_based_on_distances = 0.0
			if (try_jump_dispersal_based_on_dist == true)
				for anc_area in ancstate
					for left_area in lstate
						 jweight_for_cell_based_on_distances += dmat[anc_area,left_area]
					end
				end
				# Normalize by number of possible jump dispersals
				if (normalize_by_number_of_dispersal_events == true)
					jweight_for_cell_based_on_distances = jweight_for_cell_based_on_distances / (ancsize * lsize)
				end
			else
				# 
				jweight_for_cell_based_on_distances = 1.0
			end # end if (try_jump_dispersal_based_on_dist == true)
	
			# Calculate the final weight of this jump dispersal
			tmp_weightval = j_wt * maxent01jump[ancsize, lsize] * 1.0 * 1.0 * jweight_for_cell_based_on_distances
			weights[z] = tmp_weightval
		end
		Cijk_weights[TF] = weights
	end # End update of j event weights

#	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
#	row_weightvals_df = by(df1, :i, :weight => sum)
#	row_weightvals = row_weightvals_df[!,2]
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
#	if (in(1, row_weightvals_df[!,1]) == false)
#		row_weightvals = repeat([0], 1+length(row_weightvals_df[!,2]))
#		row_weightvals[1] = 1
#		row_weightvals[2:length(row_weightvals)] = row_weightvals_df[!,2]
#	end
#	row_weightvals
	
	# Convert the weights to conditional event probabilities
#	for i in 1:length(states_list)
#		TF = Carray_ivals .== i
#		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
#	end

#	df2 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
#	row_weightvals_df = by(df2, :i, :weight => sum)

	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
	groups = groupby(df1,:i)
	row_weightvals = collect(repeat([0.0], length(groups)))
	for g in 1:length(groups)
		row_weightvals[g] = sum(groups[g].weight)
	end
	row_weightvals
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
	if in(1, unique(df1.i)) == false
		prepend!(row_weightvals, 0)
	end
	
	# Convert the weights to conditional event probabilities
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_probs[TF] = Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
		Cijk_rates[TF] = Cijk_probs[TF] .* birthRate # by default, the birthRate is 1.0; change manually afterwards
	end

	
	# Finally, return updated Carray:
	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Carray_pair=Carray_pair, Cijk_weights=Cijk_weights, Cijk_probs=Cijk_probs, Cijk_rates=Cijk_rates, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)
	
	"""
	# Extract the values
	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Carray_pair = Carray. Carray_pair;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_probs = Carray.Cijk_probs;
	Cijk_rates = Carray.Cijk_rates;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals;
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, wt=Cijk_weights, prob=Cijk_probs, rate=Cijk_rates, val=Cijk_vals)
	row_weightvals
	"""

	return Carray
end # end update_Cijk_vals()




#######################################################
# Update the update_Cijk_vals2_noUpdate
# maxent01 = list of tables
#######################################################
function update_Cijk_vals2_noUpdate(p_Ds_v5, areas_list, states_list, bmo, maxent01, jmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))))

	Carray_ivals = p_Ds_v5.p_indices.Carray_ivals
	Carray_jvals = p_Ds_v5.p_indices.Carray_jvals
	Carray_kvals = p_Ds_v5.p_indices.Carray_kvals
	Carray_pair = p_Ds_v5.p_indices.Carray_pair
	Cijk_weights = p_Ds_v5.params.Cijk_weights
	Cijk_probs = p_Ds_v5.params.Cijk_probs;
	Cijk_rates = p_Ds_v5.params.Cijk_rates;
	Cijk_vals = p_Ds_v5.params.Cijk_vals
	Carray_event_types = p_Ds_v5.p_indices.Carray_event_types
	row_weightvals = p_Ds_v5.params.row_weightvals

	numstates = length(states_list)
	total_numareas = length(areas_list)
	
	# Get the parameters
	# Weights
	y_wt = bmo.est[bmo.rownames .== "y"][1]
	s_wt = bmo.est[bmo.rownames .== "s"][1]
	v_wt = bmo.est[bmo.rownames .== "v"][1]
	j_wt = bmo.est[bmo.rownames .== "j"][1]
	mx01y = bmo.est[bmo.rownames .== "mx01y"][1]
	mx01s = bmo.est[bmo.rownames .== "mx01s"][1]
	mx01v = bmo.est[bmo.rownames .== "mx01v"][1]
	mx01j = bmo.est[bmo.rownames .== "mx01j"][1]

	# Print y, s, v, j
	print("\n")
	txt = paste0(["update_Cijk_vals2: y=", round(y_wt, digits=4), ", s=", round(s_wt, digits=4), ", v=", round(v_wt, digits=4), ", j=", round(j_wt, digits=4)])
	print(txt)
	
	# Rates
	birthRate = bmo.est[bmo.rownames .== "birthRate"][1]
	
	# Extract the maxent01 tables, which control
	# the relative weight of the rangesizes of
	# the smaller daughter
	# (Update these tables OUTSIDE of this
	#  function!)
	maxent01symp = maxent01.maxent01symp
	maxent01sub = maxent01.maxent01sub
	maxent01vic = maxent01.maxent01vic
	maxent01jump = maxent01.maxent01jump
	
	
	# Update the "y" events (narrow sympatry / range-copying)
	TF = Carray_event_types .== "y"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = y_wt * maxent01symp[ancsize, smaller_range_size] * pair_vals[z]
			#weights[z] = y_wt * maxent01symp[ancsize, smaller_range_size] * 1.0 * 1.0
			Carray_pair
		end
		Cijk_weights[TF] = weights
	end # End update of y event weights


	# Update the "s" events (subset sympatry)
	TF = Carray_event_types .== "s"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = s_wt * maxent01sub[ancsize, smaller_range_size] * pair_vals[z]
			#weights[z] = s_wt * maxent01sub[ancsize, smaller_range_size] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	# Update the "v" events (vicariance)
	TF = Carray_event_types .== "v"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = v_wt * maxent01vic[ancsize,smaller_range_size] * pair_vals[z]
			#weights[z] = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of v event weights
	
	# Update the "j" events (jump dispersal)
	TF = Carray_event_types .== "j"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			
			# j events, modified by distance / multipliers (via input "jmat") if needed
			try_jump_dispersal_based_on_dist = true
			normalize_by_number_of_dispersal_events = true
			jweight_for_cell_based_on_distances = 0.0
			if (try_jump_dispersal_based_on_dist == true)
				for anc_area in ancstate
					for left_area in lstate
						 jweight_for_cell_based_on_distances += jmat[anc_area,left_area]
					end
				end
				# Normalize by number of possible jump dispersals
				if (normalize_by_number_of_dispersal_events == true)
					jweight_for_cell_based_on_distances = jweight_for_cell_based_on_distances / (ancsize * lsize)
				end
			else
				# 
				jweight_for_cell_based_on_distances = 1.0
			end # end if (try_jump_dispersal_based_on_dist == true)
	
			# Calculate the final weight of this jump dispersal
			smaller_range_size = min(lsize, rsize)
			tmp_weightval = j_wt * maxent01jump[ancsize, smaller_range_size] * jweight_for_cell_based_on_distances * pair_vals[z]
			#tmp_weightval = j_wt * maxent01jump[ancsize, smaller_range_size] * 1.0 * 1.0 * jweight_for_cell_based_on_distances
			
			weights[z] = tmp_weightval
		end
		Cijk_weights[TF] = weights
	end # End update of j event weights

	#df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
#	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
#	row_weightvals_df = by(df1, :i, :weight => sum)
#	row_weightvals = row_weightvals_df[!,2]
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
#	if (in(1, row_weightvals_df[!,1]) == false)
#		row_weightvals = repeat([0], 1+length(row_weightvals_df[!,2]))
#		row_weightvals[1] = 1
#		row_weightvals[2:length(row_weightvals)] = row_weightvals_df[!,2]
#	end
#	row_weightvals
	
	# Convert the weights to conditional event probabilities
#	for i in 1:length(states_list)
#		TF = Carray_ivals .== i
#		Cijk_vals[TF] = Cijk_weights[TF] .* birthRate ./ row_weightvals[i]
#	end

	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
	groups = groupby(df1,:i)
	row_weightvals = collect(repeat([0.0], length(groups)))
	for g in 1:length(groups)
		row_weightvals[g] = sum(groups[g].weight)
	end
	row_weightvals
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
	if in(1, unique(df1.i)) == false
		prepend!(row_weightvals, 0)
	end
	
	# Convert the weights to conditional event probabilities
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_probs[TF] = Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
		Cijk_rates[TF] = Cijk_probs[TF] .* birthRate # by default, the birthRate is 1.0; change manually afterwards
	end
	
	# Update
	p_Ds_v5.params.Cijk_weights[:] .= Cijk_weights
	p_Ds_v5.params.Cijk_vals[:] .= Cijk_vals
	p_Ds_v5.params.Cijk_probs[:] .= Cijk_probs
	p_Ds_v5.params.Cijk_rates[:] .= Cijk_rates
	p_Ds_v5.params.row_weightvals[:] .= row_weightvals


	"""
	# Extract the values
	Carray_event_types = p_Ds_v5.p_indices.Carray_event_types;
	Carray_ivals = p_Ds_v5.p_indices.Carray_ivals
	Carray_jvals = p_Ds_v5.p_indices.Carray_jvals
	Carray_kvals = p_Ds_v5.p_indices.Carray_kvals
	Cijk_weights = p_Ds_v5.params.Cijk_weights;
	Cijk_probs = Carray.Cijk_probs;
	Cijk_rates = Carray.Cijk_rates;
	Cijk_vals = p_Ds_v5.params.Cijk_vals;
	row_weightvals = p_Ds_v5.params.row_weightvals
		DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, wt=Cijk_weights, prob=Cijk_probs, rate=Cijk_rates, val=Cijk_vals)
	"""
	
	output = (Cijk_weights=Cijk_weights, Cijk_probs=Cijk_probs, Cijk_rates=Cijk_rates, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)
	return output
	#return p_Ds_v5
end # end update_Cijk_vals2_noUpdate()









#######################################################
# Update the Cijk_vals2!
# maxent01 = list of tables
#######################################################
function update_Cijk_vals2!(p_Ds_v5, areas_list, states_list, bmo, maxent01, jmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); printlevel=0)
	"""
	# prelim defaults:
	root_age_mult=1.5; max_range_size=NaN; include_null_range=false; bmo=bmo
	manual_states_list=NaN; area_names=LETTERS(1:numareas)

	check_if_free_params_in_mat=true;
	printlevel=0


	# defaults:
	areas_list = inputs.setup.areas_list
	states_list = inputs.setup.states_list
	bmo = inputs.bmo
	maxent01 = inputs.setup.maxent01
	jmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list)));
	jmat = inputs.setup.jmat
	printlevel=0
	"""
	
	
	Carray_event_types = p_Ds_v5.p_indices.Carray_event_types
	Carray_ivals = p_Ds_v5.p_indices.Carray_ivals
	Carray_jvals = p_Ds_v5.p_indices.Carray_jvals
	Carray_kvals = p_Ds_v5.p_indices.Carray_kvals
	Carray_pair = p_Ds_v5.p_indices.Carray_pair
	Cijk_weights = p_Ds_v5.params.Cijk_weights
	Cijk_probs = p_Ds_v5.params.Cijk_probs
	Cijk_rates = p_Ds_v5.params.Cijk_rates
	Cijk_vals = p_Ds_v5.params.Cijk_vals
	row_weightvals = p_Ds_v5.params.row_weightvals

	numstates = length(states_list)
	total_numareas = length(areas_list)
	
	# Make sure bmo is up-to-date before running this
	


	# Get the parameters
	# Weights
	y_wt = bmo.est[bmo.rownames .== "y"][1]
	s_wt = bmo.est[bmo.rownames .== "s"][1]
	v_wt = bmo.est[bmo.rownames .== "v"][1]
	j_wt = bmo.est[bmo.rownames .== "j"][1]
	mx01y = bmo.est[bmo.rownames .== "mx01y"][1]
	mx01s = bmo.est[bmo.rownames .== "mx01s"][1]
	mx01v = bmo.est[bmo.rownames .== "mx01v"][1]
	mx01j = bmo.est[bmo.rownames .== "mx01j"][1]

	# Print y, s, v, j
	if printlevel >= 2
		print("\n")
		txt = paste0(["update_Cijk_vals2: y=", round(y_wt, digits=3), ", s=", round(s_wt, digits=3), ", v=", round(v_wt, digits=3), ", j=", round(j_wt, digits=3)])
		print(txt)
	end
	
	# Rates
	birthRate = bmo.est[bmo.rownames .== "birthRate"][1]
	
	# Extract the maxent01 tables, which control
	# the relative weight of the rangesizes of
	# the smaller daughter
	# (Update these tables OUTSIDE of this
	#  function!)
	maxent01symp = maxent01.maxent01symp
	maxent01sub = maxent01.maxent01sub
	maxent01vic = maxent01.maxent01vic
	maxent01jump = maxent01.maxent01jump
	
	
	# Update the "y" events (narrow sympatry / range-copying)
	TF = Carray_event_types .== "y"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			if ancsize == 0 # null range ancestor (weird, but allowed in some simulations)
				weights[z] = y_wt * 1.0 * 1.0 * pair_vals[z]
			else
				weights[z] = y_wt * maxent01symp[ancsize, smaller_range_size] * 1.0 * 1.0 * pair_vals[z]
			end
		end
		Cijk_weights[TF] = weights
	end # End update of y event weights


	# Update the "s" events (subset sympatry)
	TF = Carray_event_types .== "s"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = s_wt * maxent01sub[ancsize, smaller_range_size] * 1.0 * 1.0 * pair_vals[z]
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	# Update the "v" events (vicariance)
	TF = Carray_event_types .== "v"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0 * pair_vals[z]
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights
	
	# Update the "j" events (jump dispersal)
	TF = Carray_event_types .== "j"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			
			# The best way to do distances in this weights framework, is to set the minimum distance equal to 1.0
			# That way any distances > 1.0 have dispersal penalties (e.g. if x=-1.0)
			# 
			# This way, if jmat_t changes, it can be used as long as jmat=1.0 or is otherwise fixed.
			#
			# (Because we have a modifier on the original jmat, which is user-set.)
			#
			# j events, modified by distance / multipliers (via input "jmat") if needed
			try_jump_dispersal_based_on_dist = true				 # irrelevant when the jmats are 1.0
			normalize_by_number_of_dispersal_events = true # irrelevant, because number of events 
																										 # is the same for each ancestral range
			jweight_for_cell_based_on_distances = 0.0
			if (try_jump_dispersal_based_on_dist == true)
				for anc_area in ancstate
					for left_area in lstate
						 jweight_for_cell_based_on_distances += jmat[anc_area,left_area]
					end
				end
				# Normalize by number of possible jump dispersals
				if (normalize_by_number_of_dispersal_events == true)
					jweight_for_cell_based_on_distances = jweight_for_cell_based_on_distances / (ancsize * lsize)
				end
			else
				# 
				jweight_for_cell_based_on_distances = 1.0
			end # end if (try_jump_dispersal_based_on_dist == true)
	
			# Calculate the final weight of this jump dispersal
			smaller_range_size = min(lsize, rsize)
			tmp_weightval = j_wt * maxent01jump[ancsize, smaller_range_size] * 1.0 * 1.0 * jweight_for_cell_based_on_distances * pair_vals[z]
			weights[z] = tmp_weightval
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	#df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
	
	# Add up the weights of the rows for each ancestral state i
	# And, divide by the sum of the weights
#	for i in 1:numstates
#		iTF = Carray_ivals .== i
#		row_weightvals[i] = sum(Cijk_weights[iTF])
#		# Multiply by birthRate here!!
#		Cijk_vals[iTF] .= Cijk_weights[iTF] .* birthRate ./ row_weightvals[i]
#	end # END for i in 1:numstates
#	Cijk_vals = p_Ds_v5.params.Cijk_vals
	
	# Convert the weights to conditional event probabilities
#	for i in 1:length(states_list)
#		TF = Carray_ivals .== i
#		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
#	end

#	df2 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
#	row_weightvals_df = by(df2, :i, :weight => sum)
	
	# Finally, return updated Carray:
#	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)

	# 2022: weird groups thing (replaced, 2023-02-06)

#	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
#	
#	groups = groupby(df1,:i)
#	row_weightvals = collect(repeat([0.0], length(groups)))
#	for g in 1:length(groups)
#		row_weightvals[g] = sum(groups[g].weight)
#	end
#	row_weightvals
#	
#	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
#	if in(1, unique(df1.i)) == false
#		prepend!(row_weightvals, 0)
#	end
	
	# Convert the weights to conditional event probabilities
#	for i in 1:length(states_list)
#		TF = Carray_ivals .== i
#		Cijk_probs[TF] = Cijk_weights[TF] ./ row_weightvals[i]
#		Cijk_rates[TF] = Cijk_vals[TF] = Cijk_probs[TF] .* birthRate # by default, the birthRate is 1.0; change manually afterwards
#	end




	# 2023-02-06 Simpler, more robust strategy than the weird "groups" thing!
	# Sum the weights for each ancestral i, divide by the sum of the weights
	for i in 1:length(row_weightvals)
		#TF = Carray_ivals .== i
		TF = p_Ds_v5.p_TFs.Ci_eq_i[i]  # This is crucial!
		row_weightvals[i] = sum(Cijk_weights[TF])
		Cijk_probs[TF] = Cijk_weights[TF] ./ row_weightvals[i]
		Cijk_rates[TF] = Cijk_vals[TF] = Cijk_probs[TF] .* birthRate # by default, the birthRate is 1.0; 
	end

	p_Ds_v5.params.Cijk_weights[:] .= Cijk_weights
	p_Ds_v5.params.Cijk_probs[:] .= Cijk_probs
	p_Ds_v5.params.Cijk_rates[:] .= Cijk_rates
	p_Ds_v5.params.Cijk_vals[:] .= Cijk_vals
	p_Ds_v5.params.row_weightvals[:] .= row_weightvals

	# Update the Cijk_rates_sub_i (where anc==i)
	for i in 1:length(states_list)
		p_Ds_v5.p_TFs.Cijk_rates_sub_i[i] .= Cijk_rates[p_Ds_v5.p_TFs.Ci_eq_i[i]]
		# p_Ds_v5.p_TFs.Cjik_rates_sub_j[i] .= Cjik_rates[p_Ds_v5.p_TFs.Ci_eq_j[i]]
		#p_Ds_v5.p_TFs.Cjik_rates_sub_j[i] .= Cijk_rates[p_Ds_v5.p_indices.Carray_jvals .== i] # uppass
		p_Ds_v5.p_TFs.Cjik_rates_sub_j[i] .= Cijk_rates[p_Ds_v5.p_TFs.Cj_eq_j[i]] # for uppass
	end
	

	"""
	# Extract the values
	Carray_event_types = p_Ds_v5.p_indices.Carray_event_types;
	Carray_ivals = p_Ds_v5.p_indices.Carray_ivals
	Carray_jvals = p_Ds_v5.p_indices.Carray_jvals
	Carray_kvals = p_Ds_v5.p_indices.Carray_kvals
	Cijk_weights = p_Ds_v5.params.Cijk_weights;
	Cijk_probs = p_Ds_v5.params.Carray.Cijk_probs;
	Cijk_rates = p_Ds_v5.params.Carray.Cijk_rates;
	Cijk_vals = p_Ds_v5.params.Cijk_vals;
	Cijk_rates_sub_i = p_Ds_v5.p_TFs.Cijk_rates_sub_i;
	Cjik_rates_sub_i = p_Ds_v5.p_TFs.Cjik_rates_sub_j;
	row_weightvals = p_Ds_v5.params.row_weightvals
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, wt=Cijk_weights, prob=Cijk_probs, rate=Cijk_rates, val=Cijk_vals)
	"""

	return p_Ds_v5
end # end update_Cijk_vals2!()


"""
# Update the BioGeoBEARS_inputs_object after changing e.g. "j"
"""
function inputs_updater_v1!(inputs)
	#global inputs
	# bmo_updater_v1! doesn't update u etc.
	bmo_updater_v1!(inputs.bmo)
	return inputs
end # END function inputs_updater_v1!(inputs)

"""
# Update the BioGeoBEARS_inputs_object after changing e.g. "j"
"""
function inputs_updater_v2!(inputs)
	#global inputs
	# bmo_updater_v1! doesn't update u etc.
	#bmo_updater_v1!(inputs.bmo)
	# bmo_updater_v2 DOES update u etc.
	inputs.bmo.est = bmo_updater_v2(inputs.bmo, inputs.setup.bmo_rows);
	return inputs
end # END function inputs_updater_v1!(inputs)



"""
# Update the BioGeoBEARS_model_object after changing e.g. "j"
#
# bmo_updater_v1! *doesn't* update u, etc. (speed advantage; use v2 for full version)
#
# bmo_updater_v1_SLOW -- full updater, slow due to lots of == 
# bmo_updater_v1! -- very fast, but requires inputs.setup.bmo_rows, and doesn't update "u" etc.
# bmo_updater_v2 -- 2x slower, requires inputs.setup.bmo_rows, but *DOES* update "u" etc.


"""
function bmo_updater_v1!(bmo)
	#global bmo
	j_wt = bmo.est[bmo.rownames .== "j"][1]
	ysv = bmo.est[bmo.rownames .== "ysv"][1]
	ys = bmo.est[bmo.rownames .== "ys"][1]
	y = bmo.est[bmo.rownames .== "y"][1]
	s = bmo.est[bmo.rownames .== "s"][1]
	v = bmo.est[bmo.rownames .== "v"][1]
	
	# Update
	ysv_func = bmo.type[bmo.rownames .== "ysv"][1]
	if ysv_func == "3-j"
		ysv = 3.0-j_wt
		ys = ysv*2/3
		y = ysv*1/3
		s = ysv*1/3
		v = ysv*1/3
	end
	if ysv_func == "2-j" # assumes "s" is 0.0
		ysv = 2.0-j_wt
		ys = ysv*1/2
		y = ysv*1/2
		v = ysv*1/2
	end
	if ysv_func == "1-j"
		ysv = 1.0-j_wt
		ys = ysv*1.0
		y_val = ysv*1.0
	end
	
	
	bmo.est[bmo.rownames .== "ysv"] .= ysv
	bmo.est[bmo.rownames .== "ys"] .= ys
	bmo.est[bmo.rownames .== "y"] .= y
	bmo.est[bmo.rownames .== "s"] .= s
	bmo.est[bmo.rownames .== "v"] .= v
	
	return bmo
end # END



"""
# p_Ds_v5 updater #1 - Updates p_Ds_v5 **FROM AN INPUTS OBJECT**
# 
# p_Ds_v5 contains the actual rates that go straight into the lnL calculations

pars = [0.03505038, 0.02832370]
parnames = ["d", "e"]

# Assuming parameters are in the order of the bmo list of "free" parameters,
# 1. update the bmo
# 2. then update the Qij and Cijk arrays in p_Ds_v5 (the rate parameters arrays)
inputs.bmo.rownames[inputs.bmo.type.=="free"]

pars = [0.03505038, 0.02832370]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs)
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)

pars = [0.01, 0.001]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs)
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)

pars = [0.1, 0.2]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs)
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)


lnL = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func(pars)
"""
function p_Ds_v5_updater_v1!(p_Ds_v5, inputs; check_if_free_params_in_mat=true, printlevel=0)
	defaults="""
	root_age_mult=1.5; max_range_size=NaN; include_null_range=false; bmo=bmo
	manual_states_list=NaN; area_names=LETTERS(1:numareas)

	check_if_free_params_in_mat=true;
	printlevel=0

	"""
	
	if check_if_free_params_in_mat == true
		free_param_names = inputs.bmo.rownames[inputs.bmo.type .== "free"]
		# Are all of the free param names %in% the p_Ds_v5
		#TF = in([unique(p_Ds_v5.p_indices.Qarray_event_types); unique(p_Ds_v5.p_indices.Carray_event_types)]).(free_param_names)
		TF1 = in(free_param_names).(["j"])[1]
		TF2 = in(unique(p_Ds_v5.p_indices.Carray_event_types)).(["j"])[1]
		if TF1 != TF2
			txt = paste0(["ERROR in p_Ds_v5_updater_v1!(): j's must be in both places.\nj_in_bmo:", TF1, "\nj_in_params:", TF2, ".\n(Check e.g. prtCp(p_Ds_v5), column 'events'.)"])
			print("\n")
			print(txt)
			print("\n")
			throw(txt)
		end # END if TF1 != TF2
	end # END if check_if_free_params_in_mat == true
	
	# Extract the parameters
	d = inputs.bmo.est[inputs.bmo.rownames .== "d"][1]
	a = inputs.bmo.est[inputs.bmo.rownames .== "a"][1]
	e = inputs.bmo.est[inputs.bmo.rownames .== "e"][1]
	x = inputs.bmo.est[inputs.bmo.rownames .== "x"][1]
	w = inputs.bmo.est[inputs.bmo.rownames .== "w"][1]
	n = inputs.bmo.est[inputs.bmo.rownames .== "n"][1]
	x2 = inputs.bmo.est[inputs.bmo.rownames .== "x2"][1]
	x3 = inputs.bmo.est[inputs.bmo.rownames .== "x3"][1]
	u = inputs.bmo.est[inputs.bmo.rownames .== "u"][1]
	
	# Update the dmat and elist
	inputs.setup.dmat_base .= d
	inputs.setup.amat_base .= a
	inputs.setup.elist_base .= e
	
	# Create the final dmat (dispersal rate multipliers) by multiplying the various input dmats
	inputs.setup.dmat .= inputs.setup.dmat_base .* inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	inputs.setup.amat .= inputs.setup.amat_base .* inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	inputs.setup.elist .= inputs.setup.elist_base .* inputs.setup.area_of_areas.^u
	
	# jmat does not include the "d" or "j" parameter
	inputs.setup.jmat .= inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	
	
	# Now update the p_Ds_v5 (the rates) for Q matrix
	p_Ds_v5 = update_Qij_vals2!(p_Ds_v5, inputs.setup.areas_list, inputs.setup.states_list, inputs.setup.dmat, inputs.setup.elist, inputs.setup.amat; return_df=false);

	# Update the mus & psi
	#p_Ds_v5.params.mu_vals[:] .= inputs.bmo.est[inputs.setup.bmo_rows.deathRate]
	
	# multi_area_ranges_have_zero_mu: set all multi-area ranges to 0.0 mu
  # 
	# mu_func: the extinction rate for a range could be:
	# mu:   bmo.deathRate, perhaps multiplied by area^u_mu
	# mu+e: bmo.deathRate, plus the rate of e^u_e for single-area ranges
	# If you set mu+e, and set mu=0.0, then the single-area deathRate is e
	area_of_ranges = repeat([0.0], inputs.setup.numstates)
	single_area_range_TF = repeat([false], inputs.setup.numstates)
	for i in 1:inputs.setup.numstates
		area_of_ranges[i] = sum(inputs.setup.area_of_areas[inputs.setup.states_list[i]])
		single_area_range_TF[i] = length(inputs.setup.states_list[i]) == 1
	end
	area_of_ranges
	single_area_range_TF
	
	u_mu = inputs.bmo.est[inputs.setup.bmo_rows.u_mu]
	u_e = inputs.bmo.est[inputs.setup.bmo_rows.u_e]
	if (inputs.setup.mu_func[1] == "mu")
		p_Ds_v5.params.mu_vals[:] .= inputs.bmo.est[inputs.setup.bmo_rows.deathRate] .* (area_of_ranges .^ u_mu)
	elseif ((p.setup.mu_func[1] == "mu+e") || (p.setup.mu_func[1] == "e+mu"))
		p_Ds_v5.params.mu_vals[:] .= single_area_range_TF .* (p_Ds_v5.params.mu_vals .* (area_of_ranges .^ u_mu)) .+ (inputs.bmo.est[inputs.setup.bmo_rows.e] .* (area_of_ranges .^ u_e))	
	end # END if (inputs.setup.mu_func == "mu")
	
	# Psi: constant for now
	p_Ds_v5.params.psi_vals[:] .= inputs.bmo.est[inputs.setup.bmo_rows.psiRate]

  # Correct "Inf" max_extinction_rates
 	p_Ds_v5.params.mu_vals[p_Ds_v5.params.mu_vals .> inputs.setup.max_extinction_rate] .= inputs.setup.max_extinction_rate
 	
 	# Alternative way to set extinction to 0.0 for multi-area ranges
 	if (inputs.setup.multi_area_ranges_have_zero_mu[1] == true)
	 	p_Ds_v5.params.mu_vals[length.(inputs.setup.states_list) .> 1] .= 0.0
	end


	
	# (The mu_vals_t might change with time; if so, see update_mus_time_t
	# p.params.mu_vals_t[i] = p.params.mu_vals[i] * get_area_of_range(t, p.states_as_areas_lists[i], p.setup.area_of_areas)^p.bmo.est[p.setup.bmo_rows.u_mu]

	#print("\n")
	#print(prtQp(p_Ds_v5))
	
	# Now update the p_Ds_v5 (the rates) for C matrix
	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v5)
		print("\np_Ds_v5_updater_v1, before:")
		print(round.(pdf_before.wt[1:5], digits=3))
	end
	update_Cijk_vals2!(p_Ds_v5, inputs.setup.areas_list, inputs.setup.states_list, inputs.bmo, inputs.setup.maxent01, inputs.setup.jmat);
	if printlevel >= 2
		pdf_after = prtCp(p_Ds_v5)
		print("\np_Ds_v5_updater_v1, after:")
		print(round.(pdf_before.wt[1:5], digits=3))
	end
	
	return p_Ds_v5	
end # END function p_Ds_v5_updater_v1!()













end # END module Optimizers
