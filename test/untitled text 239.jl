function func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="lnL", printlevel=1)

	# Get the results object (with old p_Ds_v5)
	# and treetable (trdf)
	res = inputs.res
	trdf = inputs.trdf
	
	# Set the result if parameters out of bounds
	nan_lnL = -Inf
	
	# Update the parameters
	inputs.bmo.est[inputs.bmo.type.=="free"] .= pars;
	
	# Get parameter names and lower/upper bounds
	parnames = inputs.bmo.rownames[inputs.bmo.type.=="free"];
	lower = inputs.bmo.min[bmo.type .== "free"]
	upper = inputs.bmo.max[bmo.type .== "free"]
	#bmo = inputs.bmo
	
	# Update
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est before bmo_updater_v1(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est[inputs.bmo.rownames .== "j"] .= 1.6
	#inputs.bmo.est[:] = bmo_updater_v1(inputs.bmo);
	bmo_updater_v1(inputs.bmo);
	if printlevel >= 2
		print("\nfunc_to_optimize, inputs.bmo.est after bmo_updater_v1(): ")
		print(round.(inputs.bmo.est[[1,2,9,12,13,14]], digits=4))
	end
	#inputs.bmo.est .= bmo.est

	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v5)
		print("\nfunc_to_optimize: clado weights before p_Ds_v5_updater_v1:")
		print(round.(pdf_before.wt[1:5], digits=4))
	end
#	output = p_Ds_v5_updater_v1(p_Ds_v5, inputs);  # WORKS 2022-03-10
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
		prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5)
		# This solution is an interpolator
		sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
		p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
		
		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
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
		txt = paste0([txt, "Julia_sum_lq=", round(Julia_sum_lq; digits=4), ", rootstates_lnL=", round(rootstates_lnL; digits=4), ",	Julia_total_lnLs1=", Julia_total_lnLs1, ", bgb_lnL=", round(bgb_lnL, digits=4)])
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
		inputs.p_Ds_v5.params.Qij_vals[:] .= p_Ds_v5.params.Qij_vals
		inputs.p_Ds_v5.params.Cijk_weights[:] .= p_Ds_v5.params.Cijk_weights
		inputs.p_Ds_v5.params.Cijk_vals[:] .= p_Ds_v5.params.Cijk_vals
		inputs.p_Ds_v5.params.row_weightvals[:] .= p_Ds_v5.params.row_weightvals
		
		# results (res)
		inputs.res.regime[:] .= res.regime
		inputs.res.node_state[:] .= res.node_state
		inputs.res.node_Lparent_state[:] .= res.node_Lparent_state
		inputs.res.node_Rparent_state[:] .= res.node_Rparent_state
		inputs.res.root_nodeIndex[:] .= res.root_nodeIndex
		inputs.res.numNodes[:] .= res.numNodes
		inputs.res.uppass_edgematrix[:] .= res.uppass_edgematrix
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
