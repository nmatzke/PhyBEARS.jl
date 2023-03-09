module Flow

__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/


print("PhyBEARS: loading Flow dependencies...")
using Distributed			# for Distributed.nprocs()
using LinearAlgebra  # for mul! (matrix multiplication)
#using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA
using DifferentialEquations
using Base.Threads  # <-- this is good for @spawn, Distributed.@spawn is BAD, it produces Futures that have to be scheduled etc]
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloBits
using Sundials # for CVODE_BDF(linear_solver=:GMRES)
#using Plots						# for plot
using DataFrames          # for DataFrame()

using PhyloBits.TreeTable  # for get_nonrootnodes_trdf
using PhyloBits.TrUtils # for flat2() (similar to unlist), and 
                          # convert_is_js_to_single_index(),
                          # get_a_most_common_value()
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
print("...done.\n")

export iterative_downpass_nonparallel_ClaSSE_v6_fromFlow!, iterative_downpass_parallel_ClaSSE_v6_fromFlow!, parameterized_ClaSSE_As_v6, parameterized_ClaSSE_As_v6_parallel, parameterized_ClaSSE_As_v6_sub_i, parameterized_ClaSSE_As_v5, check_linearDynamics_of_As, calc_Gs_SSE_condnums!, calc_Gs_SSE, calc_Gs_SSE_v7simd, calc_Gs_SSE_parallel, calc_Gs_SSE_sub_i, calc_Gs_SSE!, branchOp_ClaSSE_Gs_v5, iterative_downpass_nonparallel_FlowClaSSE_v5!, iterative_downpass_wGflowArray_nonparallel!, run_Gs, parameterized_ClaSSE_As_v7!, parameterized_ClaSSE_As_v7_simd!, sum_Qij_vals_inbounds_simd_A!, sum_Cijk_vals_inbounds_simd_A!




"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_nonparallel_ClaSSE_v6_fromFlow!(res; trdf, p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true)
	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()
	
	# Get basic tree info
	numTips = sum(trdf.nodeType .== "tip")
	numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")
	
	# Re-set the node states for new downpass
	TF = trdf[:,:nodeType] .== "tip"
	res.node_state[TF] .= "ready_for_branchOp"
	res.node_state[TF .== false] .= "not_ready"
	res.node_Lparent_state[TF] .= "NA"
	res.node_Lparent_state[TF .== false] .= "not_ready"
	res.node_Rparent_state[TF] .= "NA"
	res.node_Rparent_state[TF .== false] .= "not_ready"
		
	# Setup
	current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		#print("\niteration_number: ")
		#print(iteration_number)
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "done"
				return()
			end
	
			# Retrieve the inputs for the calculation down the branch
			
			# Use the RAW likelihoods (don't normalize)
			#u0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
			#u0 = u0 ./ (sum(u0))
			# You can use the normalized likelihoods, see correction at bottom
			u0 = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
			
			# Use the NORMALIZED (rescaled to sum to 1) likelihoods
			# Doesn't work -- claims an interpolation error for going beyond range
			# branchOp on current_nodeIndex=4ERROR: LoadError: Solution interpolation 
			# cannot extrapolate past the final timepoint. Either solve on a longer 
			# timespan or use the local extrapolation from the integrator interface.
			#u0 = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
			
			brlen = trdf[current_nodeIndex, :brlen]
			age_branchtop = trdf[current_nodeIndex, :node_age]
			age_branchbot = age_branchtop + brlen
			tspan = (age_branchtop, age_branchbot)
			#p_Ds_v5 = inputs.p_Ds_v5

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			#print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
#			if (parallel_TF == true)
#				push!(tasks, Base.Threads.@spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			tmp_results = branchOp_ClaSSE_Ds_v6(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			#tmp_results = branchOp(current_nodeIndex, res, num_iterations)
			push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			#print("\ni: ")
			#print(i)
			if (tasks_fetched_TF[i] == false)
				#if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = tasks[i]
					nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
# 					end
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
# 					print("\n\n12345\n\n")
# 					print("res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex]:\n")
# 					print(res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex])
# 					print("\n\nnodeData_at_bottom:\n")
# 					print(nodeData_at_bottom)
# 					print("\n\n12345\n\n")

					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
					# Can cause an underflow error 2022-03-10: "log will only return a complex result if called with a complex argument. Try log(Complex(x))."
					if sum_nodeData_at_bottom < 0.0
						sum_nodeData_at_bottom = 1e-10000
					end
					res.lq_at_branchBot[spawned_nodeIndex] = log(sum_nodeData_at_bottom) 
					res.like_at_branchBot[spawned_nodeIndex] = sum_nodeData_at_bottom

					# Get the ancestor nodeIndex
					uppass_edgematrix = res.uppass_edgematrix
					TF = uppass_edgematrix[:,2] .== spawned_nodeIndex
					parent_nodeIndex = uppass_edgematrix[TF,1][1]

					# Get the left daughter nodeIndex (1st in the uppass_edgematrix)
					edge_rows_TF = uppass_edgematrix[:,1] .== parent_nodeIndex
	
					# Error trap
					if (sum(edge_rows_TF) > 3)
						txt = "STOP error in: iterative_downpass_Gflow_nonparallel_v1() - node has more than 2 edges"
					end
					
					# Standard bifurcating node
					if (sum(edge_rows_TF) == 2)
						left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
						right_nodeIndex = uppass_edgematrix[edge_rows_TF,2][2]

						# Update the state of the parent_node's daughters
						if (spawned_nodeIndex == left_nodeIndex)
							res.node_Lparent_state[parent_nodeIndex] = "ready"
						end
						if (spawned_nodeIndex == right_nodeIndex)
							res.node_Rparent_state[parent_nodeIndex] = "ready"
						end
					end 

					# Singleton node (assumes direct-ancestor nodes are always "left"
					if (sum(edge_rows_TF) == 1)
						left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
						res.node_Lparent_state[parent_nodeIndex] = "ready"
					end
	
					# Update the state of the current node
					res.node_state[spawned_nodeIndex] = "done"
				#end
			end # END if (tasks_fetched_TF[i] == false)
		end # END for i in 1:num_tasks

		# Update which nodes are SINGLETONS and are complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = trdf.nodeType .== "direct"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"

		# Update nodes when the singletons above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, Base.Threads.@spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end

	
		# Update which nodes have had both parents complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = res.node_Rparent_state .== "ready"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"
	
		# Update nodes when the branches above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, Base.Threads.@spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			res = nodeOp_ClaSSE_v6!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass_nonparallel(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
			error(txt)
		end
		
		# Test for concluding the while loop
		are_we_done && break
	end # END while(are_we_done == false)
	
	# This breaks it for some reason:
	# ERROR: setfield! immutable struct of type Res cannot be changed
	#global res.number_of_whileLoop_iterations = iteration_number

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	Julia_sum_lq_old = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
	nonroot_nodes = get_nonrootnodes_trdf(trdf)
	sum_likes_internal_branch_tops = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
	Julia_sum_lq = Julia_sum_lq_old + sum_likes_internal_branch_tops

	# Add the root probabilities

	# Assuming diversitree options:
	# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
	# i.e., the root state probs are just the root_Ds/sum(root_Ds)
	d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
	root_stateprobs = d_root_orig/sum(d_root_orig)
	rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
	Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL

	# Without normalization
	#Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL
	# *With* normalization
	#Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL + log(sum(res.sumLikes_at_node_at_branchTop[1:(length(res.sumLikes_at_node_at_branchTop)-1)]))
	# Redundant with: sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
	
	# Consider the pure-birth log-likelihoods
	# The Yule-process ML birthrate is just (# internal nodes - 1)/total_tree_length
	ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
	yuleBirthRate = (numInternal-1) / ttl_tree_length
	yuleDeathRate = 0.0					# Yule process has 0 extinction
	bd = bd_liks_trdf(trdf, yuleBirthRate, yuleDeathRate)
	bd_lnL_noTopo = bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times
	
	# Convert to BioGeoBEARS lnL under Yule process assumption
	# Check if the first state/geographic range is null
	#if res.inputs.setup.states_list[1] == []
	#	include_null_range = true
	#end
	numstates = length(res.normlikes_at_each_nodeIndex_branchTop[1])
	equal_root_prob2 = log(1/(numstates-include_null_range)) 
	bgb_root_lnL = log(sum(d_root_orig)) + 1.0
	
	# res5t match
	res5t = Julia_sum_lq + equal_root_prob2 + bgb_root_lnL - (1-log(1/yuleBirthRate))
	# ...matches these in R: compare_BGB_diversitree_DEC_v1.R
	# bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/(birthRate))
	# bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 
	# bgb1 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 + bgb_root_lnL
	# (bgb1 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

	# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
	bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
	bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2
	bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2
	bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)
	bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo
	
	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_nonparallel_ClaSSE_v6_fromFlow!








"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
#
# Launch multi-core Julia with e.g.: 
# julia --procs auto  # (uses as many cores as are available)
# julia -p auto       # (same)
# julia --procs 10    # (for 10 cores)
#
# Once you are inside Julia:
# length(Sys.cpu_info())
# nprocs
"""
function iterative_downpass_parallel_ClaSSE_v6_fromFlow!(res; trdf, p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true)

	# Check the number of threads; this function will hang unless
	# there are multiple threads
	numthreads = Base.Threads.nthreads()
	num_processes = Distributed.nprocs()

	if ((numthreads > 1) == false) || ((num_processes > 1) == false)
		txt = paste0(["\n\nSTOP ERROR in iterative_downpass_parallel_ClaSSE_v6_fromFlow!(): this function runs ODEs down branches in *parallel*. Therefore it needs multiple threads and multiple processes/workers to function. Without that, the function will hang. You have:\n\n     Base.Threads.nthreads() = ", numthreads, "\n     Distributed.nprocs()    = ", num_processes, "\n\nYou need to re-start Julia with multiple threads, using this command (run at the command line, not inside julia):\n\n'julia -t auto -p auto' (or replace 'auto' with 4 for 4 threads & processes)\n\nYou can add processes, but not threads, from inside julia with 'tmptask = @async Distributed.addprocs(7)'.\n\n"]);
		print(txt)
	end # END if (numthreads > 1) == false


	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()
	
	# Get basic tree info
	numTips = sum(trdf.nodeType .== "tip")
	numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")
	
	# Re-set the node states for new downpass
	TF = trdf[:,:nodeType] .== "tip"
	res.node_state[TF] .= "ready_for_branchOp"
	res.node_state[TF .== false] .= "not_ready"
	res.node_Lparent_state[TF] .= "NA"
	res.node_Lparent_state[TF .== false] .= "not_ready"
	res.node_Rparent_state[TF] .= "NA"
	res.node_Rparent_state[TF .== false] .= "not_ready"
		
	# Setup
	current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		#print("\niteration_number: ")
		#print(iteration_number)
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "done"
				return()
			end
	
			# Retrieve the inputs for the calculation down the branch
			
			# Use the RAW likelihoods (don't normalize)
			#u0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
			#u0 = u0 ./ (sum(u0))
			# You can use the normalized likelihoods, see correction at bottom
			u0 = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
			
			# Use the NORMALIZED (rescaled to sum to 1) likelihoods
			# Doesn't work -- claims an interpolation error for going beyond range
			# branchOp on current_nodeIndex=4ERROR: LoadError: Solution interpolation 
			# cannot extrapolate past the final timepoint. Either solve on a longer 
			# timespan or use the local extrapolation from the integrator interface.
			#u0 = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
			
			brlen = trdf[current_nodeIndex, :brlen]
			age_branchtop = trdf[current_nodeIndex, :node_age]
			age_branchbot = age_branchtop + brlen
			tspan = (age_branchtop, age_branchbot)
			#p_Ds_v5 = inputs.p_Ds_v5

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			#print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
#			if (parallel_TF == true)
#				push!(tasks, Base.Threads.@spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			#tmp_results = branchOp_ClaSSE_Ds_v6(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			push!(tasks, Base.Threads.@spawn branchOp_ClaSSE_Ds_v6(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options))
			# MethodError: no method matching istaskstarted(::Future)
			#push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			#print("\ni: ")
			#print(i)
			if (tasks_fetched_TF[i] == false)
				if (istaskstarted(tasks[i]) == true) && (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					#(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = tasks[i]
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
					nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
# 					end
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
# 					print("\n\n12345\n\n")
# 					print("res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex]:\n")
# 					print(res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex])
# 					print("\n\nnodeData_at_bottom:\n")
# 					print(nodeData_at_bottom)
# 					print("\n\n12345\n\n")

					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
					# Can cause an underflow error 2022-03-10: "log will only return a complex result if called with a complex argument. Try log(Complex(x))."
					if sum_nodeData_at_bottom < 0.0
						sum_nodeData_at_bottom = 1e-10000
					end
					res.lq_at_branchBot[spawned_nodeIndex] = log(sum_nodeData_at_bottom) 
					res.like_at_branchBot[spawned_nodeIndex] = sum_nodeData_at_bottom

					# Get the ancestor nodeIndex
					uppass_edgematrix = res.uppass_edgematrix
					TF = uppass_edgematrix[:,2] .== spawned_nodeIndex
					parent_nodeIndex = uppass_edgematrix[TF,1][1]

					# Get the left daughter nodeIndex (1st in the uppass_edgematrix)
					edge_rows_TF = uppass_edgematrix[:,1] .== parent_nodeIndex
	
					# Error trap
					if (sum(edge_rows_TF) > 3)
						txt = "STOP error in: iterative_downpass_Gflow_nonparallel_v1() - node has more than 2 edges"
					end
					
					# Standard bifurcating node
					if (sum(edge_rows_TF) == 2)
						left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
						right_nodeIndex = uppass_edgematrix[edge_rows_TF,2][2]

						# Update the state of the parent_node's daughters
						if (spawned_nodeIndex == left_nodeIndex)
							res.node_Lparent_state[parent_nodeIndex] = "ready"
						end
						if (spawned_nodeIndex == right_nodeIndex)
							res.node_Rparent_state[parent_nodeIndex] = "ready"
						end
					end 

					# Singleton node (assumes direct-ancestor nodes are always "left"
					if (sum(edge_rows_TF) == 1)
						left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
						res.node_Lparent_state[parent_nodeIndex] = "ready"
					end
	
					# Update the state of the current node
					res.node_state[spawned_nodeIndex] = "done"
				end # END if (istaskdone(tasks[i]) == true)
			end # END if (tasks_fetched_TF[i] == false)
		end # END for i in 1:num_tasks

		# Update which nodes are SINGLETONS and are complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = trdf.nodeType .== "direct"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"

		# Update nodes when the singletons above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, Base.Threads.@spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end

	
		# Update which nodes have had both parents complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = res.node_Rparent_state .== "ready"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"
	
		# Update nodes when the branches above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, Base.Threads.@spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			res = nodeOp_ClaSSE_v6!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass_nonparallel(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
			error(txt)
		end
		
		# Test for concluding the while loop
		are_we_done && break
	end # END while(are_we_done == false)
	
	# This breaks it for some reason:
	# ERROR: setfield! immutable struct of type Res cannot be changed
	#global res.number_of_whileLoop_iterations = iteration_number

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	Julia_sum_lq_old = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
	nonroot_nodes = get_nonrootnodes_trdf(trdf)
	sum_likes_internal_branch_tops = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
	Julia_sum_lq = Julia_sum_lq_old + sum_likes_internal_branch_tops

	# Add the root probabilities

	# Assuming diversitree options:
	# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
	# i.e., the root state probs are just the root_Ds/sum(root_Ds)
	d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
	root_stateprobs = d_root_orig/sum(d_root_orig)
	rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
	Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL

	# Without normalization
	#Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL
	# *With* normalization
	#Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL + log(sum(res.sumLikes_at_node_at_branchTop[1:(length(res.sumLikes_at_node_at_branchTop)-1)]))
	# Redundant with: sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
	
	# Consider the pure-birth log-likelihoods
	# The Yule-process ML birthrate is just (# internal nodes - 1)/total_tree_length
	ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
	yuleBirthRate = (numInternal-1) / ttl_tree_length
	yuleDeathRate = 0.0					# Yule process has 0 extinction
	bd = bd_liks_trdf(trdf, yuleBirthRate, yuleDeathRate)
	bd_lnL_noTopo = bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times
	
	# Convert to BioGeoBEARS lnL under Yule process assumption
	# Check if the first state/geographic range is null
	#if res.inputs.setup.states_list[1] == []
	#	include_null_range = true
	#end
	numstates = length(res.normlikes_at_each_nodeIndex_branchTop[1])
	equal_root_prob2 = log(1/(numstates-include_null_range)) 
	bgb_root_lnL = log(sum(d_root_orig)) + 1.0
	
	# res5t match
	res5t = Julia_sum_lq + equal_root_prob2 + bgb_root_lnL - (1-log(1/yuleBirthRate))
	# ...matches these in R: compare_BGB_diversitree_DEC_v1.R
	# bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/(birthRate))
	# bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 
	# bgb1 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 + bgb_root_lnL
	# (bgb1 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))

	# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
	bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
	bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2
	bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2
	bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)
	bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo
	
	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_parallel_ClaSSE_v6_fromFlow!





# Construct interpolation function for calculating linear dynamics A, 
# at any timepoint t

#function get_LinearDynamics_A(inputs)

# 	modelD.getLinearDynamics(ages[a], dynamics); // get linear dynamics at this age
# 		// calculate largest singular value (sigma1) of the dynamics at this age
# 		// then kappa_rate <= 2*sigma1, since sigma1(exp(t*A))/sigma2(exp(t*A)) <= exp(t*2*sigma1(A)) [So and Thompson (2000) Singular values of Matrix Exponentials. Theorem 2.1]
#
# NOTE: In Julia, 
# When p=2, the operator norm is the spectral norm, equal to the largest singular value of A.
# 
# 
#
# This version excludes Xi (and Xj), just like castor's get_LinearDynamics_A
# calculation of A
parameterized_ClaSSE_As_v6 = (t, A, p; max_condition_number=1e8, print_warnings=false) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	# Terms for calculating floats
	# p.terms
	# Normally, the claSSE equation boils down to this:
	# du[i] = -(terms[1] + terms[2] + mu[i] + psi[i])*u[i] + terms[3] + terms[4]
	# ...but for A, we just need terms[1] and terms[2]
	
	two = 1.0
	# Iterate through the ancestral states
  @inbounds for i in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		#A[i,i] = A[i,i]  + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i]) # *u[i]  
		A[i,i] = A[i,i] - sum(Qij_vals[Qi_eq_i]) - sum(Cijk_vals[Ci_eq_i]) - mu[i] #+ 2*sum(Cijk_vals[Ci_sub_i])*uE[i]
		
		# case 2: anagenetic change (non-diagonal)
		@inbounds for m in 1:length(Qi_sub_i)
			A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_eq_i][m] #* u[Qj_sub_i[m]])
		end
		
		# case 34: change + eventual extinction (non-diagonal)
		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			# excluding the u[], i.e. the Ds, i.e. the Xs, just as is done in 
			# 2*speciation_rates[r]*current_E[r]
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i] * ((u[Ck_sub_i] * uE[Cj_sub_i]) + (u[Cj_sub_i] * uE[Ck_sub_i]))
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * ((uE[Ck_sub_i[m]] * uE[Cj_sub_i][m]) + (uE[Cj_sub_i[m]] * uE[Ck_sub_i[m]]))
			
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * (uE[Cj_sub_i[m]] + uE[Ck_sub_i[m]])
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * uE[i]#(uE[Cj_sub_i[m]] + uE[Ck_sub_i[m]])


			# These are several equal ways, & close to the diversitree truth
			# 1.
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * 0.5*((uE[Cj_sub_i[m]]) + (uE[Ck_sub_i[m]]))
			# 2.
			#rate_sp_then_ex1 = Cijk_vals[Ci_sub_i[m]] * 0.5*(Cijk_vals[Cj_sub_i[m]] * uE[Ck_sub_i[m]])
			#rate_sp_then_ex2 = Cijk_vals[Ci_sub_i[m]] * 0.5*(Cijk_vals[Ck_sub_i[m]] * uE[Cj_sub_i[m]])
			# 3. 
			#rate_sp_then_ex1 = Cijk_vals[Ci_sub_i[m]] * 0.5*((uE[Cj_sub_i[m]]) + (uE[Ck_sub_i[m]]))
			#rate_sp_then_ex2 = Cijk_vals[Ci_sub_i[m]] * 0.5*((uE[Ck_sub_i[m]]) + (uE[Cj_sub_i[m]]))
			#rate_sp_then_ex1 = Cijk_vals[Ci_eq_i][m] * uE[Ck_sub_i[m]] # rate of going to j, times k going extinct
			#rate_sp_then_ex2 = Cijk_vals[Ci_eq_i][m] * uE[Cj_sub_i[m]] # rate of going to k, times j going extinct
			rate_sp_then_ex1 = Cijk_vals[Ci_eq_i][m] * uE[Ck_sub_i[m]] # rate of going to j, times k going extinct
			rate_sp_then_ex2 = Cijk_vals[Ci_eq_i][m] * uE[Cj_sub_i[m]] # rate of going to k, times j going extinct
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex1
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex2

			#A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex
			#A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex

			
			# WORKS for DEC-SSE, 2020-06-25
			# julia> mul!(Xc_from_flow, Gflow_to_01l(tc), X0)
			# 3-element Array{Float64,1}:
			#  0.0010385611333310657
			#  0.5288136426527802   
			#  0.03387280014939652  
			# 
			# julia> Xc_from_flow2 = Gflow_to_01g(tc) * X0
			# 3-element Array{Float64,1}:
			#  0.0010385159443680523
			#  0.5287977443578      
			#  0.033874388008434016 
			# 
			# julia> ground_truth_Ds_interpolator(tc)
			# 3-element Array{Float64,1}:
			#  0.0010385159929410246
			#  0.528797722695622    
			#  0.03387438254182548  			
		end
		
# 		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
# 			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
# 			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
# 				 (u[Ck_sub_i].*uE[Cj_sub_i] 
# 			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))

  end # End @inbounds for i in 1:n
	
	# Error check; large growth rate in the condition number
	# (indicated by the norm of the "A" matrix giving the linear dynamics)
	# suggests that the G matrix is approaching singularity and numerical
	# errors will accumulate.
	# Check the maximum condition number of A; if it is exceeded, either raise the max cond number,
	# or divide the tree into smaller chunks
	# (after Louca & Pennell)
	if (print_warnings == true)
		sigma1_of_A = opnorm(A,2)  # the 1-norm should be adequate here (fastest calc.)
		if (2*sigma1_of_A > log(max_condition_number))
			warning_txt = join(["WARNING in parameterized_ClaSSE_As_v6 at t=", string(t), ": 2*opnorm(A,1)>log(max_condition_number) (", string(round(2*sigma1_of_A; digits=2)), " > ", string(round(log(max_condition_number); digits=2)), ")\n"], "")
			display(warning_txt)
		end # END if (2*sigma1_of_A > log(max_condition_number))
	end # END if (print_warnings == true)
 	return(A)
end # END parameterized_ClaSSE_As_v6




parameterized_ClaSSE_As_v6_sub_i = (t, A, p; max_condition_number=1e8, print_warnings=false) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  #Qij_vals = p.params.Qij_vals
  #Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
	# Pre-calculated solution of the Es
	#sol_Es = p.sol_Es_v5
	uE = p.uE
	#uE = sol_Es(t)
	uE = p.sol_Es_v5(t)
	
	
	two = 1.0
	# Iterate through the ancestral states
  @inbounds for i in 1:n
	#Base.Threads.@threads for i in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qij_singleNum_sub_i = p.p_TFs.Qij_singleNum_sub_i[i]
		Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i[i]
		
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i[i]
		Cij_singleNum_sub_i = p.p_TFs.Cij_singleNum_sub_i[i]
		Cik_singleNum_sub_i = p.p_TFs.Cik_singleNum_sub_i[i]

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		#A[i,i] = A[i,i]  + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i]) # *u[i]  
		#A[i,i] = A[i,i] - sum(Qij_vals[Qi_eq_i]) - sum(Cijk_vals[Ci_eq_i]) - mu[i] #+ 2*sum(Cijk_vals[Ci_sub_i])*uE[i]
		A[i,i] += -sum(Qij_vals_sub_i) - sum(Cijk_rates_sub_i) - mu[i] #+ 2*sum(Cijk_vals[Ci_sub_i])*uE[i]
		
		# case 2: anagenetic change (non-diagonal)
		#@inbounds for m in 1:length(Qi_sub_i)
			#A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_eq_i][m] #* u[Qj_sub_i[m]])
		#	A[Qi_sub_i[m],Qj_sub_i[m]] += Qij_vals_sub_i[m]
		#end
		
		# 200 seconds
#		A[Qi_sub_i,Qj_sub_i] .+= Qij_vals_sub_i
#		A[Ci_sub_i,Cj_sub_i] .+= Cijk_rates_sub_i .* uE[Ck_sub_i]
#		A[Ci_sub_i,Ck_sub_i] .+= Cijk_rates_sub_i .* uE[Cj_sub_i]

		A[Qij_singleNum_sub_i] .+= Qij_vals_sub_i
		A[Cij_singleNum_sub_i] .+= Cijk_rates_sub_i .* uE[Ck_sub_i]
		A[Cik_singleNum_sub_i] .+= Cijk_rates_sub_i .* uE[Cj_sub_i]
		

		
		# case 34: change + eventual extinction (non-diagonal)
#		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			# excluding the u[], i.e. the Ds, i.e. the Xs, just as is done in 
			# 2*speciation_rates[r]*current_E[r]
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i] * ((u[Ck_sub_i] * uE[Cj_sub_i]) + (u[Cj_sub_i] * uE[Ck_sub_i]))
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * ((uE[Ck_sub_i[m]] * uE[Cj_sub_i][m]) + (uE[Cj_sub_i[m]] * uE[Ck_sub_i[m]]))
			
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * (uE[Cj_sub_i[m]] + uE[Ck_sub_i[m]])
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * uE[i]#(uE[Cj_sub_i[m]] + uE[Ck_sub_i[m]])


			# These are several equal ways, & close to the diversitree truth
			# 1.
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * 0.5*((uE[Cj_sub_i[m]]) + (uE[Ck_sub_i[m]]))
			# 2.
			#rate_sp_then_ex1 = Cijk_vals[Ci_sub_i[m]] * 0.5*(Cijk_vals[Cj_sub_i[m]] * uE[Ck_sub_i[m]])
			#rate_sp_then_ex2 = Cijk_vals[Ci_sub_i[m]] * 0.5*(Cijk_vals[Ck_sub_i[m]] * uE[Cj_sub_i[m]])
			# 3. 
			#rate_sp_then_ex1 = Cijk_vals[Ci_sub_i[m]] * 0.5*((uE[Cj_sub_i[m]]) + (uE[Ck_sub_i[m]]))
			#rate_sp_then_ex2 = Cijk_vals[Ci_sub_i[m]] * 0.5*((uE[Ck_sub_i[m]]) + (uE[Cj_sub_i[m]]))
			#rate_sp_then_ex1 = Cijk_vals[Ci_eq_i][m] * uE[Ck_sub_i[m]] # rate of going to j, times k going extinct
			#rate_sp_then_ex2 = Cijk_vals[Ci_eq_i][m] * uE[Cj_sub_i[m]] # rate of going to k, times j going extinct
			#rate_sp_then_ex1 = Cijk_vals[Ci_eq_i][m] * uE[Ck_sub_i[m]] # rate of going to j, times k going extinct
			#rate_sp_then_ex2 = Cijk_vals[Ci_eq_i][m] * uE[Cj_sub_i[m]] # rate of going to k, times j going extinct
			#A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex1
			#A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex2
			#A[Ci_sub_i[m],Cj_sub_i[m]] += Cijk_rates_sub_i[m] * uE[Ck_sub_i[m]]
			#A[Ci_sub_i[m],Ck_sub_i[m]] += Cijk_rates_sub_i[m] * uE[Cj_sub_i[m]]

			#A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex
			#A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex

			
			# WORKS for DEC-SSE, 2020-06-25
			# julia> mul!(Xc_from_flow, Gflow_to_01l(tc), X0)
			# 3-element Array{Float64,1}:
			#  0.0010385611333310657
			#  0.5288136426527802   
			#  0.03387280014939652  
			# 
			# julia> Xc_from_flow2 = Gflow_to_01g(tc) * X0
			# 3-element Array{Float64,1}:
			#  0.0010385159443680523
			#  0.5287977443578      
			#  0.033874388008434016 
			# 
			# julia> ground_truth_Ds_interpolator(tc)
			# 3-element Array{Float64,1}:
			#  0.0010385159929410246
			#  0.528797722695622    
			#  0.03387438254182548  			
#		end
		
# 		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
# 			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
# 			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
# 				 (u[Ck_sub_i].*uE[Cj_sub_i] 
# 			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))

  end # End @inbounds for i in 1:n
	
	# Error check; large growth rate in the condition number
	# (indicated by the norm of the "A" matrix giving the linear dynamics)
	# suggests that the G matrix is approaching singularity and numerical
	# errors will accumulate.
	# Check the maximum condition number of A; if it is exceeded, either raise the max cond number,
	# or divide the tree into smaller chunks
	# (after Louca & Pennell)
	if (print_warnings == true)
		sigma1_of_A = opnorm(A,2)  # the 1-norm should be adequate here (fastest calc.)
		if (2*sigma1_of_A > log(max_condition_number))
			warning_txt = join(["WARNING in parameterized_ClaSSE_As_v6 at t=", string(t), ": 2*opnorm(A,1)>log(max_condition_number) (", string(round(2*sigma1_of_A; digits=2)), " > ", string(round(log(max_condition_number); digits=2)), ")\n"], "")
			display(warning_txt)
		end # END if (2*sigma1_of_A > log(max_condition_number))
	end # END if (print_warnings == true)
 	return(A)
end # END parameterized_ClaSSE_As_v6_sub_i



  



parameterized_ClaSSE_As_v6_parallel = (t, A, p; max_condition_number=1e8, print_warnings=false) -> begin
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  #Qij_vals = p.params.Qij_vals
  #Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
	# Pre-calculated solution of the Es
	#sol_Es = p.sol_Es_v5
	uE = p.uE
	#uE = sol_Es(t)
	uE = p.sol_Es_v5(t)
	
	
	two = 1.0
	# Iterate through the ancestral states
	Base.Threads.@threads for i in 1:n
	#Base.Threads.@threads for i in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qij_singleNum_sub_i = p.p_TFs.Qij_singleNum_sub_i[i]
		Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i[i]
		
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i[i]
		Cij_singleNum_sub_i = p.p_TFs.Cij_singleNum_sub_i[i]
		Cik_singleNum_sub_i = p.p_TFs.Cik_singleNum_sub_i[i]

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		#A[i,i] = A[i,i]  + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i]) # *u[i]  
		#A[i,i] = A[i,i] - sum(Qij_vals[Qi_eq_i]) - sum(Cijk_vals[Ci_eq_i]) - mu[i] #+ 2*sum(Cijk_vals[Ci_sub_i])*uE[i]
		A[i,i] += -sum(Qij_vals_sub_i) - sum(Cijk_rates_sub_i) - mu[i] #+ 2*sum(Cijk_vals[Ci_sub_i])*uE[i]
		
		# case 2: anagenetic change (non-diagonal)
		#@inbounds for m in 1:length(Qi_sub_i)
			#A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_eq_i][m] #* u[Qj_sub_i[m]])
		#	A[Qi_sub_i[m],Qj_sub_i[m]] += Qij_vals_sub_i[m]
		#end
		
		# 200 seconds
#		A[Qi_sub_i,Qj_sub_i] .+= Qij_vals_sub_i
#		A[Ci_sub_i,Cj_sub_i] .+= Cijk_rates_sub_i .* uE[Ck_sub_i]
#		A[Ci_sub_i,Ck_sub_i] .+= Cijk_rates_sub_i .* uE[Cj_sub_i]

		A[Qij_singleNum_sub_i] .+= Qij_vals_sub_i
		A[Cij_singleNum_sub_i] .+= Cijk_rates_sub_i .* uE[Ck_sub_i]
		A[Cik_singleNum_sub_i] .+= Cijk_rates_sub_i .* uE[Cj_sub_i]
		

		
		# case 34: change + eventual extinction (non-diagonal)
#		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			# excluding the u[], i.e. the Ds, i.e. the Xs, just as is done in 
			# 2*speciation_rates[r]*current_E[r]
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i] * ((u[Ck_sub_i] * uE[Cj_sub_i]) + (u[Cj_sub_i] * uE[Ck_sub_i]))
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * ((uE[Ck_sub_i[m]] * uE[Cj_sub_i][m]) + (uE[Cj_sub_i[m]] * uE[Ck_sub_i[m]]))
			
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * (uE[Cj_sub_i[m]] + uE[Ck_sub_i[m]])
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * uE[i]#(uE[Cj_sub_i[m]] + uE[Ck_sub_i[m]])


			# These are several equal ways, & close to the diversitree truth
			# 1.
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * 0.5*((uE[Cj_sub_i[m]]) + (uE[Ck_sub_i[m]]))
			# 2.
			#rate_sp_then_ex1 = Cijk_vals[Ci_sub_i[m]] * 0.5*(Cijk_vals[Cj_sub_i[m]] * uE[Ck_sub_i[m]])
			#rate_sp_then_ex2 = Cijk_vals[Ci_sub_i[m]] * 0.5*(Cijk_vals[Ck_sub_i[m]] * uE[Cj_sub_i[m]])
			# 3. 
			#rate_sp_then_ex1 = Cijk_vals[Ci_sub_i[m]] * 0.5*((uE[Cj_sub_i[m]]) + (uE[Ck_sub_i[m]]))
			#rate_sp_then_ex2 = Cijk_vals[Ci_sub_i[m]] * 0.5*((uE[Ck_sub_i[m]]) + (uE[Cj_sub_i[m]]))
			#rate_sp_then_ex1 = Cijk_vals[Ci_eq_i][m] * uE[Ck_sub_i[m]] # rate of going to j, times k going extinct
			#rate_sp_then_ex2 = Cijk_vals[Ci_eq_i][m] * uE[Cj_sub_i[m]] # rate of going to k, times j going extinct
			#rate_sp_then_ex1 = Cijk_vals[Ci_eq_i][m] * uE[Ck_sub_i[m]] # rate of going to j, times k going extinct
			#rate_sp_then_ex2 = Cijk_vals[Ci_eq_i][m] * uE[Cj_sub_i[m]] # rate of going to k, times j going extinct
			#A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex1
			#A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex2
			#A[Ci_sub_i[m],Cj_sub_i[m]] += Cijk_rates_sub_i[m] * uE[Ck_sub_i[m]]
			#A[Ci_sub_i[m],Ck_sub_i[m]] += Cijk_rates_sub_i[m] * uE[Cj_sub_i[m]]

			#A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex
			#A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex

			
			# WORKS for DEC-SSE, 2020-06-25
			# julia> mul!(Xc_from_flow, Gflow_to_01l(tc), X0)
			# 3-element Array{Float64,1}:
			#  0.0010385611333310657
			#  0.5288136426527802   
			#  0.03387280014939652  
			# 
			# julia> Xc_from_flow2 = Gflow_to_01g(tc) * X0
			# 3-element Array{Float64,1}:
			#  0.0010385159443680523
			#  0.5287977443578      
			#  0.033874388008434016 
			# 
			# julia> ground_truth_Ds_interpolator(tc)
			# 3-element Array{Float64,1}:
			#  0.0010385159929410246
			#  0.528797722695622    
			#  0.03387438254182548  			
#		end
		
# 		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
# 			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
# 			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
# 				 (u[Ck_sub_i].*uE[Cj_sub_i] 
# 			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))

  end # End @inbounds for i in 1:n
	
	# Error check; large growth rate in the condition number
	# (indicated by the norm of the "A" matrix giving the linear dynamics)
	# suggests that the G matrix is approaching singularity and numerical
	# errors will accumulate.
	# Check the maximum condition number of A; if it is exceeded, either raise the max cond number,
	# or divide the tree into smaller chunks
	# (after Louca & Pennell)
	if (print_warnings == true)
		sigma1_of_A = opnorm(A,2)  # the 1-norm should be adequate here (fastest calc.)
		if (2*sigma1_of_A > log(max_condition_number))
			warning_txt = join(["WARNING in parameterized_ClaSSE_As_v6 at t=", string(t), ": 2*opnorm(A,1)>log(max_condition_number) (", string(round(2*sigma1_of_A; digits=2)), " > ", string(round(log(max_condition_number); digits=2)), ")\n"], "")
			display(warning_txt)
		end # END if (2*sigma1_of_A > log(max_condition_number))
	end # END if (print_warnings == true)
 	return(A)
end # END parameterized_ClaSSE_As_v6






# Construct interpolation function for calculating linear dynamics A, 
# at any timepoint t

# From Castor:
#	// provide matrix encoding linear rates of change of the current state X, i.e. return A(t), where:
#	//   dX/dt = A(t)*X(t) (if inverse==false)
#	// or:
#	//   dX/dt = X(t)*A(t) (if inverse==true)
#	// note that, in principle, A may also depend on the current state X, i.e. A=A(t,X(t))
#	// The returned A must be in row-major format
#	void getLinearDynamics(double age, std::vector<double> &A) const{
#		const MuSSEstateE current_E = E(age);
#		// The mapping A is essentially the transition_matrix, plus some additional terms on the diagonal
#		A = transition_rates;
#		for(long r=0; r<Nstates; ++r){
#			A[r*Nstates+r] += - (speciation_rates[r]+extinction_rates[r]+sampling_rates[r]) + #2*speciation_rates[r]*current_E[r]; // add more terms to diagonal
#		}
#		if(inverse) A *= -1;
#	}

#function get_LinearDynamics_A(inputs)
# 	modelD.getLinearDynamics(ages[a], dynamics); // get linear dynamics at this age
# 		// calculate largest singular value (sigma1) of the dynamics at this age
# 		// then kappa_rate <= 2*sigma1, since sigma1(exp(t*A))/sigma2(exp(t*A)) <= exp(t*2*sigma1(A)) [So and Thompson (2000) Singular values of Matrix Exponentials. Theorem 2.1]
#
# NOTE: In Julia, 
# When p=2, the operator norm is the spectral norm, equal to the largest singular value of A.
# 
# This version excludes Xi (and Xj), just like castor's get_LinearDynamics_A
# calculation of A
parameterized_ClaSSE_As_v5 = (t, A, p; max_condition_number=1e8, print_warnings=false) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	two = 1.0
	# Iterate through the ancestral states
  @inbounds for i in 1:n
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		#A[i,i] = A[i,i]  + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i]) # *u[i]  
		A[i,i] = A[i,i] - sum(Qij_vals[Qi_sub_i]) - sum(Cijk_vals[Ci_sub_i]) - mu[i] #+ 2*sum(Cijk_vals[Ci_sub_i])*uE[i]
		
		# case 2: anagenetic change (non-diagonal)
		@inbounds for m in 1:length(Qi_sub_i)
			A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_sub_i[m]] #* u[Qj_sub_i[m]])
		end
		
		# case 34: change + eventual extinction (non-diagonal)
		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			# excluding the u[], i.e. the Ds, i.e. the Xs, just as is done in 
			# 2*speciation_rates[r]*current_E[r]
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i] * ((u[Ck_sub_i] * uE[Cj_sub_i]) + (u[Cj_sub_i] * uE[Ck_sub_i]))
			rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * uE[i]#(uE[Cj_sub_i[m]] + uE[Ck_sub_i[m]])
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex
			
			# WORKS for DEC-SSE, 2020-06-25
			# julia> mul!(Xc_from_flow, Gflow_to_01l(tc), X0)
			# 3-element Array{Float64,1}:
			#  0.0010385611333310657
			#  0.5288136426527802   
			#  0.03387280014939652  
			# 
			# julia> Xc_from_flow2 = Gflow_to_01g(tc) * X0
			# 3-element Array{Float64,1}:
			#  0.0010385159443680523
			#  0.5287977443578      
			#  0.033874388008434016 
			# 
			# julia> ground_truth_Ds_interpolator(tc)
			# 3-element Array{Float64,1}:
			#  0.0010385159929410246
			#  0.528797722695622    
			#  0.03387438254182548  			
		end
		
# 		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
# 			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
# 			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
# 				 (u[Ck_sub_i].*uE[Cj_sub_i] 
# 			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))

  end # End @inbounds for i in 1:n
	
	# Error check; large growth rate in the condition number
	# (indicated by the norm of the "A" matrix giving the linear dynamics)
	# suggests that the G matrix is approaching singularity and numerical
	# errors will accumulate.
	# Check the maximum condition number of A; if it is exceeded, either raise the max cond number,
	# or divide the tree into smaller chunks
	# (after Louca & Pennell)
	if (print_warnings == true)
		sigma1_of_A = opnorm(A,2)  # the 1-norm should be adequate here (fastest calc.)
		if (2*sigma1_of_A > log(max_condition_number))
			warning_txt = join(["WARNING in parameterized_ClaSSE_As_v5 at t=", string(t), ": 2*opnorm(A,1)>log(max_condition_number) (", string(round(sigma1_of_A; digits=2)), " > ", string(round(log(max_condition_number); digits=2)), ")\n"], "")
			display(warning_txt)
		end # END if (2*sigma1_of_A > log(max_condition_number))
	end # END if (print_warnings == true)
 	return(A)
end # END parameterized_ClaSSE_As_v5


# Check the linear dynamics (matrix A) for timespans where the kappa rate
# exceeds log(max_condition_number). max_condition_number is usually between
# 1e4 (slower) and 1e8 (faster)
# 
# "exponential growth rate of the condition number ("kappa_rate") based on the 
#  largest singular value of the dynamics A" 
#
# Using  A[:,:] is required, so that "A" doesn't mutate
# 
function check_linearDynamics_of_As(tvals, p_Ds_v5; max_condition_number=1e8)

	# build an A matrix to hold the linear dynamics
	n = p_Ds_v5.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))
	
	# build arrays to hold the output for each t
	upper_bound_kappa_rates_A = collect(repeat([0.0], length(tvals)))
	condbigTF = collect(repeat([false], length(tvals)))
	
	for i in 1:length(tvals)
		t = tvals[i]
		# A(t) is just the instantaneous rates matrix at time t;
		# It doesn't update.
		A_at_t = Flow.parameterized_ClaSSE_As_v6(t, A[:,:], p_Ds_v5)
		sigma1_of_A = opnorm(A_at_t,1)  # the 1-norm should be adequate here (fastest calc.)
		upper_bound_kappa_rates_A[i] = sqrt(n)*2*sigma1_of_A
		if (upper_bound_kappa_rates_A[i] > log(max_condition_number))
			condbigTF[i] = true
		end
	end
	
	
	# Create dataframe
	# condbigTF = Is the growth rathe of condition number too big, i.e. bigger than log(max_condition_number)
	kappa_Arates_df = DataFrames.DataFrame(tvals=tvals, upper_bound_kappa_rates_A=upper_bound_kappa_rates_A, ln_max_condnum=log(max_condition_number), condbigTF=condbigTF)

	return(kappa_Arates_df)
end # END function check_linearDynamics_of_As








# Calculate Ds down a branch, using Louca & Pennell "Flow" (G) algorithm
#
# Modifies branchOp to do Ds calculation down a branch, using the G matrix
#
# The fakeX0s are also stored, for re-use if needed
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
"""
cd("/GitHub/PhyBEARS.jl/notes/")
#include("tst_Flow.jl")

include("ModelLikes.jl")
import .ModelLikes
#using .Tmp

include("Flow.jl")
import .Flow

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Profile     # for @profile
using DataFrames  # for DataFrame
using PhyloBits
using PhyloBits.TrUtils # for flat2() (similar to unlist)
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
using PhyloBits

inputs = ModelLikes.setup_DEC_SSE(2, readTopology("((chimp:10,human:10):10,gorilla:20);"))
# inputs = ModelLikes.setup_MuSSE(2, readTopology("((chimp:10,human:10):10,gorilla:20);"))
res = inputs.res
trdf = inputs.trdf
n = inputs.p_Ds_v5.n
solver_options = inputs.solver_options
solver_options.save_everystep
p_Ds_v5 = inputs.p_Ds_v5  # contains model parameters, and the "Es" solver/interpolator
trdf = inputs.trdf
root_age = maximum(trdf[!, :node_age])

# Ground truth with standard ClaSSE integration
u0 = collect(repeat([0.0], n)) # likelihoods at the branch top
u0[2] = 1.0
tspan = (0.0, 10.0)   # age at the branch top and branch bottom
prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
ground_truth_Ds_interpolator = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorTsit5 = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorCvode = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
	
# Set up the G flow
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2)
A = reshape(tmpzero, (n,n))
G0 = Matrix{Float64}(I, n, n) # initialize flow matrix G

pG = (n=n, p_Ds_v5=p_Ds_v5, A=A)
tspan_for_G = (0.0, 1.1*root_age) # Extend well beyond the root to avoid weirdness at the end
prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan_for_G, pG)


tspan = (0.0, 10.0)   # age at the branch top and branch bottom
u0 = collect(repeat([0.0], n)) # likelihoods at the branch top
u0[2] = 1.0

# Set up the interpolator for G
Gflow_to_01_Cvode  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Tsit5  = solve(prob_Gs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Lsoda  = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Get node ages > 0.0
tmp_node_ages = trdf[!,:node_age]
node_ages = sort(tmp_node_ages[tmp_node_ages .> 0.0])

Gflow = Gflow_to_01_Cvode

current_nodeIndex = 1

(tmp_threadID, sol_Ds, fakeX0s, spawned_nodeIndex, calc_start_time) = Flow.branchOp_ClaSSE_Gs_v5(current_nodeIndex, trdf; u0=u0, tspan=tspan, Gflow=Gflow, solver_options=solver_options)

# Calculate ClaSSE likelihood with "using the Flow, man"
# Modifies inputs into object "res"
(total_calctime_in_sec, iteration_number) = Flow.iterative_downpass_nonparallel_FlowClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10)
resFlow = deepcopy(res)
resFlow.likes_at_each_nodeIndex_branchTop
resFlow.logsum_likes_at_nodes
sum(resFlow.logsum_likes_at_nodes)

# Compare to "standard" ClaSSE
(total_calctime_in_sec, iteration_number) = Flow.iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10)
resClaSSE = deepcopy(res)
resClaSSE.likes_at_each_nodeIndex_branchTop
resClaSSE.logsum_likes_at_nodes
sum(resClaSSE.logsum_likes_at_nodes)

# Compare
sum(resFlow.logsum_likes_at_nodes)
sum(resClaSSE.logsum_likes_at_nodes)
sum(resFlow.logsum_likes_at_nodes) - sum(resClaSSE.logsum_likes_at_nodes)

"""

function branchOp_ClaSSE_Gs_v5(current_nodeIndex, trdf; u0, tspan, Gflow, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Is it a tip?
	parent_node_indices = trdf[!,:leftNodeIndex]
	tipTF = parent_node_indices[current_nodeIndex] == -999
	
	if (tipTF == true)
		branchBot_vals = Gflow(tspan[2]) * u0
		fakeX0s = u0
	else
		# Reverse Gflow from current_nodeIndex up to t=0
		# u0 = starting likelihoods at the top of this branch
		fakeX0s = factorize(Gflow(tspan[1])) \ u0
		branchBot_vals = Gflow(tspan[2]) * fakeX0s
	end
	
	# Example slow operation
	#y = countloop(num_iterations, current_nodeIndex)
	#prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
	#sol_Ds = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol)
	sol_Ds = branchBot_vals

	#nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	#nodeData_at_bottom = nodeData_at_top / 2.0
	#nodeData_at_bottom = sol_Ds.u
	
	return(tmp_threadID, sol_Ds, fakeX0s, spawned_nodeIndex, calc_start_time)
end



# Assuming input 
function branchOp_ClaSSE_Gs_v6(current_nodeIndex, trdf; u0, tspan, Gflow, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Is it a tip?
	parent_node_indices = trdf[!,:leftNodeIndex]
	tipTF = parent_node_indices[current_nodeIndex] == -999
	
	if (tipTF == true)
		branchBot_vals = Gflow(tspan[2]) * u0
		fakeX0s = u0
	else
		# Reverse Gflow from current_nodeIndex up to t=0
		# u0 = starting likelihoods at the top of this branch
		fakeX0s = factorize(Gflow(tspan[1])) \ u0
		branchBot_vals = Gflow(tspan[2]) * fakeX0s
	end
	
	# Example slow operation
	#y = countloop(num_iterations, current_nodeIndex)
	#prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
	#sol_Ds = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol)
	sol_Ds = branchBot_vals

	#nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	#nodeData_at_bottom = nodeData_at_top / 2.0
	#nodeData_at_bottom = sol_Ds.u
	
	return(tmp_threadID, sol_Ds, fakeX0s, spawned_nodeIndex, calc_start_time)
end






"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_nonparallel_FlowClaSSE_v5!(res; trdf, p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10)
	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()
	
	# Setup
	current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "done"
				return()
			end
	
			# Retrieve the inputs for the calculation down the branch
			
			# Use the RAW likelihoods (don't normalize)
			u0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
			u0 = u0 ./ (sum(u0))
			
			# Use the NORMALIZED (rescaled to sum to 1) likelihoods
			# Doesn't work -- claims an interpolation error for going beyond range
			# branchOp on current_nodeIndex=4ERROR: LoadError: Solution interpolation 
			# cannot extrapolate past the final timepoint. Either solve on a longer 
			# timespan or use the local extrapolation from the integrator interface.
			#u0 = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
			
			brlen = trdf[current_nodeIndex, :brlen]
			age_branchtop = trdf[current_nodeIndex, :node_age]
			age_branchbot = age_branchtop + brlen
			tspan = (age_branchtop, age_branchbot)
			#p_Ds_v5 = inputs.p_Ds_v5

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
# 			if (parallel_TF == true)
# 				push!(tasks, Base.Threads.@spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
# 			else
				tmp_results = Flow.branchOp_ClaSSE_Gs_v5(current_nodeIndex, trdf; u0=u0, tspan=tspan, Gflow=Gflow, solver_options=solver_options)
				#tmp_results = branchOp(current_nodeIndex, res, num_iterations)
				push!(tasks, tmp_results)
# 			end
			push!(tasks_fetched_TF, false)
		end
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			if (tasks_fetched_TF[i] == false)
				#if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					(tmp_threadID, sol_Ds, fakeX0s, spawned_nodeIndex, calc_start_time) = tasks[i]
					#nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
					nodeData_at_bottom = sol_Ds .+ 0.0  # Gflow just outputs the vector
# 					end
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
# 					print("\n\n12345\n\n")
# 					print("res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex]:\n")
# 					print(res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex])
# 					print("\n\nnodeData_at_bottom:\n")
# 					print(nodeData_at_bottom)
# 					print("\n\n12345\n\n")
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					# Get the ancestor nodeIndex
					uppass_edgematrix = res.uppass_edgematrix
					TF = uppass_edgematrix[:,2] .== spawned_nodeIndex
					parent_nodeIndex = uppass_edgematrix[TF,1][1]

					# Get the left daughter nodeIndex (1st in the uppass_edgematrix)
					edge_rows_TF = uppass_edgematrix[:,1] .== parent_nodeIndex
					left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
					right_nodeIndex = uppass_edgematrix[edge_rows_TF,2][2]

					# Update the state of the parent_node's daughters
					if (spawned_nodeIndex == left_nodeIndex)
						res.node_Lparent_state[parent_nodeIndex] = "ready"
					end
					if (spawned_nodeIndex == right_nodeIndex)
						res.node_Rparent_state[parent_nodeIndex] = "ready"
					end

					# Update the state of the current node
					res.node_state[spawned_nodeIndex] = "done"
				#end
			end
		end
	
		# Update which nodes have had both parents complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = res.node_Rparent_state .== "ready"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"
	
		# Update nodes when the branches above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, Base.Threads.@spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			nodeOp_ClaSSE_v5(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass_nonparallel(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
			error(txt)
		end
		
		# Test for concluding the while loop
		are_we_done && break
	end
	
	# This breaks it for some reason:
	# ERROR: setfield! immutable struct of type Res cannot be changed
	#global res.number_of_whileLoop_iterations = iteration_number

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	return(total_calctime_in_sec, iteration_number)
end # END iterative_downpass_nonparallel_FlowClaSSE_v5!





#######################################################
# Notes on matrix norms from Louca & Pennell, and 
# how to implement in PhyBEARS. NJM 2022-03-23
#######################################################
# Louca & Pennell:
# "A standard measure for how close the matrix G(t) is to singularity for numerical
# purposes is the “condition number”, denoted kappa(t) (Cline et al., 1979, Turing, 1948); a greater kappa(t) generally
# means that G(t) is harder to invert numerically. The condition number is given by the ratio of the largest
# over the smallest singular value, s1/sn, where s1, .., sn are the singular values in decreasing size (Watkins,
# 2010)."
# 
# An upper bound on kappa_Gt can be found:
# obs_condnum_kappa_Gt_upper_bound = exp(2*t*sigma1_A)   # upper bound on kappa, the norm of Gt
# ln_obs_condnum_kappa_Gt_upper_bound = 2*t*sigma1_A     # upper kappa rate, the log of the norm of Gt
#
#
# castor code:
# // calculate largest singular value (sigma1) of the dynamics at this age
# // then kappa_rate <= 2*sigma1, since sigma1(exp(t*A))/sigma2(exp(t*A)) <= exp(t*2*sigma1(A)) [So and Thompson (2000)
# Singular values of Matrix Exponentials. Theorem 2.1]
#		double kappa_rate;
#		const double sigma1 = get_largest_singular_value(Nstates, Nstates, dynamics, 1000, 1e-3);
#		if(!std::isnan(sigma1)){
#			kappa_rate = 2*sigma1;
####################################################### 
# Julia Linear Algebra:
# https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/
# When p=2, the operator norm is the spectral norm, equal to the largest singular value of A.
#######################################################
# NJM says: 
#
# The instantaneous rate matrix A only changes as a result of the Es changing (increasing back in time).
# 
# Now, sigma1(A) is itself the 2-norm of A. But the 2-norm is expensive.  These notes:
#
# https://www.math.usm.edu/lambers/mat610/sum10/lecture2.pdf
# 
# ...say (pp. 3-4):
#
# m = #rows
# n = #columns
# 
# That is, the (l)1-norm of a matrix is its maximum column sum.
# That is, the (l)Inf-norm of a matrix is its maximum row sum.
# 
# "That is, the (l)2-norm of a matrix is the square root of the largest eigenvalue of ATA, 
# which is guaranteed to be nonnegative, as can be shown using the vector 2-norm. We see 
# that unlike the vector (l)2-norm, the matrix (l)2-norm is much more difficult to 
# compute than the matrix (l)1-norm or (l)Inf-norm."
#
# opnorm(A,2) <= norm(A) (Frobenius norm) <= sqrt(n)*opnorm(A,2)
# 
# 1/sqrt(n)*opnorm(A,Inf) <= opnorm(A,2) <= sqrt(m)*opnorm(A,Inf)
#
# We need a cheap upper-bound on sigma1A (i.e., "How big could our opnorm(A,2) possibly be?"), this could be: 
# 
# sqrt(m)*opnorm(A,Inf)
#
# Or, alternatively, sqrt(m)*opnorm(transpose(A),1)
# 
# (Matrix A in our usage is square, so m=n)
#
# On the other hand, we only have to calculate the opnorm(A,2) for a series of 
# tvalues once per likelihood calculation, probably not rate-limiting!
# 
# Calculate G flow matrix down time series, outputting various
# condition numbers and kappa rates
calc_Gs_SSE_condnums! = (dG, G, pG, t) -> begin
	#A = pG.A             # Initial A
	p_Ds_v5 = pG.p_Ds_v5 # Calculate Ds down a timespan
	max_allowed_condition_number = 1e4 # (slower, more accurate)
	max_allowed_condition_number = 1e8 # (faster, less accurate)
	ln_max_allowed_condition_number = log(max_allowed_condition_number)
	
	# A as a function of time t
	A = parameterized_ClaSSE_As_v5(t, pG.A[:,:], pG.p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)
#	mul!(dG, A, G)
	mul!(dG, parameterized_ClaSSE_As_v5(t, pG.A[:,:], p_Ds_v5), G)
	
	# Condition numbers of G according to Julia
	# cond(M, p::Real=2)
	# "Condition number of the matrix M, computed using the operator p-norm. Valid values for p are 1, 2 (default), or Inf."
	condG1 = cond(G,1)
	condG2 = cond(G,2)
	condGInf = cond(G,Inf)
	tmpstr = paste0(["\nAt time t=", string(round(t, digits=6)), ", Julia's cond() gives condition numbers of G=", string(round(condG1, digits=6)), ", ", string(round(condG2, digits=6)), ", ", string(round(condGInf, digits=6))])
	print(tmpstr)
	#display(cond(G))
	
	print("\nopnorms(A) with p=1, p=2, p=Inf, opnorm(transpose(A),1), sqrt(nrow(A))*||A||Inf:\n")
	print(paste0(["opnorm(A,1) = ", opnorm(A,1)]))
	print(paste0(["\nopnorm(A,2) = ", opnorm(A,2)]))
	print(paste0(["\nopnorm(A,Inf) = ", opnorm(A,Inf)]))
	print(paste0(["\nopnorm(transpose(A),1) = ", opnorm(transpose(A),1)]))
	print(paste0(["\nFast min of opnorm(A,2): 1/(sqrt(Rncol(A)))*opnorm(A,Inf) = ", 1/(sqrt(Rncol(A)))*opnorm(A,Inf)]))
	print(paste0(["\nFast max of opnorm(A,2): sqrt(Rnrow(A))*opnorm(A,Inf) = ", sqrt(Rnrow(A))*opnorm(A,Inf)]))
	#print("\n")
	#display(opnorm(A,Inf))

	true_sigma1_of_A = opnorm(A,2)
	upper_bound_condition_number = exp(2*t*opnorm(A,2))
	upper_bound_kappa_growth_rate = 2*true_sigma1_of_A
	cheap_upper_bound_kappa_growth_rate_p1 = 2*sqrt(Rnrow(A)) * opnorm(A,1)
	cheap_upper_bound_kappa_growth_rate_pInf = 2*sqrt(Rnrow(A)) * opnorm(A,Inf)
	print(paste0(["\ntrue 2*sigma1_of_A (using opnorm(A,2)), true upper_bound_condition_number of A=exp(2*t*sigma1_of_A) = ", string(round(upper_bound_condition_number,digits=2))]))
	print(paste0(["\ntrue upper_bound_kappa_growth_rate=2*sigma1_of_A = ", string(round(upper_bound_kappa_growth_rate,digits=2))]))
	print(paste0(["\ncheap upper_bound_kappa_growth_rate=2*sqrt(Rnrow(A)) * opnorm(A,1) = ", string(round(cheap_upper_bound_kappa_growth_rate_p1,digits=2))]))
	print(paste0(["\ncheap cheap_upper_bound_kappa_growth_rate_pInf=2*sqrt(Rnrow(A)) * opnorm(A,Inf) = ", string(round(cheap_upper_bound_kappa_growth_rate_pInf,digits=2))]))
	print(paste0(["\nln_max_allowed_condition_number=", string(round(ln_max_allowed_condition_number,digits=2))]))
	print("\n")
	
	# No return needed, what is returned is G (as .u)
	
	#display(dG)
	#return(dG)
end # End calc_Gs_SSE_condnums



# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
calc_Gs_SSE! = (dG, G, pG, t) -> begin
	# Have to use pG.A[:,:] to avoid overwriting A, which screws everything up!
	#A = parameterized_ClaSSE_As_v5(t, pG.A[:,:], pG.p_Ds_v5)
	#display(A)
	#dG = A * G
	# Have to use pG.A[:,:] to avoid overwriting A, which screws everything up!
	#mul!(dG, A(t), G(t))  # Equation 20
	mul!(dG, parameterized_ClaSSE_As_v6(t, pG.A[:,:], pG.p_Ds_v5), G)
	#display(dG)
	#return(dG)
end # End calc_Gs_SSE



# Calculate G flow matrix down time series
# (no printed outputs)
calc_Gs_SSE = (dG, G, pG, t; max_condition_number=1e8) -> begin
	n = pG.n
	#tmpzero = repeat([0.0], n^2)
	#A = reshape(tmpzero, (n,n))

	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v6(t, pG.A .+ 0.0, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)

	# The new dG is A %*% G
	mul!(dG, A, G)

	# No return needed, what is returned is G (as .u)
	return(dG)
end # End calc_Gs_SSE


# Calculate G flow matrix down time series
# (no printed outputs)
calc_Gs_SSE_v7simd = (dG, G, pG, t) -> begin
	n = pG.n
	#tmpzero = repeat([0.0], n^2)
	#A = reshape(tmpzero, (n,n))
	A = pG.A

	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v7_simd!(A, t, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)

	# The new dG is A %*% G
	mul!(dG, A, G)

	# No return needed, what is returned is G (as .u)
	return(dG)
end # End calc_Gs_SSE




# Calculate G flow matrix down time series
# (no printed outputs)
calc_Gs_SSE_parallel = (dG, G, pG, t; max_condition_number=1e8) -> begin
	n = pG.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v6_parallel(t, A, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)

	# The new dG is A %*% G
	mul!(dG, A, G)

	# No return needed, what is returned is G (as .u)
	return(dG)
end # End calc_Gs_SSE_parallel


calc_Gs_SSE_sub_i = (dG, G, pG, t; max_condition_number=1e8) -> begin
	n = pG.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v6_sub_i(t, A, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)

	# The new dG is A %*% G
	mul!(dG, A, G)

	# No return needed, what is returned is G (as .u)
	return(dG)
end # End Gmaps.calc_Gs_SSE_sub_i




run_Gs = (p_Ds_v5) -> begin
	n = p_Ds_v5.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	G0 = reshape(tmpzero, (n,n))
	for i in 1:n
		G0[i,i] = 1.0
	end
	G = G0
	
	for i in 1:100
		t = 0.01
		A = parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
		G = A * G
		#Base.showarray(G)
		display(G)
	end

end


# This version includes Xi, Xj in the A equation (but this can't be more efficient I don't think,
# since these will be different on every branch)
parameterized_ClaSSE_A_v5xx = (du,u,A,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	two = 1.0
	# Iterate through the ancestral states
  @inbounds for i in 1:n
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		A[i,i] = A[i,i] + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])#*u[i]  
		
		# case 2: anagenetic change
		@inbounds for m in 1:length(Qi_sub_i)
			A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_sub_i[m]]# * u[Qj_sub_i[m]]
		end
		
		# case 34: change + eventual extinction
		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			rate_sp_then_ex = (u[Ck_sub_i] * uE[Cj_sub_i]) + ( uE[Ck_sub_i])
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex
		end
		
# 		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
# 			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
# 			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
# 				 (u[Ck_sub_i].*uE[Cj_sub_i] 
# 			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))

  end # End @inbounds for i in 1:n
end





"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_wGflowArray_nonparallel!(res; trdf, p_Ds_v5, Gflow, Gseg_times, Gflows_array, Gflows_array_totals, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true)
	#######################################################
	# Use Gflow on a downpass
	#######################################################
	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()

	# Re-set the node states for new downpass
	TF = trdf[:,:nodeType] .== "tip"
	res.node_state[TF] .= "ready_for_branchOp"
	res.node_state[TF .== false] .= "not_ready"
	res.node_Lparent_state[TF] .= "NA"
	res.node_Lparent_state[TF .== false] .= "not_ready"
	res.node_Rparent_state[TF] .= "NA"
	res.node_Rparent_state[TF .== false] .= "not_ready"
		
	# Setup
	current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false
	
	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "done"
				return()
			end
	
			# Retrieve the inputs for the calculation down the branch
			
			# Use the RAW likelihoods (don't normalize)
			u0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
			#u0 = u0 ./ (sum(u0))
			
			# Use the NORMALIZED (rescaled to sum to 1) likelihoods
			# If you get this error: claims an interpolation error for going beyond range
			# branchOp on current_nodeIndex=4ERROR: LoadError: Solution interpolation 
			# cannot extrapolate past the final timepoint. Either solve on a longer 
			# timespan or use the local extrapolation from the integrator interface.
			#u0 = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
			
			brlen = trdf[current_nodeIndex, :brlen]
			age_branchtop = trdf[current_nodeIndex, :node_age]
			age_branchbot = age_branchtop + brlen
			tspan = (age_branchtop, age_branchbot)
			#p_Ds_v5 = inputs.p_Ds_v5

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			#print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
#			if (parallel_TF == true)
#				push!(tasks, Base.Threads.@spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			tmp_results = branchOp_ClaSSE_wGflowArray(current_nodeIndex, res, Gflow, Gseg_times, Gflows_array, Gflows_array_totals, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			#tmp_results = branchOp(current_nodeIndex, res, num_iterations)
			push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			if (tasks_fetched_TF[i] == false)
				#if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = tasks[i]
					
					# Differences for Gflow version
					#nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
					nodeData_at_bottom = sol_Ds .+ 0.0
# 					end
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
# 					print("\n\n12345\n\n")
# 					print("res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex]:\n")
# 					print(res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex])
# 					print("\n\nnodeData_at_bottom:\n")
# 					print(nodeData_at_bottom)
# 					print("\n\n12345\n\n")

					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
					res.lq_at_branchBot[spawned_nodeIndex] = log(sum_nodeData_at_bottom)
					res.like_at_branchBot[spawned_nodeIndex] = sum_nodeData_at_bottom

					# Get the ancestor nodeIndex
					uppass_edgematrix = res.uppass_edgematrix
					TF = uppass_edgematrix[:,2] .== spawned_nodeIndex
					parent_nodeIndex = uppass_edgematrix[TF,1][1]

					# Get the left daughter nodeIndex (1st in the uppass_edgematrix)
					edge_rows_TF = uppass_edgematrix[:,1] .== parent_nodeIndex
					
					# Error trap
					if (sum(edge_rows_TF) > 3)
						txt = "STOP error in: iterative_downpass_wGflowArray_nonparallel() - node has more than 2 edges"
					end
					
					# Standard bifurcating node
					if (sum(edge_rows_TF) == 2)
						left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
						right_nodeIndex = uppass_edgematrix[edge_rows_TF,2][2]

						# Update the state of the parent_node's daughters
						if (spawned_nodeIndex == left_nodeIndex)
							res.node_Lparent_state[parent_nodeIndex] = "ready"
						end
						if (spawned_nodeIndex == right_nodeIndex)
							res.node_Rparent_state[parent_nodeIndex] = "ready"
						end
					end 

					# Singleton node (assumes direct-ancestor nodes are always "left"
					if (sum(edge_rows_TF) == 1)
						left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
						res.node_Lparent_state[parent_nodeIndex] = "ready"
					end
					
					# Update the state of the current node
					res.node_state[spawned_nodeIndex] = "done"
				#end
			end
		end
	
		# Update which nodes are SINGLETONS and are complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = trdf.nodeType .== "direct"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"

		# Update nodes when the singletons above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, Base.Threads.@spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end
	
		
		# Update which nodes have had BOTH parents complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = res.node_Rparent_state .== "ready"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"
		
		# Update nodes when the branches above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")

		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, Base.Threads.@spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			
			#res = nodeOp_ClaSSE_v5!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			res = nodeOp_ClaSSE_v6!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass_nonparallel(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
			error(txt)
		end
		
		# Test for concluding the while loop
		are_we_done && break
	end # END while(are_we_done == false)
		
	# This breaks it for some reason:
	# ERROR: setfield! immutable struct of type Res cannot be changed
	#global res.number_of_whileLoop_iterations = iteration_number

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	# 2022-03-17_ update to new way of calculating
	#Julia_sum_lq = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
	Julia_sum_lq_old = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
	nonroot_nodes = get_nonrootnodes_trdf(trdf)
	sum_likes_internal_branch_tops = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])
	Julia_sum_lq = Julia_sum_lq_old + sum_likes_internal_branch_tops


	# Add the root probabilities

	# Assuming diversitree options:
	# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
	# i.e., the root state probs are just the root_Ds/sum(root_Ds)
	d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
	root_stateprobs = d_root_orig/sum(d_root_orig)
	rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
	Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL
	
	# Without normalization
	#Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL
	# *With* normalization
	#Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL + log(sum(res.sumLikes_at_node_at_branchTop[1:(length(res.sumLikes_at_node_at_branchTop)-1)]))
	# Redundant with: sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))[nonroot_nodes])


	
	# Consider the pure-birth log-likelihoods
	# Get basic tree info
	numTips = sum(trdf.nodeType .== "tip")
	numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")
	
	# The Yule-process ML birthrate is just (# internal nodes - 1)/total_tree_length
	ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
	yuleBirthRate = (numInternal-1) / ttl_tree_length
	yuleDeathRate = 0.0					# Yule process has 0 extinction
	bd = bd_liks_trdf(trdf, yuleBirthRate, yuleDeathRate)
	bd_lnL_noTopo = bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times

	# Convert to BioGeoBEARS lnL under Yule process assumption
	# Check if the first state/geographic range is null
	#if res.inputs.setup.states_list[1] == []
	#	include_null_range = true
	#end
	numstates = length(res.normlikes_at_each_nodeIndex_branchTop[1])
	equal_root_prob2 = log(1/(numstates-include_null_range)) 
	bgb_root_lnL = log(sum(d_root_orig)) + 1.0

	# res5t match
	res5t = Julia_sum_lq + equal_root_prob2 + bgb_root_lnL - (1-log(1/yuleBirthRate))

	# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
	bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
	bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2
	bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2
	bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)
	bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo

	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_wGflowArray_nonparallel!



"""
# Gflow with Gmap, i.e. step through a series of segments

Gseg_times = seq(0.0, 2.5, 0.001);  # Times for Gsegments
pG = (n=n, p_Ds_v5=p_inputsG, A=A);
(Gflows_array, Gflows_array_totals) = Flow.calc_Gmaps_SSE(pG, Gseg_times);

# Showing the Gmap strategy. But, it seems to require ultra-small increments in order to get close in lnLs
Gflow_tmp = Matrix{Float64}(I, p_Ds_v5.n, p_Ds_v5.n) ;
for i in 1:2500
	mul!(Gflow_tmp, Gflows_array[:,:,i], Gflow_tmp[:,:])
end
Gflow_tmp
Gflow(2.5)

# This is a very close match, but does not produce quite the same answers...

# The Gmap strategy works OK, with many increments...
res_Gflow_v6a = Flow.iterative_downpass_wGflowArray_nonparallel!(res; trdf, p_Ds_v5, Gflow, Gseg_times, Gflows_array, Gflows_array_totals, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true, include_null_range=true)
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6a, rootstates_lnL_GFv6a, Julia_total_lnLs1_GFv6a) = res_Gflow_v6a
# (0.208, 2, -5.917627381773604, -2.882205463015512, -8.799832844789115, -3.9520899056719028)

# ...but seems to break down with fewer...
Gseg_times = seq(0.0, 2.5, 0.1);  # Times for Gsegments
pG = (n=n, p_Ds_v5=p_inputsG, A=A);
(Gflows_array, Gflows_array_totals) = Flow.calc_Gmaps_SSE(pG, Gseg_times);
res_Gflow_v6b = Flow.iterative_downpass_wGflowArray_nonparallel!(res; trdf, p_Ds_v5, Gflow, Gseg_times, Gflows_array, Gflows_array_totals, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true, include_null_range=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6b, rootstates_lnL_GFv6b, Julia_total_lnLs1_GFv6b) = res_Gflow_v6b

"""

# Calculate Ds down a branch, using Gflow
#
# Modifies branchOp to do Ds calculation down a branch, using Gflows
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_ClaSSE_wGflowArray(current_nodeIndex, res, Gflow, Gseg_times, Gflows_array, Gflows_array_totals; tspan, p_Ds_v5, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()

	Gflow_total = Matrix{Float64}(I, p_Ds_v5.n, p_Ds_v5.n) ;
	Gflow_total_old = Matrix{Float64}(I, p_Ds_v5.n, p_Ds_v5.n) ;
	
	# We are calculating from the current_nodeIndex node to the bottom of the branch below it
	
	# Get the fakeX0 at current_nodeIndex
	if (tspan[1] == 0.0)
		fakeX0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	else
		# Given res, update the fakeX0 in res at this node
		#age_branchtop = trdf[current_nodeIndex, :node_age]
		age_branchtop = tspan[1]
		
		# Calculate the fakeX0 for use calculating Xp at the rootward end of this node's branch
		tmpX0s = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] # assumes normalized likelihoods

		# Key trick from Louca & Pennell

		# Multiply Gflow from the particular segments of interest
		start_seg = 1  # (segment 1 starts at 0.0 mybp)
		end_seg = find_time_in_time_segments(age_branchtop, Gseg_times)
#		Gflow_to_current_branchtop = Matrix{Float64}(I, p_Ds_v5.n, p_Ds_v5.n) ;
#		for index in start_seg:end_seg 
#			mul!(Gflow_to_current_branchtop, Gflows_array[:,:,index], Gflow_to_current_branchtop[:,:])
			#Alternate method
			#Gflow_to_current_branchtop = Gflows_array[:,:,index] * Gflow_to_current_branchtop[:,:]
#		end

#		fakeX0 = factorize(Gflow(age_branchtop)) \ tmpX0s
#		fakeX0 = factorize(Gflow_to_current_branchtop) \ tmpX0s
		fakeX0 = factorize(Gflows_array_totals[:,:,end_seg]) \ tmpX0s
		res.fakeX0s_at_each_nodeIndex_branchTop[current_nodeIndex] = fakeX0
	end
		
	# Solve for the likelihoods at the rootward end of the branch
	# (instantaneously after the speciation event at Xp, typically)
	#sol_Ds = Gflow(tspan[2]) * fakeX0

	# If the high condition number requires sub-segments and re-normalizing
	# Could replace with:
	# sol_Ds = Gflow(seg1) * Gflow(seg2) * Gflow(seg3) * fakeX0
	
	# Multiply Gflow from the particular segments of interest
	start_seg = find_time_in_time_segments(tspan[1], Gseg_times)
	end_seg = find_time_in_time_segments(tspan[2], Gseg_times)
	
	# Iterate through the segments, multiply the matrices by fakeX0
#	Xparent_tmp = deepcopy(fakeX0)
#	Gflow_to_current_branchtop = Matrix{Float64}(I, p_Ds_v5.n, p_Ds_v5.n) ;
#	mul!(Xparent_tmp, Gflows_array[:,:,start_seg], fakeX0)
#	for index in 1:end_seg 
#		mul!(Gflow_to_current_branchtop, Gflows_array[:,:,index], Gflow_to_current_branchtop[:,:])
		#mul!(Xparent_tmp, Gflows_array[:,:,index], Xparent_tmp[:])
		######## mul! works fine!!  Xparent_tmp = Gflows_array[:,:,index] * Xparent_tmp[:]
#	end
	sol_Ds = Gflows_array_totals[:,:,end_seg] * fakeX0
#	sol_Ds = Gflow_to_current_branchtop * fakeX0
#	fakeX0 = factorize(Gflows_array_totals[:,:,end_seg] \ tmpX0s

#	sol_Ds = Xparent_tmp
	
	# Normal Gflow
	#sol_Ds = Gflow(tspan[2]) * fakeX0

#	print("\n")	
#	print(Gflows_array_totals[:,:,end_seg] )

#	print("\n")	
#	print(Gflow(tspan[2]) )
	
	
	# Return a tuple of the results for updating res in the main function
	return(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)
end # END branchOp_ClaSSE_wGflowArray


"""
Gseg_times = seq(0.0, 2.4, 0.1)  # Times for Gsegments
t = 0.1
find_time_in_time_segments(t, Gseg_times)
t = 0.2
find_time_in_time_segments(t, Gseg_times)

"""
function find_time_in_time_segments(t, Gseg_times)
	diffs = abs.(Gseg_times .- t)
	minval = minimum(diffs)
	TF = diffs .== minval
	index = (1:length(TF))[TF][1] # Returns scalar (and, warning, first hit if identicals)
	return index
end






"""
include("/Users/nmat471/HD/GitHub/PhyBEARS.jl/notes/parameterized_ClaSSE_As_v7.jl")
"""



parameterized_ClaSSE_As_v7! = (A, t, p) -> begin
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  psi = p.params.psi_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals

	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	# Zero-out A
	A .= 0.0;
	
	# loop through i=states
	@inbounds for i in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		
		Qi_eq_i_index = p.p_TFs.Qi_eq_i_index[i]
		Ci_eq_i_index = p.p_TFs.Ci_eq_i_index[i]
		# For a particular i/state, loop through all of the transitions from that state,
		# for the Q matrix and the C matrix
		# Q matrix
		@inbounds @simd for mi in 1:length(Qi_eq_i_index)
			# looping through mi, with +=, does a sum
			A[i,i] -= Qij_vals[Qi_eq_i_index[mi]] # term2 / part of case1
			# 
			A[Qi_sub_i[mi],Qj_sub_i[mi]] += Qij_vals[Qi_eq_i_index[mi]] # term3 / case2
		end
		
		# C matrix
		@inbounds @simd for mi in 1:length(Ci_eq_i_index)
			# term1, part of case1
			A[i,i] -= Cijk_vals[Ci_eq_i_index[mi]]
			
			# term4, Case 3/4
			A[Ci_sub_i[mi],Cj_sub_i[mi]] += Cijk_vals[Ci_eq_i_index[mi]] * uE[Ck_sub_i[mi]]
			A[Ci_sub_i[mi],Ck_sub_i[mi]] += Cijk_vals[Ci_eq_i_index[mi]] * uE[Cj_sub_i[mi]]
		end		
		
		# Additional parts of term1
		A[i,i] -= p.params.mu_vals[i]
		A[i,i] -= p.params.psi_vals[i]
	end
	return(A)
end # END parameterized_ClaSSE_As_v7! = (A, t, p)


parameterized_ClaSSE_As_v7_simd! = (A, t, p) -> begin
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  psi = p.params.psi_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals

	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	# Zero-out A
	A .= 0.0;
	
	# loop through i=states
	@inbounds for i in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		
		Qi_eq_i_index = p.p_TFs.Qi_eq_i_index[i]
		Ci_eq_i_index = p.p_TFs.Ci_eq_i_index[i]
		# For a particular i/state, loop through all of the transitions from that state,
		# for the Q matrix and the C matrix
		# Q matrix
		sum_Qij_vals_inbounds_simd_A!(A, i, Qij_vals, Qi_sub_i, Qj_sub_i, Qi_eq_i_index)
		
		# C matrix
		sum_Cijk_vals_inbounds_simd_A!(A, i, Cijk_vals, Ci_eq_i_index, Ci_sub_i, Cj_sub_i, Ck_sub_i, uE)		
		
		# Additional parts of term1
		A[i,i] -= p.params.mu_vals[i]
		A[i,i] -= p.params.psi_vals[i]
	end
	return(A)
end # END parameterized_ClaSSE_As_v7_simd! = (A, t, p)





function sum_Qij_vals_inbounds_simd_A!(A, i, Qij_vals, Qi_sub_i, Qj_sub_i, Qi_eq_i_index)
	@inbounds @simd for mi in 1:length(Qi_eq_i_index)
		# looping through mi, with +=, does a sum
		A[i,i] -= Qij_vals[Qi_eq_i_index[mi]] # term2 / part of case1
		# 
		A[Qi_sub_i[mi],Qj_sub_i[mi]] += Qij_vals[Qi_eq_i_index[mi]] # term3 / case2
	end
	return A
end;



function sum_Cijk_vals_inbounds_simd_A!(A, i, Cijk_vals, Ci_eq_i_index, Ci_sub_i, Cj_sub_i, Ck_sub_i, uE)
	# C matrix
	@inbounds @simd for mi in 1:length(Ci_eq_i_index)
		# term1, part of case1
		A[i,i] -= Cijk_vals[Ci_eq_i_index[mi]]
		
		# term4, Case 3/4
		A[Ci_sub_i[mi],Cj_sub_i[mi]] += Cijk_vals[Ci_eq_i_index[mi]] * uE[Ck_sub_i[mi]]
		A[Ci_sub_i[mi],Ck_sub_i[mi]] += Cijk_vals[Ci_eq_i_index[mi]] * uE[Cj_sub_i[mi]]
	end		
	return A
end;



end # end module Flow

