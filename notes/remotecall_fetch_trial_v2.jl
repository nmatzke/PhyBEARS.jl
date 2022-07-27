# Using Distributed.spawnat
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
#using TreeTable 
@everywhere using PhyBEARS.SSEs
@everywhere using PhyloBits.TreeTable # for get_nonrootnodes_trdf
@everywhere using PhyloBits.TrUtils
#@everywhere include("/GitHub/PhyBEARS.jl/notes/remotecall_fetch_trial_v2.jl")


# v7c: reducing GC
# v7: This will load multiple jobs onto threads on main processor; x2 improvement
# Using @spawn - so adds to main thread, not other cores
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
function iterative_downpass_parallel_ClaSSE_v7cc!(res; trdf, p_Ds_v7, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true, printlevel=0)
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
	#current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false
	sol_Ds_alg = Any[]
	spawned_nodeIndices = Any[]

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
				res.node_state[current_nodeIndex] = "not_ready"
				break
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
			push!(spawned_nodeIndices, current_nodeIndex)
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			#print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
#			if (parallel_TF == true)
#				push!(tasks, Base.Threads.@spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			#tmp_results = branchOp_ClaSSE_Ds_v6(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			# This will load multiple jobs onto threads on main processor; x2 improvement
			push!(tasks, Base.Threads.@spawn branchOp_ClaSSE_Ds_v7_retDs(u0, tspan))
			# branchOp_ClaSSE_Ds_v7(current_nodeIndex, res, u0=deepcopy(u0), tspan=tspan, p_Ds_v7=deepcopy(p_Ds_v7), solver_options=solver_options))
			# MethodError: no method matching istaskstarted(::Future)
			#push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
		# Touch all the tasks to kick them off?
		# Print only every 100th iteration
		sleep(0.001) # A tiny wait time keeps it from hanging
		#if (mod(iteration_number, 50) == 0)
		if printlevel >= 1
			if (mod(iteration_number, 50) == 0)
				num_tasks_running = sum(istaskstarted.(tasks)) - sum(istaskdone.(tasks))
				print(paste0(["\n#tasks running=", num_tasks_running, ","]))
				end
		end
		
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
					#(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
					(tmp_threadID, sol_Ds, calc_start_time) = fetch(tasks[i])
					spawned_nodeIndex = spawned_nodeIndices[i]
					#tmp_threadID = worker_for_task[manual_iter]
					#calc_start_time = calc_start_times[manual_iter]
					#sol_Ds = fetch(tasks[manual_iter])
					finalize(tasks[i])
					#finalize(jobs[i])
					
					# Return worker to pool
					#old_worker = worker_for_task[manual_iter]
					#push!(curr_free_workers, old_worker)
					#print("\ncurr_free_workers:")
					#print(curr_free_workers)
					#print("\n")
										#(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = fetch(tasks[manual_iter])
#					push!(sol_Ds_alg, string(sol_Ds.alg))
#					nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
					push!(sol_Ds_alg, "manual")
					nodeData_at_bottom = sol_Ds .+ 0.0


					#push!(sol_Ds_alg, string(sol_Ds.alg))
					#nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
# 					end

					# Error checks for conditional likelihoods below 0.0 or above 1.0:
					msg = ""
					TF = nodeData_at_bottom .< 0.0
					if sum(TF) > 0
						nan_val = 1e-15
						#txt = paste0(["\nUnderflow issue in iterative_downpass_Gflow_nonparallel_v2!(): at nodeIndex ", string(spawned_nodeIndex), ", nodeData_at_bottom had negative values. Probably this is due to very bad parameter values. Correcting these values to nan_val=", string(nan_val), ". Printing nodeData_at_bottom to screen:\n"])
						#print(txt)
						#print(nodeData_at_bottom)
						nodeData_at_bottom[TF] .= nan_val
						msg = "u" # Underflow error
					end # END if sum(TF) > 0
					# Error check:
					TF = nodeData_at_bottom .> 1.0
					if sum(TF) > 0
						correction_val = 1.0
						#txt = paste0(["\Overflow issue in iterative_downpass_Gflow_nonparallel_v2!(): at nodeIndex ", string(spawned_nodeIndex), ", nodeData_at_bottom had values > 1.0. Probably this is due to very bad parameter values. Correcting these values so that the conditional likelihoods sum to 1.0=", string(nan_val), ". Printing nodeData_at_bottom to screen:\n"])
						#print(txt)
						#print(nodeData_at_bottom)
						nodeData_at_bottom[TF] .= correction_val
						nodeData_at_bottom .= nodeData_at_bottom ./ sum(nodeData_at_bottom)
						msg = "o" # Overflow error
					end # END if sum(TF) > 0


					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					
					# Task finished; update, and print, if desired
					tasks_fetched_TF[i] = true
					#if (mod(iteration_number, 100) == 0)
					if printlevel >= 1
						txt = paste0([i, msg, ","])
						print(txt)
					end


					
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
					#if sum_nodeData_at_bottom < 0.0
					#	sum_nodeData_at_bottom = 1e-10000
					#end
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
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v7)
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
			res = nodeOp_ClaSSE_v6!(current_nodeIndex, res, p_Ds_v5=p_Ds_v7)
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

	if printlevel >= 1
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
	nonroot_nodes = TreePass.get_nonrootnodes_trdf(trdf)
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
	
	# Extract the most common solver:
	most_common_solver = get_a_most_common_value(sol_Ds_alg)
	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, most_common_solver)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_parallel_ClaSSE_v7cc!




#######################################################
# Spawnat with workers assigned explicitly
#######################################################

# Using Distributed.spawnat, workers mannually assigned
"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Parallel version using manually-assigned workers
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
function iterative_downpass_parallel_ClaSSE_v7dd!(res; trdf, p_Ds_v7, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true, printlevel=0)
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
	#current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	jobs = Any[]
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false
	sol_Ds_alg = Any[]

	list_of_workers = Distributed.workers()
	curr_free_workers = Distributed.workers()
	worker_for_task = Any[]


	#res.node_state[res.node_state .== "calculating_branchOp"] .= "ready_for_branchOp"
	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		#print("\niteration_number: ")
		#print(iteration_number)
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		#print(paste0(["\nNumber of indexes_ready=", length(indexes_ready)]))
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			#res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "not_ready"
				break
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
			#print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
#			if (parallel_TF == true)
#				push!(tasks, Base.Threads.@spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			#tmp_results = branchOp_ClaSSE_Ds_v6(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			
			manual_check = """
			x = TreePass.branchOp_ClaSSE_Ds_v7(current_nodeIndex, res, u0=tmp_u0, tspan=tspan, p_Ds_v7=p_Ds_v7, solver_options=solver_options);
			(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = x;
			nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)]
			sum(nodeData_at_bottom)
			"""
			
			if (length(curr_free_workers) > 0)
				new_worker = curr_free_workers[1]
				deleteat!(curr_free_workers, 1)
				push!(worker_for_task, new_worker) # record which worker was used
				
				#print(["\nAdding task current_nodeIndex=", current_nodeIndex, " on worker ", new_worker])

				res.node_state[current_nodeIndex] = "calculating_branchOp"
				res.calc_spawn_start[current_nodeIndex] = Dates.now()
				
				u0_tmp=deepcopy(u0)
				#p_tmp=deepcopy(p_Ds_v7)
				
				#push!(jobs, Distributed.remotecall_fetch(TreePass.branchOp_ClaSSE_Ds_v7, new_worker, current_nodeIndex, res, u0=u0_tmp, tspan=tspan, p_Ds_v7=p_tmp, solver_options=solver_options))
				push!(jobs, Distributed.remotecall(branchOp_ClaSSE_Ds_v7, new_worker, current_nodeIndex, res, u0=u0_tmp, tspan=tspan, p_Ds_v7=p_Ds_v7, solver_options=solver_options))
				push!(tasks, @async fetch(jobs[length(jobs)]))  # fetch works secondarily on tasks
			# MethodError: no method matching istaskstarted(::Future)
			#push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
				push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
			end
		end # END for current_nodeIndex in indexes_ready
		#print(paste0(["\nlength(tasks)=", length(tasks)]))
		#print(paste0(["\nsum(istaskdone.(tasks))=", sum(istaskdone.(tasks))]))

		if printlevel >= 1
			if (mod(iteration_number, 50) == 0)
				num_tasks_running = sum(istaskstarted.(tasks)) - sum(istaskdone.(tasks))
				print(paste0(["\n#tasks running=", num_tasks_running, ","]))
			end
		end	

	
		# Check which jobs are done, fetch them, and update status of that node
		#print("\nCheck for finished jobs:")
		num_tasks = length(jobs)
		manual_iter = 0 
		for i in 1:num_tasks
			#print("\ni: ")
			#print(i)
			manual_iter = manual_iter+1
			if (tasks_fetched_TF[manual_iter] == false)
				if (istaskstarted(tasks[manual_iter]) == true) && (istaskdone(tasks[manual_iter]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = fetch(jobs[manual_iter])
					finalize(tasks[manual_iter])
					finalize(jobs[manual_iter])
					
					# Return worker to pool
					old_worker = worker_for_task[manual_iter]
					push!(curr_free_workers, old_worker)
					#print("\ncurr_free_workers:")
					#print(curr_free_workers)
					#print("\n")
					
					#(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = fetch(tasks[manual_iter])
					push!(sol_Ds_alg, string(sol_Ds.alg))
					nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
# 					end

					# Error checks for conditional likelihoods below 0.0 or above 1.0:
					msg = ""
					TF = nodeData_at_bottom .< 0.0
					if sum(TF) > 0
						nan_val = 1e-15
						#txt = paste0(["\nUnderflow issue in iterative_downpass_Gflow_nonparallel_v2!(): at nodeIndex ", string(spawned_nodeIndex), ", nodeData_at_bottom had negative values. Probably this is due to very bad parameter values. Correcting these values to nan_val=", string(nan_val), ". Printing nodeData_at_bottom to screen:\n"])
						#print(txt)
						#print(nodeData_at_bottom)
						nodeData_at_bottom[TF] .= nan_val
						msg = "u" # Underflow error
					end # END if sum(TF) > 0
					# Error check:
					TF = nodeData_at_bottom .> 1.0
					if sum(TF) > 0
						correction_val = 1.0
						#txt = paste0(["\Overflow issue in iterative_downpass_Gflow_nonparallel_v2!(): at nodeIndex ", string(spawned_nodeIndex), ", nodeData_at_bottom had values > 1.0. Probably this is due to very bad parameter values. Correcting these values so that the conditional likelihoods sum to 1.0=", string(nan_val), ". Printing nodeData_at_bottom to screen:\n"])
						#print(txt)
						#print(nodeData_at_bottom)
						nodeData_at_bottom[TF] .= correction_val
						nodeData_at_bottom .= nodeData_at_bottom ./ sum(nodeData_at_bottom)
						msg = "o" # Overflow error
					end # END if sum(TF) > 0


					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					
					# Task finished; update, and print, if desired
					tasks_fetched_TF[manual_iter] = true
					if printlevel >= 1
						txt = paste0([i, msg, ","])
						print(txt)
					end
					
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
					#if sum_nodeData_at_bottom < 0.0
					#	sum_nodeData_at_bottom = 1e-10000
					#end
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
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v7)
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
			res = nodeOp_ClaSSE_v6!(current_nodeIndex, res, p_Ds_v5=p_Ds_v7)
			# (updates res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		count_nodes_finished(res.node_state)
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

	if printlevel >= 1
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
	
	# Extract the most common solver:
	most_common_solver = get_a_most_common_value(sol_Ds_alg)
	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, most_common_solver)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_parallel_ClaSSE_v7dd!




function branchOp_ClaSSE_Ds_v7_retDs(u0, tspan)
	calc_start_time = Dates.now()
	tmp_threadID = Hwloc.getinfo()[:Core]
	
	prob_Ds_v7 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v7_simd_sums, u0, tspan, p_Ds_v7)
	sol_Ds = solve(prob_Ds_v7, solver_options.solver, dense=false, save_start=false, save_end=true, save_everystep=false, abstol=solver_options.abstol, reltol=solver_options.reltol)
	
	Ds = sol_Ds[length(sol_Ds)]
	return(tmp_threadID, Ds, calc_start_time)
end


#######################################################
# Spawnat with workers assigned via @spawnat :any ;
# Functions pre-allocated
#######################################################
function iterative_downpass_parallel_ClaSSE_v7ee!(res; trdf, p_Ds_v7, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true, printlevel=0)
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
	#current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	jobs = Any[]
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false
	sol_Ds_alg = Any[]
	spawned_nodeIndices = Any[]
	calc_start_times = Any[]
	list_of_workers = Distributed.workers()
	curr_free_workers = Distributed.workers()
	worker_for_task = Any[]


	#res.node_state[res.node_state .== "calculating_branchOp"] .= "ready_for_branchOp"
	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		#print("\niteration_number: ")
		#print(iteration_number)
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		#print(paste0(["\nNumber of indexes_ready=", length(indexes_ready)]))
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			#res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "not_ready"
				break
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
			#print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
#			if (parallel_TF == true)
#				push!(tasks, Base.Threads.@spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			#tmp_results = branchOp_ClaSSE_Ds_v6(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			
			manual_check = """
			x = TreePass.branchOp_ClaSSE_Ds_v7(current_nodeIndex, res, u0=tmp_u0, tspan=tspan, p_Ds_v7=p_Ds_v7, solver_options=solver_options);
			(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = x;
			nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)]
			sum(nodeData_at_bottom)
			"""
			
			#if (length(curr_free_workers) > 0)
				#new_worker = curr_free_workers[1]
				#deleteat!(curr_free_workers, 1)
				#push!(worker_for_task, new_worker) # record which worker was used
				
				#print(["\nAdding task current_nodeIndex=", current_nodeIndex, " on worker ", new_worker])

				res.node_state[current_nodeIndex] = "calculating_branchOp"
				res.calc_spawn_start[current_nodeIndex] = Dates.now()
				push!(spawned_nodeIndices, current_nodeIndex)
				push!(calc_start_times, Dates.now())
				u0_tmp=deepcopy(u0)
				#p_tmp=deepcopy(p_Ds_v7)
				
				#push!(jobs, Distributed.remotecall_fetch(TreePass.branchOp_ClaSSE_Ds_v7, new_worker, current_nodeIndex, res, u0=u0_tmp, tspan=tspan, p_Ds_v7=p_tmp, solver_options=solver_options))
				#push!(jobs, @spawnat :any branchOp_ClaSSE_Ds_v7_retDs(u0_tmp, tspan))
				push!(jobs, @spawnat :any branchOp_ClaSSE_Ds_v7_retDs(u0_tmp, tspan))
				push!(tasks, @async fetch(jobs[length(jobs)]))  # fetch works secondarily on tasks
			# MethodError: no method matching istaskstarted(::Future)
			#push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
				push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
			#end
		end # END for current_nodeIndex in indexes_ready
		#print(paste0(["\nlength(tasks)=", length(tasks)]))
		#print(paste0(["\nsum(istaskdone.(tasks))=", sum(istaskdone.(tasks))]))
		
		# Touch all the tasks to kick them off?
		# Print only every 100th iteration
		sleep(0.001) # A tiny wait time keeps it from hanging
		#if (mod(iteration_number, 50) == 0)
		if printlevel >= 1
			if (mod(iteration_number, 50) == 0)
				num_tasks_running = sum(istaskstarted.(tasks)) - sum(istaskdone.(tasks))
				print(paste0(["\n#tasks running=", num_tasks_running, ","]))
			end
		end
	
		# Check which jobs are done, fetch them, and update status of that node
		#print("\nCheck for finished jobs:")
		num_tasks = length(jobs)
		manual_iter = 0 
		for i in 1:num_tasks
			#print("\ni: ")
			#print(i)
			#manual_iter = manual_iter+1
			if (tasks_fetched_TF[i] == false)
				if (istaskstarted(tasks[i]) == true) && (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					#(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = fetch(tasks[manual_iter])
					(tmp_threadID, sol_Ds, calc_start_time) = fetch(tasks[i])
					spawned_nodeIndex = spawned_nodeIndices[i]
					#tmp_threadID = worker_for_task[manual_iter]
					#calc_start_time = calc_start_times[manual_iter]
					#sol_Ds = fetch(tasks[manual_iter])
					finalize(tasks[i])
					finalize(jobs[i])
					
					# Return worker to pool
					#old_worker = worker_for_task[manual_iter]
					#push!(curr_free_workers, old_worker)
					#print("\ncurr_free_workers:")
					#print(curr_free_workers)
					#print("\n")
										#(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = fetch(tasks[manual_iter])
#					push!(sol_Ds_alg, string(sol_Ds.alg))
#					nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
					push!(sol_Ds_alg, "manual")
					nodeData_at_bottom = sol_Ds .+ 0.0
# 					end

					# Error checks for conditional likelihoods below 0.0 or above 1.0:
					msg = ""
					TF = nodeData_at_bottom .< 0.0
					if sum(TF) > 0
						nan_val = 1e-15
						#txt = paste0(["\nUnderflow issue in iterative_downpass_Gflow_nonparallel_v2!(): at nodeIndex ", string(spawned_nodeIndex), ", nodeData_at_bottom had negative values. Probably this is due to very bad parameter values. Correcting these values to nan_val=", string(nan_val), ". Printing nodeData_at_bottom to screen:\n"])
						#print(txt)
						#print(nodeData_at_bottom)
						nodeData_at_bottom[TF] .= nan_val
						msg = "u" # Underflow error
					end # END if sum(TF) > 0
					# Error check:
					TF = nodeData_at_bottom .> 1.0
					if sum(TF) > 0
						correction_val = 1.0
						#txt = paste0(["\Overflow issue in iterative_downpass_Gflow_nonparallel_v2!(): at nodeIndex ", string(spawned_nodeIndex), ", nodeData_at_bottom had values > 1.0. Probably this is due to very bad parameter values. Correcting these values so that the conditional likelihoods sum to 1.0=", string(nan_val), ". Printing nodeData_at_bottom to screen:\n"])
						#print(txt)
						#print(nodeData_at_bottom)
						nodeData_at_bottom[TF] .= correction_val
						nodeData_at_bottom .= nodeData_at_bottom ./ sum(nodeData_at_bottom)
						msg = "o" # Overflow error
					end # END if sum(TF) > 0


					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					
					# Task finished; update, and print, if desired
					tasks_fetched_TF[i] = true
					#if (mod(iteration_number, 100) == 0)
					if printlevel >= 1
						txt = paste0([i, msg, ","])
						print(txt)
					end
					
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
					#if sum_nodeData_at_bottom < 0.0
					#	sum_nodeData_at_bottom = 1e-10000
					#end
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
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v7)
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
			res = nodeOp_ClaSSE_v6!(current_nodeIndex, res, p_Ds_v5=p_Ds_v7)
			# (updates res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		count_nodes_finished(res.node_state)
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

	if printlevel >= 1
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
	
	# Extract the most common solver:
	most_common_solver = get_a_most_common_value(sol_Ds_alg)
	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, most_common_solver)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_parallel_ClaSSE_v7ee!




