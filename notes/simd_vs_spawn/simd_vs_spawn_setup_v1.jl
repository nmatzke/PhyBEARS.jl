
using LoopVectorization

function sum_Cijk_rates_Ds_inbounds_simd_tmp(Cijk_rates_sub_i, tmp_u, tmp_uE, Cj_sub_i, Ck_sub_i)
    term1 = 0.0
    term4 = 0.0
    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
    	term4 += Cijk_rates_sub_i[it] * (tmp_u[Ck_sub_i[it]] * tmp_uE[Cj_sub_i[it]] + tmp_u[Cj_sub_i[it]] * tmp_uE[Ck_sub_i[it]])
    end
    return term1, term4
end;

function sum_Cijk_rates_Es_inbounds_simd_tmp(Cijk_rates_sub_i, tmp_u, Cj_sub_i, Ck_sub_i)
    term1 = 0.0
    term4 = 0.0
    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
    	term4 += Cijk_rates_sub_i[it] * tmp_u[Cj_sub_i[it]] * tmp_u[Ck_sub_i[it]]
    end
    return term1, term4
end;

function sum_Qij_vals_inbounds_simd_tmp(Qij_vals_sub_i, tmp_u, Qj_sub_i)
    term2 = 0.0
    term3 = 0.0
    @inbounds @simd for it=1:length(Qij_vals_sub_i)
    	term2 += Qij_vals_sub_i[it]
    	term3 += Qij_vals_sub_i[it] * tmp_u[Qj_sub_i[it]]
    end
    return term2, term3
end;




#Ds_v7_simd_sums_tmp = (du,u,p,t) -> begin
function Ds_v7_simd_sums_tmp(du,u,p,t)

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
	uE = p.sol_Es_v5(t)
	terms = p.terms
  for i in 1:n
		terms[1] = 0.0
		terms[2] = 0.0
		terms[3] = 0.0
		terms[4] = 0.0

		terms[1], terms[4] = sum_Cijk_rates_Ds_inbounds_simd_tmp(p.p_TFs.Cijk_rates_sub_i[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i])
	
		terms[2], terms[3] = sum_Qij_vals_inbounds_simd_tmp(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i])
		
		du[i] = -(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4]
  end
end




Ds_v8_simd_sums_tmp = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
#  Qij_vals = p.params.Qij_vals
#  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
#	Qarray_ivals = p.p_indices.Qarray_ivals
#	Qarray_jvals = p.p_indices.Qarray_jvals
#	Carray_ivals = p.p_indices.Carray_ivals
#	Carray_jvals = p.p_indices.Carray_jvals
#	Carray_kvals = p.p_indices.Carray_kvals

#	Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i
#	Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i
	
	# Pre-calculated solution of the Es
#	sol_Es = p.sol_Es_v5
#	uE = p.uE
	uE = p.sol_Es_v5(t)
	terms = p.terms
  for i in 1:n
		#Qi_sub_i = p.p_TFs.Qi_sub_i[i]
#		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		#Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		#Ci_sub_i = p.p_TFs.Ci_sub_i[i]
#		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
#		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		#Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "D" (likelihood of tip data)
#		du[i] = -(sum(Cijk_rates_sub_i[i]) + sum(Qij_vals_sub_i[i]) + mu[i])*u[i] +  # case 1: no event
#			(sum(Qij_vals_sub_i[i] .* u[Qj_sub_i])) + 	# case 2	
#			(sum(Cijk_rates_sub_i[i] .*                                               # case 3/4: change + eventual extinction
#				 (u[Ck_sub_i].*uE[Cj_sub_i] 
#			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
		
#		Cijk_rates_temp = @view p.p_TFs.Cijk_rates_sub_i[i]
#		Qij_rates_temp = @view p.p_TFs.Qij_vals_sub_i[i]
		
		terms[1] = 0.0
		terms[2] = 0.0
		terms[3] = 0.0
		terms[4] = 0.0

		terms[1], terms[4] = sum_Cijk_rates_Ds_tturbo(p.p_TFs.Cijk_rates_sub_i[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i])
	
		terms[2], terms[3] = sum_Qij_vals_tturbo(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i])
		
		du[i] = -(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4]
  end
end






function branchOp_ClaSSE_Ds_v5_tmp(current_nodeIndex, res; u0, tspan, p_Ds_v5, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Example slow operation
	#y = countloop(num_iterations, current_nodeIndex)
# 	print("\n")
# 	print("branchOp_ClaSSE_Ds_v5: d: ")
# 	print(p_Ds_v5.params.Qij_vals[1])
# 	print("branchOp_ClaSSE_Ds_v5: e: ")
# 	print(p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)])
# 	print("\n")
	
	prob_Ds_v5 = DifferentialEquations.ODEProblem(Ds_v5_tmp, u0, tspan, p_Ds_v5)

		sol_Ds = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol)
	
	# Error catch seems to slow it down!!
	"""
	try
		sol_Ds = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol)
	catch e
		txt = paste0(["\nError in branchOp_ClaSSE_Ds_v5(): solve(prob_Ds_v5, ", SubString(string(solver_options.solver), 1:15), "..., save_everystep=", string(solver_options.save_everystep), ", abstol=", string(solver_options.abstol), ", reltol=", string(solver_options.reltol), "): this error may indicate an impossible parameter combination (e.g. all rates are 0.0). Returning 1e-15 for all Ds."])
		print(txt)
		# Return a fake sol_Ds function
		nan_val = 1e-15
		function sol_Ds_tmp(t, tmp_n)
			res = collect(repeat([nan_val], tmp_n))
		end
		sol_Ds = t -> sol_Ds_tmp(t, length(u0)) # returns a array of nan_val, regardless of time t input
	end
	"""
# 	print(sol_Ds[length(sol_Ds)])
# 	print("\n")
	
	#nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	#nodeData_at_bottom = nodeData_at_top / 2.0
	#nodeData_at_bottom = sol_Ds.u
	
	return(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)
end



# Calculate Ds down a branch
#
# Modifies branchOp to do Ds calculation down a branch
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_ClaSSE_Ds_v6_tmp(current_nodeIndex, res; u0, tspan, p_Ds_v5, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Example slow operation
	#y = countloop(num_iterations, current_nodeIndex)
# 	print("\n")
# 	print("branchOp_ClaSSE_Ds_v5: d: ")
# 	print(p_Ds_v5.params.Qij_vals[1])
# 	print("branchOp_ClaSSE_Ds_v5: e: ")
# 	print(p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)])
# 	print("\n")
	
	prob_Ds_v5 = DifferentialEquations.ODEProblem(Ds_v6_tmp, u0, tspan, p_Ds_v5)
	#sol_Ds = solve(prob_Ds_v5, solver_options.solver, saveat=solver_options.saveat, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol, dense=solver_options.dense)
	# Better: when we are solving on a particular branch, all we need is the 
	# end-point (bottom of the branch); so, dense=false and saveat=[] (empty).
	sol_Ds = solve(prob_Ds_v5, solver_options.solver, saveat=solver_options.saveat, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol)

	
# 	print(sol_Ds[length(sol_Ds)])
# 	print("\n")
	
	#nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	#nodeData_at_bottom = nodeData_at_top / 2.0
	#nodeData_at_bottom = sol_Ds.u
	
	return(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)
end # END branchOp_ClaSSE_Ds_v6






# Calculate Ds down a branch
#
# Modifies branchOp to do Ds calculation down a branch
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_ClaSSE_Ds_v6_solverFree(current_nodeIndex, res; u0, tspan, p_Ds_v5, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Example slow operation
	#y = countloop(num_iterations, current_nodeIndex)
# 	print("\n")
# 	print("branchOp_ClaSSE_Ds_v5: d: ")
# 	print(p_Ds_v5.params.Qij_vals[1])
# 	print("branchOp_ClaSSE_Ds_v5: e: ")
# 	print(p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)])



# Calculate Ds down a branch
function branchOp_ClaSSE_Ds_v7_tmp(current_nodeIndex, res; u0, tspan, p_Ds_v7, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	prob_Ds_v7 = DifferentialEquations.ODEProblem(Ds_v7_simd_sums_tmp2, u0, tspan, p_Ds_v7)

	sol_Ds = solve(prob_Ds_v7, solver_options.solver, dense=false, save_start=false, save_end=true, save_everystep=false, abstol=solver_options.abstol, reltol=solver_options.reltol)
	
	return(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)
end






function parallel_down_tree!(res; trdf, p_Ds_v7, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true)
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
	sol_Ds_alg = Any[]

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
			push!(tasks, Base.Threads.@spawn branchOp_ClaSSE_Ds_v7_tmp(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v7=p_Ds_v7, solver_options=solver_options))
			# MethodError: no method matching istaskstarted(::Future)
			#push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			#print("\nTask i: ")
			#print(i)
			#print(": ")
			#print(tasks[i])
			if (tasks_fetched_TF[i] == false)
				tasks[i]
				if (istaskstarted(tasks[i]) == true) && (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					#(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = tasks[i]
					#print("\nFetching #")
					#print(i)
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
					#print("...fetched.")
					push!(sol_Ds_alg, string(sol_Ds.alg))
					nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
# 					end
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
# 					#print("\n\n12345\n\n")
# 					#print("res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex]:\n")
# 					#print(res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex])
# 					#print("\n\nnodeData_at_bottom:\n")
# 					#print(nodeData_at_bottom)
# 					#print("\n\n12345\n\n")

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

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		#print(txt)
		#print("\n")
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
end # END parallel_down_tree!





#######################################################
# Data
#######################################################

function load_ps_511()
	tmpstr = "(n = 511, params = (mu_vals = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, " ⋯ 18856777 bytes ⋯ " 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], terms = [0.10000000000000003, 0.18, 0.018000000001105538, 0.0010000000001136261])"
	
	# https://discourse.julialang.org/t/reading-repr-string-to-object/28548
	p_Ds_v5 = eval(Meta.parse(tmpstr))
	return p_Ds_v5
end


function load_geog_511()

tmpstr="""
tipnames	A	B	C	D	E	F	G	H	I
Agathis_australis	0	0	0	0	0	0	0	0	1
Agathis_atropurpurea	0	0	0	0	0	0	1	0	0
Agathis_borneensis	0	0	0	1	0	0	0	0	0
Agathis_dammara	0	0	0	1	0	0	0	0	0
Agathis_lanceolata	0	0	0	0	0	0	0	1	0
Agathis_macrophylla	0	0	0	0	1	1	0	0	0
Agathis_microstachya	0	0	0	0	0	0	1	0	0
Agathis_montana	0	0	0	0	0	0	0	1	0
Agathis_moorei	0	0	0	0	0	0	0	1	0
Agathis_ovata	0	0	0	0	0	0	0	1	0
Agathis_robusta	0	0	0	0	1	0	1	0	0
Agathis_silbae	0	0	0	0	1	0	0	0	0
Agathis_alba	0	0	0	1	0	0	0	0	0
Agathis_corbassonii	0	0	0	0	0	0	0	1	0
Agathis_endertii	0	0	0	1	0	0	0	0	0
Agathis_obtusa	0	0	0	0	1	0	0	0	0
Agathis_palmerstonii	0	0	0	0	1	0	1	0	0
Agathis_philippinensis	0	0	0	1	0	0	0	0	0
Agathis_vitiensis	0	0	0	0	1	1	0	0	0
Araucaria_angustifolia	1	0	0	0	0	0	0	0	0
Araucaria_araucana	1	0	0	0	0	0	0	0	0
Araucaria_bidwillii	0	0	0	0	0	0	1	0	0
Araucaria_bernieri	0	0	0	0	0	0	0	1	0
Araucaria_biramulata	0	0	0	0	0	0	0	1	0
Araucaria_columnaris	0	0	0	0	0	0	0	1	0
Araucaria_cunninghamii	0	0	0	0	0	0	1	0	0
Araucaria_excelsa	0	0	0	0	0	0	1	0	0
Araucaria_heterophylla	0	0	0	0	0	0	1	0	0
Araucaria_humboldtensis	0	0	0	0	0	0	0	1	0
Araucaria_hunsteinii	0	0	0	0	1	0	0	0	0
Araucaria_laubenfelsii	0	0	0	0	0	0	0	1	0
Araucaria_luxurians	0	0	0	0	0	0	0	1	0
Araucaria_montana	0	0	0	0	0	0	0	1	0
Araucaria_muelleri	0	0	0	0	0	0	0	1	0
Araucaria_nemorosa	0	0	0	0	0	0	0	1	0
Araucaria_rulei	0	0	0	0	0	0	0	1	0
Araucaria_schmidii	0	0	0	0	0	0	0	1	0
Araucaria_scopulorum	0	0	0	0	0	0	0	1	0
Araucaria_subulata	0	0	0	0	0	0	0	1	0
Wollemia_nobilis	0	0	0	0	0	0	1	0	0
Lepidothamnus_laxifolius	0	0	0	0	0	0	0	0	1
Lepidothamnus_fonkii	1	0	0	0	0	0	0	0	0
Lepidothamnus_intermedius	0	0	0	0	0	0	0	0	1
Halocarpus_biformis	0	0	0	0	0	0	0	0	1
Halocarpus_bidwillii	0	0	0	0	0	0	0	0	1
Halocarpus_kirkii	0	0	0	0	0	0	0	0	1
Phyllocladus_alpinus	0	0	0	0	0	0	0	0	1
Phyllocladus_trichomanoides_var_trichomanoides	0	0	0	0	0	0	0	0	1
Phyllocladus_aspleniifolius	0	0	0	0	0	0	1	0	0
Phyllocladus_hypophyllus	0	0	0	1	1	0	0	0	0
Phyllocladus_toatoa	0	0	0	0	0	0	0	0	1
Parasitaxus_usta	0	0	0	0	0	0	0	1	0
Lagarostrobos_franklinii	0	0	0	0	0	0	1	0	0
Manoao_colensoi	0	0	0	0	0	0	0	0	1
Prumnopitys_ladei	0	0	0	0	0	0	1	0	0
Prumnopitys_ferruginea	0	0	0	0	0	0	0	0	1
Prumnopitys_ferruginoides	0	0	0	0	0	0	0	1	0
Prumnopitys_amara	0	0	0	1	1	0	1	0	0
Prumnopitys_andina	1	0	0	0	0	0	0	0	0
Prumnopitys_exigua	1	0	0	0	0	0	0	0	0
Prumnopitys_montana	1	0	0	0	0	0	0	0	0
Prumnopitys_taxifolia	0	0	0	0	0	0	0	0	1
Saxegothaea_conspicua	1	0	0	0	0	0	0	0	0
Microcachrys_tetragona	0	0	0	0	0	0	1	0	0
Pherosphaera_fitzgeraldii	0	0	0	0	0	0	1	0	0
Pherosphaera_hookeriana	0	0	0	0	0	0	1	0	0
Acmopyle_pancheri	0	0	0	0	0	0	0	1	0
Acmopyle_sahniana	0	0	0	0	0	1	0	0	0
Dacrycarpus_vieillardii	0	0	0	0	0	0	0	1	0
Dacrycarpus_dacrydioides	0	0	0	0	0	0	0	0	1
Dacrycarpus_compactus	0	0	0	0	1	0	0	0	0
Dacrycarpus_expansus	0	0	0	0	1	0	0	0	0
Dacrycarpus_cinctus	0	0	0	1	1	0	0	0	0
Dacrycarpus_kinabaluensis	0	0	0	1	0	0	0	0	0
Falcatifolium_falciforme	0	0	0	1	0	0	0	0	0
Falcatifolium_gruezoi	0	0	0	1	0	0	0	0	0
Falcatifolium_papuanum	0	0	0	0	1	0	0	0	0
Falcatifolium_taxoides	0	0	0	0	0	0	0	1	0
Dacrydium_gibbsiae	0	0	0	1	0	0	0	0	0
Dacrydium_elatum	0	0	1	1	0	0	0	0	0
Dacrydium_pectinatum	0	0	1	1	0	0	0	0	0
Dacrydium_comosum	0	0	0	1	0	0	0	0	0
Dacrydium_magnum	0	0	0	1	1	0	0	0	0
Dacrydium_gracile	0	0	0	1	0	0	0	0	0
Dacrydium_xanthandrum	0	0	0	1	1	0	0	0	0
Dacrydium_beccarii	0	0	0	1	1	0	0	0	0
Dacrydium_cupressinum	0	0	0	0	0	0	0	0	1
Dacrydium_x_suprinii	0	0	0	0	0	0	0	1	0
Dacrydium_guillauminii	0	0	0	0	0	0	0	1	0
Dacrydium_lycopodioides	0	0	0	0	0	0	0	1	0
Dacrydium_araucarioides	0	0	0	0	0	0	0	1	0
Dacrydium_balansae	0	0	0	0	0	0	0	1	0
Dacrydium_nausoriense	0	0	0	0	0	1	0	0	0
Dacrydium_nidulum	0	0	0	1	1	1	0	0	0
Retrophyllum_rospigliosii	1	0	0	0	0	0	0	0	0
Retrophyllum_vitiense	0	0	0	1	1	1	0	0	0
Retrophyllum_comptonii	0	0	0	0	0	0	0	1	0
Retrophyllum_minus	0	0	0	0	0	0	0	1	0
Afrocarpus_gaussenii	0	1	0	0	0	0	0	0	0
Afrocarpus_gracilior	0	1	0	0	0	0	0	0	0
Afrocarpus_mannii	0	1	0	0	0	0	0	0	0
Afrocarpus_falcatus	0	1	0	0	0	0	0	0	0
Afrocarpus_dawei	0	1	0	0	0	0	0	0	0
Afrocarpus_usambarensis	0	1	0	0	0	0	0	0	0
Nageia_motleyi	0	0	1	1	0	0	0	0	0
Nageia_wallichiana	0	0	1	1	1	0	0	0	0
Nageia_fleuryi	0	0	1	0	0	0	0	0	0
Nageia_formosensis	0	0	1	0	0	0	0	0	0
Nageia_nagi	0	0	1	0	0	0	0	0	0
Nageia_nankoensis	0	0	1	0	0	0	0	0	0
Podocarpus_smithii	0	0	0	0	0	0	1	0	0
Podocarpus_alpinus	0	0	0	0	0	0	1	0	0
Podocarpus_nivalis	0	0	0	0	0	0	0	0	1
Podocarpus_gnidioides	0	0	0	0	0	0	0	1	0
Podocarpus_lawrencei	0	0	0	0	0	0	1	0	0
Podocarpus_acutifolius	0	0	0	0	0	0	0	0	1
Podocarpus_totara	0	0	0	0	0	0	0	0	1
Podocarpus_cunninghamii	0	0	0	0	0	0	0	0	1
Podocarpus_hallii	0	0	0	0	0	0	0	0	1
Podocarpus_atjehensis	0	0	0	1	1	0	0	0	0
Podocarpus_angustifolius	1	0	0	0	0	0	0	0	0
Podocarpus_salignus	1	0	0	0	0	0	0	0	0
Podocarpus_capuronii	0	1	0	0	0	0	0	0	0
Podocarpus_madagascariensis_var_madagascariensis	0	1	0	0	0	0	0	0	0
Podocarpus_henkelii	0	1	0	0	0	0	0	0	0
Podocarpus_elongatus	0	1	0	0	0	0	0	0	0
Podocarpus_latifolius	0	1	0	0	0	0	0	0	0
Podocarpus_milanjianus	0	1	0	0	0	0	0	0	0
Podocarpus_parlatorei	1	0	0	0	0	0	0	0	0
Podocarpus_lambertii	1	0	0	0	0	0	0	0	0
Podocarpus_sprucei	1	0	0	0	0	0	0	0	0
Podocarpus_transiens	1	0	0	0	0	0	0	0	0
Podocarpus_hispaniolensis	1	0	0	0	0	0	0	0	0
Podocarpus_matudae	1	0	0	0	0	0	0	0	0
Podocarpus_matudae_var_reichei	1	0	0	0	0	0	0	0	0
Podocarpus_barretoi	1	0	0	0	0	0	0	0	0
Podocarpus_urbanii	1	0	0	0	0	0	0	0	0
Podocarpus_purdieanus	1	0	0	0	0	0	0	0	0
Podocarpus_aristulatus	1	0	0	0	0	0	0	0	0
Podocarpus_ekmanii	1	0	0	0	0	0	0	0	0
Podocarpus_magnifolius	1	0	0	0	0	0	0	0	0
Podocarpus_coriaceus	1	0	0	0	0	0	0	0	0
Podocarpus_trinitensis	1	0	0	0	0	0	0	0	0
Podocarpus_celatus	1	0	0	0	0	0	0	0	0
Podocarpus_oleifolius	1	0	0	0	0	0	0	0	0
Podocarpus_guatemalensis	1	0	0	0	0	0	0	0	0
Podocarpus_rusbyi	1	0	0	0	0	0	0	0	0
Podocarpus_sellowii	1	0	0	0	0	0	0	0	0
Podocarpus_brasiliensis	1	0	0	0	0	0	0	0	0
Podocarpus_costaricensis	1	0	0	0	0	0	0	0	0
Podocarpus_drouynianus	0	0	0	0	0	0	1	0	0
Podocarpus_glaucus	0	0	0	1	1	0	0	0	0
Podocarpus_spinulosus	0	0	0	0	0	0	1	0	0
Podocarpus_rumphii	0	0	1	1	1	0	0	0	0
Podocarpus_teysmannii	0	0	0	1	0	0	0	0	0
Podocarpus_longifoliolatus	0	0	0	0	0	0	0	1	0
Podocarpus_polyspermus	0	0	0	0	0	0	0	1	0
Podocarpus_elatus	0	0	0	0	0	0	1	0	0
Podocarpus_grayae	0	0	0	0	0	0	1	0	0
Podocarpus_decumbens	0	0	0	0	0	0	0	1	0
Podocarpus_beecherae	0	0	0	0	0	0	0	1	0
Podocarpus_novae_caledoniae	0	0	0	0	0	0	0	1	0
Podocarpus_lucienii	0	0	0	0	0	0	0	1	0
Podocarpus_sylvestris	0	0	0	0	0	0	0	1	0
Podocarpus_dispermus	0	0	0	0	0	0	1	0	0
Podocarpus_salomonensis	0	0	0	0	1	0	0	0	0
Podocarpus_brevifolius	0	0	0	1	0	0	0	0	0
Podocarpus_ledermannii	0	0	0	0	1	0	0	0	0
Podocarpus_archboldii	0	0	0	0	1	0	0	0	0
Podocarpus_brassii_var_brassii	0	0	0	0	1	0	0	0	0
Podocarpus_brassii_var_humilis	0	0	0	0	1	0	0	0	0
Podocarpus_spathoides	0	0	0	1	0	0	0	0	0
Podocarpus_crassigemmis	0	0	0	0	1	0	0	0	0
Podocarpus_rubens	0	0	0	1	1	0	0	0	0
Podocarpus_affinis	0	0	0	0	0	1	0	0	0
Podocarpus_decipiens	0	0	0	0	0	1	0	0	0
Podocarpus_pallidus	0	0	0	0	0	1	0	0	0
Podocarpus_degeneri	0	0	0	0	0	1	0	0	0
Podocarpus_insularis	0	0	0	0	1	0	0	0	0
Podocarpus_ramosii	0	0	0	1	0	0	0	0	0
Podocarpus_gibbsiae	0	0	0	1	0	0	0	0	0
Podocarpus_thailandensis	0	0	1	0	0	0	0	0	0
Podocarpus_deflexus	0	0	0	1	0	0	0	0	0
Podocarpus_assamica	0	0	1	0	0	0	0	0	0
Podocarpus_subtropicalis	0	0	1	0	0	0	0	0	0
Podocarpus_philippinensis	0	0	0	1	0	0	0	0	0
Podocarpus_polystachyus	0	0	1	1	1	0	0	0	0
Podocarpus_costalis	0	0	1	1	0	0	0	0	0
Podocarpus_pilgeri	0	0	1	1	1	0	0	0	0
Podocarpus_annamiensis	0	0	1	0	0	0	0	0	0
Podocarpus_pseudobracteatus	0	0	0	0	1	0	0	0	0
Podocarpus_fasciculus	0	0	1	0	0	0	0	0	0
Podocarpus_nakaii	0	0	1	0	0	0	0	0	0
Podocarpus_macrophyllus_var_maki	0	0	1	0	0	0	0	0	0
Podocarpus_macrophyllus_var_macrophyllus	0	0	1	0	0	0	0	0	0
Podocarpus_chingianus	0	0	1	0	0	0	0	0	0
Podocarpus_forrestii	0	0	1	0	0	0	0	0	0
"""

df = DataFrame(CSV.File(IOBuffer(tmpstr), delim=" ", ignorerepeated=true));
size(df)

end


function get_tr197()
	trstr = "(((Wollemia_nobilis:39.804869076837576,(Agathis_australis:32.13251618806202,((Agathis_lanceolata:6.278503190706726,(Agathis_ovata:5.311712544649567,(Agathis_moorei:4.171941865111998,(Agathis_corbassonii:2.1433691556215186,Agathis_montana:2.1433691556215186):2.0285727094904797):1.1397706795375688):0.9667906460571585):4.104304213330659,(Agathis_endertii:8.084991541042209,(((Agathis_macrophylla:1.6356521502720767,Agathis_silbae:1.6356521502720767):3.2377046799015643,(Agathis_vitiensis:2.8112612227663973,(Agathis_obtusa:0.9918591913604268,Agathis_robusta:0.991859191360427):1.8194020314059705):2.0620956074072434):2.582786680343774,((Agathis_borneensis:1.8836906690726878,Agathis_philippinensis:1.8836906690726878):2.948714442878276,(Agathis_atropurpurea:2.9391012341025213,(Agathis_microstachya:2.2043437440147313,(Agathis_dammara:1.66965184996365,(Agathis_alba:1.2289279625586755,Agathis_palmerstonii:1.2289279625586758):0.44072388740497415):0.5346918940510821):0.7347574900877896):1.8933038778484423):2.6237383985664504):0.6288480305247948):2.2978158629951757):21.749708784024634):7.672352888775556):14.566850939750992,(((Araucaria_angustifolia:6.666745371262853,Araucaria_araucana:6.666745371262853):9.349389887070462,(Araucaria_bidwillii:3.730936557063434,Araucaria_hunsteinii:3.730936557063434):12.28519870126988):24.237434494579738,((Araucaria_cunninghamii:4.598552536499169,(Araucaria_excelsa:1.3989005143682764,Araucaria_heterophylla:1.3989005143682767):3.1996520221308926):7.298386976571458,((Araucaria_columnaris:1.7739014085399278,(Araucaria_luxurians:0.6165857881680196,Araucaria_nemorosa:0.6165857881680195):1.1573156203719082):2.1141890876363383,((Araucaria_schmidii:1.584808753688203,(Araucaria_bernieri:0.8365902972877224,Araucaria_subulata:0.8365902972877223):0.7482184564004806):1.4878480306336304,((Araucaria_humboldtensis:0.7249708025583026,Araucaria_scopulorum:0.7249708025583026):1.9725557859306475,((Araucaria_laubenfelsii:1.4906771753887456,Araucaria_rulei:1.4906771753887456):1.08342370858036,(Araucaria_muelleri:1.407894071671472,(Araucaria_biramulata:0.8164272524114736,Araucaria_montana:0.8164272524114734):0.5914668192599988):1.1662068122976341):0.12342570451984436):0.37513019583288276):0.8154337118544328):8.008849016894361):28.356630239842424):14.118150263675517):130.94940932529911,((((Halocarpus_biformis:4.0294895143773735,(Halocarpus_bidwillii:1.8091540923352885,Halocarpus_kirkii:1.8091540923352887):2.220335422042085):84.67683893467317,((Lepidothamnus_fonkii:14.039734462064818,(Lepidothamnus_intermedius:6.0207777603501444,Lepidothamnus_laxifolius:6.0207777603501444):8.018956701714675):50.19706367688423,(Phyllocladus_alpinus:6.862712482934043,(Phyllocladus_trichomanoides_var_trichomanoides:5.063598488504128,(Phyllocladus_hypophyllus:3.480022310727679,(Phyllocladus_aspleniifolius:1.031447446744104,Phyllocladus_toatoa:1.031447446744104):2.448574863983575):1.583576177776448):1.7991139944299155):57.374085656015005):24.46953031010149):14.03565872623679,((Parasitaxus_usta:80.02579488014688,(Lagarostrobos_franklinii:60.57037341086034,Manoao_colensoi:60.57037341086034):19.455421469286527):13.144523134571827,((Prumnopitys_ladei:25.906388660724787,(Prumnopitys_ferruginea:1.84209811771863,Prumnopitys_ferruginoides:1.84209811771863):24.06429054300616):30.857362708817433,(Prumnopitys_amara:38.98315352605804,(Prumnopitys_taxifolia:15.343957729248567,(Prumnopitys_montana:9.544321018743553,(Prumnopitys_andina:2.096967423056695,Prumnopitys_exigua:2.096967423056695):7.447353595686858):5.799636710505014):23.639195796809474):17.78059784348418):36.406566645176476):9.571669160568632):27.502500427044097,(Saxegothaea_conspicua:111.60678006791831,(Microcachrys_tetragona:101.42757861730149,((Pherosphaera_fitzgeraldii:4.590771198025255,Pherosphaera_hookeriana:4.590771198025255):87.69287594542895,((Acmopyle_pancheri:26.923489035931894,Acmopyle_sahniana:26.923489035931894):55.482655557235695,(((Dacrycarpus_vieillardii:8.143100848481694,(Dacrycarpus_dacrydioides:5.645563935670678,((Dacrycarpus_cinctus:1.414112215247252,Dacrycarpus_kinabaluensis:1.414112215247252):1.9507920125022162,(Dacrycarpus_compactus:0.4191268941396414,Dacrycarpus_expansus:0.4191268941396414):2.945777333609827):2.2806597079212105):2.4975369128110145):46.66190922017411,((Falcatifolium_taxoides:15.096334862001745,(Falcatifolium_papuanum:12.051454592248692,(Falcatifolium_falciforme:2.1510933011045257,Falcatifolium_gruezoi:2.1510933011045252):9.900361291144169):3.044880269753053):9.548661564996502,((Dacrydium_elatum:9.548977974561561,(Dacrydium_comosum:4.813244318448209,(Dacrydium_magnum:2.1804241752231373,(Dacrydium_gracile:0.661456712904113,Dacrydium_xanthandrum:0.6614567129041129):1.5189674623190248):2.6328201432250715):4.735733656113351):3.171839635866119,(((Dacrydium_beccarii:1.3092078035637533,Dacrydium_gibbsiae:1.3092078035637533):4.8941531112822165,(Dacrydium_cupressinum:1.9744537903902224,Dacrydium_pectinatum:1.9744537903902226):4.228907124455747):1.7911232706705622,(Dacrydium_guillauminii:3.1911309103540737,(Dacrydium_nidulum:2.502697699291244,((Dacrydium_lycopodioides:0.6817736157305194,Dacrydium_x_suprinii:0.6817736157305195):1.594005959308822,(Dacrydium_nausoriense:1.8332945704495787,(Dacrydium_araucarioides:1.248501291475857,Dacrydium_balansae:1.248501291475857):0.5847932789737214):0.44248500458976303):0.22691812425190339):0.6884332110628284):4.803353275162459):4.7263334249111475):11.924178816570569):30.160013641657557):9.254668501527583,(((Retrophyllum_rospigliosii:11.017967635142545,(Retrophyllum_vitiense:7.957284981001082,(Retrophyllum_comptonii:5.448017018359398,Retrophyllum_minus:5.448017018359398):2.509267962641684):3.0606826541414645):17.552814344047064,((Afrocarpus_gaussenii:3.4505349761560846,(Afrocarpus_mannii:1.4804853358520997,(Afrocarpus_gracilior:1.1694899781419938,(Afrocarpus_falcatus:0.7798104003081242,(Afrocarpus_dawei:0.212686568403859,Afrocarpus_usambarensis:0.21268656840385905):0.5671238319042654):0.38967957783386986):0.3109953577101059):1.9700496403039849):13.35313466866873,((Nageia_motleyi:1.5233215289274986,Nageia_wallichiana:1.5233215289274986):5.683086333675103,(Nageia_fleuryi:4.840744888319199,(Nageia_formosensis:2.8885037492906194,(Nageia_nagi:1.1915682042127238,Nageia_nankoensis:1.1915682042127236):1.6969355450778956):1.9522411390285788):2.3656629742834037):9.597261782222214):11.767112334364793):11.243189792939841,(((Podocarpus_smithii:14.835578963442016,(((Podocarpus_acutifolius:1.4477871107944957,Podocarpus_totara:1.4477871107944957):0.6496827460277856,(Podocarpus_cunninghamii:0.42926111421608115,Podocarpus_hallii:0.42926111421608115):1.6682087426062004):2.1475790139201076,((Podocarpus_alpinus:0.4781174700864218,Podocarpus_nivalis:0.4781174700864218):2.465167437767522,(Podocarpus_gnidioides:0.8527006147995211,Podocarpus_lawrencei:0.8527006147995211):2.0905842930544223):1.3017639628884456):10.590530092699627):2.704261839545893,(((Podocarpus_atjehensis:11.115039408045929,(Podocarpus_angustifolius:1.4830592450284013,Podocarpus_salignus:1.4830592450284013):9.631980163017527):3.016847486204359,(Podocarpus_capuronii:7.770956900100764,(Podocarpus_madagascariensis_var_madagascariensis:4.966147230953731,(Podocarpus_henkelii:3.1666334519972152,(Podocarpus_elongatus:1.3384292541229066,(Podocarpus_latifolius:0.3750559068177822,Podocarpus_milanjianus:0.3750559068177822):0.9633733473051243):1.8282041978743084):1.7995137789565159):2.8048096691470326):6.360929994149523):0.7835924475116389,((Podocarpus_parlatorei:2.856474005980407,(Podocarpus_lambertii:1.9986493608359106,(Podocarpus_sprucei:0.8578094217056789,Podocarpus_transiens:0.8578094217056784):1.140839939130232):0.8578246451444962):8.521688353631982,(((Podocarpus_urbanii:2.664947037098949,(Podocarpus_purdieanus:1.4739400206020359,(Podocarpus_aristulatus:0.7134022441852732,Podocarpus_ekmanii:0.7134022441852735):0.7605377764167625):1.1910070164969138):2.946668385442674,(Podocarpus_magnifolius:3.33263920467052,(Podocarpus_sellowii:3.067323031103573,(Podocarpus_costaricensis:2.053794203504945,(Podocarpus_coriaceus:0.3602320911220435,Podocarpus_trinitensis:0.3602320911220435):1.6935621123829003):1.0135288275986287):0.26531617356694737):2.278976217871103):0.0033759548746914447,((Podocarpus_celatus:0.7196231019318928,Podocarpus_oleifolius:0.7196231019318928):4.7959065416720605,((Podocarpus_rusbyi:3.630716994902208,(Podocarpus_matudae:1.2567372883830905,Podocarpus_matudae_var_reichei:1.2567372883830905):2.3739797065191173):0.88241493991484,((Podocarpus_barretoi:1.1843727299569595,Podocarpus_guatemalensis:1.1843727299569597):2.2408715864321773,(Podocarpus_brasiliensis:2.5339637430965185,Podocarpus_hispaniolensis:2.5339637430965194):0.8912805732926183):1.08788761842791):1.0023977087869058):0.09946173381236179):5.7631709821960735):3.5373169821495374):2.6243614612259822):7.950599873921014,((Podocarpus_drouynianus:7.632306168370535,(Podocarpus_glaucus:2.5429429229470597,Podocarpus_spinulosus:2.5429429229470597):5.089363245423476):3.605910910663421,(((Podocarpus_rumphii:1.0857386068370574,Podocarpus_teysmannii:1.0857386068370574):6.139444821831144,(((Podocarpus_elatus:2.4507209240045307,Podocarpus_grayae:2.4507209240045307):2.5526642562762234,(Podocarpus_longifoliolatus:4.198400970726024,Podocarpus_polyspermus:4.198400970726024):0.8049842095547293):1.0787419672303091,(Podocarpus_decumbens:2.7991956693996807,(Podocarpus_lucienii:2.4385546221564045,(Podocarpus_sylvestris:1.7024698872677737,(Podocarpus_beecherae:0.47678827647144295,Podocarpus_novae_caledoniae:0.47678827647144295):1.2256816107963306):0.736084734888631):0.36064104724327617):3.282931478111382):1.143056281157138):1.8274086022951996,(Podocarpus_dispermus:7.816256070726183,(((Podocarpus_ramosii:2.2964722901958665,((Podocarpus_degeneri:0.7927113827072496,Podocarpus_insularis:0.7927113827072495):1.456766346998542,(Podocarpus_decipiens:1.4888814830567987,(Podocarpus_affinis:0.937915052463093,Podocarpus_pallidus:0.9379150524630928):0.5509664305937055):0.7605962466489928):0.04699456049007544):1.6641066509734164,((Podocarpus_crassigemmis:1.0883299535891224,Podocarpus_rubens:1.0883299535891224):2.0676977608517495,(Podocarpus_ledermannii:2.4951716123338,(Podocarpus_archboldii:0.9150969890222431,(Podocarpus_brassii_var_brassii:0.37542706726578773,Podocarpus_brassii_var_humilis:0.37542706726578784):0.5396699217564553):1.5800746233115572):0.6608561021070711):0.8045512267284116):1.7380991901378438,(((Podocarpus_annamiensis:0.9755549789981713,Podocarpus_pseudobracteatus:0.9755549789981713):2.2688364196335504,(Podocarpus_pilgeri:2.381429825282206,((Podocarpus_fasciculus:0.4923587675910192,Podocarpus_nakaii:0.4923587675910192):1.4925586687900965,(Podocarpus_costalis:1.6602063975980985,(Podocarpus_macrophyllus_var_maki:1.0980358104941321,(Podocarpus_forrestii:0.7767637011911765,(Podocarpus_chingianus:0.6752880280250302,Podocarpus_macrophyllus_var_macrophyllus:0.6752880280250297):0.10147567316614647):0.3212721093029556):0.5621705871039664):0.324711038783017):0.39651238890109086):0.8629615733495157):2.1075578643256785,(Podocarpus_brevifolius:4.741207160586413,(Podocarpus_salomonensis:3.5238674713644174,(Podocarpus_thailandensis:2.3676403215517974,((Podocarpus_gibbsiae:0.8909889136376613,Podocarpus_spathoides:0.8909889136376613):0.8718950640640806,(Podocarpus_philippinensis:0.827631296592605,(Podocarpus_subtropicalis:0.636792884387208,(Podocarpus_assamica:0.43011288729918074,(Podocarpus_deflexus:0.2929552506039599,Podocarpus_polystachyus:0.29295525060396005):0.13715763669522074):0.2066799970880268):0.19083841220539743):0.9352526811091367):0.6047563438500554):1.15622714981262):1.2173396892219954):0.6107421023709882):0.34672886834972605):2.117577939419056):1.2363359602372173):2.1856250480705555):14.252223597874966):14.323531095220527):24.245706798053938):18.346466022984202):9.877502550286621):9.143931473847275):10.179201450616816):18.637707534413124):55.076641739556266);"
	
	tr = readTopology(trstr)
	return tr
end

