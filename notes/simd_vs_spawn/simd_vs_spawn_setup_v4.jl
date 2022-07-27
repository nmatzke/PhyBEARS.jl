# Plain version
Ds_v5_tmp = (du,u,p,t) -> begin
#	let p=p;

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
	
  @inbounds for i in 1:n
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		du[i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 	# case 2	
			(sum(Cijk_vals[Ci_eq_i] .*                                               
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
 # end

	end # END let p=p;
end


# @simd-enhanced version; 10 times faster on a single core
Ds_v7_simd_sums = (du,u,p,t) -> begin
#	let p=p;
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
#	end # END let p=p;
end

# @simd-enhanced version; 10 times faster on a single core
# Manually allocate terms here (minor GC cost)
Ds_v8_simd_sums = (du,u,p,t) -> begin
#	let p=p;
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
	uE = p.sol_Es_v5(t)
	terms = Vector{Float64}(undef, 4)
	
  for i in 1:n
		terms .= 0.0
		terms[1], terms[4] = sum_Cijk_rates_Ds_inbounds_simd_tmp(p.p_TFs.Cijk_rates_sub_i[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i])
	
		terms[2], terms[3] = sum_Qij_vals_inbounds_simd_tmp(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i])
		
		du[i] = -(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4]
  end
#	end # END let p=p;
end




# Inner-loop functions
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


# Only needs to be calculated once, then interpolator re-used by Ds_v7_simd_sums via p.sol_Es_v5(t)
Es_v7_simd_sums = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
	terms = p.terms
  for i in 1:n
		terms[1] = 0.0
		terms[2] = 0.0
		terms[3] = 0.0
		terms[4] = 0.0

		terms[1], terms[4] = sum_Cijk_rates_Es_inbounds_simd_tmp(p.p_TFs.Cijk_rates_sub_i[i], u, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i])
	
		terms[2], terms[3] = sum_Qij_vals_inbounds_simd_tmp(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i])
		
		du[i] = mu[i] -(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4]
  end
end






# This is the core operation; plain version (no @simd)
function core_op_plain(u, tspan, p_Ds_v7)
	prob_Ds_v5 = DifferentialEquations.ODEProblem(Ds_v5_tmp, u, tspan, p_Ds_v7);

	sol_Ds_v5 = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol=1e-12, reltol=1e-9);
	return sol_Ds_v5
end


# This is the core operation; @simd version
function core_op_simd(u, tspan, p_Ds_v7)
	prob_Ds_v7 = DifferentialEquations.ODEProblem(Ds_v7_simd_sums, u, tspan, p_Ds_v7);

	sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol=1e-12, reltol=1e-9);
	return sol_Ds_v7
end

function core_op_simd_v8(u, tspan, p_Ds_v7)
	prob_Ds_v7 = DifferentialEquations.ODEProblem(Ds_v8_simd_sums, u, tspan, p_Ds_v7);

	sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol=1e-12, reltol=1e-9);
	return sol_Ds_v7
end



function serial_with_plain_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=10)
	start_time = Dates.now()
	for i in 1:number_of_solves
		# Temporary u
		solve_results1[i,:] .= 0.0
		
		# Change the ith state from 0.0 to 1.0
		solve_results1[i,i] = 1.0
		solve_results1

		sol_Ds_v7 = core_op_plain(solve_results1[i,:], tspan, p_Ds_v7)
		solve_results1[i,:] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
	#	print("\n")
	#	print(round.(sol_Ds_v7[length(sol_Ds_v7)], digits=3))
	end
	
	end_time = Dates.now()
	duration = (end_time - start_time).value / 1000.0
	sum_of_solutions = sum(sum.(solve_results1))
	return (duration, sum_of_solutions)
end


function serial_with_simd_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=10)
	start_time = Dates.now()
	for i in 1:number_of_solves
		# Temporary u
		solve_results1[i,:] .= 0.0
		
		# Change the ith state from 0.0 to 1.0
		solve_results1[i,i] = 1.0
		solve_results1

		sol_Ds_v7 = core_op_simd(solve_results1[i,:], tspan, p_Ds_v7)
		solve_results1[i,:] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
	#	print("\n")
	#	print(round.(sol_Ds_v7[length(sol_Ds_v7)], digits=3))
	end
	
	end_time = Dates.now()
	duration = (end_time - start_time).value / 1000.0
	sum_of_solutions = sum(sum.(solve_results1))
	return (duration, sum_of_solutions)
end

function serial_with_simd_v8(tspan, p_Ds_v7, solve_results1; number_of_solves=10)
	start_time = Dates.now()
	for i in 1:number_of_solves
		# Temporary u
		solve_results1[i,:] .= 0.0
		
		# Change the ith state from 0.0 to 1.0
		solve_results1[i,i] = 1.0
		solve_results1

		sol_Ds_v7 = core_op_simd_v8(solve_results1[i,:], tspan, p_Ds_v7)
		solve_results1[i,:] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
	#	print("\n")
	#	print(round.(sol_Ds_v7[length(sol_Ds_v7)], digits=3))
	end
	
	end_time = Dates.now()
	duration = (end_time - start_time).value / 1000.0
	sum_of_solutions = sum(sum.(solve_results1))
	return (duration, sum_of_solutions)
end

#######################################################
# This version STALLS unless many slow print/wait statements
# are inserted. DO NOT USE!
#######################################################
function parallel_with_plain_v5(tspan, p_Ds_v7, solve_results2; number_of_solves=10)
	# Check that threads and workers/cores are available!
	numthreads = Base.Threads.nthreads()
	num_workers = length(Distributed.workers())
	if (numthreads == 1)
		txt = join(["STOP ERROR in parallel_with_plain_v5(): Threads.nthreads() must be higher than 1. Re-start Julia with e.g. 'JULIA_NUM_THREADS=8 julia' to fix this. Also make sure to use Distributed.addprocs to add N-1 processors."])
		print("\n")
		print(txt)
		print("\n")
		error(txt)
	end
	if (num_workers > (numthreads-1))
		txt = join(["STOP ERROR in parallel_with_plain_v5(): You have Distributed.workers()=", string(num_workers), " workers/processors enabled with Distributed.addprocs(), so you should have had at least ", string(num_workers+1), " threads available. Instead, Threads.nthreads()=", string(Threads.nthreads()), ". Re-start Julia with e.g. JULIA_NUM_THREADS=", string(num_workers+1), " julia to fix this. Also make sure to use Distributed.addprocs to add N-1 processors."])
		print("\n")
		print(txt)
		print("\n")
		error(txt)
	end
	
	start_time = Dates.now()
	list_of_workers = Distributed.workers()
	curr_free_workers = Distributed.workers()
		
	# Individual ODE solutions will occur over different timeperiods,
	# initial values, and parameters.  We'd just like to load up the 
	# cores for the first jobs in the list, then add jobs as earlier
	# jobs finish.
	jobs = Any[]
	tasks = Any[]
	task_waits = Any[]
	tasks_started_TF = Bool[]
	tasks_fetched_TF = Bool[]
	task_numbers = Any[]
	worker_for_task = Any[]
	task_inc = 0
	are_we_done = false
	current_running_tasks = Any[]
	
	# List the tasks
	for i in 1:number_of_solves
		# Temporary u
		solve_results2[i,:] .= 0.0
		
		# Change the ith state from 0.0 to 1.0
		solve_results2[i,i] = 1.0

		task_inc = task_inc + 1
		push!(tasks_started_TF, false) # Add a "false" to tasks_started_TF
		push!(tasks_fetched_TF, false) # Add a "false" to tasks_started_TF
		push!(worker_for_task, 0)
		push!(task_numbers, task_inc)
	end
	
	# Total number of tasks
	num_tasks = length(tasks_fetched_TF)

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1

		print("\nHere2")		

		# Launch tasks when thread (core) is available
		for j in 1:number_of_solves
			if (tasks_fetched_TF[j] == false)
				if (tasks_started_TF[j] == false) && (length(curr_free_workers) > 0)
					# Start a task
					new_worker = curr_free_workers[1]
					worker_for_task[j] = new_worker
					push!(jobs, Distributed.remotecall(core_op_plain, new_worker, solve_results2[j,:].+0.0, tspan, p_Ds_v7));
					txt = join(["\nLoaded task #", task_numbers[j], " on worker #", new_worker, ". Task: "])
					sleep(0.2)
					print(txt)
					
					
					push!(tasks, @async fetch(jobs[length(jobs)]))
					
					#print(tasks[j])
					#schedule(tasks[length(tasks)])
					deleteat!(curr_free_workers, 1)
					tasks_started_TF[j] = true;
					push!(current_running_tasks, task_numbers[j])
				end
			end
		end
		
		#print("\n")
		# Try this as an updater
		#fetch.(tasks);
		print("\nHere3")		
		
		# Check for finished tasks
		tasks_to_check_TF = ((tasks_started_TF.==true) .+ (tasks_fetched_TF.==false)).==2
		if sum(tasks_to_check_TF .== true) > 0
			for k in 1:sum(tasks_to_check_TF)
				if (tasks_fetched_TF[current_running_tasks[k]] == false)
					if (tasks_started_TF[current_running_tasks[k]] == true) && (istaskdone(tasks[k]) == true)
						print("\nHere4")		
						sol_Ds_v7 = fetch(tasks[k]);
						print("\nHere5")		
						#sol_Ds_v7 = jobs[k].v;
						finalize(tasks[k]);
						finalize(jobs[k]);
						print("\nHere6")		
						solve_results2[current_running_tasks[k],:] .= sol_Ds_v7.u[length(sol_Ds_v7.u)].+0.0
						tasks_fetched_TF[current_running_tasks[k]] = true
						current_tasknum = current_running_tasks[k]
						finished_worker = worker_for_task[current_running_tasks[k]]
						push!(curr_free_workers, finished_worker)
						deleteat!(jobs, k)
						deleteat!(tasks, k)
						deleteat!(current_running_tasks, k)
						print("\nHere7")		
						print("\nFinished task #")
						print(current_tasknum)
						print(", current task k=")
						print(k)
						break # break out of this loop, since you have modified current_running_tasks
					end
				end
			end
		end

		are_we_done = sum(tasks_fetched_TF) == length(tasks_fetched_TF)
		# Test for concluding the while loop
		are_we_done && break
	end # END while(are_we_done == false)

	end_time = Dates.now()
	duration = (end_time - start_time).value / 1000.0
	sum_of_solutions = sum(sum.(solve_results2))
	print("\n")
	return (duration, sum_of_solutions)
end
