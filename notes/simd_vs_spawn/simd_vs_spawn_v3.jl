
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Statistics 			# for mean(), max()
using DataFrames  # for e.g. DataFrame()
using Dates						# for e.g. DateTime, Dates.now()
using DifferentialEquations # for ODEProblem
using BenchmarkTools	# for @benchmark

"""
# start Julia with multiple cores using:
JULIA_NUM_THREADS=10 julia --startup-file=yes
"""


# Check that you have multiple threads
numthreads = Base.Threads.nthreads()

# Load the model structure/rates (all precalculated for speed)
include("/GitHub/PhyBEARS.jl/test/model_p_object.jl")
p_Es_v5 = load_ps_127();


# Load the ODE functions
include("/GitHub/PhyBEARS.jl/test/simd_vs_spawn_setup_v2.jl")

# Set up output object
numstates = 127
number_of_solves = 10
solve_results1 = collect(repeat([collect(repeat([0.0], numstates))], number_of_solves));
solve_results2 = collect(repeat([collect(repeat([0.0], numstates))], number_of_solves));
length(solve_results1)
length(solve_results1[1])
sum(sum.(solve_results1))


# Precalculate the Es for use in the Ds
Es_tspan = (0.0, 60.0)
prob_Es_v7 = DifferentialEquations.ODEProblem(Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
sol_Es_v7 = solve(prob_Es_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, 
abstol=1e-12, reltol=1e-9);

p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v7);


# Set up ODE inputs
u = collect(repeat([0.0], numstates));
u[2] = 1.0
du = similar(u)
du .= 0.0
p = p_Ds_v7;
t = 1.0

# ODE functions to integrate (single-step; ODE solvers will run this many many times)
@time Ds_v5_tmp(du,u,p,t)
@time Ds_v5_tmp(du,u,p,t)
@time Ds_v7_simd_sums(du,u,p,t)
@time Ds_v7_simd_sums(du,u,p,t)

#@btime Ds_v5_tmp(du,u,p,t)
# 7.819 ms (15847 allocations: 1.09 MiB)

#@btime Ds_v7_simd_sums(du,u,p,t)
# 155.858 μs (3075 allocations: 68.66 KiB)



include("/GitHub/PhyBEARS.jl/test/simd_vs_spawn_setup_v2.jl")

tspan = (0.0, 1.0)
prob_Ds_v7 = DifferentialEquations.ODEProblem(Ds_v7_simd_sums, p_Ds_v7.uE, tspan, p_Ds_v7);

sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol=1e-12, reltol=1e-9);

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

@time core_op_plain(u, tspan, p_Ds_v7);
@time core_op_plain(u, tspan, p_Ds_v7);
@time core_op_simd(u, tspan, p_Ds_v7);
@time core_op_simd(u, tspan, p_Ds_v7);


function serial_with_plain_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=10)
	start_time = Dates.now()
	tspan = (0.0, 1.0)
	for i in 1:number_of_solves
		u .= 0.0
		u[i] = 1.0

		sol_Ds_v7 = core_op_plain(u, tspan, p_Ds_v7)
		solve_results1[i] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
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
	tspan = (0.0, 1.0)
	for i in 1:number_of_solves
		u .= 0.0
		u[i] = 1.0

		sol_Ds_v7 = core_op_simd(u, tspan, p_Ds_v7)
		solve_results1[i] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
	#	print("\n")
	#	print(round.(sol_Ds_v7[length(sol_Ds_v7)], digits=3))
	end
	
	end_time = Dates.now()
	duration = (end_time - start_time).value / 1000.0
	sum_of_solutions = sum(sum.(solve_results1))
	return (duration, sum_of_solutions)
end

# Output is (runtime, sum_of_solutions)
serial_with_plain_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
serial_with_plain_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
# (1.136, 8.129628626179963)

serial_with_simd_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
serial_with_simd_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
# (0.065, 8.129628626179963)





function parallel_with_plain_v5(tspan, p_Ds_v7, solve_results2; number_of_solves=10)
	start_time = Dates.now()
	tspan = (0.0, 1.0)
	
	# Individual ODE solutions will occur over different timeperiods,
	# initial values, and parameters.  We'd just like to load up the 
	# cores for the first jobs in the list, then add jobs as earlier
	# jobs finish.
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	task_numbers = Any[]
	task_inc = 0
	are_we_done = false
	
	# List the tasks
	for i in 1:number_of_solves
		u .= 0.0
		u[i] = 1.0

		task_inc = task_inc + 1
		# THE ONLY DIFFERENCE IS HERE
		push!(tasks, Base.Threads.@spawn core_op_plain(u, tspan, p_Ds_v7))
		push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		push!(task_numbers, task_inc)
	end

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1

		num_tasks = length(tasks)
		for j in 1:num_tasks
			if (tasks_fetched_TF[j] == false)
				if (istaskstarted(tasks[j]) == true) && (istaskdone(tasks[j]) == true)
					sol_Ds_v7 = fetch(tasks[j])
					solve_results2[task_numbers[j]] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
					print("\nFinished task j=")
					print(j)
					tasks_fetched_TF[j] = true
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
	return (duration, sum_of_solutions)
end


function parallel_with_simd_v7(tspan, p_Ds_v7, solve_results2; number_of_solves=10)
	start_time = Dates.now()
	tspan = (0.0, 1.0)
	
	# Individual ODE solutions will occur over different timeperiods,
	# initial values, and parameters.  We'd just like to load up the 
	# cores for the first jobs in the list, then add jobs as earlier
	# jobs finish.
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	task_numbers = Any[]
	task_inc = 0
	are_we_done = false
	
	# List the tasks
	for i in 1:number_of_solves
		u .= 0.0
		u[i] = 1.0

		task_inc = task_inc + 1
		# THE ONLY DIFFERENCE IS HERE
		push!(tasks, Base.Threads.@spawn core_op_simd(u, tspan, p_Ds_v7))
		push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		push!(task_numbers, task_inc)
	end

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1

		num_tasks = length(tasks)
		for j in 1:num_tasks
			if (tasks_fetched_TF[j] == false)
				if (istaskstarted(tasks[j]) == true) && (istaskdone(tasks[j]) == true)
					sol_Ds_v7 = fetch(tasks[j])
					solve_results2[task_numbers[j]] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
					print("\nFinished task j=")
					print(j)
					tasks_fetched_TF[j] = true
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
	return (duration, sum_of_solutions)
end


parallel_with_plain_v5(tspan, p_Ds_v7, solve_results2; number_of_solves=number_of_solves)
# Faster!
# (runtime, sum_of_solutions)
# (0.193, 8.129628626179963)


parallel_with_plain_v5(tspan, p_Ds_v7, solve_results2; number_of_solves=number_of_solves)
# Hangs sometimes, on the 2nd go!

parallel_with_simd_v7(tspan, p_Ds_v7, solve_results2; number_of_solves=number_of_solves)
# This just HANGS



