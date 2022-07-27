
versioninfo()
notes="""
My setup is:
Julia Version 1.7.3
Commit 742b9abb4d (2022-05-06 12:58 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin21.4.0)
  CPU: Intel(R) Xeon(R) CPU E5-2697 v2 @ 2.70GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, ivybridge)
"""


# Startup notes
notes="""
# "If $JULIA_NUM_THREADS is set to auto, then the number of threads will be set to the number of CPU threads."
JULIA_NUM_THREADS=8 julia --startup-file=no
Threads.nthreads(): 8 # Number of CPU threads
"""
using Distributed

function Rnames(obj)
	flat2(fieldnames(typeof(obj)))
end


# And multiple cores/workers to match
t = @async Distributed.addprocs(7)
Distributed.nprocs()
# Check that you have same number of processors and threads
Distributed.nprocs()
numthreads = Base.Threads.nthreads()

using Dates						# for e.g. DateTime, Dates.now()
using DifferentialEquations # for ODEProblem
using Sundials				# for CVODE_BDF
using SciMLBase

Distributed.workers()
# Check that you have same number of processors and threads
Distributed.nprocs()
numthreads = Base.Threads.nthreads()

@everywhere using Distributed
@everywhere using Dates						# for e.g. DateTime, Dates.now()
@everywhere using DifferentialEquations # for ODEProblem
@everywhere using Sundials				# for CVODE_BDF
@everywhere using SciMLBase


# Download & include the pre-saved model structure/rates (all precalculated for speed; 1.8 MB)
#include("/GitHub/PhyBEARS.jl/test/model_p_object.jl")
url = "https://gist.githubusercontent.com/nmatzke/ed99ab8f5047794eb25e1fdbd5c43b37/raw/b3e6ddff784bd3521d089642092ba1e3830699c0/model_p_object.jl"
download(url,  "model_p_object.jl")
@everywhere include("model_p_object.jl")

# Load the ODE functions
#url = "https://gist.githubusercontent.com/nmatzke/f116258c78bd43ab7a448f07c4290516/raw/24a210261fd2e090b8ed27bc64a59a1ff9ec62cd/simd_vs_spawn_setup_v2.jl"
#download(url,  "simd_vs_spawn_setup_v2.jl")
#@everywhere include("simd_vs_spawn_setup_v2.jl")

#include("/GitHub/PhyBEARS.jl/test/simd_vs_spawn_setup_v2.jl")
@everywhere include("/GitHub/PhyBEARS.jl/test/simd_vs_spawn_setup_v4.jl")

# Load the pre-saved model structure/rates (all precalculated for speed; 1.8 MB)
@everywhere p_Es_v5 = load_ps_127();



# Set up output object
numstates = 127
number_of_solves = 10

solve_results1 = Array{Float64, 2}(undef, number_of_solves, numstates)
solve_results1 .= 0.0
solve_results2 = Array{Float64, 2}(undef, number_of_solves, numstates)
solve_results2 .= 0.0
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
du = similar(u);
du .= 0.0;
p = p_Ds_v7;
t = 1.0;

# ODE functions to integrate (single-step; ODE solvers will run this many many times)
@time Ds_v5_tmp(du,u,p,t)
@time Ds_v5_tmp(du,u,p,t)
@time Ds_v7_simd_sums(du,u,p,t)
@time Ds_v7_simd_sums(du,u,p,t)

#@btime Ds_v5_tmp(du,u,p,t)
# 7.819 ms (15847 allocations: 1.09 MiB)

#@btime Ds_v7_simd_sums(du,u,p,t)
# 155.858 μs (3075 allocations: 68.66 KiB)



tspan = (0.0, 1.0)
prob_Ds_v7 = DifferentialEquations.ODEProblem(Ds_v7_simd_sums, p_Ds_v7.uE, tspan, p_Ds_v7);

sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol=1e-12, reltol=1e-9);


@time core_op_plain(u, tspan, p_Ds_v7);
@time core_op_plain(u, tspan, p_Ds_v7);
@time core_op_simd(u, tspan, p_Ds_v7);
@time core_op_simd(u, tspan, p_Ds_v7);


# Output is (runtime, sum_of_solutions)
serial_with_plain_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
serial_with_plain_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
serial_with_plain_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
# (duration, sum_of_solutions)
# (1.1, 8.731365050398926)
# (0.878, 8.731365050398926)
# (0.898, 8.731365050398926)

serial_with_simd_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
serial_with_simd_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
serial_with_simd_v7(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
# (duration, sum_of_solutions)
# (0.046, 8.731365050398928)
# (0.042, 8.731365050398928)
# (0.046, 8.731365050398928)

@everywhere include("/GitHub/PhyBEARS.jl/test/simd_vs_spawn_setup_v4.jl")

j=2
new_worker = j
solve_results2[j,j] = 1.0
x=Distributed.remotecall(core_op_plain, new_worker, solve_results2[j,:], tspan, p_Ds_v7)
@async fetch(x)

function istaskdone2(x)
	doneTF = x.v != nothing;
	return doneTF
end



parallel_with_plain_v5(tspan, p_Ds_v7, solve_results2; number_of_solves=10)
x = Dates.now()



function parallel_with_simd_v7(tspan, p_Ds_v7, solve_results2; number_of_solves=10)
	start_time = Dates.now()
	number_of_threads = Base.Threads.nthreads()
	curr_numthreads = Base.Threads.nthreads()
		
	# Individual ODE solutions will occur over different timeperiods,
	# initial values, and parameters.  We'd just like to load up the 
	# cores for the first jobs in the list, then add jobs as earlier
	# jobs finish.
	tasks = Any[]
	tasks_started_TF = Bool[]
	tasks_fetched_TF = Bool[]
	task_numbers = Any[]
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
		push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		push!(task_numbers, task_inc)
	end
	
	# Total number of tasks
	num_tasks = length(tasks_fetched_TF)

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		
		# Launch tasks when thread (core) is available
		for j in 1:num_tasks
			if (tasks_fetched_TF[j] == false)
				if (tasks_started_TF[j] == false) && (curr_numthreads > 0)
					# Start a task
					push!(tasks, Base.Threads.@spawn core_op_simd(solve_results2[j,:], tspan, p_Ds_v7))
					curr_numthreads = curr_numthreads-1;
					tasks_started_TF[j] = true;
					push!(current_running_tasks, task_numbers[j])
				end
			end
		end
		
		# Check for finished tasks
		tasks_to_check_TF = ((tasks_started_TF.==true) .+ (tasks_fetched_TF.==false)).==2
		if sum(tasks_to_check_TF .== true) > 0
			for k in 1:sum(tasks_to_check_TF)
				if (tasks_fetched_TF[current_running_tasks[k]] == false)
					if (istaskstarted(tasks[k]) == true) && (istaskdone(tasks[k]) == true)
						sol_Ds_v7 = fetch(tasks[k]);
						solve_results2[current_running_tasks[k],:] .= sol_Ds_v7.u[length(sol_Ds_v7.u)].+0.0
						tasks_fetched_TF[current_running_tasks[k]] = true
						current_tasknum = current_running_tasks[k]
						deleteat!(tasks, k)
						deleteat!(current_running_tasks, k)
						curr_numthreads = curr_numthreads+1;
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

tspan = (0.0, 1.0)
parallel_with_plain_v5(tspan, p_Ds_v7, solve_results2; number_of_solves=number_of_solves)
# Faster than serial plain version
# (duration, sum_of_solutions)
# (0.351, 8.731365050398926)
# (0.343, 8.731365050398926)
# (0.366, 8.731365050398926)

parallel_with_simd_v7(tspan, p_Ds_v7, solve_results2; number_of_solves=number_of_solves)
# Dramatically slower than serial simd version
# (duration, sum_of_solutions)
# (136.966, 9.61313614002137)
# (141.843, 9.616688089683372)


