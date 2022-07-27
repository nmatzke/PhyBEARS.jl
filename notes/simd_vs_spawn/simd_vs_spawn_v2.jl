
# Check that you have multiple threads
numthreads = Base.Threads.nthreads()

using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using Statistics 			# for mean(), max()
using PhyBEARS.TreeTable # for prt()
using PhyBEARS.TrUtils # for flat2() (similar to unlist)
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.ModelLikes
using PhyBEARS.Flow
using PhyBEARS.Parsers
using PhyBEARS.Gmaps
using PhyBEARS.Optimizers

using Profile					# for @profile
using BenchmarkTools	# for @benchmark
#using PProf						# for pprof()


# Pre-allocated transition rates for a complex model with 511 states
p_Es_v5 = load_ps_511()



parallel_down_tree


# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(Ds_v7_simd_sums_tmp, p_Es_v5.uE, Es_tspan, p_Es_v5);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol);

p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v7);



# Set up output object
global numstates = n
number_of_solves = 100
solve_results1 = collect(repeat([collect(repeat([0.0], numstates))], number_of_solves))
solve_results2 = collect(repeat([collect(repeat([0.0], numstates))], number_of_solves))
length(solve_results1)
length(solve_results1[1])
sum(sum.(solve_results1))


# Set up ODE inputs
u = collect(repeat([0.0], numstates))
u[2] = 1.0
du = similar(u)
du .= 0.0
p = p_Ds_v7;
t = 1.0

# ODE functions to integrate
@time parameterized_ClaSSE_Ds_v7_simd_sums(du,u,p,t)

tspan = (0.0, 1.0)
prob_Ds_v7 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v7_simd_sums, p_Ds_v7.uE, tspan, p_Ds_v7);

sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);


function core_op(u, tspan, p_Ds_v7)
	prob_Ds_v7 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v7_simd_sums, u, tspan, p_Ds_v7);

	sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
	return sol_Ds_v7
end




tspan = (0.0, 1.0);
t_start = Dates.now();
for i in 1:20
	u .= 0.0
	u[i] = 1.0

	sol_Ds_v7 = core_op(u, tspan, p_Ds_v7)
	solve_results1[i] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
#	print("\n")
#	print(round.(sol_Ds_v7[length(sol_Ds_v7)], digits=3))
end
t_end = Dates.now();
sum(sum.(solve_results1))
t_end - t_start



tasks = Any[]
tasks_fetched_TF = Bool[]
are_we_done = false

iteration_number = 0
t_start = Dates.now()
while(are_we_done == false)
	iteration_number = iteration_number+1

	for i in 1:20
		u .= 0.0
		u[i] = 1.0

		sol_Ds_v7 = core_op(u, tspan, p_Ds_v7)	

		push!(tasks, Base.Threads.@spawn core_op(u, tspan, p_Ds_v7))
		push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
	end

	num_tasks = length(tasks)
	for i in 1:num_tasks
		if (tasks_fetched_TF[i] == false)
			if (istaskstarted(tasks[i]) == true) && (istaskdone(tasks[i]) == true)
				results = fetch(tasks[i])
				solve_results2[i] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
				tasks_fetched_TF[i] = true
			end
		end
	end

	are_we_done = sum(tasks_fetched_TF) == length(tasks_fetched_TF)
	# Test for concluding the while loop
	are_we_done && break
end # END while(are_we_done == false)
t_end = Dates.now()
sum(sum.(solve_results2))
t_end - t_start

