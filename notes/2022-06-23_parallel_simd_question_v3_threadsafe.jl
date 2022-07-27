# Problem setup:
# "p_Es_v5" contains a large list of pre-allocated indices and rates for the model. 
# These have to be summed over all indices i,j,k. 
# This means that the inner-loop functions can be optimized with @inbounds @simd


using Dates						# for e.g. DateTime, Dates.now()
using DifferentialEquations # for ODEProblem
using Sundials				# for CVODE_BDF

# Check that you have multiple threads
numthreads = Base.Threads.nthreads()

# Download & include() the pre-saved model structure/rates (all precalculated for speed; 1.8 MB)
#include("/GitHub/PhyBEARS.jl/test/model_p_object.jl")
url = "https://gist.githubusercontent.com/nmatzke/ed99ab8f5047794eb25e1fdbd5c43b37/raw/b3e6ddff784bd3521d089642092ba1e3830699c0/model_p_object.jl"
download(url,  "model_p_object.jl")
include("model_p_object.jl")
p_Es_v5 = load_ps_127();


# ODE to solve many times; @simd-enhanced version; 10+ times faster on a single core
ODE_to_solve_Ds_plain = (du,u,p,t) -> begin
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
	# Pre-calculated solution of the Es
	uE = collect(repeat([0.0], 127)) # zero out the uE for now
	
  @inbounds for i in 1:n
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]
		du[i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 
			(sum(Cijk_vals[Ci_eq_i] .*                                               
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
	end
end

# This is the core operation to solve "Ds"; plain version (no @simd)
function core_op_plain(u, tspan, p_Ds_v7)
	prob_Ds_v5 = DifferentialEquations.ODEProblem(ODE_to_solve_Ds_plain, u.+0.0, tspan, p_Ds_v7);
	sol_Ds_v5 = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol=1e-12, reltol=1e-9);
	return sol_Ds_v5
end

# Do 8 solves in serial
function serial_solves_plain(tspan, p_Ds_v7, solve_results1; number_of_solves=8)
	start_time = Dates.now()
	for i in 1:number_of_solves
		solve_results1[i,:] .= 0.0
		solve_results1[i,i] = 1.0
		sol_Ds_v7 = core_op_plain(solve_results1[i,:], tspan, p_Ds_v7)
		solve_results1[i,:] .= 	sol_Ds_v7.u[length(sol_Ds_v7.u)]
	end
	duration = (Dates.now() - start_time).value / 1000.0
	sum_of_solutions = sum(sum.(solve_results1))
	return (duration, sum_of_solutions)
end

# Do 8 solves in parallel
function parallel_solves_plain(tspan, p_Ds_v7, solve_results2; number_of_solves=8)
	print("\nSolved task #: ")
	start_time = Dates.now()
	# Make tasks, then run them
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	for i in 1:number_of_solves
		solve_results2[i,i] = 1.0 
		# USING "PLAIN", NON-SIMD OPERATION HERE
		push!(tasks, Base.Threads.@spawn core_op_plain(solve_results2[i,:].+0.0, tspan, p_Ds_v7));
		push!(tasks_fetched_TF, false)
	end

	are_we_done = false;
	done_count = 0;
	while(are_we_done == false)
		for k in 1:number_of_solves
			if (istaskstarted(tasks[k]) == true) && (istaskdone(tasks[k]) == true) && (tasks_fetched_TF[k] == false)
				sol_Ds_v7 = fetch(tasks[k]);
				solve_results2[k,:] .= sol_Ds_v7.u[length(sol_Ds_v7.u)].+0.0
				done_count = done_count + 1;
				tasks_fetched_TF[k] = true
				print(k)
				print(" ")
				break
			end
		end
		if (done_count >= number_of_solves)
			are_we_done = true
			are_we_done && break
		end
	end
	duration = (Dates.now() - start_time).value / 1000;
	sum_of_solutions = sum(sum.(solve_results2));
	print("\n")
	return (duration, sum_of_solutions)
end

# Set up input; output objects
numstates = 127
number_of_solves = 8

solve_results1 = Array{Float64, 2}(undef, number_of_solves, numstates);
solve_results1 .= 0.0;
solve_results2 = Array{Float64, 2}(undef, number_of_solves, numstates);
solve_results2 .= 0.0;
size(solve_results1)

# Set up ODE Ds inputs
tspan = (0.0, 1.0)
p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms);
p = p_Ds_v7;
tspan = (0.0, 1.0)

# Single solve operations for Ds; but solve() will be run many times
@time core_op_plain(u, tspan, p_Ds_v7);
@time core_op_plain(u, tspan, p_Ds_v7);

# Multiple ODE solves, serial versions
serial_solves_plain(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
serial_solves_plain(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
serial_solves_plain(tspan, p_Ds_v7, solve_results1; number_of_solves=number_of_solves)
#	(1.042, 7.048516354666927)
# (1.04, 7.048516354666927)
# (1.028, 7.048516354666927)

# Multiple ODE solves, parallel versions; answers differ
parallel_solves_plain(tspan, p_Ds_v7, solve_results2; number_of_solves=8)
parallel_solves_plain(tspan, p_Ds_v7, solve_results2; number_of_solves=8)
# Solved task #: 8 1 6 2 5 4 7 3 
# (0.331, 7.048516354666927)
# Solved task #: 8 1 6 2 5 4 7 3 
# (0.215, 7.680340772075006)


# Usually hangs on the 3rd try...
parallel_solves_plain(tspan, p_Ds_v7, solve_results2; number_of_solves=8)
# Hangs on:
# Solved task #: 




