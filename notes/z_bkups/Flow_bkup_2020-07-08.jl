module Flow
print("\n\nStarting module 'Flow'...loading dependencies...\n")
using LinearAlgebra  # for mul! (matrix multiplication)
using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA
using DifferentialEquations
using Distributed
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloNetworks
#using Plots						# for plot
using DataFrames          # for DataFrame()
using PhyBEARS.TrUtils # for flat2() (similar to unlist)
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs

export parameterized_ClaSSE_As_v5, check_linearDynamics_of_As, calc_Gs_SSE_condnums!, calc_Gs_SSE, calc_Gs_SSE!, branchOp_ClaSSE_Gs_v5, iterative_downpass_nonparallel_FlowClaSSE_v5!, run_Gs


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
parameterized_ClaSSE_As_v5 = (t, A, p; max_condition_number=1e8, print_warnings=true) -> begin

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
		sigma1_of_A = opnorm(A,1)  # the 1-norm should be adequate here (fastest calc.)
		if (2*sigma1_of_A > log(max_condition_number))
			warning_txt = join(["WARNING in parameterized_ClaSSE_As_v5 at t=", string(t), ": 2*opnorm(A,1)>log(max_condition_number) (", string(round(sigma1_of_A; digits=2)), " > ", string(round(log(max_condition_number); digits=2)), ")\n"], "")
			display(warning_txt)
		end # END if (2*sigma1_of_A > log(max_condition_number))
	end # END if (print_warnings == true)
 	return(A)
end


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
		A_at_t = Flow.parameterized_ClaSSE_As_v5(t, A[:,:], p_Ds_v5)
		sigma1_of_A = opnorm(A,1)  # the 1-norm should be adequate here (fastest calc.)
		upper_bound_kappa_rates_A[i] = 2*sigma1_of_A
		if (upper_bound_kappa_rates_A[i] > log(max_condition_number))
			condbigTF[i] = true
		end
	end
	
	
	# Creat dataframe
	# condbigTF = Is the growth rathe of condition number too big, i.e. bigger than log(max_condition_number)
	kappa_Arates_df = DataFrames.DataFrame(tvals=tvals, ub_kappa_ratesA=upper_bound_kappa_rates_A, condbigTF=condbigTF)

	return(kappa_Arates_df)
end



# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
calc_Gs_SSE! = (dG, G, pG, t) -> begin
	# Have to use pG.A[:,:] to avoid overwriting A, which screws everything up!
	#A = parameterized_ClaSSE_As_v5(t, pG.A[:,:], pG.p_Ds_v5)
	#display(A)
	#dG = A * G
	# Have to use pG.A[:,:] to avoid overwriting A, which screws everything up!
	#mul!(dG, A(t), G(t))  # Equation 20
	mul!(dG, parameterized_ClaSSE_As_v5(t, pG.A[:,:], pG.p_Ds_v5), G)
	#display(dG)
	#return(dG)
end # End calc_Gs_SSE



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
using PhyloNetworks
using PhyBEARS.TrUtils # for flat2() (similar to unlist)
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
using PhyloNetworks

<<<<<<< HEAD
inputs = ModelLikes.setup_DEC_SSE(2, readTopology("((chimp:1,human:1):1,gorilla:2);"))
#inputs = ModelLikes.setup_MuSSE_biogeo(2, readTopology("((chimp:10,human:10):10,gorilla:20);"))
=======
inputs = ModelLikes.setup_DEC_SSE(2, readTopology("((chimp:10,human:10):10,gorilla:20);"))
# inputs = ModelLikes.setup_MuSSE(2, readTopology("((chimp:10,human:10):10,gorilla:20);"))
>>>>>>> 8719c19396bb1ead82f32ce725d8a52c371756cf
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
# 				push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
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
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
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





# Calculate G flow matrix down time series, outputting various
# condition numbers and kappa rates
calc_Gs_SSE_condnums! = (dG, G, pG, t) -> begin
	#A = pG.A             # Initial A
	p_Ds_v5 = pG.p_Ds_v5 # Calculate Ds down a timespan
	
	# A as a function of time t
	A = parameterized_ClaSSE_As_v6(t, pG.A[:,:], p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)
#	mul!(dG, A, G)
	mul!(dG, parameterized_ClaSSE_As_v6(t, pG.A[:,:], p_Ds_v5), G)

	condG1 = cond(G,1)
	condG2 = cond(G,2)
	condGInf = cond(G,Inf)
	tmpstr = paste0(["\nAt time t=", string(round(t, digits=6)), ", Condition number of G=", string(round(condG1, digits=6)), ", ", string(round(condG2, digits=6)), ", ", string(round(condGInf, digits=6))])
	print(tmpstr)
	#display(cond(G))
	
	# From Louca & Pennell code:
	# phylogenetics_cpp_routines.cpp
	# max_condition_number,			// (INPUT) unitless number, 
	# the maximum acceptable condition number for the Gmap 
	# (as estimated from the linearized dynamics), when choosing 
	# the integration interval size. A larger max_condition number 
	# leads to fewer age-splits, thus faster computation but also 
	# lower accuracy. Hence, this number controls the trade-off 
	# between speed and accuracy. Typical values are 
	# 1e4 (slower, more accurate) up to 
	# 1e8 (faster, less accurate).
	
	# The condition number is kappa:
	# https://julia.quantecon.org/tools_and_techniques/iterative_methods_sparsity.html
	# The cond() operation can be expensive, so the inequality in Louca, Supp. Mat. Eq. 3 is useful
	# It looks *extremely* conservative, it blows up at e.g.
	# time = 0.00015
	#
	# upper_bound_condition_number: 18425.777249466213
	# 
	# At time t=0.00015, Condition number=2.175764, 1.658813, 1.70285
	# opnorms p=1 & p=2:
	# 32759.807937823098
	# 30274.76762619003
	# 30274.767626190034
	# 17479.145238634184
	#
	# vs.
	# 
	# upper_bound_condition_number: 1.0971658583069032e6
	# 
	# At time t=0.0002, Condition number=3.261165, 2.255482, 2.334215
	# opnorms p=1 & p=2:
	# 34805.07811454515
	# 32164.890247616975
	# 32164.89024761698
	# 18570.408042916435
	# 
	# Probably we could just use opnorm(A,1), or even opnorm(A,1)^2
	
	
	# Matrix norms
	# See: 
	# Lambers, Jim (2009). Vector Norms & Matrix Norms. pp. 1-16.
	# https://www.math.usm.edu/lambers/mat610/sum10/lecture2.pdf
	# These notes have the inequalities between norm forms
	
	# https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/
	# Note: operator norm = matrix norm
	# When p=1, ||A||1, much faster, seems to always be bigger
	# When p=2, ||A||2, the operator norm is the spectral norm, equal to the largest singular value of A
	# When p=Inf, ||A||Inf, this is just opnorm(t(A),1), and 
	#
	# ||A||2 <= sqrt(ncol(A))*||A||Inf
	# ||A||2 <= sqrt(nrow(A))*||A||Inf
	# 
	# (doesn't seem to be true, actually. But cond of opnorm(1) does seem conservative
	
	# Actually, instead of tracking the condition number kappa, they are tracking the 
	# *growth rate* of kappa:
	#
	# // determine linear dynamics (matrix form) of D at various representative 
	# ages, and keep track of the largest estimated exponential growth rate of 
	# the condition number of Gmap ("kappa_rate")
	# 
	# // calculate largest singular value (sigma1) of the dynamics at this age
	#	// then kappa_rate <= 2*sigma1, since 
	# // sigma1(exp(t*A))/sigma2(exp(t*A)) <= exp(t*2*sigma1(A))
	# [So and Thompson (2000) Singular values of Matrix Exponentials. Theorem 2.1]
	# 
	
	# From Louca & Pennell:
	# max_condition_number,	
	# // (INPUT) unitless number, the maximum acceptable condition number for the Gmap 
	# (as estimated from the linearized dynamics), when choosing the integration 
	# interval size. A larger max_condition number leads to fewer age-splits, thus 
	# faster computation but also lower accuracy. Hence, this number controls the 
	# trade-off between speed and accuracy. Typical values are 1e4 (slower, more accurate) 
	# up to 1e8 (faster, less accurate).
	#
	# The comparison is done in:
	# const double DeltaT 		= max(root_age/max_Nintervals, 
	# min(1.0000001*root_age,log(max_condition_number)/max_kappa_rate));
	# 
	# If the maximum observed kappa was 20, and the max condition number
	# was 1e8, log(1e8) would be 18.4
	# 

# Rescaling after Delta_int time interval:
#
# Just take mean of the vector:
#
# inline double vector_mean(const std::vector<double> &values){
# 	double S = 0;
# 	for(long i=0; i<values.size(); ++i) S += values[i];
# 	return (S/values.size());
# }
# 
# 		const double initial_mean = vector_mean(initial);
# 		if(initial_mean<=0) return false;
# 		scale = log(initial_mean);
# 		shape = initial/initial_mean;
# 		return true; 
# 	}
# 
# 	// record a new time series point, provided by the numerical solver
# 	void registerState(double age, const MuSSEstateD &state){
# 		trajectory.push_back(state); 
# 		ages.push_back(age); 
# 
# 		// make sure entries are in [0,1]
# 		const long i = trajectory.size()-1;
# 		for(long s=0; s<trajectory[i].size(); ++s) trajectory[i][s] = max(0.0, min(1.0, trajectory[i][s]));
# 	}
# 	
# 	
# 	// record a new trajectory point, provided by the numerical solver
# 	// The state X to be recorded is provided in rescaled format, i.e. state = exp(scale) * shape
# 	// You can either record shape and scale separately, or combine them to obtain the actual state
# 	void registerScaledState(double age, const MuSSEstateD &shape, const double scale){
# 		trajectory_shape.push_back(shape);
# 		trajectory_scale.push_back(scale);
# 		ages.push_back(age); 






	
	print("\nopnorms(A, p=1, p=2, p=Inf, sqrt(nrow(A))*||A||Inf:\n")
	display(opnorm(A,1))
	display(opnorm(A,2))
	display(opnorm(transpose(A),2))
	display(1/(sqrt(size(A)[1]))*opnorm(transpose(A),2))
	print("\n")
	#display(opnorm(A,Inf))

	sigma1_of_A = opnorm(A,1)
	upper_bound_condition_number = exp(2*t*sigma1_of_A)
	upper_bound_kappa_growth_rate = 2*sigma1_of_A
	tmpstr = paste0(["\nupper_bound_condition_number of A=", string(upper_bound_condition_number), "\nupper_bound_kappa_growth_rate=", string(upper_bound_kappa_growth_rate)])
	print(tmpstr)
	print("\n")
	


	# No return needed, what is returned is G (as .u)
	
	#display(dG)
	#return(dG)
end # End calc_Gs_SSE_condnums



# Calculate G flow matrix down time series
# (no printed outputs)
calc_Gs_SSE = (dG, G, pG, t; max_condition_number=1e8) -> begin
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)

	# The new dG is A %*% G
	mul!(dG, A, G)

	# No return needed, what is returned is G (as .u)
	return(dG)
end # End calc_Gs_SSE




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


end # end module Flow

