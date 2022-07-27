#######################################################
# Construct Gmaps: A series of Gflow objects, with
# customizable time-steps
#######################################################
module Gmaps

print("PhyBEARS: loading Gmaps dependencies...")

using DifferentialEquations
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using DoubleFloats		# for Double64
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

using PhyloBits.TreeTable  # for get_nonrootnodes_trdf
using PhyloBits.TrUtils # for flat2() (similar to unlist)
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.Flow

print("...done.\n")

export find_Gflow_increment_size, get_timebin, interp_from_Gmap, construct_Gmap_interpolator, construct_Gmap_interpolator_float64, construct_Gmap_interpolator_float64_parallel, construct_Gmap_interpolator_double64, construct_Gmap_interpolator_double64_parallel

"""
include("/GitHub/PhyBEARS.jl/src/Gmaps.jl")
"""


# Check the linear dynamics (matrix A) for timespans where the kappa rate
# exceeds log(max_condition_number). max_condition_number is usually between
# 1e4 (slower) and 1e8 (faster)
# 
# "exponential growth rate of the condition number ("kappa_rate") based on the 
#  largest singular value of the dynamics A" 
#
# Using  A[:,:] is required, so that "A" doesn't mutate
# 
function find_Gflow_increment_size(tvals, p_Ds_v5; max_condition_number=1e8)

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
		A_at_t = Flow.parameterized_ClaSSE_As_v6(t, A[:,:], p_Ds_v5; max_condition_number=max_condition_number, print_warnings=false)
		sigma1_of_A = opnorm(A_at_t,2)  # use the 2-norm here (this will save us steps)
		#upper_bound_kappa_rates_A[i] = sqrt(n)*2*sigma1_of_A # if using opnorm(A_at_t, 1)
		upper_bound_kappa_rates_A[i] = 2*sigma1_of_A
		if (upper_bound_kappa_rates_A[i] > log(max_condition_number))
			return tvals[i] 
		end
	end
	
	return tvals[length(tvals)]
end # END function find_Gflow_increment_size(tvals, p_Ds_v5; max_condition_number=1e8)




"""
#######################################################
# Script for Gmaps.jl
#######################################################

# Construct a series of interpolators at given time increments
using DifferentialEquations
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using PhyBEARS.Flow
using PhyBEARS.Parsers

include("/GitHub/PhyBEARS.jl/src/Gmaps.jl")
using .Gmaps

tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")

lgdata_fn = "/GitHub/PhyBEARS.jl/Rsrc/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.est[bmo.rownames .== "birthRate"] .= 0.32881638319078066
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
bmo.est[:] = bmo_updater_v1(bmo);
numareas = 4
n = 16            # 4 areas, 16 states

# Set up the model
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-6;
solver_options.reltol = 1e-6;

p_Ds_v5 = inputs.p_Ds_v5;
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

# Solve the Es
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);

res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6

tspan = (0.0, 1.1 * maximum(trdf.node_age))
prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG);
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-6, reltol = 1e-6);

root_age = trdf[tr.root,:node_age]
Gseg_times = seq(0.0, root_age, 0.1);
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap = construct_Gmap_interpolator(pG, Gseg_times; abstol=1e-6, reltol=1e-6);

Gflows_array_totals[:,:,1]
interp_from_Gmap(0.1, Gmap)
Gflow_to_01_GMRES(0.1)


Gflows_array_totals[:,:,2]
interp_from_Gmap(0.2, Gmap)
Gflow_to_01_GMRES(0.2)


Gflows_array_totals[:,:,3]
interp_from_Gmap(0.3, Gmap)
Gflow_to_01_GMRES(0.3)

Gflows_array_totals[:,:,50]
interp_from_Gmap(5.0, Gmap)
Gflow_to_01_GMRES(5.0)



Gflow = t -> interp_from_Gmap(t, Gmap)

# The Gmap strategy works OK, with many increments...
res_Gflow_v6a = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6a, rootstates_lnL_GFv6a, Julia_total_lnLs1_GFv6a, bgb_lnl_GFv6a) = res_Gflow_v6a

abstol=1e-9
reltol=1e-9
inc = 1
"""


"""
# Get the time-bin a segment falls into
# (assumes first number is NOT 0!)
t = 1.5
timebins = [0.5, 1.0, 1.5, 2.0]
get_timebin(1.5, timebins)
get_timebin(0.4, timebins)
"""
function get_timebin(t, timebins)
	#print(paste0([string(t), ", "]))
	if issorted(timebins) == false
		txt = "STOP ERROR in get_timebin(). The input 'timebins' must be sorted to be monotonic increasing. Use sort!() to fix."
		throw(txt)
	end
	if timebins[1] < 1.0e-10
		txt = "STOP ERROR in get_timebin(). The input 'timebins' must not start with 0. Please use popfirst!(timebins) to remove."
		throw(txt)
	end

	TF2 = round(t, digits=8) .≤ timebins   # it seems like .<= produces error
	
	# Take the last hit 
	indices = collect(1:length(timebins))
	index = indices[TF2][1]
		
	return index
end # END function get_timebin(t, timebins)

# Get the Ds at time t
function interp_from_Gmap(t, Gmap)
	n = dim(Gmap.Gflows_array_totals[:,:,1])[1]
	Gflow_output = Matrix{Float64}(I, n, n);
	index = get_timebin(t, Gmap.Gseg_times)
	if index == 1
		Gflow_total_old = Matrix{Float64}(I, n, n);
		old_time = 0.0
	else
		Gflow_total_old = Gmap.Gflows_array_totals[:,:,(index-1)]  # Get the starting Gflow
		old_time = Gmap.Gseg_times[index-1]
	end

	Gflow_res = Gmap.Gflows_dict[index](t)
	mul!(Gflow_output, Gflow_total_old, Gflow_res)
	return Gflow_output
end # END function interp_from_Gmap(t, Gmap)






# Function that can switch between float64 and double 64 versions (requires "using DoubleFloats")
function construct_Gmap_interpolator(pG, Gseg_times; abstol=1e-6, reltol=1e-6, use_double=false)
	if (use_double == false)
		Gmap = construct_Gmap_interpolator_float64(pG, Gseg_times; abstol=abstol, reltol=reltol)
	else
		Gmap = construct_Gmap_interpolator_double64(pG, Gseg_times; abstol=abstol, reltol=reltol)
	end
	"""
	(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap
	"""
	return Gmap
end # END construct_Gmap_interpolator



# Standard calculation with float64 improves matrix multiplication -- parallelized
function construct_Gmap_interpolator_float64_parallel(pG, Gseg_times; abstol=1e-6, reltol=1e-6)
	# Check number of threads, create tasks list:
	numthreads = Base.Threads.nthreads()
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	# Check for 0.0 in first slot; remove
	if (Gseg_times[1] < 1.0e-10)
		popfirst!(Gseg_times) # remove first element from array
	end
	num_segments = length(Gseg_times)
	Gflows_array = zeros(Float64,pG.n,pG.n,num_segments)
	Gflow_total = Matrix{Float64}(I, pG.n, pG.n);
	Gflows_array_totals = zeros(Float64,pG.n,pG.n,num_segments)
	Gflows_dict = Dict()
	
	# Create tasks list
	old_time = 0.0
	for inc in 1:num_segments
		if inc == 1
			Gflow_total_old = Matrix{Float64}(I, pG.n, pG.n);
		else
			Gflow_total_old = Gflows_array_totals[:,:,(inc-1)]
		end
		G0 = Matrix{Float64}(I, pG.n, pG.n) ;
		# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
		pG.A[:,:] .= 0.0
		new_time = Gseg_times[inc]
		end_tspan = new_time * 1.01
		tspan = (old_time, end_tspan)
		prob_Gs_v5 = DifferentialEquations.ODEProblem(calc_Gs_SSE_sub_i, G0, tspan, pG);
		Gflows_dict[inc] = Base.Threads.@spawn solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=abstol, reltol=reltol);
		Gres = Gflows_dict[inc](new_time)
		Gflows_array[:,:,inc] = Gres
		mul!(Gflow_total, Gflow_total_old, Gres)
		Gflows_array_totals[:,:,inc] = Gflow_total
		old_time = new_time
	end
	
	# Return results
	Gmap = (Gseg_times=Gseg_times, Gflows_array=Gflows_array, Gflows_array_totals=Gflows_array_totals, Gflows_dict=Gflows_dict);
	"""
	(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap
	"""
	return Gmap
end # END function construct_Gmap_interpolator_float64



# Standard calculation with float64 improves matrix multiplication
function construct_Gmap_interpolator_float64(pG, Gseg_times; abstol=1e-6, reltol=1e-6)
	# Check for 0.0 in first slot; remove
	if (Gseg_times[1] < 1.0e-10)
		popfirst!(Gseg_times) # remove first element from array
	end
	num_segments = length(Gseg_times)
	Gflows_array = zeros(Float64,pG.n,pG.n,num_segments)
	Gflow_total = Matrix{Float64}(I, pG.n, pG.n);
	Gflows_array_totals = zeros(Float64,pG.n,pG.n,num_segments)
	Gflows_dict = Dict()
	
	old_time = 0.0
	for inc in 1:num_segments
		if inc == 1
			Gflow_total_old = Matrix{Float64}(I, pG.n, pG.n);
		else
			Gflow_total_old = Gflows_array_totals[:,:,(inc-1)]
		end
		G0 = Matrix{Float64}(I, pG.n, pG.n) ;
		# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
		pG.A[:,:] .= 0.0
		new_time = Gseg_times[inc]
		end_tspan = new_time * 1.01
		tspan = (old_time, end_tspan)
		prob_Gs_v5 = DifferentialEquations.ODEProblem(calc_Gs_SSE_sub_i, G0, tspan, pG);
		Gflows_dict[inc] = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=abstol, reltol=reltol);
		Gres = Gflows_dict[inc](new_time)
		Gflows_array[:,:,inc] = Gres
		mul!(Gflow_total, Gflow_total_old, Gres)
		Gflows_array_totals[:,:,inc] = Gflow_total
		old_time = new_time
	end
	
	# Return results
	Gmap = (Gseg_times=Gseg_times, Gflows_array=Gflows_array, Gflows_array_totals=Gflows_array_totals, Gflows_dict=Gflows_dict);
	"""
	(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap
	"""
	return Gmap
end # END function construct_Gmap_interpolator_float64

# using Double64 might improve matrix multiplication
function construct_Gmap_interpolator_double64_parallel(pG, Gseg_times; abstol=1e-6, reltol=1e-6)
	# Check number of threads, create tasks list:
	numthreads = Base.Threads.nthreads()
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	# Check for 0.0 in first slot; remove
	if (Gseg_times[1] < 1.0e-10)
		popfirst!(Gseg_times) # remove first element from array
	end
	num_segments = length(Gseg_times)
	Gflows_array = zeros(Double64,pG.n,pG.n,num_segments)
	Gflow_total = Matrix{Double64}(I, pG.n, pG.n);
	Gflows_array_totals = zeros(Double64,pG.n,pG.n,num_segments)
	Gflows_dict = Dict()
	
	old_time = 0.0
	for inc in 1:num_segments
		if inc == 1
			#Gflow_total_old = Matrix{Float64}(I, pG.n, pG.n);
			Gflow_total_old = Matrix{Double64}(I, pG.n, pG.n);
		else
			Gflow_total_old = Gflows_array_totals[:,:,(inc-1)]
		end
		G0 = Matrix{Float64}(I, pG.n, pG.n) ;  # Can't be double going into DifferentialEquations
		# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
		pG.A[:,:] .= 0.0
		new_time = Gseg_times[inc]
		end_tspan = new_time * 1.01
		tspan = (old_time, end_tspan)
		#txt = paste0(["\ninc=", string(inc), ": ", string(old_time), ", ", string(new_time)])
		#print(txt)
		prob_Gs_v5 = DifferentialEquations.ODEProblem(calc_Gs_SSE_sub_i, G0, tspan, pG);
		Gflows_dict[inc] = Base.Threads.@spawn solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=abstol, reltol=reltol);
		Gres = Gflows_dict[inc](new_time)
		Gflows_array[:,:,inc] = Double64.(Gres)
		mul!(Gflow_total, Gflow_total_old, Gres)
		Gflows_array_totals[:,:,inc] = Gflow_total
		old_time = new_time
	end
	
	# Return results
	Gmap = (Gseg_times=Gseg_times, Gflows_array=Gflows_array, Gflows_array_totals=Gflows_array_totals, Gflows_dict=Gflows_dict);
	"""
	(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap
	"""
	return Gmap
end # END function construct_Gmap_interpolator_double64


# using Double64 might improve matrix multiplication
function construct_Gmap_interpolator_double64(pG, Gseg_times; abstol=1e-6, reltol=1e-6)
	# Check for 0.0 in first slot; remove
	if (Gseg_times[1] < 1.0e-10)
		popfirst!(Gseg_times) # remove first element from array
	end
	num_segments = length(Gseg_times)
	Gflows_array = zeros(Double64,pG.n,pG.n,num_segments)
	Gflow_total = Matrix{Double64}(I, pG.n, pG.n);
	Gflows_array_totals = zeros(Double64,pG.n,pG.n,num_segments)
	Gflows_dict = Dict()
	
	old_time = 0.0
	for inc in 1:num_segments
		if inc == 1
			#Gflow_total_old = Matrix{Float64}(I, pG.n, pG.n);
			Gflow_total_old = Matrix{Double64}(I, pG.n, pG.n);
		else
			Gflow_total_old = Gflows_array_totals[:,:,(inc-1)]
		end
		G0 = Matrix{Float64}(I, pG.n, pG.n) ;  # Can't be double going into DifferentialEquations
		# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
		pG.A[:,:] .= 0.0
		new_time = Gseg_times[inc]
		end_tspan = new_time * 1.01
		tspan = (old_time, end_tspan)
		#txt = paste0(["\ninc=", string(inc), ": ", string(old_time), ", ", string(new_time)])
		#print(txt)
		prob_Gs_v5 = DifferentialEquations.ODEProblem(calc_Gs_SSE_sub_i, G0, tspan, pG);
		Gflows_dict[inc] = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=abstol, reltol=reltol);
		Gres = Gflows_dict[inc](new_time)
		Gflows_array[:,:,inc] = Double64.(Gres)
		mul!(Gflow_total, Gflow_total_old, Gres)
		Gflows_array_totals[:,:,inc] = Gflow_total
		old_time = new_time
	end
	
	# Return results
	Gmap = (Gseg_times=Gseg_times, Gflows_array=Gflows_array, Gflows_array_totals=Gflows_array_totals, Gflows_dict=Gflows_dict);
	"""
	(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap
	"""
	return Gmap
end # END function construct_Gmap_interpolator_double64










end # END module Gmaps
