#######################################################
# Likelihood calculations: input tree, tipdata, model,
#                          get likelihood back
#######################################################

module ModelLikes

print("\n\nStarting module 'ModelLikes'...loading dependencies...\n")
using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA           # for lsoda()
using Sundials        # for CVODE_BDF(linear_solver=:GMRES)
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

# (1) List all function names here:
export say_hello2, setup_MuSSE_biogeo, workprecision, setup_DEC_SSE, calclike_DEC_SSE, setup_DEC_Cmat2

#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/PhyBEARS.jl/notes/")
include("tst2.jl")
"""
#######################################################


#######################################################
# (2) write the functions here
#######################################################

say_hello2() = println("Hello dude2!")
say_hello2()
"""
# Set up a DEC-like model for numareas areas
# (numstates = (2^numareas)-1
# Will calculate Es over 120% of root depth.

include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
inputs = ModelLikes.setup_MuSSE_biogeo()
prtQi(inputs)
prtCi(inputs)


"""

function setup_MuSSE_biogeo(numstates=2, tr=readTopology("((chimp:1,human:1):1,gorilla:2);"); root_age_mult=1.5, in_params=NaN)
	#numareas=2
	#tr=readTopology("((chimp:1,human:1):1,gorilla:2);")
	areas_list = collect(1:numstates)
	total_numareas = length(areas_list)
	numareas = length(areas_list)
	
	type_string = string(typeof(in_params))
	if (startswith(type_string, "NamedTuple") == false) && (isnan(in_params) == true)
		in_params = (birthRate=0.2, deathRate=0.0, d_val=0.0, e_val=0.0, a_val=0.0, j_val=0.0)
	end
	
	states_list = collect(1:numstates)
	n = length(states_list)

	res = construct_Res(tr, n)
	rootnodenum = tr.root
	trdf = prt(tr, rootnodenum)
	tipnodes = trdf[!,1][trdf[!,10].=="tip"]
	
	birthRate = in_params.birthRate
	deathRate = in_params.deathRate
	
	d_val = in_params.d_val
	e_val = in_params.e_val
	a_val = in_params.a_val
	j_val = in_params.j_val
	
	dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list)))
	amat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list)))
	elist = repeat([1.0], length(areas_list))
	
	# Set up a pure-"a" Qmat
	# (a = anagenetic range-switching between single areas)
	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["a"])
	Qarray_event_types = Qmat.Qarray_event_types
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	
	# Update Qij parameters, manually
	dTF = Qarray_event_types .== "d"
	Qmat.Qij_vals[dTF] = Qmat.Qij_vals[dTF] .* in_params.d_val
	eTF = Qarray_event_types .== "e"
	Qmat.Qij_vals[eTF] = Qmat.Qij_vals[eTF] .* in_params.e_val
	aTF = Qarray_event_types .== "a"
	Qmat.Qij_vals[aTF] = Qmat.Qij_vals[aTF] .* in_params.a_val
	
	Qmat = (Qarray_event_types=Qarray_event_types, Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals)
	
	# Default values of y, s, v, and j
	Cparams = default_Cparams()
# 	maxent_constraint_01 = 0.0
# 	maxent01symp = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
# 	maxent01sub = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
# 	maxent01jump = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
# 	maxent_constraint_01 = 0.0
# 	maxent01vic = relative_probabilities_of_vicariants(total_numareas, maxent_constraint_01)
# 	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
# 	Carray = setup_DEC_Cmat(areas_list, states_list, maxent01, Cparams)
	
	# Set up a sympatry-only Cmat, manually
	Carray_ivals = collect(1:n)
	Carray_jvals = collect(1:n)
	Carray_kvals = collect(1:n)
	Cijk_weights = repeat([1.0], n)
	Cijk_vals = repeat([birthRate], n)
	Carray_event_types = repeat(["y"], n) # y=sYmpatric speciation (for MuSSE)
	
	Carray = (Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, Carray_event_types=Carray_event_types)
	
	prtQ(Qmat)
	prtC(Carray)
	
	# Set up mu (extinction) rates, manually
	mu_vals = repeat([deathRate], n)

	params = (mu_vals=mu_vals, Qij_vals=Qmat.Qij_vals, Cijk_weights=Cijk_weights, Cijk_vals=Carray.Cijk_vals)

	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	p_indices = (Qarray_ivals=Qmat.Qarray_ivals, Qarray_jvals=Qmat.Qarray_jvals, Qarray_event_types=Qmat.Qarray_event_types, Carray_ivals=Carray.Carray_ivals, Carray_jvals=Carray.Carray_jvals, Carray_kvals=Carray.Carray_kvals, Carray_event_types=Carray.Carray_event_types)

	# True/False statements by index
	# The calculation of dEi and dDi for state i involves many
	# ==i and !=i operations across Q and C. These only need to be 
	# done once per problem (may or may not save time to 
	# pre-calculate).
	# 
	# Pre-allocating the Carray_ivals .== i, Qarray_jvals[Qarray_ivals .== i
	# Reduces GC (Garbage Collection) from 40% to ~5%
	# 10+ times speed improvement (!)
	Qi_eq_i = Any[]
	Ci_eq_i = Any[]

	Qi_sub_i = Any[]
	Qj_sub_i = Any[]

	# These are the (e.g.) j state-indices (left descendant) when the ancestor==state i
	Ci_sub_i = Any[]
	Cj_sub_i = Any[]
	Ck_sub_i = Any[]

	# The push! operation may get slow at huge n
	# This will have to change for non-Mk models
	for i in 1:n
		push!(Qi_eq_i, Qmat.Qarray_ivals .== i)
		push!(Qi_sub_i, Qmat.Qarray_ivals[Qarray_ivals .== i])
		push!(Qj_sub_i, Qmat.Qarray_jvals[Qarray_ivals .== i])

		push!(Ci_eq_i, Carray.Carray_ivals .== i)
		push!(Ci_sub_i, Carray.Carray_ivals[Carray.Carray_ivals .== i])
		push!(Cj_sub_i, Carray.Carray_jvals[Carray.Carray_ivals .== i])
		push!(Ck_sub_i, Carray.Carray_kvals[Carray.Carray_ivals .== i])
	end

	# Inputs to the Es calculation
	p_TFs = (Qi_eq_i=Qi_eq_i, Ci_eq_i=Ci_eq_i, Qi_sub_i=Qi_sub_i, Qj_sub_i=Qj_sub_i, Ci_sub_i=Ci_sub_i, Cj_sub_i=Cj_sub_i, Ck_sub_i=Ck_sub_i)
	p_orig = (n=n, params=params, p_indices=p_indices)
	p = p_orig
	p_Es_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs)
	
	# Solutions to the E vector
	uE = repeat([0.0], n)
	max_t = root_age_mult*trdf[tr.root,:node_age]
	Es_tspan = (0.0, max_t) # 110% of tree root age
	#by_t = max_t / 10.0
	#Es_tspan = collect(0.0:by_t:max_t)

	#######################################################
	# Downpass with ClaSSE
	#######################################################
	# Solve for the Ds
	du = repeat([0.0], n)
	u0 = repeat([0.0], n)
	u0[2] = 1.0  # Starting likelihood
	#tspan = (0.0, 2.0*trdf[tr.root,:node_age]) # 110% of tree root age
	current_nodeIndex = 1
	res = construct_Res(tr, n)
	
	res.likes_at_each_nodeIndex_branchTop
	for i in 1:tr.numTaxa
		# [:] avoids creating a linked reference
		res.likes_at_each_nodeIndex_branchTop[tipnodes[i]] = u0[:];
		res.normlikes_at_each_nodeIndex_branchTop[tipnodes[i]] = u0[:] / sum(u0[:]);
		res.Es_at_each_nodeIndex_branchTop = uE[:];
		res.Es_at_each_nodeIndex_branchBot = uE[:];
	end
	#res.likes_at_each_nodeIndex_branchTop[6] = u0;
	res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]

	# Updates res
	res_orig = res
	res_orig.likes_at_each_nodeIndex_branchTop

	solver_options = construct_SolverOpt()
	solver_options.solver=CVODE_BDF(linear_solver=:GMRES) # requires Sundials
	solver_options.save_everystep = true		# false produced pathological interpolation 
																					# for lsoda, Tsit5; GMRES worked either way
	solver_options.abstol = 1.0e-6
	solver_options.reltol = 1.0e-6
	
	p_Ds_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs, prob=prob_Es_v5, sol_Es_v5=sol_Es_v5, uE=uE)
	"""
	res = inputs.res
	trdf = inputs.trdf
	solver_options = inputs.solver_options
	p_Ds_v5 = inputs.p_Ds_v5
	"""
	
	inputs = (res=res, trdf=trdf, solver_options=solver_options, p_Ds_v5=p_Ds_v5, p_Es_v5=p_Es_v5, Es_tspan=Es_tspan)
	
	return inputs
end # End function setup_MuSSE_biogeo


# Makes sol_Es from the stored function parameterized_ClaSSE_Es_v5
"""
	print("\nSolving the Es once, for the whole tree timespan...")
	
	# Solve the Es
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, uE, Es_tspan, p_Es_v5)

	print(Es_tspan)
	
	# This solution is a linear interpolator
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
	
	print("...solving Es has finished, creating interpolator 'sol_Es_v5'.\n")
"""





# Set up a DEC-like model for numareas areas
# (numstates = (2^numareas)-1
# Will calculate Es over 120% of root depth.
function setup_DEC_SSE(numareas=2, tr=readTopology("((chimp:1,human:1):1,gorilla:2);"); root_age_mult=1.5, max_range_size=NaN, include_null_range=false, in_params=NaN)
	#numareas=2
	#tr=readTopology("((chimp:1,human:1):1,gorilla:2);")
	areas_list = collect(1:numareas)
	total_numareas = length(areas_list)
	
	if (isnan(in_params) == true)
		in_params = (birthRate=0.2, deathRate=0.0, d_val=0.0, e_val=0.0, a_val=0.0, j_val=0.0)
	end
	
	# Check if max_range_size=NaN
	if (isnan(max_range_size) == true)
		max_range_size = numareas
	end
	
	states_list = areas_list_to_states_list(areas_list, max_range_size, include_null_range)
	n = length(states_list)

	
	res = construct_Res(tr, n)
	rootnodenum = tr.root
	trdf = prt(tr, rootnodenum)
	tipnodes = trdf[!,1][trdf[!,10].=="tip"]
	
	birthRate = 0.2
	deathRate = 0.1
	
	d_val = 0.0
	e_val = 0.0
	j_val = 0.2
	
	dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list)))
	amat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list)))
	elist = repeat([1.0], length(areas_list))
	
	dmat = d_val .* dmat
	elist = e_val .* elist
	
	
	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	print("\n")
	print("elist")
	print(elist)
	print("\n")
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	Qarray_event_types = Qmat.Qarray_event_types


	Cparams = default_Cparams()
	
	Cparams.j = j_val
	wt = (3.0 - j_val) / 3.0
	Cparams.y = wt
	Cparams.s = wt
	Cparams.v = wt
	Cparams
	
	maxent_constraint_01 = 0.0
	maxent01symp = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
	maxent01sub = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
	maxent01jump = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
	maxent_constraint_01 = 0.0
	maxent01vic = relative_probabilities_of_vicariants(total_numareas, maxent_constraint_01)
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
	Carray = setup_DEC_Cmat(areas_list, states_list, maxent01, Cparams)

	# Possibly varying parameters
	mu_vals = repeat([deathRate], n)

	params = (mu_vals=mu_vals, Qij_vals=Qmat.Qij_vals, Qarray_event_types=Qarray_event_types, Cijk_vals=Carray.Cijk_vals)
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	p_indices = (Qarray_ivals=Qmat.Qarray_ivals, Qarray_jvals=Qmat.Qarray_jvals, Carray_ivals=Carray.Carray_ivals, Carray_jvals=Carray.Carray_jvals, Carray_kvals=Carray.Carray_kvals)

	# True/False statements by index
	# The calculation of dEi and dDi for state i involves many
	# ==i and !=i operations across Q and C. These only need to be 
	# done once per problem (may or may not save time to 
	# pre-calculate).
	# 
	# Pre-allocating the Carray_ivals .== i, Qarray_jvals[Qarray_ivals .== i
	# Reduces GC (Garbage Collection) from 40% to ~5%
	# 10+ times speed improvement (!)
	Qi_eq_i = Any[]
	Ci_eq_i = Any[]

	Qi_sub_i = Any[]
	Qj_sub_i = Any[]

	# These are the (e.g.) j state-indices (left descendant) when the ancestor==state i
	Ci_sub_i = Any[]
	Cj_sub_i = Any[]
	Ck_sub_i = Any[]

	# The push! operation may get slow at huge n
	# This will have to change for non-Mk models
	for i in 1:n
		push!(Qi_eq_i, Qmat.Qarray_ivals .== i)
		push!(Qi_sub_i, Qmat.Qarray_ivals[Qarray_ivals .== i])
		push!(Qj_sub_i, Qmat.Qarray_jvals[Qarray_ivals .== i])

		push!(Ci_eq_i, Carray.Carray_ivals .== i)
		push!(Ci_sub_i, Carray.Carray_ivals[Carray.Carray_ivals .== i])
		push!(Cj_sub_i, Carray.Carray_jvals[Carray.Carray_ivals .== i])
		push!(Ck_sub_i, Carray.Carray_kvals[Carray.Carray_ivals .== i])
	end

	# Inputs to the Es calculation
	p_TFs = (Qi_eq_i=Qi_eq_i, Ci_eq_i=Ci_eq_i, Qi_sub_i=Qi_sub_i, Qj_sub_i=Qj_sub_i, Ci_sub_i=Ci_sub_i, Cj_sub_i=Cj_sub_i, Ck_sub_i=Ck_sub_i)
	p_orig = (n=n, params=params, p_indices=p_indices)
	p = p_orig
	p_Es_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs)
	
	# Solutions to the E vector
	uE = repeat([0.0], n)
	max_t = root_age_mult*trdf[tr.root,:node_age]
	Es_tspan = (0.0, max_t) # 110% of tree root age
	#by_t = max_t / 10.0
	#Es_tspan = collect(0.0:by_t:max_t)

	#######################################################
	# Downpass with ClaSSE
	#######################################################
	# Solve for the Ds
	du = repeat([0.0], n)
	u0 = repeat([0.0], n)
	u0[2] = 1.0  # Starting likelihood
	#tspan = (0.0, 2.0*trdf[tr.root,:node_age]) # 110% of tree root age
	current_nodeIndex = 1
	res = construct_Res(tr, n)
	
	res.likes_at_each_nodeIndex_branchTop
	for i in 1:tr.numTaxa
		# [:] avoids creating a linked reference
		res.likes_at_each_nodeIndex_branchTop[tipnodes[i]] = u0[:];
		res.normlikes_at_each_nodeIndex_branchTop[tipnodes[i]] = u0[:] / sum(u0[:]);
	end
	#res.likes_at_each_nodeIndex_branchTop[6] = u0;
	res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]

	# Updates res
	res_orig = res
	res_orig.likes_at_each_nodeIndex_branchTop

	solver_options = construct_SolverOpt()
	solver_options.solver=CVODE_BDF(linear_solver=:GMRES) # requires Sundials
	solver_options.save_everystep = true		# false produced pathological interpolation 
																					# for lsoda, Tsit5; GMRES worked either way
	solver_options.abstol = 1.0e-6
	solver_options.reltol = 1.0e-6
	
	print("\nSolving the Es once, for the whole tree timespan...")
	
	# Solve the Es
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, uE, Es_tspan, p_Es_v5)
	
	# This solution is a linear interpolator
	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
	
	print("...solving Es has finished, creating interpolator 'sol_Es_v5'.\n")
	
	p_Ds_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs, prob=prob_Es_v5, sol_Es_v5=sol_Es_v5, uE=uE)
	"""
	res = inputs.res
	trdf = inputs.trdf
	solver_options = inputs.solver_options
	p_Ds_v5 = inputs.p_Ds_v5
	"""
	
	inputs = (res=res, trdf=trdf, solver_options=solver_options, p_Ds_v5=p_Ds_v5, p_Es_v5=p_Es_v5, Es_tspan=Es_tspan)
	
	return inputs
end # End function setup_DEC_SSE










"""
areas_list = [1,2,3]
states_list = areas_list_to_states_list(areas_list, 3, true)
predeclare_array_length=10000000
Carray = setup_DEC_Cmat(areas_list, states_list)

i = 8
j = i
k = i
ancstate = states_list[i]
ancsize = length(ancstate)
lstate = states_list[j]
lsize = length(lstate)
rstate = states_list[k]
rsize = length(rstate)	

ancstate == lstate == rstate

areas_list = [1,2,3]
states_list = areas_list_to_states_list(areas_list, 3, true)
Cparams=(y=1.0,s=1.0,v=1.0,j=0.0)
total_numareas = length(areas_list)
maxent_constraint_01 = 0.0
maxent01symp = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
maxent01sub = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
maxent01jump = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
maxent_constraint_01 = 0.5
maxent01vic = relative_probabilities_of_vicariants(total_numareas, maxent_constraint_01)
maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
predeclare_array_length=10000000
Carray = setup_DEC_Cmat2(areas_list, states_list, Cparams)

"""

function setup_DEC_Cmat2(areas_list, states_list, maxent01, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); predeclare_array_length=Integer(min(length(states_list)*length(states_list)*round((length(states_list)/2)), 10000000)), min_precision=1e-9)
	numareas = length(areas_list)
	numstates = length(states_list)
	
	maxent01symp = maxent01.maxent01symp
	maxent01sub = maxent01.maxent01sub
	maxent01vic = maxent01.maxent01vic
	maxent01jump = maxent01.maxent01jump
	
	# Weights
	y_wt = Cparams.y
	s_wt = Cparams.s
	v_wt = Cparams.v
	j_wt = Cparams.j
	
	# MAKE SURE the states_list is sorted by size
	range_sizes = length.(states_list)
	rangesize_diffs = collect(repeat([0.0], numstates-1))
	for i in 2:length(range_sizes)
	 rangesize_diffs[i-1] = length(states_list[i]) - length(states_list[i-1]) + 0.0
	end
	TF = rangesize_diffs .< 0
	if sum(TF) > 0
		txt = "\nSTOP ERROR in setup_DEC_Cmat(): Your states_list is not ordered in ascending rangesize. Printing states_list:\n"
		print(txt)
		print("\nstates_list")
		print(states_list)
		error(txt)
	end
	
	
	# Pre-declare arrays
	# This contains the sum of the weights for each ancestor row
	row_weightvals = collect(repeat([0.0], numstates))
	
	# Int32 ranges from -32768 to 32767, this suggests a maximum number of states
	Carray_ivals = collect(repeat([0], predeclare_array_length))
	Carray_jvals = collect(repeat([0], predeclare_array_length))
	Carray_kvals = collect(repeat([0], predeclare_array_length))
	Cijk_weights = collect(repeat([0.0], predeclare_array_length))
	Carray_event_types = collect(repeat([""], predeclare_array_length))
	numC = 0 # counter of the number of allow cladogenesis events
	
	
	# Go through:
	# i = ancestor state index
	# j = left state index
	# k = right state index
	
	# Preallocate this vector ONCE, size = numareas * 2
	#tmp_merged_vec = repeat([0], 2*numareas)
	
	for i in 1:numstates
		ancstate = states_list[i]
		ancsize = length(ancstate)
		# Exclude the null range as an ancestor for cladogenetic range-changing events
		# (ancestor of range NULL give 0 likelihood to the observed data)
		if (ancsize == 0)
			continue
		end
		
		# Loop through left and right descendant states
		# Left states
		for j in 1:numstates
			lstate = states_list[j]
			lsize = length(lstate)
			if (lsize == 0)
				continue # go to the next j
			end
		
			# Right states
			for k in 1:numstates # We only have to do half of the possible right events;
													 # reverse each to make a pair
				rstate = states_list[k]
				rsize = length(rstate)
				if (rsize == 0)
					continue # go to the next k
				end
				
				if (y_wt > min_precision)
					# Sympatry (range-copying)
# 204:
					if (all([lsize == ancsize, rsize==ancsize, lstate==ancstate, rstate==ancstate]))

#2244					if (all([lstate==ancstate, rstate==ancstate]))
						# Check if the weight > 0.0
						# ancsize, lsize, rsize are the same so we don't have to 
						# choose the smaller daugher
						tmp_weightval = y_wt * maxent01symp[ancsize, lsize] * 1.0 * 1.0
						if tmp_weightval > min_precision
							# Record the range-copying event
							numC += 1
							Carray_event_types[numC] = "y"
							Carray_ivals[numC] = i
							Carray_jvals[numC] = j
							Carray_kvals[numC] = k
							Cijk_weights[numC] = tmp_weightval
							row_weightvals[i] += tmp_weightval
							continue
						end # end if tmp_weightval > 0.0
					end # end if (all([lsize == ancsize, rsize==ancsize, lstate==ancstate, rstate==ancstate])
				end # end if (y_wt > min_precision)
				
				# If one of the descendants is identical to the ancestor, 
				# (after we've excluded sympatry)
				# we can have jump dispersal or subset speciation
				if ( all([ancsize==rsize, ancstate==rstate]) )
					# Subset sympatry
					if (s_wt > min_precision)
						# Check for subset sympatry: lstate smaller than rstate, lstate inside rstate
						if ((lsize < rsize) && (array_in_array(lstate, rstate) == true) )
							# Check if the weight > 0.0
							# lsize is smaller by definition
							# choose the smaller daughter
							tmp_weightval = s_wt * maxent01sub[ancsize, lsize] * 1.0 * 1.0
							if tmp_weightval > min_precision
								# Record the range-copying event
								numC += 1
								Carray_event_types[numC] = "s"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = j
								Carray_kvals[numC] = k
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval

								# Same event, flip left/right descendant states
								numC += 1
								Carray_event_types[numC] = "s"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = k
								Carray_kvals[numC] = j
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval
								continue
							end # end if tmp_weightval > 0.0
						end # end if ((array_in_array(lstate, rstate) == true) && (lsize < rsize))
					end # end if (s_wt > min_precision)
					
					# Jump dispersal
					if (j_wt > min_precision)
						# If the left descendant is of size 1, & different from right, then
						# jump dispersal
						if ( (lsize == 1) && (array_in_array(lstate, rstate) == false) )
							# Historically, on analogy with other DEC cladogenesis events
							# the weight of each j event was not influenced by the number of
							# ranges of the ancestor. That is followed here. Obviously, 
							# other models could be imagined! (e.g., d events *are* influenced
							# by ancestor rangesize).
							try_jump_dispersal_based_on_dist = true
							normalize_by_number_of_dispersal_events = true
							jweight_for_cell_based_on_distances = 0.0
							if (try_jump_dispersal_based_on_dist == true)
								for anc_area in ancstate
									for left_area in lstate
										 jweight_for_cell_based_on_distances += dmat[anc_area,left_area]
									end
								end
								# Normalize by number of possible jump dispersals
								if (normalize_by_number_of_dispersal_events == true)
									jweight_for_cell_based_on_distances = jweight_for_cell_based_on_distances / (ancsize * lsize)
								end
							else
								# 
								jweight_for_cell_based_on_distances = 1.0
							end # end if (try_jump_dispersal_based_on_dist == true)
					
							# Calculate the final weight of this jump dispersal
							tmp_weightval = j_wt * maxent01jump[ancsize, lsize] * 1.0 * 1.0 * jweight_for_cell_based_on_distances
					
							print("\n")
							print([i, j, k, tmp_weightval])
							if (tmp_weightval > min_precision)
								print("yes")
								# Record the jump-dispersal event
								numC += 1
								Carray_event_types[numC] = "j"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = j
								Carray_kvals[numC] = k
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval

								# Same event, flip left/right descendant states
								numC += 1
								Carray_event_types[numC] = "j"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = k
								Carray_kvals[numC] = j
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval
								continue
							end # if (tmp_weightval > 0.0)
							# end of jump dispersal
						end # end if ( (lsize == 1) && (array_in_array(lstate, rstate) == false) )
					end # end if (j_wt > min_precision)
				end # end if ( (ancstate == rstate) )
			
				# Vicariance
				if (v_wt > min_precision)
					# Check if the combined vector equals the ancestor vector					
					#tmp_merged_vec .= repeat([0], 2*numareas)
					#combined_size = length(lstate)+length(rstate)
					
					#tmp_merged_vec[1:length(lstate)] = lstate
					#tmp_merged_vec[(length(lstate)+1):(length(lstate)+length(rstate))] = rstate
					#combined_vector = sort(tmp_merged_vec)
					if ( is_event_vicariance(ancstate, lstate, rstate) )
						smaller_range_size = min(lsize, rsize)
						tmp_weightval = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0
						if (tmp_weightval > min_precision)
							# Record the jump-dispersal event
							numC += 1
							Carray_event_types[numC] = "v"
							Carray_ivals[numC] = i
							Carray_jvals[numC] = j
							Carray_kvals[numC] = k
							Cijk_weights[numC] = tmp_weightval
							row_weightvals[i] += tmp_weightval
							continue
							# Same event, flip left/right descendant states
							# You won't hit it again, as k >= i
# 							numC += 1
# 							Carray_event_types[numC] = "v"
# 							Carray_ivals[numC] = i
# 							Carray_jvals[numC] = k
# 							Carray_kvals[numC] = j
# 							Cijk_weights[numC] = tmp_weightval
# 							row_weightvals[i] += tmp_weightval
						end # end if (tmp_weightval > 0.0)
					end # end if ( combined_vector == sort(ancstate) )
				end # end if ( (v>0.0) && ..
				# End of Vicariance
			end # end i (right state indices)
		end # end j (left state indices)
	end # end i (ancestor state indices)
	
	TF = Carray_event_types .!= ""
	Carray_ivals = Carray_ivals[TF]
	Carray_jvals = Carray_jvals[TF]
	Carray_kvals = Carray_kvals[TF]
	Cijk_weights = Cijk_weights[TF]
	Carray_event_types = Carray_event_types[TF]
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
	end

	
	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)
	
	"""
	# Extract the values
	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	row_weightvals
	"""
	
	return Carray
end # end setup_DEC_Cmat()










function calclike_DEC_SSE(tr, tipdata)
end


end # end module ModelLikes

