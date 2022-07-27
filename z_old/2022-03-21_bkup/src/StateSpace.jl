module StateSpace
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

using BioGeoJulia.MaxentInterp # for discrete_maxent_distrib_of_smaller_daughter_ranges
using BioGeoJulia.TrUtils # for e.g. flat2
using Combinatorics  # for e.g. combinations()
using DataFrames     # for e.g. DataFrame()
using PhyloNetworks
using Distributions  # for quantile

print("\nBioGeoJulia: loading StateSpace.jl")

export CparamsStructure, default_Cparams, sumy, sums, sumv, sumj, prtQ, prtQi, prtQp, prtC, prtCi, prtCp, numstates_from_numareas, areas_list_to_states_list, states_list_to_txt, get_default_inputs, run_model, setup_MuSSE, setup_DEC_DEmat, construct_BioGeoBEARS_model_object, update_Qij_vals, update_Qij_vals2!, p_Ds_v5_updater_v1, bmo_updater_v1, array_in_array, is_event_vicariance, setup_DEC_Cmat, setup_DEC_Cmat2, totals_prtC, setup_DEC_Cmat3, update_maxent01, update_Cijk_vals, update_Cijk_vals2_noUpdate, update_Cijk_vals2!, p_Ds_v5_updater_v1!, bmo_updater_v1!





# Cladogenetic parameter weights structure
# "mutable" means you can change the values referred to by the keys
mutable struct CparamsStructure
	y::Float64
	s::Float64
	v::Float64
	j::Float64
end

"""
Cparams = default_Cparams()
"""
function default_Cparams()
	y = 1.0
	s = 1.0
	v = 1.0
	j = 0.0
	Cparams = CparamsStructure(y, s, v, j)
end



function sumy(x)
	sum(x .== "y")
end

function sums(x)
	sum(x .== "s")
end

function sumv(x)
	sum(x .== "v")
end

function sumj(x)
	sum(x .== "j")
end



"""
# Print a Qarray to a data.frame
"""
function prtQ(Qarray)
	Qdf = DataFrame(event=Qarray.Qarray_event_types, i=Qarray.Qarray_ivals, j=Qarray.Qarray_jvals, val=Qarray.Qij_vals)
	return Qdf
end

"""
# Print a Qarray from inputs to a data.frame
"""
function prtQi(inputs)
	Qdf = DataFrame(event=inputs.p_Ds_v5.p_indices.Qarray_event_types, i=inputs.p_Ds_v5.p_indices.Qarray_ivals, j=inputs.p_Ds_v5.p_indices.Qarray_jvals, val=inputs.p_Ds_v5.params.Qij_vals)
	return Qdf
end


"""
# Print a p_Ds_v5 Qarray to a data.frame
"""
function prtQp(p_Ds_v5)
	Qdf = DataFrame(event=p_Ds_v5.p_indices.Qarray_event_types, i=p_Ds_v5.p_indices.Qarray_ivals, j=p_Ds_v5.p_indices.Qarray_jvals, val=p_Ds_v5.params.Qij_vals)
	return Qdf
end




"""
# Print a Carray to a data.frame
"""
function prtC(Carray)
	Cdf = DataFrame(event=Carray.Carray_event_types, i=Carray.Carray_ivals, j=Carray.Carray_jvals, k=Carray.Carray_kvals, 	pair=Carray.Carray_pair, wt=Carray.Cijk_weights, val=Carray.Cijk_vals)
	return Cdf
end

"""
# Print a Carray from inputs to a data.frame
"""
function prtCi(inputs)
	Cdf = DataFrame(event=inputs.p_Ds_v5.p_indices.Carray_event_types, i=inputs.p_Ds_v5.p_indices.Carray_ivals, j=inputs.p_Ds_v5.p_indices.Carray_jvals, k=inputs.p_Ds_v5.p_indices.Carray_kvals, pair=inputs.p_Ds_v5.p_indices.Carray_pair, wt=inputs.p_Ds_v5.params.Cijk_weights, val=inputs.p_Ds_v5.params.Cijk_vals)
	return Cdf
end

"""
# Print a p_Ds_v5 Carray to a data.frame
"""
function prtCp(p_Ds_v5)
	Cdf = DataFrame(event=p_Ds_v5.p_indices.Carray_event_types, i=p_Ds_v5.p_indices.Carray_ivals, j=p_Ds_v5.p_indices.Carray_jvals, k=p_Ds_v5.p_indices.Carray_kvals, pair=p_Ds_v5.p_indices.Carray_pair, wt=p_Ds_v5.params.Cijk_weights, val=p_Ds_v5.params.Cijk_vals)
	return Cdf
end




"""
# Print downpass "res" object likelihoods to an R-like table

"""


"""
numstates_from_numareas(3,3,false)
numstates_from_numareas(3,3,true)
numstates_from_numareas(10,1,false)
numstates_from_numareas(10,2,false)
numstates_from_numareas(10,3,false)
numstates_from_numareas(10,10,false)
numstates_from_numareas(10,10,true)
"""


"""
	numstates_from_numareas(numareas[, maxareas, include_null_range])

Calculates the number of possible states in the state space (e.g., the number of possible geographic 
ranges). The inputs are the number of areas (`numareas`), and the maximum range size (`maxareas`).

If `numareas` = `maxareas`, then the size of the state space is 2^`numareas`. This is
because, for each area, a species could be absent (0) or present (1). Therefore, if
there is one area, there are 2x1=2 possible ranges (0 and 1). If there are two areas, there 
are 2x2=4 possible ranges (00, 01, 10, 11). If there are three areas, there are 2x2x2=8 
possible ranges (000, 001, 010, 011, 100, 101, 110, 111), etc.

The `include_null_range` input, if set to `true` allows the "all absent" range (e.g. 000). 
If set to `false`, this range is disallowed, decreasing the size of the state space by 1.

*NOTE:* The size of the state space is a fundamental constraint on the computational speed
of likelihood calculations for biogeographical models using matrix exponentiation, ODEs 
(ordinary differential equations), etc. Researchers must think carefully about how large 
of a state space they require to test their hypotheses, and whether their analysis will 
run fast enough (or at all!).  If `numareas`=20 and `maxareas`=20, the size of the state
space is 1,048,576. Trying to calculate likelihoods for such a huge state space will likely
just fill up computer memory and then crash the program.

Researchers (and reviewers and editors!) should clearly recognize that any computational
inference is inherently a compromise between a very complex reality, and the memory and 
speed limitations of computers. Different people might reach different conclusions about
where exactly this compromise should end up. I tend to think that running somewhat simpler and 
faster models, thus allowing more time to run variant models, model checks, etc., is more 
valuable than setting up the most complex model you can think of, devoting months of 
computing time to it, and then either (a) losing all of that work when it crashes, or (b)
treating the output as gospel because you don't have the time or money available to do
anything else.


# Examples
```julia-repl
julia> numstates_from_numareas(3,3,false)
7

julia> numstates_from_numareas(3,3,true)
8

julia> numstates_from_numareas(10,1,false)
10

julia> numstates_from_numareas(10,2,false)
55

julia> numstates_from_numareas(10,3,false)
175

julia> numstates_from_numareas(10,10,false)
1023

julia> numstates_from_numareas(10,10,true)
1024

julia> numstates_from_numareas(20,20,true)
1048576
```
"""
function numstates_from_numareas(numareas=3, maxareas=1, include_null_range=false)
	# The formula for the number of geographic states, based on the number of areas,
	# is:
	#
	# sum_from_k=1_to_m (N choose k)
	# 
	numstates = 0
	
	# to avoid "does not accept keyword arguments"
	n = numareas + 0
	
	for k in 1:maxareas
		tmp_numstates = binomial(n, k)   # see also Rchoose
		numstates = numstates + tmp_numstates
	end
	
	if include_null_range == true
		numstates += 1
	end
	return numstates
end



# areas_list (1-based) to states_list (1-based)
"""
areas_list = collect(1:3)
states_list = areas_list_to_states_list(areas_list, 1, false)
states_list = areas_list_to_states_list(areas_list, 1, true)
states_list = areas_list_to_states_list(areas_list, 3, false)
states_list = areas_list_to_states_list(areas_list, 3, true)
areas_list_to_states_list()
"""


"""
	areas_list_to_states_list(areas_list[, maxareas, include_null_range])

Provides the list of possible states (e.g., geographic ranges) in the state space. The 
inputs are:

* `areas_list` - The list of areas. Each area is described with a number. This is done
with the `collect` function, e.g. collect(1:3).

* `maxareas` - The maximum number of areas occupied per geographic range. See
`numstates_from_numareas` for a discussion of how the state space grows (rapidly!) 
with `numareas` and `maxareas`.

* `include_null_range` - if set to `true`, allows the "all absent" range (e.g. 000). 
If set to `false`, this range is disallowed, decreasing the size of the state space by 1.

NOTE: The size of the state space is a fundamental constraint on the computational speed
of likelihood calculations for biogeographical models using matrix exponentiation, ODEs 
(ordinary differential equations), etc. Researchers must think carefully about how large 
of a state space they require to test their hypotheses, and whether their analysis will 
run fast enough (or at all!).  If `numareas`=20 and `maxareas`=20, the size of the state
space is 1,048,576. Trying to calculate likelihoods for such a huge state space will likely
just fill up computer memory and then crash the program.

Researchers (and reviewers and editors!) should clearly recognize that any computational
inference is inherently a compromise between a very complex reality, and the memory and 
speed limitations of computers. Different people might reach different conclusions about
where exactly this compromise should end up. I tend to think that running somewhat simpler and 
faster models, thus allowing more time to run variant models, model checks, etc., is more 
valuable than setting up the most complex model you can think of, devoting months of 
computing time to it, and then either (a) losing all of that work when it crashes, or (b)
treating the output as gospel because you don't have the time or money available to do
anything else.


# Examples
```julia-repl
julia> areas_list = collect(1:3)
3-element Array{Int64,1}:
 1
 2
 3

julia> states_list = areas_list_to_states_list(areas_list, 1, false)
3-element Array{Array{Any,1},1}:
 [1]
 [2]
 [3]

julia> states_list = areas_list_to_states_list(areas_list, 1, true)
4-element Array{Array{Any,1},1}:
 [] 
 [1]
 [2]
 [3]

julia> states_list = areas_list_to_states_list(areas_list, 3, false)
7-element Array{Array{Any,1},1}:
 [1]      
 [2]      
 [3]      
 [1, 2]   
 [1, 3]   
 [2, 3]   
 [1, 2, 3]

julia> states_list = areas_list_to_states_list(areas_list, 3, true)
8-element Array{Array{Any,1},1}:
 []       
 [1]      
 [2]      
 [3]      
 [1, 2]   
 [1, 3]   
 [2, 3]   
 [1, 2, 3]
```
"""
function areas_list_to_states_list(areas_list=collect(1:3), maxareas=3, include_null_range=false)
	
	# Initialize the states_list to the correct size
	numareas = length(areas_list)
	if maxareas > numareas
		maxareas = numareas
	end
	numstates = numstates_from_numareas(numareas, maxareas, include_null_range)
	# Empty list with numstates states
	states_list = repeat([[]], numstates)
	
	# Populate the states_list
	# If include_null_range=true, the first state is NULL, i.e. [] (empty list)
	if include_null_range == true
		state_num = 1
	else
		state_num = 0
	end
	
	# Fill in the states_list
	for k in 1:maxareas
		tmp_states_list = collect(Combinatorics.combinations(areas_list, k))
		for j in 1:length(tmp_states_list)
			state_num += 1
			states_list[state_num] = tmp_states_list[j]
		end
	end
	
	return states_list
end


"""
# Get text list of ranges

areas_list = collect(1:3)
area_names = ["A", "B", "C"]
states_list = areas_list_to_states_list(areas_list, 1, false)
states_list_to_txt(states_list, area_names; delim="")
states_list = areas_list_to_states_list(areas_list, 1, true)
states_list_to_txt(states_list, area_names; delim="")
states_list = areas_list_to_states_list(areas_list, 3, false)
states_list_to_txt(states_list, area_names; delim="")
states_list = areas_list_to_states_list(areas_list, 3, true)
states_list_to_txt(states_list, area_names; delim="")
"""
function states_list_to_txt(states_list, area_names; delim="")
	ranges_list = collect(repeat([""], length(states_list)))
	for i in 1:length(ranges_list)
		areas_in_this_state = area_names[states_list[i]]
		if (length(areas_in_this_state) == 0)
			areas_in_this_state = "_"
		end
		ranges_list[i] = join(areas_in_this_state, delim)
	end # END for i in 1:length(ranges_list)
	return(ranges_list)
end # function states_list_to_txt(states_list, area_names)


#######################################################
# Run models
#######################################################
function get_default_inputs(n=2)
	# Load tree
	great_ape_newick_string = "((chimp:1,human:1):1,gorilla:2);"
	tr = readTopology(great_ape_newick_string)
	rootnodenum = tr.root
	trdf = prt(tr, rootnodenum)
	#trdf
	
	# Set up a simple MuSSE model
	p_Es_v5 = setup_MuSSE(n; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
	
	# Solutions to the E vector
	u0_Es = repeat([0.0], 1*n)
	uE = repeat([0.0], n)
	tspan = (0.0, 1.2*trdf[tr.root,:node_age]) # 110% of tree root age

	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, u0_Es, tspan, p_Es_v5)
	sol_Es_v5 = solve(prob_Es_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);

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

	# Populate the tip likelihoods -- MODIFY THIS
	res.likes_at_each_nodeIndex_branchTop[1] = u0;
	res.likes_at_each_nodeIndex_branchTop[2] = u0;
	res.likes_at_each_nodeIndex_branchTop[4] = u0;

	solver_options = construct_SolverOpt()
	solver_options.solver=Tsit5()
	solver_options.abstol = 1.0e-6
	solver_options.reltol = 1.0e-6

	inputs = (res=res, trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
	return inputs
end


function run_model(inputs=default_inputs())
	res = inputs.res
	trdf = inputs.trdf
	p_Ds_v5 = inputs.p_Ds_v5
	solver_options = inputs.solver_options
	(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res, trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10)
	return (total_calctime_in_sec, iteration_number)
	# return res?
end








#######################################################
# Set up models
#######################################################

# Set up a MuSSE (or BiSSE) model with n states
# This default model assumes:
#   - same birthRate & deathRate across all states
#   - transitions only possible to adjacent states
# 
"""
p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
p_Es_v5 = setup_MuSSE(3; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
p_Es_v5 = setup_MuSSE(4, birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)

# Anagenetic transition matrix
hcat(p_Es_v5.p_indices.Qarray_ivals, p_Es_v5.p_indices.Qarray_jvals, p_Es_v5.params.Qij_vals)
DataFrame(p_Es_v5.p_TFs.Qi_eq_i)
p_Es_v5.p_TFs.Qi_sub_i
p_Es_v5.p_TFs.Qj_sub_i

# Cladogenetic transition matrix
hcat(p_Es_v5.p_indices.Carray_ivals, p_Es_v5.p_indices.Carray_jvals, p_Es_v5.p_indices.Carray_kvals, p_Es_v5.params.Cijk_vals)
DataFrame(p_Es_v5.p_TFs.Ci_eq_i)
p_Es_v5.p_TFs.Ci_sub_i
p_Es_v5.p_TFs.Cj_sub_i
p_Es_v5.p_TFs.Ck_sub_i
"""
function setup_MuSSE(n=2; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
	# Define Qarray - zeros
	Qarray_ivals = collect(1:(n-1))
	Qarray_jvals = collect(2:n)
	Qarray_ivals = vcat(Qarray_ivals, collect(2:n))
	Qarray_jvals = vcat(Qarray_jvals, collect(1:(n-1)))
	Qij_vals = vcat(repeat([q01],(n-1)), repeat([q10],(n-1)))
	
	# Carray: Cladogenetic parameters
	# column 1: state i
	# column 2: state j
	# column 3: state k
	# column 4: lambda_ijk (a parameter that is at least possibly nonzero, under the given model)
	Carray_ivals = collect(1:n)
	Carray_jvals = collect(1:n)
	Carray_kvals = collect(1:n)
	Carray_pair = repeat([1], n)
	Cijk_vals = repeat([birthRate], n)

# 	Qij_vals[((Qarray_ivals .== 1) .+ (Qarray_jvals .!= 1)) .== 2]
# 	hcat(Qarray_ivals, Qarray_jvals, Qij_vals)
# 	hcat(Carray_ivals, Carray_jvals, Carray_kvals, Cijk_vals)
	
	# Extinction rates
	mu_vals = repeat([deathRate], n)
	
	# Assemble a "params" tuple
	
	# Possibly varying parameters
	params = (mu_vals=mu_vals, Qij_vals=Qij_vals, Cijk_vals=Cijk_vals)

	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	p_indices = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Carray_pair=Carray_pair)

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
		push!(Qi_eq_i, Qarray_ivals .== i)
		push!(Qi_sub_i, Qarray_ivals[Qarray_ivals .== i])
		push!(Qj_sub_i, Qarray_jvals[Qarray_ivals .== i])

		push!(Ci_eq_i, Carray_ivals .== i)
		push!(Ci_sub_i, Carray_ivals[Carray_ivals .== i])
		push!(Cj_sub_i, Carray_jvals[Carray_ivals .== i])
		push!(Ck_sub_i, Carray_kvals[Carray_ivals .== i])
	end
	
	# Inputs to the Es calculation
	p_TFs = (Qi_eq_i=Qi_eq_i, Ci_eq_i=Ci_eq_i, Qi_sub_i=Qi_sub_i, Qj_sub_i=Qj_sub_i, Ci_sub_i=Ci_sub_i, Cj_sub_i=Cj_sub_i, Ck_sub_i=Ck_sub_i)
	p_Es_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs)
	
	tmptxt="""
	n = p_Es_v5.n
	params = p_Es_v5.params
	p_indices = p_Es_v5.p_indices
	p_TFs = p_Es_v5.p_TFs
	"""
	
	return(p_Es_v5)
end








#######################################################
# Set up a sparse Qmat for the DEC model
#######################################################
# It will contain references to the parameters
# dmat: a numareas x numareas matrix of "d" values (range-expansion dispersal)
# elist: a numareas array of "e" values (range-contraction, extirpation, "local extinction")
# amat: a numareas x numareas matrix of "a" values
### mult_mat: a numareas x numareas matrix of dispersal multipliers
### e_mult: a numareas list of extirpation multipliers
### exclude_zeros if true, exclude from the matrix 

"""
numareas = 3
areas_list = collect(1:numareas)
states_list = areas_list_to_states_list(areas_list, 3, true)
numstates = length(states_list)
amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
elist = repeat([0.123], numstates)
allowed_event_types=["d","e"]

Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
Qarray_event_types = Qmat.Qarray_event_types
hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)


states_list = areas_list_to_states_list(areas_list, 3, false)
Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
Qarray_event_types = Qmat.Qarray_event_types
hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

# An "a" matrix for single areas - 2 areas
states_list = areas_list_to_states_list(areas_list, 1, false)
Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["a"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
Qarray_event_types = Qmat.Qarray_event_types
hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

# An "a" matrix for single areas - 5 areas

areas_list = collect(1:5)
states_list = areas_list_to_states_list(areas_list, 1, false)
numareas = length(areas_list)
numstates = length(states_list)
amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
elist = repeat([0.123], numstates)
Qmat = StateSpace.setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["a"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
Qarray_event_types = Qmat.Qarray_event_types
hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)


"""

# states_list=areas_list_to_states_list(areas_list, length(areas_list), true), dmat=reshape(repeat([0.1], numstates),(length(areas_list),length(areas_list))), elist=repeat([0.01],length(areas_list)), amat=reshape(repeat([0.0], numstates),(length(areas_list),length(areas_list))); allowed_event_types=["d","e"]

function setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	# Set up items to iterate over
	numstates = length(states_list)
	statenums = collect(1:numstates)
	range_sizes = length.(states_list)
	#areas_list = sort(unique(flat2(states_list)))
	numareas = length(areas_list)
	
	# Get the approx size of the nonzero rates (for pre-allocation)
	if (in("e", allowed_event_types))
		num_e_rates = numstates
	else
		num_e_rates = 0
	end

	
	# Count the approximate number of d events by
	# iterating counts upwards (assumes states_list is size-ordered)
	#num_d_rates = ceil((numstates^2)/2)
	if (in("d", allowed_event_types))
		num_d_rates = 0
		lengths_states_list = length.(states_list)
		for i in 1:(length(states_list)-1)
			for j in (i+1):length(states_list)
				if ((lengths_states_list[i]+1) == lengths_states_list[j])
					num_d_rates += 1
				end
			end
		end
	else
		num_d_rates = 0
	end

	# Single-area transitions
	if (in("a", allowed_event_types))
		num_a_rates = (numareas^2) - numareas
	else
		num_a_rates = 0
	end
	
	num_nonzero_rates = num_e_rates + num_d_rates + num_a_rates
	
	# Initialize empty arrays
	Qarray_ivals = repeat([0], num_nonzero_rates)
	Qarray_jvals = repeat([0], num_nonzero_rates)
	Qij_vals =  repeat([0.0], num_nonzero_rates)  # 0-element Array{Any,1}  
	                     # This is populated by calculating through the others
	#base_vals = Array{Float64, num_nonzero_rates}  # base rates -- d, e, a
	#mod_vals = Array{Float64, num_nonzero_rates}  # base modifiers -- e.g., "2" for AB -> ABC
	#mult_vals = Array{Float64, num_nonzero_rates} # multipliers from mult_mat, e_mult 
	# (ie modification by distance, area, etc.)
	Qarray_event_types = repeat([""], num_nonzero_rates)
	index = 0
	
	
	# Events of "a" type: anagenetic range-switching
	# restricted to single-area ranges (!)
	if (in("a", allowed_event_types))
		rangesize_eq1_TF = range_sizes .== 1
		statenums_of_size1 = statenums[rangesize_eq1_TF]
		
		for i in 1:(length(statenums_of_size1)-1)			# starting states
			for j in (i+1):length(statenums_of_size1)		# ending states
				statenum_ival = statenums_of_size1[i]
				statenum_jval = statenums_of_size1[j]
				starting_state = states_list[statenum_ival]
				ending_state = states_list[statenum_jval]
				size_i = length(starting_state)
				size_j = length(ending_state)
			
				# "a" events -- anagenetic transitions between single-area states
				# Only use a's, if a>0, or 0 a's are desired
				# Do "amat[i,j][]" instead of "amat[i,j]" in case it's a Ref()
				if (starting_state != ending_state) && (size_i == 1) && (size_j == 1) # single-areas, not null
					starting_areanum = starting_state[] # only works with length 1
					ending_areanum = ending_state[]     # only works with length 1
					if (amat[starting_areanum, ending_areanum][] > 0) || (exclude_zeros == false)
						# Add to arrays
						# Forward event
						index += 1
						Qarray_event_types[index] = "a"
						
						
						Qarray_ivals[index] = statenum_ival
						Qarray_jvals[index] = statenum_jval
						Qij_vals[index] = amat[starting_areanum,ending_areanum]
						
						# Reverse event
						index += 1
						Qarray_event_types[index] = "a"
						Qarray_ivals[index] = statenum_jval
						Qarray_jvals[index] = statenum_ival
						Qij_vals[index] = amat[ending_areanum, starting_areanum]

	# 					push!(Qarray_event_types, "a")
	# 					push!(Qarray_ivals, i)
	# 					push!(Qarray_jvals, j)
	# 					push!(base_vals, amat[starting_areanum,ending_areanum])
	# 					push!(mod_vals, 1)
	# 					push!(mult_vals, 1)
					end # end if (amat
				end # ending if a[] > 0
			end # ending for j in (i+1):length(statenums_of_size1
		end # ending for i in 1:(length(statenums_of_size1)-1)
	end # ending if (in("a", allowed_event_types))
	
	
	# Events of "d" type: anagenetic range-expansion dispersal
	if (in("d", allowed_event_types))
		for i in 1:(numstates-1)			# starting states
			for j in (i+1):numstates		# ending states
				starting_state = states_list[i]
				ending_state = states_list[j]
				size_i = length(starting_state)
				size_j = length(ending_state)
						
				# "d" events -- anagenetic range-expansion events
				# Is the ending range 1 area more than the starting range?
				if (starting_state != ending_state) && ((size_i+1) == size_j) && (size_i != 0) # state i is 1 smaller; not null
					starting_areanums = starting_state
					ending_areanums = ending_state
					end_areanums_not_found_in_start_areas = setdiff(ending_areanums, starting_areanums)
					if length(end_areanums_not_found_in_start_areas) == 1
						# Add to arrays
						index += 1
						Qarray_event_types[index] = "d"
						Qarray_ivals[index] = i
						Qarray_jvals[index] = j
				
						# Add up the d events
						tmp_d_sum = 0.0
						for k in 1:size_i
							# Because there is only 1 end_areanums_not_found_in_start_areas
							tmp_d_sum += dmat[starting_areanums[k], end_areanums_not_found_in_start_areas[1]][]
						end

						Qij_vals[index] = tmp_d_sum
	# 					push!(Qarray_event_types, "d")
	# 					push!(Qarray_ivals, i)
	# 					push!(Qarray_jvals, j)
	# 					push!(base_vals, tmp_d_sum)
	# 					push!(mod_vals, size_i)
	# 					push!(mult_vals, 1)
						
					end # ending if length(end_areanums_not_found_in_start_areas) == 1
				end # ending if (starting_state != ending_state)...
			end # ending j loop
		end # ending i loop
	end # ending if (in("d", allowed_event_types)
		
	# Events of "e" type: anagenetic range-loss/extirpation
	# NOTE: we could combine with "d" with some effort
	if (in("e", allowed_event_types))
		for i in 2:numstates			# starting states
			for j in 1:(i-1)		# ending states
				starting_state = states_list[i]
				ending_state = states_list[j]
				size_i = length(starting_state)
				size_j = length(ending_state)

				if (starting_state != ending_state) && ((size_i-1) == size_j) && (size_i != 0) # state i is 1 bigger; not null
					starting_areanums = starting_state
					ending_areanums = ending_state
					start_areanums_not_found_in_end_areas = setdiff(starting_areanums, ending_areanums)
					if length(start_areanums_not_found_in_end_areas) == 1
						# Add to arrays
						index += 1
						Qarray_event_types[index] = "e"
						Qarray_ivals[index] = i
						Qarray_jvals[index] = j
						# Because there is only 1 area in start_areanums_not_found_in_end_areas
						Qij_vals[index] = elist[start_areanums_not_found_in_end_areas[1]]
					end # ending if length(start_areanums_not...
				end # ending if (starting_state != ending_state)...
			end # ending j loop
		end # ending i loop
	end # ending if (in("e", allowed_event_types)
	
	keepTF = Qarray_event_types .!= ""
	Qarray_ivals = Qarray_ivals[keepTF]
	Qarray_jvals = Qarray_jvals[keepTF]
	Qij_vals = Qij_vals[keepTF]
	Qarray_event_types = Qarray_event_types[keepTF]
	
	# Return results
	Qmat = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals, Qarray_event_types=Qarray_event_types)
	
	"""
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	Qarray_event_types = Qmat.Qarray_event_types
	
	hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)
	"""
	
	return Qmat
end # end setup_DEC_DEmat()


#######################################################
# bmo: A BioGeoBEARS_model_object DataFrame
# Much like BioGeoBEARS, this will contain the 
# higher-level parameters that are then used to 
# set the individual rates.
# 
# A DataFrame.
#######################################################
# defaults
"""
bmo = construct_BioGeoBEARS_model_object()
"""
function construct_BioGeoBEARS_model_object()
	type_vec = ["free", "free", "fixed", "fixed", "fixed", "fixed", "fixed", "fixed", "fixed", "3-j", "ysv*2/3", "ysv*1/3", "ysv*1/3", "ysv*1/3", "fixed", "fixed", "fixed", "mx01", "mx01", "mx01", "mx01", "fixed", "fixed", "fixed", "fixed"]
	init_vec = [0.01, 0.01, 0, 1, 0, 0, 1, 0, 0, 2.99999, 1.99999, 1, 1, 1, 0.3288164, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.1, 1, 0]
	#init_vec = [0.01, 0.01, 0, 1, 0, 0, 1, 0, 0, 2.99999, 1.99999, 1, 1, 1, 0.3288164, 0.0, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 0.5, 0.1, 1, 0]
	min_vec = [1e-12, 1e-12, 1e-12, 1e-12, -2.5, -10, -10, -10, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 0.0, 0.0, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 0.005, 0.005, 0.005]
	max_vec = [4.999999999999, 4.999999999999, 4.999999999999, 0.999999999999, 2.5, 10, 10, 10, 2.99999, 3, 2, 1, 1, 1, 100.0, 100.0, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.995, 0.995, 0.995]
	#est_vec = [0.01, 0.01, 0, 1, 0, 0, 1, 0, 0, 2.99999, 1.99999, 1, 1, 1, 0.3288164, 0.0, 1e-04, 1e-04, 1e-04, 1e-04, 1e-04, 0.5, 0.1, 1, 0]
	est_vec = [0.01, 0.01, 0, 1, 0, 0, 1, 0, 0, 2.99999, 1.99999, 1, 1, 1, 0.3288164, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.1, 1, 0]
		 
	note_vec = ["works", "works", "works", "non-stratified only", "works", "works", "works", "works", "works", "works", "works", "works", "works", "works", "works", "works", "works", "works", "works", "works", "works", "no", "yes", "yes", "yes"]
	desc_vec = ["anagenesis: rate of 'dispersal' (range expansion)", "anagenesis: rate of 'extinction' (range contraction)", "anagenesis: rate of range-switching (i.e. for a standard char.)", "anagenesis: exponent on branch lengths", "exponent on distance (modifies d, j, a)", "exponent on environmental distance (modifies d, j, a)", "exponent on manual dispersal multipliers (modifies d, j, a)", "anagenesis: exponent on extinction risk with area (modifies e)", "cladogenesis: relative per-event weight of jump dispersal", "cladogenesis: y+s+v", "cladogenesis: y+s", "cladogenesis: relative per-event weight of sympatry (range-copying)", "cladogenesis: relative per-event weight of subset speciation", "cladogenesis: relative per-event weight of vicariant speciation", "speciation rate", "extinction rate (lineages)", "cladogenesis: controls range size of smaller daughter", "cladogenesis: controls range size of smaller daughter", "cladogenesis: controls range size of smaller daughter", "cladogenesis: controls range size of smaller daughter", "cladogenesis: controls range size of smaller daughter", "root: controls range size probabilities of root", "mean frequency of truly sampling OTU of interest", "detection probability per true sample of OTU of interest", "false detection of OTU probability per true taphonomic control sample"]
	
	rownames = ["d", "e", "a", "b", "x", "n", "w", "u", "j", "ysv", "ys", "y", "s", "v", "birthRate", "deathRate", "mx01", "mx01j", "mx01y", "mx01s", "mx01v", "mx01r", "mf", "dp", "fdp"]
	# Load into a DataFrame
	bmo = DataFrames.DataFrame(rownames=rownames, type=type_vec, init=init_vec, min=min_vec, max=max_vec, est=est_vec, note=note_vec, desc=desc_vec)
	return bmo
end




"""
# Update Qij_vals
# Takes dmat, amat, elist as inputs
# dmat is d * various dispersal multipliers
# amat is a * various dispersal multipliers
# elist is "e" for each area

# Update Qij_vals
numareas = 3
areas_list = collect(1:numareas)
states_list = areas_list_to_states_list(areas_list, 3, true)
numstates = length(states_list)
amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
elist = repeat([0.123], numstates)
allowed_event_types=["d","e"]

Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
Qarray_event_types = Qmat.Qarray_event_types
Qmat1_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

# Update!
dmat = reshape(repeat([0.5], numareas^2), (numareas,numareas))
Qmat2 = update_Qij_vals(Qmat, areas_list, states_list, dmat, elist, amat )
Qmat2

Qarray_ivals = Qmat2.Qarray_ivals
Qarray_jvals = Qmat2.Qarray_jvals
Qij_vals = Qmat2.Qij_vals
Qarray_event_types = Qmat2.Qarray_event_types
Qmat2_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

Qmat1_df
Qmat2_df

"""
function update_Qij_vals(Qmat, areas_list, states_list, dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))), elist=repeat([1.0], length(areas_list)), amat=dmat; return_df=false)

	numstates = length(states_list)
	statenums = collect(1:numstates)
	range_sizes = length.(states_list)
	#areas_list = sort(unique(flat2(states_list)))
	numareas = length(areas_list)
	
	if (Rclass(Qmat) == "DataFrame")
		Qarray_ivals = Qmat[!,:i]
		Qarray_jvals = Qmat[!,:j]
		Qij_vals = Qmat[!,:val]
		Qarray_event_types = Qmat[!,:event]
	else
		Qarray_ivals = Qmat.Qarray_ivals
		Qarray_jvals = Qmat.Qarray_jvals
		Qij_vals = Qmat.Qij_vals
		Qarray_event_types = Qmat.Qarray_event_types
	end
	
	# Update the "d" events (anagenetic range expansion)
	TF = Qarray_event_types .== "d"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			end_areanums_not_found_in_start_areas = setdiff(ending_areanums, starting_areanums)
			
			# Add up the d events
			tmp_d_sum = 0.0
			for k in 1:ancsize
				# Because there is only 1 end_areanums_not_found_in_start_areas
				tmp_d_sum += dmat[starting_areanums[k], end_areanums_not_found_in_start_areas[1]][]
			end
			# Store the result
			rates[z] = tmp_d_sum
		end
		Qij_vals[TF] = rates
	end # End update of d event weights


	# Update the "a" events (anagenetic range expansion)
	TF = Qarray_event_types .== "a"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanum = ancstate # because it has size=1 by def.
			ending_areanum = decstate   # because it has size=1 by def.

			# Store the result
			rates[z] = amat[starting_areanum,ending_areanum]
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Update the "e" events (anagenetic range extinction/
	# local extirpation/range loss)
	TF = Qarray_event_types .== "e"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			start_areanums_not_found_in_end_areas = setdiff(starting_areanums, ending_areanums)

			# Store the result
			rates[z] = elist[start_areanums_not_found_in_end_areas[1]]
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Return results
	if return_df == true
		Qmat2 = DataFrames.DataFrame(Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals, Qarray_event_types=Qarray_event_types)
	else
		Qmat2 = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals, Qarray_event_types=Qarray_event_types)
	end
	return Qmat2
end # end function update_Qij_vals




"""
# Update Qij_vals2!
# Takes dmat, amat, elist as inputs
# dmat is d * various dispersal multipliers
# amat is a * various dispersal multipliers
# elist is "e" for each area

# Update Qij_vals2
numareas = 3
areas_list = collect(1:numareas)
states_list = areas_list_to_states_list(areas_list, 3, true)
numstates = length(states_list)
amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
elist = repeat([0.123], numstates)
allowed_event_types=["d","e"]

Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
Qarray_ivals = Qmat.Qarray_ivals
Qarray_jvals = Qmat.Qarray_jvals
Qij_vals = Qmat.Qij_vals
Qarray_event_types = Qmat.Qarray_event_types
Qmat1_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

# Update!
dmat = reshape(repeat([0.5], numareas^2), (numareas,numareas))
Qmat2 = update_Qij_vals2!(Qmat, areas_list, states_list, dmat, elist, amat )
Qmat2

Qarray_ivals = Qmat2.Qarray_ivals
Qarray_jvals = Qmat2.Qarray_jvals
Qij_vals = Qmat2.Qij_vals
Qarray_event_types = Qmat2.Qarray_event_types
Qmat2_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)

Qmat1_df
Qmat2_df

"""
function update_Qij_vals2!(p_Ds_v5, areas_list, states_list, dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))), elist=repeat([1.0], length(areas_list)), amat=dmat; return_df=false)

	numstates = length(states_list)
	statenums = collect(1:numstates)
	range_sizes = length.(states_list)
	#areas_list = sort(unique(flat2(states_list)))
	numareas = length(areas_list)
	
	Qarray_ivals = p_Ds_v5.p_indices.Qarray_ivals
	Qarray_jvals = p_Ds_v5.p_indices.Qarray_jvals
	Qij_vals = p_Ds_v5.params.Qij_vals
	Qarray_event_types = p_Ds_v5.p_indices.Qarray_event_types
	
	# Update the "d" events (anagenetic range expansion)
	TF = Qarray_event_types .== "d"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			end_areanums_not_found_in_start_areas = setdiff(ending_areanums, starting_areanums)
			
			# Add up the d events
			tmp_d_sum = 0.0
			for k in 1:ancsize
				# Because there is only 1 end_areanums_not_found_in_start_areas -- 
				# (all d dispersals add just 1 area)
				tmp_d_sum += dmat[starting_areanums[k], end_areanums_not_found_in_start_areas[1]][]
			end
			# Store the result
			rates[z] = tmp_d_sum
		end
		Qij_vals[TF] = rates
	end # End update of d event weights


	# Update the "a" events (anagenetic range expansion)
	TF = Qarray_event_types .== "a"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanum = ancstate # because it has size=1 by def.
			ending_areanum = decstate   # because it has size=1 by def.

			# Store the result
			rates[z] = amat[starting_areanum,ending_areanum]
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Update the "e" events (anagenetic range extinction/
	# local extirpation/range loss)
	TF = Qarray_event_types .== "e"
	if (sum(TF) > 0)
		ivals = Qarray_ivals[TF]
		jvals = Qarray_jvals[TF]
		rates = Qij_vals[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			decstate = states_list[j]
			decsize = length(decstate)
			
			starting_areanums = ancstate
			ending_areanums = decstate
			# Because only 1 area is lost at a time under this model
			start_areanums_not_found_in_end_areas = setdiff(starting_areanums, ending_areanums)

			# Store the result
			rates[z] = elist[start_areanums_not_found_in_end_areas[1]]
		end
		Qij_vals[TF] = rates
	end # End update of a event weights

	# Return results
	if return_df == true
		Qmat2 = DataFrames.DataFrame(Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qij_vals=Qij_vals, Qarray_event_types=Qarray_event_types)
		return Qmat2
	else
		p_Ds_v5.params.Qij_vals[:] .= Qij_vals
		return p_Ds_v5
	end
	return NaN  # shouldn't get here
end # end function update_Qij_vals2!

"""
# Update the BioGeoBEARS_model_object after changing e.g. "j"
"""
function bmo_updater_v1(bmo)
	j_wt = bmo.est[bmo.rownames .== "j"][1]
	ysv = bmo.est[bmo.rownames .== "ysv"][1]
	ys = bmo.est[bmo.rownames .== "ys"][1]
	y = bmo.est[bmo.rownames .== "y"][1]
	s = bmo.est[bmo.rownames .== "s"][1]
	v = bmo.est[bmo.rownames .== "v"][1]
	
	# Update
	ysv_func = bmo.type[bmo.rownames .== "ysv"][1]
	if ysv_func == "3-j"
		ysv = 3-j_wt
	end
	if ysv_func == "2-j"
		ysv = 2-j_wt
	end
	if ysv_func == "1-j"
		ysv = 1-j_wt
	end
	
	ys = ysv*2/3
	y = ysv*1/3
	s = ysv*1/3
	v = ysv*1/3
	
	bmo.est[bmo.rownames .== "ysv"] .= ysv
	bmo.est[bmo.rownames .== "ys"] .= ys
	bmo.est[bmo.rownames .== "y"] .= y
	bmo.est[bmo.rownames .== "s"] .= s
	bmo.est[bmo.rownames .== "v"] .= v
	
	return bmo.est
end # END bmo_updater_v1



"""
# p_Ds_v5 updater #1
# p_Ds_v5 contains the actual rates that go straight into the lnL calculations

pars = [0.03505038, 0.02832370]
parnames = ["d", "e"]

# Assuming parameters are in the order of the bmo list of "free" parameters,
# 1. update the bmo
# 2. then update the Qij and Cijk arrays in p_Ds_v5 (the rate parameters arrays)
inputs.bmo.rownames[inputs.bmo.type.=="free"]

pars = [0.03505038, 0.02832370]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
output = p_Ds_v5_updater_v1(p_Ds_v5, inputs);  # WORKS 2022-03-10
p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)

pars = [0.01, 0.001]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
output = p_Ds_v5_updater_v1(p_Ds_v5, inputs);  # WORKS 2022-03-10
p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)

pars = [0.1, 0.2]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
output = p_Ds_v5_updater_v1(p_Ds_v5, inputs);  # WORKS 2022-03-10
p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)


lnL = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func(pars)
"""
function p_Ds_v5_updater_v1(p_Ds_v5, inputs; check_if_free_params_in_mat=true)
	
	if check_if_free_params_in_mat == true
		free_param_names = inputs.bmo.rownames[inputs.bmo.type .== "free"]
		# Are all of the free param names %in% the p_Ds_v5
		#TF = in([unique(p_Ds_v5.p_indices.Qarray_event_types); unique(p_Ds_v5.p_indices.Carray_event_types)]).(free_param_names)
		TF1 = in(free_param_names).(["j"])[1]
		TF2 = in(unique(p_Ds_v5.p_indices.Carray_event_types)).(["j"])[1]
		if TF1 != TF2
			txt = paste(["ERROR in p_Ds_v5_updater_v1(): j's must be in both places. j_in_bmo:", TF1, "; j_in_params:", TF2])
			print("\n")
			print(txt)
			print("\n")
			throw(txt)
		end # END if TF1 != TF2
	end # END if check_if_free_params_in_mat == true
	
	# Extract the parameters
	d = inputs.bmo.est[inputs.bmo.rownames .== "d"][1]
	a = inputs.bmo.est[inputs.bmo.rownames .== "a"][1]
	e = inputs.bmo.est[inputs.bmo.rownames .== "e"][1]
	x = inputs.bmo.est[inputs.bmo.rownames .== "x"][1]
	w = inputs.bmo.est[inputs.bmo.rownames .== "w"][1]
	n = inputs.bmo.est[inputs.bmo.rownames .== "n"][1]
	x2 = 1.0
	x3 = 1.0 
	u = inputs.bmo.est[inputs.bmo.rownames .== "u"][1]
	
	# Update the dmat and elist
	inputs.setup.dmat_base .= d
	inputs.setup.amat_base .= a
	inputs.setup.elist_base .= e
	
	# Create the final dmat (dispersal rate multipliers) by multiplying the various input dmats
	inputs.setup.dmat .= inputs.setup.dmat_base .* inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	inputs.setup.amat .= inputs.setup.amat_base .* inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	inputs.setup.elist .= inputs.setup.elist_base .* inputs.setup.area_of_areas.^u
	
	# jmat does not include the "d" parameter
	inputs.setup.jmat .= inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	
	# Update the mus
	p_Ds_v5.params.mu_vals[:] .= inputs.bmo.est[inputs.bmo.rownames .== "deathRate"][1]
	
	# Now update the p_Ds_v5 (the rates) for Q matrix
	p_Ds_v5 = update_Qij_vals2!(p_Ds_v5, inputs.setup.areas_list, inputs.setup.states_list, inputs.setup.dmat, inputs.setup.elist, inputs.setup.amat; return_df=false);
	
	#print("\n")
	#print(prtQp(p_Ds_v5))
	
	# Now update the p_Ds_v5 (the rates) for C matrix
	pdf_before = prtCp(p_Ds_v5)
	print("\np_Ds_v5_updater_v1, before:")
	print(round.(pdf_before.wt[1:5], digits=4))
	output = update_Cijk_vals2(p_Ds_v5, inputs.setup.areas_list, inputs.setup.states_list, inputs.bmo, inputs.setup.maxent01, inputs.setup.jmat);
	p_Ds_v5.params.Cijk_weights[:] .= output.Cijk_weights
	p_Ds_v5.params.Cijk_vals[:] .= output.Cijk_vals
	p_Ds_v5.params.row_weightvals[:] .= output.row_weightvals
	pdf_after = prtCp(p_Ds_v5)
	print("\np_Ds_v5_updater_v1, after:")
	print(round.(pdf_before.wt[1:5], digits=4))
		
	return output
end # END function p_Ds_v5_updater_v1()








#######################################################
# Set up a Cmat (cladogenesis event-weights matrix)
#######################################################

function array_in_array(items, collection)
	TF = collect(repeat([false], length(items)))
	for i in 1:length(items)
		TF[i] = in(items[i], collection)
	end
	if (sum(TF) == length(TF))
		return true
	else
		return false
	end
end






"""
ancstate = [1, 2, 3,4];
lstate = [1, 2];
rstate = [4];
is_event_vicariance(ancstate, lstate, rstate)


ancstate = [1, 2, 3,4];
lstate = [1, 2];
rstate = [2, 4];
is_event_vicariance(ancstate, lstate, rstate)

ancstate = [1, 2, 3,4];
lstate = [1, 2];
rstate = [3, 4];
is_event_vicariance(ancstate, lstate, rstate)

"""


function is_event_vicariance(ancstate, lstate, rstate)
	ancsize = length(ancstate)
	lsize = length(lstate)
	rsize = length(rstate)
	
	if (ancsize != lsize+rsize)
		return false
	end
	
	# If any lstate areas are in rstate, then not vicariance
	for i in 1:lsize
		if in(lstate[i], rstate)
			return false
		end
	end
	
	# If any lstate areas not in ancstate, then not vicariance
	for i in 1:lsize
		if in(lstate[i], ancstate) == false
			return false
		end
	end

	# If any rstate areas not in ancstate, then not vicariance
	for i in 1:rsize
		if in(rstate[i], ancstate) == false
			return false
		end
	end
	
	# Otherwise, vicariance
	return true
end # END function is_event_vicariance(ancstate, lstate, rstate)




"""
areas_list = [1,2,3]
states_list = areas_list_to_states_list(areas_list, 3, true)
predeclare_array_length=10000000
Carray = setup_DEC_Cmat(areas_list, states_list);
prtC(Carray)

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
Carray = setup_DEC_Cmat(areas_list, states_list, maxent01, Cparams)
prtC(Carray)

"""

function setup_DEC_Cmat(areas_list, states_list, maxent01=NaN, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); predeclare_array_length=Integer(min(length(states_list)*length(states_list)*round((length(states_list)/2)), 10000000)), min_precision=1e-9)
	numareas = length(areas_list)
	total_numareas = numareas
	numstates = length(states_list)
	
	# Check if max_range_size=NaN
	type_string = string(typeof(maxent01))
	if (startswith(type_string, "NamedTuple") == false) && (isnan(maxent01) == true)
		maxent01 = maxent01_defaults(total_numareas)
	end
	
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
	Carray_pair = collect(repeat([1], predeclare_array_length))  # "pair" indicates if there are 1 or 2 such events
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
							Carray_pair[numC] = 1
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
								Carray_pair[numC] = 1
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval

								# Same event, flip left/right descendant states
								numC += 1
								Carray_event_types[numC] = "s"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = k
								Carray_kvals[numC] = j
								Carray_pair[numC] = 1
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
					
							#print("\n")
							#print([i, j, k, tmp_weightval])
							if (tmp_weightval > min_precision)
								#print("yes")
								# Record the jump-dispersal event
								numC += 1
								Carray_event_types[numC] = "j"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = j
								Carray_kvals[numC] = k
								Carray_pair[numC] = 1
								Cijk_weights[numC] = tmp_weightval
								row_weightvals[i] += tmp_weightval

								# Same event, flip left/right descendant states
								numC += 1
								Carray_event_types[numC] = "j"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = k
								Carray_kvals[numC] = j
								Carray_pair[numC] = 1
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
							Carray_pair[numC] = 1
							Cijk_weights[numC] = tmp_weightval
							row_weightvals[i] += tmp_weightval
							continue
							# Same event, flip left/right descendant states
							# You won't hit it again, as k >= i
							# This makes sense in setup_DEC_Cmat(), as k is 1:numstates
							# But in setup_DEC_Cmat3(), k is j:numstates, so you need to flip
# 							numC += 1
# 							Carray_event_types[numC] = "v"
# 							Carray_ivals[numC] = i
# 							Carray_jvals[numC] = k
# 							Carray_kvals[numC] = j
#								Carray_pair[numC] = 1
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
	Carray_pair = Carray_pair[TF]
	Cijk_weights = Cijk_weights[TF]
	Carray_event_types = Carray_event_types[TF]
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
	end

	
	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Carray_pair=Carray_pair, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)
	
	"""
	# Extract the values
	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Carray_pair = Carray.Carray_pair;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, pair=Carray_pair, wt=Cijk_weights, prob=Cijk_vals)
	row_weightvals
	"""
	
	return Carray
end # end setup_DEC_Cmat()










"""
### setup_DEC_Cmat2 HAS ISSUES 2022-03-17, WORK ON setup_DEC_Cmat3
# In setup_DEC_Cmat2, we will merge Cevent i,j,k and i,k,j, and have a x2 multiplier

areas_list = [1,2,3]
states_list = areas_list_to_states_list(areas_list, 3, true)
predeclare_array_length=10000000
Carray1 = setup_DEC_Cmat(areas_list, states_list);
Carray2 = setup_DEC_Cmat2(areas_list, states_list);
sum(prtC(Carray1)[:,:val])
sum(prtC(Carray2)[:,:val])

areas_list = [1,2,3]
states_list = areas_list_to_states_list(areas_list, 3, true)
predeclare_array_length=10000000
Carray1 = setup_DEC_Cmat(areas_list, states_list);
Carray2 = setup_DEC_Cmat2(areas_list, states_list);
prtC(Carray1)
prtC(Carray2)
sum(prtC(Carray1)[:,:val])
sum(prtC(Carray2)[:,:val])

Cparams=CparamsStructure(0.9, 0.9, 0.9, 0.3)
Carray1b = setup_DEC_Cmat(areas_list, states_list, NaN, Cparams);
Carray2b = setup_DEC_Cmat2(areas_list, states_list, NaN, Cparams);
prtC(Carray1b)
prtC(Carray2b)
sum(prtC(Carray1b)[:,:val])
sum(prtC(Carray2b)[:,:val])

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
Carray2 = setup_DEC_Cmat2(areas_list, states_list, maxent01, Cparams)
prtC(Carray2)

"""
function setup_DEC_Cmat2(areas_list, states_list, maxent01=NaN, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); predeclare_array_length=Integer(min(length(states_list)*length(states_list)*round((length(states_list)/2)), 10000000)), min_precision=1e-9)
	numareas = length(areas_list)
	total_numareas = numareas
	numstates = length(states_list)
	
	# Check if max_range_size=NaN
	type_string = string(typeof(maxent01))
	if (startswith(type_string, "NamedTuple") == false) && (isnan(maxent01) == true)
		maxent01 = maxent01_defaults(total_numareas)
	end
	
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
		txt = "\nSTOP ERROR in setup_DEC_Cmat2(): Your states_list is not ordered in ascending rangesize. Printing states_list:\n"
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
	Carray_pair = collect(repeat([1], predeclare_array_length))  # "pair" indicates if there are 1 or 2 such events
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
			for k in j:numstates # We only have to do half of the possible right events;
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
							Carray_pair[numC] = 1
							Cijk_weights[numC] = tmp_weightval
							row_weightvals[i] += tmp_weightval
							continue
						end # end if tmp_weightval > 0.0
					end # end if (all([lsize == ancsize, rsize==ancsize, lstate==ancstate, rstate==ancstate])
				end # end if (y_wt > min_precision)
				
				# If one of the descendants is identical to the ancestor, 
				# (after we've excluded sympatry)
				# we can have jump dispersal or subset speciation
				if ( all([ancsize==rsize, ancstate==rstate]) )   # Subset sympatry (no longer j events here)
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
								Carray_pair[numC] = 2
								Cijk_weights[numC] = tmp_weightval * 2.0 
								row_weightvals[i] += tmp_weightval * 2.0

								# Same event, flip left/right descendant states
								#numC += 1
								#Carray_event_types[numC] = "s"
								#Carray_ivals[numC] = i
								#Carray_jvals[numC] = k
								#Carray_kvals[numC] = j
								#Cijk_weights[numC] = tmp_weightval
								#row_weightvals[i] += tmp_weightval
								continue
							end # end if tmp_weightval > 0.0
						end # end if ((array_in_array(lstate, rstate) == true) && (lsize < rsize))
					end # end if (s_wt > min_precision)
				end # end if ( (ancstate == rstate) )
					
				# Jump dispersal
				if (j_wt > min_precision)
					if ( all([ancsize==rsize, ancstate==rstate]) )  # jump dispersal starter
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
					
							#print("\n")
							#print([i, j, k, tmp_weightval])
							if (tmp_weightval > min_precision)
								#print("yes")
								# Record the jump-dispersal event
								numC += 1
								Carray_event_types[numC] = "j"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = j
								Carray_kvals[numC] = k
								Carray_pair[numC] = 2
								Cijk_weights[numC] = tmp_weightval * 2.0
								row_weightvals[i] += tmp_weightval * 2.0

								# Same event, flip left/right descendant states
#								numC += 1
#								Carray_event_types[numC] = "j"
#								Carray_ivals[numC] = i
#								Carray_jvals[numC] = k
#								Carray_kvals[numC] = j
#								Cijk_weights[numC] = tmp_weightval
#								row_weightvals[i] += tmp_weightval
								continue
							end # if (tmp_weightval > 0.0)
							# end of jump dispersal
						end # end if ( (lsize == 1) && (array_in_array(lstate, rstate) == false) )
					end # end if ( all([ancsize==lsize, ancstate==lstate]) )
				end # end if (j_wt > min_precision)
			
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
							Carray_pair[numC] = 2
							Cijk_weights[numC] = tmp_weightval * 2.0
							row_weightvals[i] += tmp_weightval * 2.0
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
			end # end k (right state indices)
		end # end j (left state indices)
	end # end i (ancestor state indices)
	
	TF = Carray_event_types .!= ""
	Carray_ivals = Carray_ivals[TF]
	Carray_jvals = Carray_jvals[TF]
	Carray_kvals = Carray_kvals[TF]
	Carray_pair = Carray_pair[TF]
	Cijk_weights = Cijk_weights[TF]
	Carray_event_types = Carray_event_types[TF]
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
	end

	
	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Carray_pair=Carray_pair, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)
	
	"""
	# Extract the values
	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Carray_pair = Carray.Carray_pair;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, pair=Carray_pair, wt=Cijk_weights, prob=Cijk_vals)
	row_weightvals
	"""
	
	return Carray
end # end setup_DEC_Cmat2()
### setup_DEC_Cmat2 HAS ISSUES 2022-03-17, WORK ON setup_DEC_Cmat3




"""
# Sum the totals of events in the prtC(Carray) table

areas_list = [1,2,3]
states_list = areas_list_to_states_list(areas_list, 3, true)
predeclare_array_length=1000
Carray1 = setup_DEC_Cmat(areas_list, states_list);
Carray3 = setup_DEC_Cmat3(areas_list, states_list);
trt1 = prtC(Carray1); sort!(trt1, :i)
trt3 = prtC(Carray3); sort!(trt3, :i)
trt1
trt3
sum(prtC(Carray1)[:,:val])
sum(prtC(Carray3)[:,:val])

Cparams=CparamsStructure(0.75, 0.75, 0.75, 0.75)
Carray1b = setup_DEC_Cmat(areas_list, states_list, NaN, Cparams);
Carray3b = setup_DEC_Cmat3(areas_list, states_list, NaN, Cparams);
trt1b = prtC(Carray1b); 
trt3b = prtC(Carray3b); sort!(trt3b, :k); sort!(trt3b, :j); sort!(trt3b, :i)
trt1b
trt3b
sum(prtC(Carray1b)[:,:val])
sum(prtC(Carray3b)[:,:val])

# Check that they match
ijks1b = []
for q in 1:Rnrow(trt1b)
	tmpvec = []
	tmpvec = push!(tmpvec, trt1b.i[q])
	tmpvec = vcat(tmpvec, sort!([trt1b.j[q],trt1b.k[q]]) )
	push!(ijks1b, tmpvec)
end
sorted_ijks1b = sort!(unique(ijks1b))

ijks3b = []
for q in 1:Rnrow(trt3b)
	tmpvec = []
	tmpvec = push!(tmpvec, trt3b.i[q])
	tmpvec = vcat(tmpvec, sort!([trt3b.j[q],trt3b.k[q]]) )
	push!(ijks3b, tmpvec)
end
sorted_ijks3b = sort!(unique(ijks3b))
sort!(unique(ijks1b)) == sort!(unique(ijks3b))
# They match!  Do the numbers match??

events_table1b = totals_prtC(trt1b)
events_table3b = totals_prtC(trt3b)
events_table1b .== events_table3b
events_table1b .- events_table3b
"""
# Sum the totals of events in the prtC(Carray) table
function totals_prtC(Ctable)
	obs_states = sort(unique(Ctable.i))
	num_obs_states = length(obs_states)
	event_types = ["y", "s", "v", "j"]
	events_table = []
	
	for i in 1:length(obs_states)
		TF = Ctable.i .== obs_states[i]
		tmptable = Ctable[TF,:]

		num_events = []
		val_events = []
		wt_events = []

		for event_type in event_types
			TF2 = tmptable.event .== event_type
			tmp_num_events = sum(tmptable.pair[TF2])
			tmp_wt_events = sum(tmptable.wt[TF2])
			tmp_val_events = sum(tmptable.val[TF2])
			num_events = push!(num_events, tmp_num_events)
			wt_events = push!(wt_events, tmp_wt_events)
			val_events = push!(val_events, tmp_val_events)
		end # END for event_type in event_types
		tmp_events_table = vcat(num_events, wt_events, val_events)
		tmp_events_table2 = prepend!(tmp_events_table, [obs_states[i]])
		if i == 1
			events_table = transpose(tmp_events_table2)
		else
			events_table = Rrbind(events_table, transpose(tmp_events_table2))
		end
	end # END for i in 1:length(obs_states)
	
	# Make DataFrame and add names
	events_table = DataFrame(events_table, :auto)
	colnames = ["state#", "#y", "#s", "#v", "#j", "ywt", "swt", "vwt", "jwt", "yval", "sval", "vval", "jval"]
	symbols=Array{AbstractString,1}(colnames)
	rename!(events_table, symbols)

	return events_table
end # END function totals_prtC(Ctable)




"""
# In setup_DEC_Cmat3, we will merge Cevent i,j,k and i,k,j, and have a x2 multiplier

areas_list = [1,2,3]
states_list = areas_list_to_states_list(areas_list, 3, true)
predeclare_array_length=1000
Carray1 = setup_DEC_Cmat(areas_list, states_list);
Carray3 = setup_DEC_Cmat3(areas_list, states_list);
trt1 = prtC(Carray1); sort!(trt1, :i)
trt3 = prtC(Carray3); sort!(trt3, :i)
trt1
trt3
sum(prtC(Carray1)[:,:val])
sum(prtC(Carray3)[:,:val])

Cparams=CparamsStructure(0.75, 0.75, 0.75, 0.75)
Carray1b = setup_DEC_Cmat(areas_list, states_list, NaN, Cparams);
Carray3b = setup_DEC_Cmat3(areas_list, states_list, NaN, Cparams);
trt1b = prtC(Carray1b); 
trt3b = prtC(Carray3b); sort!(trt3b, :k); sort!(trt3b, :j); sort!(trt3b, :i)
trt1b
trt3b
sum(prtC(Carray1b)[:,:val])
sum(prtC(Carray3b)[:,:val])

# Check that they match
ijks1b = []
for q in 1:Rnrow(trt1b)
	tmpvec = []
	tmpvec = push!(tmpvec, trt1b.i[q])
	tmpvec = vcat(tmpvec, sort!([trt1b.j[q],trt1b.k[q]]) )
	push!(ijks1b, tmpvec)
end
sorted_ijks1b = sort!(unique(ijks1b))

ijks3b = []
for q in 1:Rnrow(trt3b)
	tmpvec = []
	tmpvec = push!(tmpvec, trt3b.i[q])
	tmpvec = vcat(tmpvec, sort!([trt3b.j[q],trt3b.k[q]]) )
	push!(ijks3b, tmpvec)
end
sorted_ijks3b = sort!(unique(ijks3b))
sort!(unique(ijks1b)) == sort!(unique(ijks3b))
# They match!  Do the numbers match??

events_table1b = totals_prtC(trt1b)
events_table3b = totals_prtC(trt3b)
events_table1b .== events_table3b
events_table1b .- events_table3b

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
predeclare_array_length=Integer(min(length(states_list)*length(states_list)*round((length(states_list)/2)), 10000000))
Carray3 = setup_DEC_Cmat3(areas_list, states_list, maxent01, Cparams)
prtC(Carray3)

"""
function setup_DEC_Cmat3(areas_list, states_list, maxent01=NaN, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); predeclare_array_length=Integer(min(length(states_list)*length(states_list)*round((length(states_list)/2)), 10000000)), min_precision=1e-9)
	numareas = length(areas_list)
	total_numareas = numareas
	numstates = length(states_list)
	
	# Check if max_range_size=NaN
	type_string = string(typeof(maxent01))
	if (startswith(type_string, "NamedTuple") == false) && (isnan(maxent01) == true)
		maxent01 = maxent01_defaults(total_numareas)
	end
	
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
		txt = "\nSTOP ERROR in setup_DEC_Cmat3(): Your states_list is not ordered in ascending rangesize. Printing states_list:\n"
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
	Carray_pair = collect(repeat([1], predeclare_array_length))  # "pair" indicates if there are 1 or 2 such events
	Cijk_weights = collect(repeat([0.0], predeclare_array_length))
	Carray_event_types = collect(repeat([""], predeclare_array_length))
	numC = 0 # counter of the number of allow cladogenesis events
	
	
	# Go through:
	# i = ancestor state index
	# j = left state index
	# k = right state index
	
	# Preallocate this vector ONCE, size = numareas * 2
	#tmp_merged_vec = repeat([0], 2*numareas)
	
	# Define a vector of rangesizes, and indices to them
	rangesizes = length.(states_list)
	# Check that rangesizes strictly increase
	rangesizes_subtracted = rangesizes[2:length(rangesizes)] - rangesizes[1:(length(rangesizes)-1)]
	if any(rangesizes_subtracted .< 0)
		txt = paste0(["STOP error in setup_DEC_Cmat3(): the range sizes of states_list *must* be in rangesize order. Printing rangesizes:"])
		print("\n")
		print(txt)
		print("\n")
		print(rangesizes)
		print("\n")
		error(txt)
	end
	
	# Set up a dictionary for each rangesize category
	if rangesizes[1] == 0
		range_size_category_indexes_dict = Dict(0 => [1])
	end
	if rangesizes[1] == 1
		range_size_category_indexes_dict = Dict(1 => [1])
	end
	
	# Fill in the rangesize category dictionary
	oldsize = 0
	newsize = 0
	for i in 1:numstates
		newsize = rangesizes[i]
		if newsize > oldsize
			range_size_category_indexes_dict[newsize] = findall(x->x==newsize, rangesizes)
		end
		oldsize = newsize
	end 
	range_size_category_indexes_dict

	############################################
	# Sympatry (range-copying)
	############################################
	if (y_wt > min_precision)
		for i in minimum(range_size_category_indexes_dict[1]):numstates
			ancstate = states_list[i]
			ancsize = length(ancstate)
			# Exclude the null range as an ancestor for cladogenetic range-changing events
			# (ancestor of range NULL give 0 likelihood to the observed data)
			# Also exclude if daughter rangesize not allowed by maxent01symp

			possible_min_rangesizes_indices = collect(1:Rncol(maxent01symp))
			TF = maxent01symp[ancsize,:] .> min_precision
			max_min_rangesize = maximum(possible_min_rangesizes_indices[TF])
			if (ancsize > max_min_rangesize) || (ancsize == 0)
				continue # go to next ancestral state i
			end

			# Store the sympatry event, if it has positive weightį
			tmp_weightval = y_wt * maxent01symp[ancsize, ancsize] * 1.0 * 1.0
			if tmp_weightval > min_precision
				# Record the range-copying event
				# (no reverse event to record)
				numC += 1
				Carray_event_types[numC] = "y"
				Carray_ivals[numC] = i
				Carray_jvals[numC] = i
				Carray_kvals[numC] = i
				Carray_pair[numC] = 1
				Cijk_weights[numC] = tmp_weightval
				row_weightvals[i] += tmp_weightval
				continue
			end # end if tmp_weightval > min_precision
		end # end for i in 1:numstates
	end # END if (y_wt > min_precision)  (ending sympatry)
	
	############################################
	# Subset sympatry (range-copying)
	############################################
	if (s_wt > min_precision)
		# Start at rangesize 2
		for i in minimum(range_size_category_indexes_dict[2]):numstates
			ancstate = states_list[i]
			ancsize = length(ancstate)
			
			# Subset events can only happen if range size >= 2 areas
			if ancsize < 2
				continue
			end # END if ancsize < 2

			# Exclude the null range as an ancestor for cladogenetic range-changing events
			# (ancestor of range NULL give 0 likelihood to the observed data)
			# Also exclude if daughter rangesize not allowed by maxent01symp
			possible_min_rangesizes_indices = collect(1:Rncol(maxent01sub))
			TF = maxent01sub[ancsize,:] .> min_precision
			max_min_rangesize = maximum(possible_min_rangesizes_indices[TF])
	
			# For every range size between 2 and the maximum minrange, list every combination
			# WITHIN the anc areas--the tough thing here is that you don't have the indexes
			
			# Go through the states of the possible smaller descendants
			for daughter_size in 1:(max_min_rangesize)
				for state_index in range_size_category_indexes_dict[daughter_size]
					rstate = states_list[state_index]
					rsize = length(rstate)
					# Include only if ALL the descendant right states are in the ancestor
					# (and exclude sympatry)
					if (array_in_array(rstate, ancstate) == true) && (ancstate != rstate)
						tmp_weightval = s_wt * maxent01sub[ancsize, rsize] * 1.0 * 1.0
						if tmp_weightval > min_precision
							# Record the range-copying event
							numC += 1
							Carray_event_types[numC] = "s"
							Carray_ivals[numC] = i
							Carray_jvals[numC] = i	# because one daughter is same as anc in subsets
							Carray_kvals[numC] = state_index
							Carray_pair[numC] = 2
							Cijk_weights[numC] = tmp_weightval * 2.0 
							row_weightvals[i] += tmp_weightval * 2.0
							continue
						end # END if tmp_weightval > min_precision
					end # END if array_in_array(rstate, ancstate) == true)
				end # END for state_index in range_size_category_indexes_dict[daughter_size]
			end # END for daughter_size in 1:(max_min_rangesize)
		end # END for i in range_size_category_indexes_dict[2]:numstates
	end # END if (s_wt > min_precision)
	
	############################################
	# Vicariance (range-splitting)
	############################################
	if (v_wt > min_precision)
		# Start at rangesize 2
		for i in minimum(range_size_category_indexes_dict[2]):numstates
			ancstate = states_list[i]
			ancsize = length(ancstate)
			
			# Vicariance events can only happen if range size >= 2 areas
			if ancsize < 2
				continue
			end # END if ancsize < 2

			# Exclude the null range as an ancestor for cladogenetic range-changing events
			# (ancestor of range NULL give 0 likelihood to the observed data)
			# Also exclude if daughter rangesize not allowed by maxent01symp
			possible_min_rangesizes_indices = collect(1:Rncol(maxent01vic))
			TF = maxent01vic[ancsize,:] .> min_precision
			max_min_rangesize = maximum(possible_min_rangesizes_indices[TF])
	
			# For every range size between 2 and the maximum minrange, list every combination
			# WITHIN the anc areas--the tough thing here is that you don't have the indexes
			
			# Go through the states of the possible smaller descendants
			for daughter_size in 1:(max_min_rangesize)
				for state_index in range_size_category_indexes_dict[daughter_size]
					rstate = states_list[state_index]
					rstate_index = state_index
					rsize = length(rstate)
					# Include only if ALL the descendant right states are in the ancestor
					# (and exclude sympatry)
					if (array_in_array(rstate, ancstate) == true) && (rstate != ancstate)
					#if ( is_event_vicariance(ancstate, lstate, rstate) )
						# Find the left state
						lstate = ancstate[(!in).(ancstate,Ref(rstate))]
						sort!(lstate)
						lstate_index = (1:numstates)[[lstate] .== states_list][1]
						lsize = length(lstate)
						smaller_range_size = min(lsize, rsize)
						tmp_weightval = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0
						if (tmp_weightval > min_precision)
							if lsize != rsize
								# Record the jump-dispersal event
								numC += 1
								Carray_event_types[numC] = "v"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = lstate_index
								Carray_kvals[numC] = rstate_index
								Carray_pair[numC] = 2
								Cijk_weights[numC] = tmp_weightval * 2.0
								row_weightvals[i] += tmp_weightval * 2.0
								continue
							end # END if lsize != rsize
							if (lsize == rsize) && (lstate_index > rstate_index)
								# These get double-recorded, so keep only half
								numC += 1
								Carray_event_types[numC] = "v"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = lstate_index
								Carray_kvals[numC] = rstate_index
								Carray_pair[numC] = 2
								Cijk_weights[numC] = tmp_weightval * 2.0
								row_weightvals[i] += tmp_weightval * 2.0
								continue
							end # END if lsize == rsize
						end # END if tmp_weightval > min_precision
					end # END if array_in_array(rstate, ancstate) == true)
				end # END for state_index in range_size_category_indexes_dict[daughter_size]
			end # END for daughter_size in 1:(max_min_rangesize)
		end # END for i in range_size_category_indexes_dict[2]:numstates
	end # END if (v_wt > min_precision)
	
	############################################
	# Jump dispersal (founder-event speciation)
	############################################
	if (j_wt > min_precision)
		# Start at rangesize 1
		for i in minimum(range_size_category_indexes_dict[1]):numstates
			ancstate = states_list[i]
			ancsize = length(ancstate)
			
			# Exclude the null range as an ancestor for cladogenetic range-changing events
			# (ancestor of range NULL give 0 likelihood to the observed data)
			# Also exclude if daughter rangesize not allowed by maxent01symp
			possible_min_rangesizes_indices = collect(1:Rncol(maxent01jump))
			TF = maxent01jump[ancsize,:] .> min_precision
			max_min_rangesize = maximum(possible_min_rangesizes_indices[TF])
	
			# For every range size between 2 and the maximum minrange, list every combination
			# WITHIN the anc areas--the tough thing here is that you don't have the indexes
			
			# Go through the states of the possible smaller descendants
			for daughter_size in 1:(max_min_rangesize)
				for state_index in range_size_category_indexes_dict[daughter_size]
					rstate = states_list[state_index]
					rstate_index = state_index
					rsize = length(rstate)
					lstate = ancstate
					lstate_index = i
					lsize = length(lstate)
					# Include only if ALL the descendant right states are in the ancestor
					# (and exclude sympatry)
					if (array_in_array(rstate, ancstate) == false)
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
								for right_area in rstate
									 jweight_for_cell_based_on_distances += dmat[anc_area,right_area]
								end
							end
							# Normalize by number of possible jump dispersals
							if (normalize_by_number_of_dispersal_events == true)
								jweight_for_cell_based_on_distances = jweight_for_cell_based_on_distances / (ancsize * rsize)
							end
						else
							# 
							jweight_for_cell_based_on_distances = 1.0
						end # end if (try_jump_dispersal_based_on_dist == true)
			
						# Calculate the final weight of this jump dispersal
						tmp_weightval = j_wt * maxent01jump[ancsize, rsize] * 1.0 * 1.0 * jweight_for_cell_based_on_distances

						if tmp_weightval > min_precision
							# Record the jump dispersal event
#							if (lsize != rsize)
								numC += 1
								Carray_event_types[numC] = "j"
								Carray_ivals[numC] = i
								Carray_jvals[numC] = i	# because one daughter is same as anc in jumps
								Carray_kvals[numC] = rstate_index	# rstate for jump
								Carray_pair[numC] = 2
								Cijk_weights[numC] = tmp_weightval * 2.0 
								row_weightvals[i] += tmp_weightval * 2.0
								continue
#							end # END if lsize != rsize
#							if (lsize == rsize) && (lstate_index > rstate_index)
#								# These get double-recorded, so keep only half
#								numC += 1
#								Carray_event_types[numC] = "j"
#								Carray_ivals[numC] = i
#								Carray_jvals[numC] = lstate_index
#								Carray_kvals[numC] = rstate_index
#								Carray_pair[numC] = 2
#								Cijk_weights[numC] = tmp_weightval * 2.0
#								row_weightvals[i] += tmp_weightval * 2.0
#								continue
#							end # END if (lsize == rsize) && (lstate_index > rstate_index)

						end # END if tmp_weightval > min_precision
					end # END if array_in_array(rstate, ancstate) == false)
				end # END for state_index in range_size_category_indexes_dict[daughter_size]
			end # END for daughter_size in 1:(max_min_rangesize)
		end # END for i in range_size_category_indexes_dict[1]:numstates
	end # END if (j_wt > min_precision)
	
	# Package the results
	TF = Carray_event_types .!= ""
	Carray_ivals = Carray_ivals[TF]
	Carray_jvals = Carray_jvals[TF]
	Carray_kvals = Carray_kvals[TF]
	Carray_pair = Carray_pair[TF]
	Cijk_weights = Cijk_weights[TF]
	Carray_event_types = Carray_event_types[TF]
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
	end

	
	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Carray_pair=Carray_pair, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)
	
	"""
	# Extract the values
	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Carray_pair = Carray.Carray_pair;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, pair=Carray_pair, wt=Cijk_weights, prob=Cijk_vals)
	row_weightvals
	"""
	
	return Carray
end # end setup_DEC_Cmat3()





"""
######################################
# update_maxent01(bmo)
######################################
#
# The "maxent01" object contains a series of tables specifying the relative weights
# of daughters of different range sizes.  This was constructed in BioGeoBEARS to 
# allow a single parameter to control the decisions about daughter range sizes.
# 
# The "maxent" part only refers to the fact that a Maximum Entropy function applied
# to discrete numerical data (e.g. numbers 1-6 on a die) can put flat or skewed 
# weights/probabilities on the various possible outcomes. 
#
# Basically, there is a mx01 parameter for each cladogenesis process: y,s,v, and j.
#
# If the mx01 parameter = 0.0, then the smaller daughter with a rangesize of 1
# has 100% probability. (Thus, DEC has mx01y, mx01s, mx01v = 0.0.)
#
# If the mx01 parameter = 0.5, then the smaller daughter with a rangesize of 1
# has 100% probability. (Thus, DIVALIKE has mx01v = 0.5.)
#
# If the mx01 parameter = 1.0 then the smaller daughter will have the 
# largest possible rangesize, given the ancestor range. (Thus, BAYAREALIKE has mx01y=0.999.)
# 
# update_maxent01(bmo) re-calculates the maxent01 list of tables, given the mx01 parameters
# in bmo (the BioGeoBEARS_model_object)
"""
function update_maxent01(bmo)
	# Generate the tables (weights for smaller daughter rangesize, for each ancestral rangesize)
	maxent01symp = relative_probabilities_of_subsets(total_numareas, bmo.est[bmo.rownames .== "mx01y"][1])
	maxent01sub = relative_probabilities_of_subsets(total_numareas, bmo.est[bmo.rownames .== "mx01s"][1])
	maxent01jump = relative_probabilities_of_subsets(total_numareas, bmo.est[bmo.rownames .== "mx01j"][1])
	maxent01vic = relative_probabilities_of_vicariants(total_numareas, bmo.est[bmo.rownames .== "mx01v"][1])
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
	return maxent01
end # END function update_maxent01(bmo)



#######################################################
# Update the Cijk_vals
#######################################################
function update_Cijk_vals(Carray, areas_list, states_list, maxent01, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); row_weightvals=null )

	if (Rclass(Carray) == "DataFrame")
		if (row_weightvals == null)
			txt = "STOP ERROR in update_Cijk_vals(): Input 'Carray' is a DataFrame, so row_weightvals cannot be null."
			println("\n")
			println(txt)
			println("\n")
			error(txt)
		end
		Carray_ivals = Carray[!,:i]
		Carray_jvals = Carray[!,:j]
		Carray_kvals = Carray[!,:k]
		Cijk_weights = Carray[!,:wt]
		Cijk_vals = Carray[!,:val]
		row_weightvals = row_weightvals
	else
		Carray_ivals = Carray.Carray_ivals;
		Carray_jvals = Carray.Carray_jvals;
		Carray_kvals = Carray.Carray_kvals;
		Carray_event_types = Carray.Carray_event_types;
		Cijk_weights = Carray.Cijk_weights;
		Cijk_vals = Carray.Cijk_vals;
		row_weightvals = Carray.row_weightvals
	end # END if (Rclass(Carray) == "DataFrame")



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
	
	# Update the "y" events (narrow sympatry / range-copying)
	TF = Carray_event_types .== "y"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			weights[z] = y_wt * maxent01symp[ancsize, lsize] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of y event weights


	# Update the "s" events (subset sympatry)
	TF = Carray_event_types .== "s"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			weights[z] = s_wt * maxent01sub[ancsize, lsize] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	# Update the "v" events (vicariance)
	TF = Carray_event_types .== "v"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of v event weights
	
	# Update the "j" events (jump dispersal)
	TF = Carray_event_types .== "j"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			
			# j events, modified by distance / multipliers (via input "dmat") if needed
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
			weights[z] = tmp_weightval
		end
		Cijk_weights[TF] = weights
	end # End update of j event weights

#	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
#	row_weightvals_df = by(df1, :i, :weight => sum)
#	row_weightvals = row_weightvals_df[!,2]
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
#	if (in(1, row_weightvals_df[!,1]) == false)
#		row_weightvals = repeat([0], 1+length(row_weightvals_df[!,2]))
#		row_weightvals[1] = 1
#		row_weightvals[2:length(row_weightvals)] = row_weightvals_df[!,2]
#	end
#	row_weightvals
	
	# Convert the weights to conditional event probabilities
#	num_clado_events = length(Cijk_weights)
#	Cijk_vals = collect(repeat([0.0], num_clado_events))
#	for i in 1:length(states_list)
#		TF = Carray_ivals .== i
#		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
#	end

#	df2 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
#	row_weightvals_df = by(df2, :i, :weight => sum)

	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
	groups = groupby(df1,:i)
	row_weightvals = collect(repeat([0.0], length(groups)))
	for g in 1:length(groups)
		row_weightvals[g] = sum(groups[g].weight)
	end
	row_weightvals
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
	if in(1, unique(df1.i)) == false
		prepend!(row_weightvals, 0)
	end
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] .* birthRate ./ row_weightvals[i]
	end


	
	# Finally, return updated Carray:
	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)

	"""
	# Extract the values
	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals;
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	row_weightvals
	"""

	return Carray
end # end update_Cijk_vals()




#######################################################
# Update the update_Cijk_vals2_noUpdate
# maxent01 = list of tables
#######################################################
function update_Cijk_vals2_noUpdate(p_Ds_v5, areas_list, states_list, bmo, maxent01, jmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))))

	Carray_ivals = p_Ds_v5.p_indices.Carray_ivals
	Carray_jvals = p_Ds_v5.p_indices.Carray_jvals
	Carray_kvals = p_Ds_v5.p_indices.Carray_kvals
	Carray_pair = p_Ds_v5.p_indices.Carray_pair
	Cijk_weights = p_Ds_v5.params.Cijk_weights
	Cijk_vals = p_Ds_v5.params.Cijk_vals
	Carray_event_types = p_Ds_v5.p_indices.Carray_event_types
	row_weightvals = p_Ds_v5.params.row_weightvals

	numstates = length(states_list)
	total_numareas = length(areas_list)
	
	# Get the parameters
	# Weights
	y_wt = bmo.est[bmo.rownames .== "y"][1]
	s_wt = bmo.est[bmo.rownames .== "s"][1]
	v_wt = bmo.est[bmo.rownames .== "v"][1]
	j_wt = bmo.est[bmo.rownames .== "j"][1]
	mx01y = bmo.est[bmo.rownames .== "mx01y"][1]
	mx01s = bmo.est[bmo.rownames .== "mx01s"][1]
	mx01v = bmo.est[bmo.rownames .== "mx01v"][1]
	mx01j = bmo.est[bmo.rownames .== "mx01j"][1]

	# Print y, s, v, j
	print("\n")
	txt = paste0(["update_Cijk_vals2: y=", round(y_wt, digits=4), ", s=", round(s_wt, digits=4), ", v=", round(v_wt, digits=4), ", j=", round(j_wt, digits=4)])
	print(txt)
	
	# Rates
	birthRate = bmo.est[bmo.rownames .== "birthRate"][1]
	
	# Extract the maxent01 tables, which control
	# the relative weight of the rangesizes of
	# the smaller daughter
	# (Update these tables OUTSIDE of this
	#  function!)
	maxent01symp = maxent01.maxent01symp
	maxent01sub = maxent01.maxent01sub
	maxent01vic = maxent01.maxent01vic
	maxent01jump = maxent01.maxent01jump
	
	
	# Update the "y" events (narrow sympatry / range-copying)
	TF = Carray_event_types .== "y"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = y_wt * maxent01symp[ancsize, smaller_range_size] * pair_vals[z]
			#weights[z] = y_wt * maxent01symp[ancsize, smaller_range_size] * 1.0 * 1.0
			Carray_pair
		end
		Cijk_weights[TF] = weights
	end # End update of y event weights


	# Update the "s" events (subset sympatry)
	TF = Carray_event_types .== "s"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = s_wt * maxent01sub[ancsize, smaller_range_size] * pair_vals[z]
			#weights[z] = s_wt * maxent01sub[ancsize, smaller_range_size] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	# Update the "v" events (vicariance)
	TF = Carray_event_types .== "v"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = v_wt * maxent01vic[ancsize,smaller_range_size] * pair_vals[z]
			#weights[z] = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0
		end
		Cijk_weights[TF] = weights
	end # End update of v event weights
	
	# Update the "j" events (jump dispersal)
	TF = Carray_event_types .== "j"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			
			# j events, modified by distance / multipliers (via input "jmat") if needed
			try_jump_dispersal_based_on_dist = true
			normalize_by_number_of_dispersal_events = true
			jweight_for_cell_based_on_distances = 0.0
			if (try_jump_dispersal_based_on_dist == true)
				for anc_area in ancstate
					for left_area in lstate
						 jweight_for_cell_based_on_distances += jmat[anc_area,left_area]
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
			smaller_range_size = min(lsize, rsize)
			tmp_weightval = j_wt * maxent01jump[ancsize, smaller_range_size] * jweight_for_cell_based_on_distances * pair_vals[z]
			#tmp_weightval = j_wt * maxent01jump[ancsize, smaller_range_size] * 1.0 * 1.0 * jweight_for_cell_based_on_distances
			
			weights[z] = tmp_weightval
		end
		Cijk_weights[TF] = weights
	end # End update of j event weights

	#df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
#	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
#	row_weightvals_df = by(df1, :i, :weight => sum)
#	row_weightvals = row_weightvals_df[!,2]
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
#	if (in(1, row_weightvals_df[!,1]) == false)
#		row_weightvals = repeat([0], 1+length(row_weightvals_df[!,2]))
#		row_weightvals[1] = 1
#		row_weightvals[2:length(row_weightvals)] = row_weightvals_df[!,2]
#	end
#	row_weightvals
	
	# Convert the weights to conditional event probabilities
#	num_clado_events = length(Cijk_weights)
#	Cijk_vals = collect(repeat([0.0], num_clado_events))
#	for i in 1:length(states_list)
#		TF = Carray_ivals .== i
#		Cijk_vals[TF] = Cijk_weights[TF] .* birthRate ./ row_weightvals[i]
#	end

	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
	groups = groupby(df1,:i)
	row_weightvals = collect(repeat([0.0], length(groups)))
	for g in 1:length(groups)
		row_weightvals[g] = sum(groups[g].weight)
	end
	row_weightvals
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
	if in(1, unique(df1.i)) == false
		prepend!(row_weightvals, 0)
	end
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] .* birthRate ./ row_weightvals[i]
	end

	p_Ds_v5.params.Cijk_weights[:] .= Cijk_weights
	p_Ds_v5.params.Cijk_vals[:] .= Cijk_vals
	p_Ds_v5.params.row_weightvals[:] .= row_weightvals


	"""
	# Extract the values
	Carray_event_types = p_Ds_v5.p_indices.Carray_event_types;
	Carray_ivals = p_Ds_v5.p_indices.Carray_ivals
	Carray_jvals = p_Ds_v5.p_indices.Carray_jvals
	Carray_kvals = p_Ds_v5.p_indices.Carray_kvals
	Cijk_weights = p_Ds_v5.params.Cijk_weights;
	Cijk_vals = p_Ds_v5.params.Cijk_vals;
	row_weightvals = p_Ds_v5.params.row_weightvals
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	"""
	
	output = (Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)
	return output
	#return p_Ds_v5
end # end update_Cijk_vals2_noUpdate()









#######################################################
# Update the Cijk_vals2!
# maxent01 = list of tables
#######################################################
function update_Cijk_vals2!(p_Ds_v5, areas_list, states_list, bmo, maxent01, jmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); printlevel=0)

	Carray_ivals = p_Ds_v5.p_indices.Carray_ivals
	Carray_jvals = p_Ds_v5.p_indices.Carray_jvals
	Carray_kvals = p_Ds_v5.p_indices.Carray_kvals
	Carray_pair = p_Ds_v5.p_indices.Carray_pair
	Cijk_weights = p_Ds_v5.params.Cijk_weights
	Cijk_vals = p_Ds_v5.params.Cijk_vals
	Carray_event_types = p_Ds_v5.p_indices.Carray_event_types
	row_weightvals = p_Ds_v5.params.row_weightvals

	numstates = length(states_list)
	total_numareas = length(areas_list)
	
	# Get the parameters
	# Weights
	y_wt = bmo.est[bmo.rownames .== "y"][1]
	s_wt = bmo.est[bmo.rownames .== "s"][1]
	v_wt = bmo.est[bmo.rownames .== "v"][1]
	j_wt = bmo.est[bmo.rownames .== "j"][1]
	mx01y = bmo.est[bmo.rownames .== "mx01y"][1]
	mx01s = bmo.est[bmo.rownames .== "mx01s"][1]
	mx01v = bmo.est[bmo.rownames .== "mx01v"][1]
	mx01j = bmo.est[bmo.rownames .== "mx01j"][1]

	# Print y, s, v, j
	if printlevel >= 2
		print("\n")
		txt = paste0(["update_Cijk_vals2: y=", round(y_wt, digits=3), ", s=", round(s_wt, digits=3), ", v=", round(v_wt, digits=3), ", j=", round(j_wt, digits=3)])
		print(txt)
	end
	
	# Rates
	birthRate = bmo.est[bmo.rownames .== "birthRate"][1]
	
	# Extract the maxent01 tables, which control
	# the relative weight of the rangesizes of
	# the smaller daughter
	# (Update these tables OUTSIDE of this
	#  function!)
	maxent01symp = maxent01.maxent01symp
	maxent01sub = maxent01.maxent01sub
	maxent01vic = maxent01.maxent01vic
	maxent01jump = maxent01.maxent01jump
	
	
	# Update the "y" events (narrow sympatry / range-copying)
	TF = Carray_event_types .== "y"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = y_wt * maxent01symp[ancsize, smaller_range_size] * 1.0 * 1.0 * pair_vals[z]
		end
		Cijk_weights[TF] = weights
	end # End update of y event weights


	# Update the "s" events (subset sympatry)
	TF = Carray_event_types .== "s"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = s_wt * maxent01sub[ancsize, smaller_range_size] * 1.0 * 1.0 * pair_vals[z]
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	# Update the "v" events (vicariance)
	TF = Carray_event_types .== "v"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			smaller_range_size = min(lsize, rsize)
			weights[z] = v_wt * maxent01vic[ancsize,smaller_range_size] * 1.0 * 1.0 * pair_vals[z]
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights
	
	# Update the "j" events (jump dispersal)
	TF = Carray_event_types .== "j"
	if (sum(TF) > 0)
		ivals = Carray_ivals[TF]
		jvals = Carray_jvals[TF]
		kvals = Carray_kvals[TF]
		pair_vals = Carray_pair[TF]
		weights = Cijk_weights[TF]
		for z in 1:sum(TF)
			i = ivals[z]
			j = jvals[z]
			k = kvals[z]
			ancstate = states_list[i]
			ancsize = length(ancstate)
			lstate = states_list[j]
			lsize = length(lstate)
			rstate = states_list[k]
			rsize = length(rstate)
			
			# j events, modified by distance / multipliers (via input "jmat") if needed
			try_jump_dispersal_based_on_dist = true
			normalize_by_number_of_dispersal_events = true
			jweight_for_cell_based_on_distances = 0.0
			if (try_jump_dispersal_based_on_dist == true)
				for anc_area in ancstate
					for left_area in lstate
						 jweight_for_cell_based_on_distances += jmat[anc_area,left_area]
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
			smaller_range_size = min(lsize, rsize)
			tmp_weightval = j_wt * maxent01jump[ancsize, smaller_range_size] * 1.0 * 1.0 * jweight_for_cell_based_on_distances * pair_vals[z]
			weights[z] = tmp_weightval
		end
		Cijk_weights[TF] = weights
	end # End update of s event weights

	#df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
	
	# Add up the weights of the rows for each ancestral state i
	# And, divide by the sum of the weights
#	for i in 1:numstates
#		iTF = Carray_ivals .== i
#		row_weightvals[i] = sum(Cijk_weights[iTF])
#		# Multiply by birthRate here!!
#		Cijk_vals[iTF] .= Cijk_weights[iTF] .* birthRate ./ row_weightvals[i]
#	end # END for i in 1:numstates
#	Cijk_vals = p_Ds_v5.params.Cijk_vals
	
	# Convert the weights to conditional event probabilities
#	num_clado_events = length(Cijk_weights)
#	Cijk_vals = collect(repeat([0.0], num_clado_events))
#	for i in 1:length(states_list)
#		TF = Carray_ivals .== i
#		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
#	end

#	df2 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
#	row_weightvals_df = by(df2, :i, :weight => sum)
	
	# Finally, return updated Carray:
#	Carray = (Carray_event_types=Carray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals, row_weightvals=row_weightvals)



	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
	groups = groupby(df1,:i)
	row_weightvals = collect(repeat([0.0], length(groups)))
	for g in 1:length(groups)
		row_weightvals[g] = sum(groups[g].weight)
	end
	row_weightvals
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
	if in(1, unique(df1.i)) == false
		prepend!(row_weightvals, 0)
	end
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] .* birthRate ./ row_weightvals[i]
	end

	p_Ds_v5.params.Cijk_weights[:] .= Cijk_weights
	p_Ds_v5.params.Cijk_vals[:] .= Cijk_vals
	p_Ds_v5.params.row_weightvals[:] .= row_weightvals


	"""
	# Extract the values
	Carray_event_types = p_Ds_v5.p_indices.Carray_event_types;
	Carray_ivals = p_Ds_v5.p_indices.Carray_ivals
	Carray_jvals = p_Ds_v5.p_indices.Carray_jvals
	Carray_kvals = p_Ds_v5.p_indices.Carray_kvals
	Cijk_weights = p_Ds_v5.params.Cijk_weights;
	Cijk_vals = p_Ds_v5.params.Cijk_vals;
	row_weightvals = p_Ds_v5.params.row_weightvals
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	"""

	return p_Ds_v5
end # end update_Cijk_vals2!()





"""
# Update the BioGeoBEARS_model_object after changing e.g. "j"
"""
function bmo_updater_v1!(bmo)
	j_wt = bmo.est[bmo.rownames .== "j"][1]
	ysv = bmo.est[bmo.rownames .== "ysv"][1]
	ys = bmo.est[bmo.rownames .== "ys"][1]
	y = bmo.est[bmo.rownames .== "y"][1]
	s = bmo.est[bmo.rownames .== "s"][1]
	v = bmo.est[bmo.rownames .== "v"][1]
	
	# Update
	ysv_func = bmo.type[bmo.rownames .== "ysv"][1]
	if ysv_func == "3-j"
		ysv = 3-j_wt
	end
	if ysv_func == "2-j"
		ysv = 2-j_wt
	end
	if ysv_func == "1-j"
		ysv = 1-j_wt
	end
	
	ys = ysv*2/3
	y = ysv*1/3
	s = ysv*1/3
	v = ysv*1/3
	
	bmo.est[bmo.rownames .== "ysv"] .= ysv
	bmo.est[bmo.rownames .== "ys"] .= ys
	bmo.est[bmo.rownames .== "y"] .= y
	bmo.est[bmo.rownames .== "s"] .= s
	bmo.est[bmo.rownames .== "v"] .= v
	
	#return bmo
end # END



"""
# p_Ds_v5 updater #1
# p_Ds_v5 contains the actual rates that go straight into the lnL calculations

pars = [0.03505038, 0.02832370]
parnames = ["d", "e"]

# Assuming parameters are in the order of the bmo list of "free" parameters,
# 1. update the bmo
# 2. then update the Qij and Cijk arrays in p_Ds_v5 (the rate parameters arrays)
inputs.bmo.rownames[inputs.bmo.type.=="free"]

pars = [0.03505038, 0.02832370]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs)
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)

pars = [0.01, 0.001]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs)
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)

pars = [0.1, 0.2]
inputs.bmo.est[inputs.bmo.type.=="free"] .= pars
p_Ds_v5 = p_Ds_v5_updater_v1!(p_Ds_v5, inputs)
prtQp(p_Ds_v5)
prtCp(p_Ds_v5)


lnL = func_to_optimize(pars, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func(pars)
"""
function p_Ds_v5_updater_v1!(p_Ds_v5, inputs; check_if_free_params_in_mat=true, printlevel=0)
	
	if check_if_free_params_in_mat == true
		free_param_names = inputs.bmo.rownames[inputs.bmo.type .== "free"]
		# Are all of the free param names %in% the p_Ds_v5
		#TF = in([unique(p_Ds_v5.p_indices.Qarray_event_types); unique(p_Ds_v5.p_indices.Carray_event_types)]).(free_param_names)
		TF1 = in(free_param_names).(["j"])[1]
		TF2 = in(unique(p_Ds_v5.p_indices.Carray_event_types)).(["j"])[1]
		if TF1 != TF2
			txt = paste0(["ERROR in p_Ds_v5_updater_v1!(): j's must be in both places. j_in_bmo:", TF1, "; j_in_params:", TF2])
			print("\n")
			print(txt)
			print("\n")
			throw(txt)
		end # END if TF1 != TF2
	end # END if check_if_free_params_in_mat == true
	
	# Extract the parameters
	d = inputs.bmo.est[inputs.bmo.rownames .== "d"][1]
	a = inputs.bmo.est[inputs.bmo.rownames .== "a"][1]
	e = inputs.bmo.est[inputs.bmo.rownames .== "e"][1]
	x = inputs.bmo.est[inputs.bmo.rownames .== "x"][1]
	w = inputs.bmo.est[inputs.bmo.rownames .== "w"][1]
	n = inputs.bmo.est[inputs.bmo.rownames .== "n"][1]
	x2 = 1.0
	x3 = 1.0 
	u = inputs.bmo.est[inputs.bmo.rownames .== "u"][1]
	
	# Update the dmat and elist
	inputs.setup.dmat_base .= d
	inputs.setup.amat_base .= a
	inputs.setup.elist_base .= e
	
	# Create the final dmat (dispersal rate multipliers) by multiplying the various input dmats
	inputs.setup.dmat .= inputs.setup.dmat_base .* inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	inputs.setup.amat .= inputs.setup.amat_base .* inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	inputs.setup.elist .= inputs.setup.elist_base .* inputs.setup.area_of_areas.^u
	
	# jmat does not include the "d" parameter
	inputs.setup.jmat .= inputs.setup.dispersal_multipliers_mat.^w .* inputs.setup.distmat.^x .* inputs.setup.envdistmat.^n .* inputs.setup.distmat2.^x2 .* inputs.setup.distmat3.^x3
	
	
	# Now update the p_Ds_v5 (the rates) for Q matrix
	p_Ds_v5 = update_Qij_vals2!(p_Ds_v5, inputs.setup.areas_list, inputs.setup.states_list, inputs.setup.dmat, inputs.setup.elist, inputs.setup.amat; return_df=false);

	# Update the mus
	p_Ds_v5.params.mu_vals[:] .= inputs.bmo.est[inputs.bmo.rownames .== "deathRate"][1]

	#print("\n")
	#print(prtQp(p_Ds_v5))
	
	# Now update the p_Ds_v5 (the rates) for C matrix
	if printlevel >= 2
		pdf_before = prtCp(p_Ds_v5)
		print("\np_Ds_v5_updater_v1, before:")
		print(round.(pdf_before.wt[1:5], digits=3))
	end
	update_Cijk_vals2!(p_Ds_v5, inputs.setup.areas_list, inputs.setup.states_list, inputs.bmo, inputs.setup.maxent01, inputs.setup.jmat);
	if printlevel >= 2
		pdf_after = prtCp(p_Ds_v5)
		print("\np_Ds_v5_updater_v1, after:")
		print(round.(pdf_before.wt[1:5], digits=3))
	end
	
	#return p_Ds_v5	
end # END function p_Ds_v5_updater_v1!()










end # Ends the module command
