module TreePass
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/
using BioGeoJulia.TrUtils # for e.g. flat2
using BioGeoJulia.TreeTable
using BioGeoJulia.SSEs 
#using BioGeoJulia.Flow 
using DataFrames			# for e.g. DataFrame()
using PhyloNetworks		# for e.g. readTopology()
using Dates						# for e.g. DateTime, Dates.now()
using Distributed			# for e.g. @spawn
using Random					# for MersenneTwister()
using DifferentialEquations # for ODEProblem
using LSODA						# for lsoda()
using LinearAlgebra		# for factorize()
using SpecialFunctions			# for e.g. logfactorial

print("\nBioGeoJulia: loading TreePass.jl")

export SolverOpt, construct_SolverOpt, Res, construct_Res, count_nodes_finished, nodeOp_average_likes, nodeOp, nodeOp_Cmat, nodeOp_Cmat2, nodeOp_singleton!, nodeOp_ClaSSE_v5!, nodeOp_ClaSSE_v6!, branchOp_example, branchOp_ClaSSE_Gflow_v1, find_time_in_time_segments, branchOp_ClaSSE_Ds_v5, branchOp, setup_inputs_branchOp_ClaSSE_Ds_v5, countloop, iterative_downpass!, iterative_downpass_Gflow_nonparallel_v1!, iterative_downpass_Gflow_nonparallel_v2!, iterative_downpass_nonparallel_ClaSSE_v5!, iterative_downpass_nonparallel_ClaSSE_v6!, iterative_downpass_nonparallel!  # branchOp_ClaSSE_Gflow_v2, 








#######################################################
# Parallel operations on binary trees
#######################################################
# Threaded downpass that spawns new processes when the 2 nodes above are done.

# Solver options structure
# "mutable" means you can change the values referred to by the keys
mutable struct SolverOpt
	solver::Any
	save_everystep::Bool
	abstol::Float64
	reltol::Float64
end

function construct_SolverOpt()
	solver = Tsit5()
	save_everystep = true
	abstol = 1e-6
	reltol = 1e-6
	solver_options = SolverOpt(solver, save_everystep, abstol, reltol)
	return solver_options
end


# Results structure
struct Res
	# The "regime" is just a number, indicating which
	# ClaSSE model is operating for this branch
	# Ideally, the regime could have 
	# - its own list of states
	# - the ClaSSE model parameters, and/or functions determining the same
	# - functions for converting states to the other regimes,
	#   forwards and backwards in time
	regime::Array{Int64}


	# The states can be 
	# "not_ready" (value at branch top NOT available)
	# "ready_for_nodeOp" (values at branches above are ready)
	# "ready_for_branchOp" (value at branch top available)
	# "calculating" (value at branch bottom being calculated)
	# "done" (value at branch bottom available)
	node_state::Array{String}
	node_Lparent_state::Array{String}
	node_Rparent_state::Array{String}
	
	# Tree structure
	root_nodeIndex::Int64
	numNodes::Int64
	uppass_edgematrix::Array{Int64}

	# Calculation timing & nodes/threads
	thread_for_each_nodeOp::Array{Int64}
	thread_for_each_branchOp::Array{Int64}
	calc_spawn_start::Array{DateTime}
	calc_start_time::Array{DateTime}
	calc_end_time::Array{DateTime}
	calc_duration::Array{Float64}
	calctime_iterations::Array{Float64}

	# Likelihood calculations
	sumLikes_at_node_at_branchTop::Array{Float64}
	lnL_at_node_at_branchTop::Array{Float64}
	lq_at_branchBot::Array{Float64}
	like_at_branchBot::Array{Float64}

	Es_at_each_nodeIndex_branchTop::Array{Array{Float64,1},1}
	Es_at_each_nodeIndex_branchBot::Array{Array{Float64,1},1}
	fakeX0s_at_each_nodeIndex_branchTop::Array{Array{Float64,1},1}
	likes_at_each_nodeIndex_branchTop::Array{Array{Float64,1},1}
	normlikes_at_each_nodeIndex_branchTop::Array{Array{Float64,1},1}
	likes_at_each_nodeIndex_branchBot::Array{Array{Float64,1},1}
	normlikes_at_each_nodeIndex_branchBot::Array{Array{Float64,1},1}
end

# Construct a default, simple results structure
# (likes_at_each_nodeIndex_branchTop is an array)
function construct_Res_old()
	# The "regime" is just a number, indicating which
	# ClaSSE model is operating for this branch
	# Ideally, the regime could have 
	# - its own list of states
	# - the ClaSSE model parameters, and/or functions determining the same
	# - functions for converting states to the other regimes,
	#   forwards and backwards in time
	regime = [1 1 1 1 1 1 1]

	n = 1 # number of states
	node_state = ["ready_for_branchOp", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready"]
	node_Lparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	node_Rparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	root_nodeIndex = 7
	numNodes = 7
	uppass_edgematrix = [7 6; 7 5; 5 4; 5 3; 3 2; 3 1]

	thread_for_each_nodeOp = collect(repeat([0], numNodes))
	thread_for_each_branchOp = collect(repeat([0], numNodes))
	calc_spawn_start = collect(repeat([Dates.now()], numNodes))
	calc_start_time = collect(repeat([Dates.now()], numNodes))
	calc_end_time = collect(repeat([Dates.now()], numNodes))
	calc_duration = collect(repeat([0.0], numNodes))
	calctime_iterations = [0.0, 0.0]

	sumLikes_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lnL_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lq_at_branchBot = collect(repeat([0.0], numNodes))
	like_at_branchBot = collect(repeat([0.0], numNodes))

	Es_at_each_nodeIndex_branchTop = collect(repeat([0.0], numNodes))
	Es_at_each_nodeIndex_branchBot = collect(repeat([0.0], numNodes))
	fakeX0s_at_each_nodeIndex_branchTop = collect(repeat([0.0], numNodes))

	likes_at_each_nodeIndex_branchTop = collect(repeat([1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0], n))
	normlikes_at_each_nodeIndex_branchTop = likes_at_each_nodeIndex_branchTop ./ sum(likes_at_each_nodeIndex_branchTop)
	likes_at_each_nodeIndex_branchBot = collect(repeat([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], n))
	normlikes_at_each_nodeIndex_branchBot = collect(repeat([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], n))

	res = Res(regime, node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, 
thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations, sumLikes_at_node_at_branchTop, lnL_at_node_at_branchTop, lq_at_branchBot, like_at_branchBot,  Es_at_each_nodeIndex_branchTop, Es_at_each_nodeIndex_branchBot, fakeX0s_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, normlikes_at_each_nodeIndex_branchBot)
	return res
end


# Construct a default, simple results structure
# (likes_at_each_nodeIndex_branchTop is an array of arrays)
function construct_Res()
	n = 1 # number of states
	numNodes = 7  # number of nodes

	# The "regime" is just a number, indicating which
	# ClaSSE model is operating for this branch (each branch is below a node)
	# Ideally, the regime could have 
	# - its own list of states
	# - the ClaSSE model parameters, and/or functions determining the same
	# - functions for converting states to the other regimes,
	#   forwards and backwards in time
	regime = collect(repeat([1], numNodes))

	node_state = ["ready_for_branchOp", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready"]
	node_Lparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	node_Rparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	root_nodeIndex = 7
	uppass_edgematrix = [7 6; 7 5; 5 4; 5 3; 3 2; 3 1]
	likes_OneNode = collect(repeat([0.0], n))
# 	likes_at_each_nodeIndex_branchTop = repeat([likes_OneNode], numNodes)
# 	likes_at_each_nodeIndex_branchBot = repeat([likes_OneNode], numNodes)
# 	normlikes_at_each_nodeIndex_branchTop = repeat([likes_OneNode], numNodes)
# 	normlikes_at_each_nodeIndex_branchBot = repeat([likes_OneNode], numNodes)
	Es_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	Es_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	fakeX0s_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	likes_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	likes_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	normlikes_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	normlikes_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	
	default_likes_at_each_nodeIndex_branchTop = [1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0]
	for i in 1:length(likes_at_each_nodeIndex_branchTop)
		likes_at_each_nodeIndex_branchTop[i][1] = default_likes_at_each_nodeIndex_branchTop[i]
	end
	#typeof(likes_at_each_nodeIndex_branchTop)
	sumLikes_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lnL_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lq_at_branchBot = collect(repeat([0.0], numNodes))
	like_at_branchBot = collect(repeat([0.0], numNodes))


	thread_for_each_nodeOp = collect(repeat([0], 7))
	thread_for_each_branchOp = collect(repeat([0], 7))

	calc_spawn_start = collect(repeat([Dates.now()], 7))
	calc_start_time = collect(repeat([Dates.now()], 7))
	calc_end_time = collect(repeat([Dates.now()], 7))
	calc_duration = collect(repeat([0.0], 7))

	calctime_iterations = [0.0, 0.0]

	res = Res(regime, node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, 
thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations, sumLikes_at_node_at_branchTop, lnL_at_node_at_branchTop, lq_at_branchBot, like_at_branchBot,  Es_at_each_nodeIndex_branchTop, Es_at_each_nodeIndex_branchBot, fakeX0s_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, normlikes_at_each_nodeIndex_branchBot)
	return res
end





"""
# Load a simple tree, see a simple list of nodeIndexes
using DataFrames
using PhyloNetworks

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr.root
# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

"""
function construct_Res(tr::HybridNetwork)
	root_nodeIndex = tr.root
	numNodes = tr.numNodes
	uppass_edgematrix = get_LR_uppass_edgematrix(tr)
	
	# Give tip nodeIndexes their nodeNodex as the "likelihood"
	indexNum_table = get_nodeIndex_PNnumber(tr)
	tipsTF = indexNum_table[:,2] .> 0
	tipLikes = indexNum_table[tipsTF,2] * 1.0

	# The "regime" is just a number, indicating which
	# ClaSSE model is operating for this branch (each branch is below a node)
	# Ideally, the regime could have 
	# - its own list of states
	# - the ClaSSE model parameters, and/or functions determining the same
	# - functions for converting states to the other regimes,
	#   forwards and backwards in time
	regime = collect(repeat([1], numNodes))

	# Set up an array of length nstates (n), to hold the likelihoods for each node
	n = 1
	blank_states = collect(repeat([0.0], n))
# 	likes_at_each_nodeIndex_branchTop = collect(repeat([blank_states], numNodes))
# 	likes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
# 	normlikes_at_each_nodeIndex_branchTop = collect(repeat([blank_states], numNodes))
# 	normlikes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
	Es_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	Es_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	fakeX0s_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)

	likes_at_each_nodeIndex_branchTop = collect(repeat([collect(repeat([0.0], n))], numNodes))
	likes_at_each_nodeIndex_branchBot = collect(repeat([collect(repeat([0.0], n))], numNodes))
	normlikes_at_each_nodeIndex_branchTop = collect(repeat([collect(repeat([0.0], n))], numNodes))
	normlikes_at_each_nodeIndex_branchBot = collect(repeat([collect(repeat([0.0], n))], numNodes))

	# Put in the tip node numbers as the fake likelihoods
	function f(numNodes, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, fakeX0s_at_each_nodeIndex_branchTop, tipsTF)
		j = 0
		for i in 1:numNodes
			if (tipsTF[i] == true)
				j = j+1
				# Transfer from 1D array to 1D array
				likes_at_each_nodeIndex_branchTop[i] = [tipLikes[j]]
				normlikes_at_each_nodeIndex_branchTop[i] = [tipLikes[j] / sum(tipLikes[j])]
				fakeX0s_at_each_nodeIndex_branchTop[i] = [tipLikes[j] / sum(tipLikes[j])]
			end
		end
	end
	# Run function f()
	f(numNodes, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, fakeX0s_at_each_nodeIndex_branchTop, tipsTF)
	
	likes_at_each_nodeIndex_branchTop
	sumLikes_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lnL_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lq_at_branchBot = collect(repeat([0.0], numNodes))
	like_at_branchBot = collect(repeat([0.0], numNodes))
	
	# Fill in the node_states
	node_state = collect(repeat(["not_ready"], numNodes))
	node_state[tipsTF] .= "ready_for_branchOp"
	node_Lparent_state = collect(repeat(["not_ready"], numNodes))
	node_Rparent_state = collect(repeat(["not_ready"], numNodes))
	node_Lparent_state[tipsTF] .= "NA"
	node_Rparent_state[tipsTF] .= "NA"
	
	# Initialize with zeros for the other items
	#likes_at_each_nodeIndex_branchBot = collect(repeat([0.0], numNodes))
	likes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
	thread_for_each_nodeOp = collect(repeat([0], numNodes))
	thread_for_each_branchOp = collect(repeat([0], numNodes))
	
	calc_spawn_start = collect(repeat([Dates.now()], numNodes))
	calc_start_time = collect(repeat([Dates.now()], numNodes))
	calc_end_time = collect(repeat([Dates.now()], numNodes))
	calc_duration = collect(repeat([0.0], numNodes))

	calctime_iterations = [0.0, 0.0]
	number_of_whileLoop_iterations = [0]	

	# Initialize res object
	res = Res(regime, node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, 
thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations, sumLikes_at_node_at_branchTop, lnL_at_node_at_branchTop, lq_at_branchBot, like_at_branchBot,  Es_at_each_nodeIndex_branchTop, Es_at_each_nodeIndex_branchBot, fakeX0s_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, normlikes_at_each_nodeIndex_branchBot)
	return res
end





"""
# Load a simple tree, see a simple list of nodeIndexes
using DataFrames
using PhyloNetworks

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr.root
# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
n=10 # number of states in the state space

"""
function construct_Res(tr::HybridNetwork, n)
	root_nodeIndex = tr.root
	numNodes = tr.numNodes
	uppass_edgematrix = get_LR_uppass_edgematrix(tr)
	

	# The "regime" is just a number, indicating which
	# ClaSSE model is operating for this branch (each branch is below a node)
	# Ideally, the regime could have 
	# - its own list of states
	# - the ClaSSE model parameters, and/or functions determining the same
	# - functions for converting states to the other regimes,
	#   forwards and backwards in time
	regime = collect(repeat([1], numNodes))

	
	# Set up an array of length nstates (n), to hold the likelihoods for each node
	blank_states = collect(repeat([0.0], n))
# 	likes_at_each_nodeIndex_branchTop = collect(repeat([blank_states], numNodes))
# 	likes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
# 	normlikes_at_each_nodeIndex_branchTop = collect(repeat([blank_states], numNodes))
# 	normlikes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
	Es_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	Es_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	fakeX0s_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	likes_at_each_nodeIndex_branchTop = collect(repeat([collect(repeat([0.0], n))], numNodes))
	likes_at_each_nodeIndex_branchBot = collect(repeat([collect(repeat([0.0], n))], numNodes))
	normlikes_at_each_nodeIndex_branchTop = collect(repeat([collect(repeat([0.0], n))], numNodes))
	normlikes_at_each_nodeIndex_branchBot = collect(repeat([collect(repeat([0.0], n))], numNodes))

	# Give tip nodeIndexes a likelihood of 1 at all states
	indexNum_table = get_nodeIndex_PNnumber(tr)
	tipsTF = indexNum_table[:,2] .> 0
	tipnums = seq(1, length(tipsTF), 1)[tipsTF]
	
# 	put_in_fake_tipLikes = false
# 	if put_in_fake_tipLikes == true
# 		for i in 1:length(tipnums)
# 			tipLikes = collect(repeat([1.0], n))
# 			likes_at_each_nodeIndex_branchTop[tipnums[i]] = tipLikes
# 			normlikes_at_each_nodeIndex_branchTop[i] = [tipLikes[j] / sum(tipLikes[j])]
# 		end
# 	end
	sumLikes_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lnL_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lq_at_branchBot = collect(repeat([0.0], numNodes))
	like_at_branchBot = collect(repeat([0.0], numNodes))
	
	# Fill in the node_states
	node_state = collect(repeat(["not_ready"], numNodes))
	node_state[tipsTF] .= "ready_for_branchOp"
	node_Lparent_state = collect(repeat(["not_ready"], numNodes))
	node_Rparent_state = collect(repeat(["not_ready"], numNodes))
	node_Lparent_state[tipsTF] .= "NA"
	node_Rparent_state[tipsTF] .= "NA"
	
	# Initialize with zeros for the other items
	thread_for_each_nodeOp = collect(repeat([0], numNodes))
	thread_for_each_branchOp = collect(repeat([0], numNodes))
	
	calc_spawn_start = collect(repeat([Dates.now()], numNodes))
	calc_start_time = collect(repeat([Dates.now()], numNodes))
	calc_end_time = collect(repeat([Dates.now()], numNodes))
	calc_duration = collect(repeat([0.0], numNodes))

	calctime_iterations = [0.0, 0.0]
	number_of_whileLoop_iterations = [0]	

	# Initialize res object
	res = Res(regime, node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, 
thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations, sumLikes_at_node_at_branchTop, lnL_at_node_at_branchTop, lq_at_branchBot, like_at_branchBot,  Es_at_each_nodeIndex_branchTop, Es_at_each_nodeIndex_branchBot, fakeX0s_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, normlikes_at_each_nodeIndex_branchBot)
	return res
end # END function construct_Res(tr::HybridNetwork, n)




# Convert this res object to a DataFrame
#function res_to_df(res)
#	
#end

function count_nodes_finished(node_state)
	sum(node_state .== "done")
end

# Average the downpass likelihoods at a node
# (default)
function nodeOp_average_likes(tmp1, tmp2)
	nodeData_at_top = (tmp1 + tmp2)/2
	return(nodeData_at_top)
end

# Use the Cmat to combine the likelihoods
nodeOp_Cmat = (tmpDs; tmp1, tmp2, p_Ds_v5) -> begin
	p = p_Ds_v5
#	hcat(p.p_indices.Carray_ivals, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals, p.params.Cijk_vals)
	
	# Go through each Ci (ancestral state index)

 # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
# 	print("\n")
# 	print("\n")
# 	print("Running nodeOp_Cmat:\n\n")
	
	# Calculate likelihoods of states just before speciation
  @inbounds for i in 1:n
  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]

# 		Calculation of "D" (likelihood of tip data)
# 		for state i=1, multiply
# 		
# 		All the speciation rates where the ancestor = i; length is length(Ci_sub_i)
# 		Cijk_vals[Ci_sub_i]
# 		
# 		All the ancestors (i's) where the ancestor = i; length is length(Ci_sub_i)
# 		Carray_ivals[Ci_sub_i]
# 
# 		All the left descendant states (j's) where the ancestor just before speciation = i; length is sum(Ci_sub_i)
# 		Carray_jvals[Ci_sub_i]
# 		
# 		All the right descendant states (k's) where the ancestor just before speciation = i; length is sum(Ci_sub_i)
# 		Carray_kvals[Ci_sub_i]
# 		
# 		Coming down from left branch, contributing to likelihood of ancestor state i;
# 		resulting length is sum(Ci_sub_i)
# 		tmp1[Carray_jvals[Ci_sub_i]]
# 
# 		Coming down from left branch, contributing to likelihood of ancestor state i;
# 		resulting length is sum(Ci_sub_i)
# 		tmp2[Carray_kvals[Ci_sub_i]]
		
		# Equivalent:
		# inputs.p_Ds_v5.params.Cijk_vals[inputs.p_Ds_v5.p_TFs.Ci_eq_i[i]]
		# 6-element Array{Float64,1}:
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		# 
		# julia> inputs.p_Ds_v5.p_TFs.Ci_sub_i[i]
		# 6-element Array{Int64,1}:
		#  3
		#  3
		#  3
		#  3
		#  3
		#  3
		# 
		# julia> inputs.p_Ds_v5.params.Cijk_vals[ inputs.p_Ds_v5.p_TFs.Ci_sub_i[i]]
		# 6-element Array{Float64,1}:
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333		
		
		# Parameter values for these events with nonzero rates
		tmpDs[i] = sum(Cijk_vals[Ci_eq_i] .* tmp1[Cj_sub_i] .* tmp2[Ck_sub_i])
# 		print(tmpDs[i])
# 		print("\n")
# 		print(Cijk_vals)
# 		print("\n")
# 		print(Ci_sub_i)
# 		print("\n")
# 		print(Cj_sub_i)
# 		print("\n")
# 		print(Ck_sub_i)
# 		print("\n")
# 		print(tmp1)
# 		print("\n")
# 		print(tmp2)
# 		print("\n")
  end
  return(tmpDs)
end




# Re-do nodeOp_Cmat for when i,j,k and i,k,j events are collapsed
nodeOp_Cmat2 = (tmpDs; tmp1, tmp2, p_Ds_v5) -> begin
	p = p_Ds_v5
#	hcat(p.p_indices.Carray_ivals, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals, p.params.Cijk_vals)
	
	# Go through each Ci (ancestral state index)

 # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	Carray_pair = p.p_indices.Carray_pair
	
# 	print("\n")
# 	print("\n")
# 	print("Running nodeOp_Cmat:\n\n")
	
	# Calculate likelihoods of states just before speciation
  @inbounds for i in 1:n
  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]

# 		Calculation of "D" (likelihood of tip data)
# 		for state i=1, multiply
# 		
# 		All the speciation rates where the ancestor = i; length is length(Ci_sub_i)
# 		Cijk_vals[Ci_sub_i]
# 		
# 		All the ancestors (i's) where the ancestor = i; length is length(Ci_sub_i)
# 		Carray_ivals[Ci_sub_i]
# 
# 		All the left descendant states (j's) where the ancestor just before speciation = i; length is sum(Ci_sub_i)
# 		Carray_jvals[Ci_sub_i]
# 		
# 		All the right descendant states (k's) where the ancestor just before speciation = i; length is sum(Ci_sub_i)
# 		Carray_kvals[Ci_sub_i]
# 		
# 		Coming down from left branch, contributing to likelihood of ancestor state i;
# 		resulting length is sum(Ci_sub_i)
# 		tmp1[Carray_jvals[Ci_sub_i]]
# 
# 		Coming down from left branch, contributing to likelihood of ancestor state i;
# 		resulting length is sum(Ci_sub_i)
# 		tmp2[Carray_kvals[Ci_sub_i]]
		
		# Equivalent:
		# inputs.p_Ds_v5.params.Cijk_vals[inputs.p_Ds_v5.p_TFs.Ci_eq_i[i]]
		# 6-element Array{Float64,1}:
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		# 
		# julia> inputs.p_Ds_v5.p_TFs.Ci_sub_i[i]
		# 6-element Array{Int64,1}:
		#  3
		#  3
		#  3
		#  3
		#  3
		#  3
		# 
		# julia> inputs.p_Ds_v5.params.Cijk_vals[ inputs.p_Ds_v5.p_TFs.Ci_sub_i[i]]
		# 6-element Array{Float64,1}:
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333		
		
		# Parameter values for these events with nonzero rates
		#tmpDs[i] = sum(Cijk_vals[Ci_eq_i] .* tmp1[Cj_sub_i] .* tmp2[Ck_sub_i])
		
		# Divide by the "pair" value (1.0 or 2.0), do calc again switching j and k
		# (Carray_pair[Ci_eq_i] .== 2)   # This gives JUST the pair=2, for the 2nd sum,
		#                                # in order to avoid double-counting pair=1 events
		yTFs = Carray_pair[Ci_eq_i] .== 1
		ysums = sum(Cijk_vals[Ci_eq_i][yTFs] .* tmp1[Cj_sub_i][yTFs] .* tmp2[Ck_sub_i][yTFs]) 
		
		#nony_TFs = Carray_pair[Ci_eq_i] .== 2
		#nony_sums1 = sum(Cijk_vals[Ci_eq_i][nony_TFs] .* tmp1[Cj_sub_i][nony_TFs] .* tmp2[Ck_sub_i][nony_TFs]) / 2.0
		#nony_sums2 = sum(Cijk_vals[Ci_eq_i][nony_TFs] .* tmp1[Ck_sub_i][nony_TFs] .* tmp2[Cj_sub_i][nony_TFs]) / 2.0
		#tmpDs[i] = ysums + nony_sums1 + nony_sums2
		tmpDs[i] = sum(Cijk_vals[Ci_eq_i]./Carray_pair[Ci_eq_i] .* tmp1[Cj_sub_i] .* tmp2[Ck_sub_i]) + sum((Carray_pair[Ci_eq_i] .== 2) .* Cijk_vals[Ci_eq_i]./Carray_pair[Ci_eq_i] .* tmp1[Ck_sub_i] .* tmp2[Cj_sub_i])

# 		print(tmpDs[i])
# 		print("\n")
# 		print(Cijk_vals)
# 		print("\n")
# 		print(Ci_sub_i)
# 		print("\n")
# 		print(Cj_sub_i)
# 		print("\n")
# 		print(Ck_sub_i)
# 		print("\n")
# 		print(tmp1)
# 		print("\n")
# 		print(tmp2)
# 		print("\n")
  end
  return(tmpDs)
end # END nodeOp_Cmat2 = (tmpDs; tmp1, tmp2, p_Ds_v5) -> begin


	




# Combine likelihoods from above
function nodeOp(current_nodeIndex, res; nodeOp_function=nodeOp_average_likes)
	res.node_state[current_nodeIndex] = "calculating_nodeOp"
	uppass_edgematrix = res.uppass_edgematrix
	
	# Record the thread, for kicks
	tmp_threadID = Threads.threadid()
	res.thread_for_each_nodeOp[current_nodeIndex] = tmp_threadID
	TF = uppass_edgematrix[:,1] .== current_nodeIndex
	if (sum(TF) == 2)
		# Get likelihoods from above (iterates up to tips)
		parent_nodeIndexes = uppass_edgematrix[TF,2]

		# Crucial
# 		tmp1 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
# 		tmp2 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		tmp1 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		tmp2 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]

		# Check that data are actually available
		if (sum(tmp1) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp1) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		if (sum(tmp2) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp2) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		#nodeData_at_top = tmp1 + tmp2
		#nodeData_at_top = (tmp1 + tmp2)/2
		nodeData_at_top = nodeOp_function(tmp1, tmp2)
		
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] = nodeData_at_top
		
		# Check if it's the root node
		if (current_nodeIndex == res.root_nodeIndex)
			res.node_state[current_nodeIndex] = "done"
		else
			res.node_state[current_nodeIndex] = "ready_for_branchOp"
		end
		return()
	elseif (sum(TF) == 0)
	  # If a tip
	  txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): shouldn't be run on a tip node."], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	else
	  txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(TF) should be 0 or 2"], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	end
	txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): shouldn't get here."], "")
	print("\n")
	print(txt)
	print("\n")
	return(error(txt))
end # END function nodeOp(current_nodeIndex, res; nodeOp_function=nodeOp_average_likes)

# Example nodeOp when the node is a singleton (one descendant)
# This one just passes the lnLs down, but others could be imagined
# (e.g. for an area appearing/disappearing)
function nodeOp_singleton!(current_nodeIndex, res; p_Ds_v5)
	res.node_state[current_nodeIndex] = "calculating_nodeOp"
	uppass_edgematrix = res.uppass_edgematrix

	# Record the thread, for kicks
	tmp_threadID = Threads.threadid()
	res.thread_for_each_nodeOp[current_nodeIndex] = tmp_threadID
	TF = uppass_edgematrix[:,1] .== current_nodeIndex
	
	# Singleton nodes
	if (sum(TF) == 1)
		# Get likelihoods from above (iterates up to tips)
		parent_nodeIndexes = uppass_edgematrix[TF,2]
		#tmp1 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		#tmp2 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		# Crucial
		# The likelihoods were taken out to produce normlikes, stored in "res.lq_at_branchBot"
		tmp1 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		#tmp2 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		#combined_branch_lnLs = res.lq_at_branchBot[parent_nodeIndexes[1]] + res.lq_at_branchBot[parent_nodeIndexes[2]]
		
		
		# Check that data are actually available
		if (sum(tmp1) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp1) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end


		#nodeData_at_top = tmp1 + tmp2
		#nodeData_at_top = (tmp1 + tmp2)/2
		#nodeData_at_top = nodeOp_function(tmp1, tmp2)
		
		# Singleton
		#nodeData_at_top = nodeOp_Cmat(nodeData_at_top, tmp1=tmp1, tmp2=tmp2, p_Ds_v5=p_Ds_v5)
		# Not a placeholder, because we are just passing the likelihoods 
		# down from the bottom of the branch above (no Cmat calculation needed)
		nodeData_at_top = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]

		# Somehow adding .+ 0.0 individualizes the assignment!
		#sum_likes_at_node = sum(nodeData_at_top)
		sum_likes_at_node = sum(nodeData_at_top)
		#sum_likes_at_node = 1.0
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] = (nodeData_at_top .+ 0.0)
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
# 		print("\n\ncurrent_nodeIndex:")
# 		print(current_nodeIndex)
# 		print("\n")
		res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] = ((nodeData_at_top .+ 0.0) ./ sum_likes_at_node)
		
		res.sumLikes_at_node_at_branchTop[current_nodeIndex] = sum_likes_at_node + 0.0
		res.lnL_at_node_at_branchTop[current_nodeIndex] = log(sum_likes_at_node) + 0.0
# 		print("\nnodeData_at_top:\n")
# 		print(nodeData_at_top)
# 
# 		print("\nres.likes_at_each_nodeIndex_branchTop[current_nodeIndex]:\n")
# 		print(res.likes_at_each_nodeIndex_branchTop[current_nodeIndex])
# 
# 		print("\nres.likes_at_each_nodeIndex_branchTop:\n")
# 		print(res.likes_at_each_nodeIndex_branchTop)

		
		# Check if it's the root node
		if (current_nodeIndex == res.root_nodeIndex)
			res.node_state[current_nodeIndex] = "done"
		else
			res.node_state[current_nodeIndex] = "ready_for_branchOp"
		end
		return(res)  # END if singleton; 2021-07-18_NJM
	end # END if (sum(TF) == 1)
end # END function nodeOp_singleton!(current_nodeIndex, res; p_Ds_v5)




# Node operation, combining probabilities from above (assumed to be fast)
function nodeOp_ClaSSE_v5!(current_nodeIndex, res; p_Ds_v5)
	res.node_state[current_nodeIndex] = "calculating_nodeOp"
	uppass_edgematrix = res.uppass_edgematrix
	
	# Record the thread, for kicks
	tmp_threadID = Threads.threadid()
	res.thread_for_each_nodeOp[current_nodeIndex] = tmp_threadID
	TF = uppass_edgematrix[:,1] .== current_nodeIndex
	
	# Binary nodes
	if (sum(TF) == 2)
		# Get likelihoods from above (iterates up to tips)
		parent_nodeIndexes = uppass_edgematrix[TF,2]
		#tmp1 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		#tmp2 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		# Crucial
		# The likelihoods were taken out to produce normlikes, stored in "res.lq_at_branchBot"
		tmp1 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		tmp2 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		#combined_branch_lnLs = res.lq_at_branchBot[parent_nodeIndexes[1]] + res.lq_at_branchBot[parent_nodeIndexes[2]]
		
		
		# Check that data are actually available
		if (sum(tmp1) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp1) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		if (sum(tmp2) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp2) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		#nodeData_at_top = tmp1 + tmp2
		#nodeData_at_top = (tmp1 + tmp2)/2
		#nodeData_at_top = nodeOp_function(tmp1, tmp2)
# 		print("\n\n")
# 		
# 		print("\ncurrent_nodeIndex:\n")
# 		print(current_nodeIndex)

		nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] .* 0.0 # Placeholder
		nodeData_at_top = nodeOp_Cmat(nodeData_at_top, tmp1=tmp1, tmp2=tmp2, p_Ds_v5=p_Ds_v5)

		# Somehow adding .+ 0.0 individualizes the assignment!
		sum_likes_at_node = sum(nodeData_at_top)
		
		# Error trap 2022-03-10
		if sum_likes_at_node < 0.0
		 sum_likes_at_node = 1e-10000
		end
		
		#sum_likes_at_node = 1.0
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] = (nodeData_at_top .+ 0.0)
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
# 		print("\n\ncurrent_nodeIndex:")
# 		print(current_nodeIndex)
# 		print("\n")
		res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] = ((nodeData_at_top .+ 0.0) ./ sum_likes_at_node)
		
		res.sumLikes_at_node_at_branchTop[current_nodeIndex] = sum_likes_at_node + 0.0
		res.lnL_at_node_at_branchTop[current_nodeIndex] = log(sum_likes_at_node) + 0.0
# 		print("\nnodeData_at_top:\n")
# 		print(nodeData_at_top)
# 
# 		print("\nres.likes_at_each_nodeIndex_branchTop[current_nodeIndex]:\n")
# 		print(res.likes_at_each_nodeIndex_branchTop[current_nodeIndex])
# 
# 		print("\nres.likes_at_each_nodeIndex_branchTop:\n")
# 		print(res.likes_at_each_nodeIndex_branchTop)

		
		# Check if it's the root node
		if (current_nodeIndex == res.root_nodeIndex)
			res.node_state[current_nodeIndex] = "done"
		else
			res.node_state[current_nodeIndex] = "ready_for_branchOp"
		end
		return(res)  # END if binary node; 2020-08-04_NJM
	elseif (sum(TF) == 1)
	  # If a singleton
	  txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): shouldn't be run on a singleton node (only 1 descendant). The function nodeOp_singleton() should be used instead."], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	elseif (sum(TF) == 0)
	  # If a tip
	  txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): shouldn't be run on a tip node."], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	else
	  txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): sum(TF) should be 0 or 2"], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	end # End if (sum(TF) == 2)
	txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): shouldn't get here."], "")
	print("\n")
	print(txt)
	print("\n")
	return(error(txt))
end # END function nodeOp_ClaSSE_v5!(current_nodeIndex, res; p_Ds_v5)





# Node operation, combining probabilities from above (assumed to be fast)
function nodeOp_ClaSSE_v6!(current_nodeIndex, res; p_Ds_v5)
	res.node_state[current_nodeIndex] = "calculating_nodeOp"
	uppass_edgematrix = res.uppass_edgematrix
	
	# Record the thread, for kicks
	tmp_threadID = Threads.threadid()
	res.thread_for_each_nodeOp[current_nodeIndex] = tmp_threadID
	TF = uppass_edgematrix[:,1] .== current_nodeIndex
	
	# Binary nodes
	if (sum(TF) == 2)
		# Get likelihoods from above (iterates up to tips)
		parent_nodeIndexes = uppass_edgematrix[TF,2]
		#tmp1 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		#tmp2 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		# Crucial
		# The likelihoods were taken out to produce normlikes, stored in "res.lq_at_branchBot"
		tmp1 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		tmp2 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		#combined_branch_lnLs = res.lq_at_branchBot[parent_nodeIndexes[1]] + res.lq_at_branchBot[parent_nodeIndexes[2]]
		
		
		# Check that data are actually available
		if (sum(tmp1) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp1) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		if (sum(tmp2) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp2) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		#nodeData_at_top = tmp1 + tmp2
		#nodeData_at_top = (tmp1 + tmp2)/2
		#nodeData_at_top = nodeOp_function(tmp1, tmp2)
# 		print("\n\n")
# 		
# 		print("\ncurrent_nodeIndex:\n")
# 		print(current_nodeIndex)

		nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] .* 0.0 # Placeholder
		nodeData_at_top = nodeOp_Cmat2(nodeData_at_top, tmp1=tmp1, tmp2=tmp2, p_Ds_v5=p_Ds_v5)

		# Somehow adding .+ 0.0 individualizes the assignment!
		sum_likes_at_node = sum(nodeData_at_top)
		
		# Error trap 2022-03-10
		if sum_likes_at_node < 0.0
		 sum_likes_at_node = 1e-10000
		end
		
		#sum_likes_at_node = 1.0
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] = (nodeData_at_top .+ 0.0)
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
# 		print("\n\ncurrent_nodeIndex:")
# 		print(current_nodeIndex)
# 		print("\n")
		res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] = ((nodeData_at_top .+ 0.0) ./ sum_likes_at_node)
		
		res.sumLikes_at_node_at_branchTop[current_nodeIndex] = sum_likes_at_node + 0.0
		res.lnL_at_node_at_branchTop[current_nodeIndex] = log(sum_likes_at_node) + 0.0
# 		print("\nnodeData_at_top:\n")
# 		print(nodeData_at_top)
# 
# 		print("\nres.likes_at_each_nodeIndex_branchTop[current_nodeIndex]:\n")
# 		print(res.likes_at_each_nodeIndex_branchTop[current_nodeIndex])
# 
# 		print("\nres.likes_at_each_nodeIndex_branchTop:\n")
# 		print(res.likes_at_each_nodeIndex_branchTop)

		
		# Check if it's the root node
		if (current_nodeIndex == res.root_nodeIndex)
			res.node_state[current_nodeIndex] = "done"
		else
			res.node_state[current_nodeIndex] = "ready_for_branchOp"
		end
		return(res)  # END if binary node; 2020-08-04_NJM
	elseif (sum(TF) == 1)
	  # If a singleton
	  txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): shouldn't be run on a singleton node (only 1 descendant). The function nodeOp_singleton() should be used instead."], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	elseif (sum(TF) == 0)
	  # If a tip
	  txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): shouldn't be run on a tip node."], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	else
	  txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): sum(TF) should be 0 or 2"], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	end # End if (sum(TF) == 2)
	txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): shouldn't get here."], "")
	print("\n")
	print(txt)
	print("\n")
	return(error(txt))
end # END function nodeOp_ClaSSE_v6!(current_nodeIndex, res; p_Ds_v5)





# Calculate down a branch
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_example(current_nodeIndex, res; num_iterations=10000000)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Example slow operation
	y = countloop(num_iterations, current_nodeIndex)

	nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	nodeData_at_bottom = nodeData_at_top / 2.0

	return(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time)
end




# Calculate Ds down a branch, using Gflow
#
# Modifies branchOp to do Ds calculation down a branch, using Gflows
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_ClaSSE_Gflow_v1(current_nodeIndex, res, Gflow; tspan, p_Ds_v5, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# We are calculating from the current_nodeIndex node to the bottom of the branch below it
	
	# Get the fakeX0 at current_nodeIndex
	if (tspan[1] == 0.0)
		fakeX0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	else
		# Given res, update the fakeX0 in res at this node
		#age_branchtop = trdf[current_nodeIndex, :node_age]
		age_branchtop = tspan[1]
		
		# Calculate the fakeX0 for use calculating Xp at the rootward end of this node's branch
		tmpX0s = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] # assumes normalized likelihoods
		
		# Key trick from Louca & Pennell
		fakeX0 = factorize(Gflow(age_branchtop)) \ tmpX0s
		res.fakeX0s_at_each_nodeIndex_branchTop[current_nodeIndex] = fakeX0
	end
		
	# Solve for the likelihoods at the rootward end of the branch
	# (instantaneously after the speciation event at Xp, typically)
	sol_Ds = Gflow(tspan[2]) * fakeX0
	# If the high condition number requires sub-segments and re-normalizing
	# Could replace with:
	# sol_Ds = Gflow(seg1) * Gflow(seg2) * Gflow(seg3) * fakeX0
	
	# Return a tuple of the results for updating res in the main function
	return(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)
end # END branchOp_ClaSSE_Gflow_v1



# Calculate Ds down a branch
#
# Modifies branchOp to do Ds calculation down a branch
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_ClaSSE_Ds_v5(current_nodeIndex, res; u0, tspan, p_Ds_v5, solver_options=solver_options)
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
	
	prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
	sol_Ds = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol)
	
# 	print(sol_Ds[length(sol_Ds)])
# 	print("\n")
	
	#nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	#nodeData_at_bottom = nodeData_at_top / 2.0
	#nodeData_at_bottom = sol_Ds.u
	
	return(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)
end








"""
	Set up inputs
"""

function setup_inputs_branchOp_ClaSSE_Ds_v5(u0, tspan, p_Ds_v5; solver="Tsit5()", 
				 save_everystep="false", abstol="1e-9", reltol="1e-9")
	
	prob_str = "prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)"
	solve_str = join(["sol_Ds = solve(prob_Ds_v5, ", solver, ", save_everystep=", save_everystep, ", abstol=", abstol, ", reltol=", reltol, ")"])
	store_str = "nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)]"
	
	# Assemble the NamedTuple of inputs
	inputs = (u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver=solver, save_everystep=save_everystep, abstol=abstol, reltol=reltol, prob_str=prob_str, solve_str=solve_str, store_str=store_str)
	return inputs
end


# Calculate Ds down a branch
#
# Modifies branchOp to do Ds calculation down a branch
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp(current_nodeIndex, res, inputs)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	
	# The old practice input was an Int64
	if (typeof(inputs) != Int64)
		u0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
		tspan = inputs.tspan
		p_Ds_v5 = inputs.p_Ds_v5
		#solver = inputs.solver
		#save_everystep = inputs.save_everystep
		#abstol = inputs.abstol
		#reltol = inputs.reltol


		# Example slow operation
		#y = countloop(num_iterations, current_nodeIndex)
		#	prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
		#	sol_Ds = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
		#print(u0)
		
		# It DOESN'T really work to set the functions by text;
		# global variable problems, etc.
		# prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
		# sol_Ds = solve(prob_Ds_v5, Tsit5(), save_everystep=false, abstol=1e-9, reltol=1e-9)
		# nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)]
		eval(Meta.parse(inputs.prob_str))
		eval(Meta.parse(inputs.solve_str))
		eval(Meta.parse(inputs.store_str))
		
	else
		nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
		nodeData_at_bottom = nodeData_at_top / 2.0
		#nodeData_at_bottom = sol_Ds
	end

	
	return(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time)
end




function countloop(num_iterations, current_nodeIndex)
	x = 0.0
	random_number_generator = MersenneTwister(current_nodeIndex);

	for i in 1:num_iterations
	   x = x + (randn(random_number_generator, 1)[1] / num_iterations)
	end
	return(x)
end




"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_Gflow_nonparallel_v1!(res; trdf, p_Ds_v5, Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false)
	#######################################################
	# Use Gflow on a downpass
	#######################################################
	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()

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
			#u0 = u0 ./ (sum(u0))
			
			# Use the NORMALIZED (rescaled to sum to 1) likelihoods
			# If you get this error: claims an interpolation error for going beyond range
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
#				push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			tmp_results = branchOp_ClaSSE_Gflow_v1(current_nodeIndex, res, Gflow, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			#tmp_results = branchOp(current_nodeIndex, res, num_iterations)
			push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
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
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = tasks[i]
					
					# Differences for Gflow version
					#nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
					nodeData_at_bottom = sol_Ds .+ 0.0
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

					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
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
				#end
			end
		end
	
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
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end
	
		
		# Update which nodes have had BOTH parents complete
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
			
			#res = nodeOp_ClaSSE_v5!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			res = nodeOp_ClaSSE_v6!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
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
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	Julia_sum_lq = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])

	# Add the root probabilities

	# Assuming diversitree options:
	# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
	# i.e., the root state probs are just the root_Ds/sum(root_Ds)
	d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
	root_stateprobs = d_root_orig/sum(d_root_orig)
	rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
	
	# Without normalization
	#Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL
	# *With* normalization
	Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL + log(sum(res.sumLikes_at_node_at_branchTop[1:(length(res.sumLikes_at_node_at_branchTop)-1)]))
	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_Gflow_nonparallel_v1!






"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5, Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true)
	#######################################################
	# Use Gflow on a downpass
	#######################################################
	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()

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
			#u0 = u0 ./ (sum(u0))
			
			# Use the NORMALIZED (rescaled to sum to 1) likelihoods
			# If you get this error: claims an interpolation error for going beyond range
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
#				push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			tmp_results = branchOp_ClaSSE_Gflow_v1(current_nodeIndex, res, Gflow, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			#tmp_results = branchOp(current_nodeIndex, res, num_iterations)
			push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
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
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = tasks[i]
					
					# Differences for Gflow version
					#nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
					nodeData_at_bottom = sol_Ds .+ 0.0
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

					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
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
						txt = "STOP error in: iterative_downpass_Gflow_nonparallel_v2() - node has more than 2 edges"
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
				#end
			end
		end
	
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
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end
	
		
		# Update which nodes have had BOTH parents complete
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
			
			#res = nodeOp_ClaSSE_v5!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			res = nodeOp_ClaSSE_v6!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
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
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	# 2022-03-17_ update to new way of calculating
	#Julia_sum_lq = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])
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
	# Get basic tree info
	numTips = sum(trdf.nodeType .== "tip")
	numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")
	
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

	# Go back to BioGeoBEARS log-likelihood, under Yule process assumptions
	bgb2 = res5t - (bd.lnL - bd.lnl_topology) - equal_root_prob2
	bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_numtips_wOneMinusDeathRate + bd.lnl_branching_times) - equal_root_prob2
	bgb2 = res5t - (bd.lnl_numBirths + bd.lnl_Births_above_root + bd.lnl_branching_times) - equal_root_prob2
	bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) - -log(1/yuleBirthRate) - (bd.lnL - bd.lnl_topology)
	bgb_lnL = Julia_sum_lq + log(sum(d_root_orig)) + log(1/yuleBirthRate) - bd_lnL_noTopo

	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_Gflow_nonparallel_v2!







"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf, p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true)
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
#				push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			tmp_results = branchOp_ClaSSE_Ds_v5(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			#tmp_results = branchOp(current_nodeIndex, res, num_iterations)
			push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			#print("\ni: ")
			#print(i)
			if (tasks_fetched_TF[i] == false)
				#if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = tasks[i]
					nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
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
				#end
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
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
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
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			res = nodeOp_ClaSSE_v5!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
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
		print(txt)
		print("\n")
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
	
	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_nonparallel_ClaSSE_v5!







"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false, include_null_range=true)
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
#				push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			tmp_results = branchOp_ClaSSE_Ds_v5(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			#tmp_results = branchOp(current_nodeIndex, res, num_iterations)
			push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			#print("\ni: ")
			#print(i)
			if (tasks_fetched_TF[i] == false)
				#if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = tasks[i]
					nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
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
				#end
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
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			
			res = nodeOp_singleton!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
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
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			res = nodeOp_ClaSSE_v6!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
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
		print(txt)
		print("\n")
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
	
	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_nonparallel_ClaSSE_v6!















"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
"""
function iterative_downpass!(res; max_iterations=10^10, num_iterations=10000000)

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	if (parallel_TF == false)
		txt = "Error in iterative_downpass!(): This function probably requires multiple threads operating to consistently compile. Try starting julia with e.g. 'JULIA_NUM_THREADS=8 julia'."
		error(txt)
	end

	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()
	
	# Setup
	current_nodeIndex = res.root_nodeIndex
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

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
			push!(tasks_fetched_TF, false)
		end
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			if (tasks_fetched_TF[i] == false)
				if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
					(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
					
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
					res.lq_at_branchBot[spawned_nodeIndex] = log(sum_nodeData_at_bottom)
					res.like_at_branchBot[spawned_nodeIndex] = sum_nodeData_at_bottom

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
				end
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
			nodeOp(current_nodeIndex, res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
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
end # END iterative_downpass!




"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_nonparallel!(res; max_iterations=10^10, num_iterations=10000000)
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

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
# 			if (parallel_TF == true)
# 				push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
# 			else
				tmp_results = branchOp(current_nodeIndex, res, num_iterations)
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
						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = tasks[i]
# 					end
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
					res.lq_at_branchBot[spawned_nodeIndex] = log(sum_nodeData_at_bottom)
					res.like_at_branchBot[spawned_nodeIndex] = sum_nodeData_at_bottom

					
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
			nodeOp(current_nodeIndex, res)
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
end # END iterative_downpass_nonparallel!













end # end of module