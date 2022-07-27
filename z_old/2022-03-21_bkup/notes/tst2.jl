#######################################################
# A test module to run Tmp.jl, not in scope Main
# 
# Development workflow recommended by:
#
# https://docs.julialang.org/en/v1/manual/workflow-tips/
#
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst2.jl")

"""

# 
#######################################################


#######################################################
# Run all of the functions here
#######################################################

module Tst2
	cd("/GitHub/BioGeoJulia.jl/notes/")
	include("ModelLikes.jl")
	import .ModelLikes
	#using .Tmp
	
	using Profile     # for @profile
	using DataFrames  # for DataFrame
	using PhyloNetworks
	using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
	using BioGeoJulia.StateSpace
	using BioGeoJulia.TreePass
	using BioGeoJulia.SSEs

	ModelLikes.say_hello2()
	
	#######################################################
	# NEXT:
	# 1. multiply conditional probabilities by lambda
	# 2. load some geographic states, and a phylogeny
	# 3. input into likelihood
	# 4. maximum likelihood (or speed tests)
	#######################################################
	
	ModelLikes.setup_DEC_SSE(2)
#	ModelLikes.setup_DEC_SSE(4)
# 	ModelLikes.setup_DEC_SSE(10)
# 	ModelLikes.setup_DEC_SSE(20)
# 	ModelLikes.setup_DEC_SSE(100)
	
	# Update Qij_vals
	numareas = 7
	areas_list = collect(1:numareas)
	states_list = areas_list_to_states_list(areas_list, numareas, true)
	numstates = length(states_list)
	numstates
	amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
	dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
	elist = repeat([0.123], numstates)
	allowed_event_types=["d","e"]

	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	event_type_vals = Qmat.event_type_vals
	Qmat1_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)

	# Update!
	dmat = reshape(repeat([0.5], numareas^2), (numareas,numareas))
	Qmat2 = ModelLikes.update_Qij_vals(Qmat, areas_list, states_list, dmat, elist, amat )
	Qmat2
	
	Qarray_ivals = Qmat2.Qarray_ivals
	Qarray_jvals = Qmat2.Qarray_jvals
	Qij_vals = Qmat2.Qij_vals
	event_type_vals = Qmat2.event_type_vals
	Qmat2_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, event_type_vals)
	
	Qmat1_df
	Qmat2_df



	tr=readTopology("((chimp:1,human:1):1,gorilla:2);")

	numareas_vec = [2,3,4,5,6,7,8,9,10,11]
#	numareas_vec = [3]
	numstates_vec = repeat([0.0], length(numareas_vec))
	calctimes = repeat([0.0], length(numareas_vec))
	ind=1
	for ind in 1:length(numareas_vec)
		numareas = numareas_vec[ind]
		areas_list = collect(1:numareas)
		
		inputs = ModelLikes.setup_DEC_SSE(numareas, tr);
		
		"""
		tmpdf = DataFrame(Ci=inputs.p_Ds_v5.p_indices.Carray_ivals, Cj=inputs.p_Ds_v5.p_indices.Carray_jvals, Ck=inputs.p_Ds_v5.p_indices.Carray_kvals, vals=inputs.p_Ds_v5.params.Cijk_vals)
		showall(tmpdf)
		"""
		
# 		@profile inputs = ModelLikes.setup_DEC_SSE(numareas, tr);
# 		Profile.print()
		states_list = areas_list_to_states_list(areas_list, length(areas_list), true);
		numstates = length(states_list);

		res = inputs.res
		trdf = inputs.trdf
		solver_options = inputs.solver_options
		solver_options.solver=Tsit5()
		p_Ds_v5 = inputs.p_Ds_v5

		(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res, trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10);

#  		@profile (total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res, trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10);
#  		Profile.print()

		res.likes_at_each_nodeIndex_branchTop
		res.sum_likes_at_nodes
		res.logsum_likes_at_nodes
		log.(res.sum_likes_at_nodes[res.sum_likes_at_nodes.!=0.0])
		sum(log.(res.sum_likes_at_nodes[res.sum_likes_at_nodes.!=0.0]))

		total_calctime_in_sec
		iteration_number


# 		print("\n")
# 		print(res.likes_at_each_nodeIndex_branchTop)
# 		print("\n")
# 		print(res.sum_likes_at_nodes)
# 		print("\n")
# 		print(res.logsum_likes_at_nodes)
# 		print("\n")
		print(log.(res.sum_likes_at_nodes[res.sum_likes_at_nodes.!=0.0]))
		print("\n")
		print(sum(log.(res.sum_likes_at_nodes[res.sum_likes_at_nodes.!=0.0])))
		print("\n")

		print("DEC-SSE total_calctime_in_sec:")
		print(total_calctime_in_sec)
		print("\n")
		numstates_vec[ind] = numstates
		calctimes[ind] = total_calctime_in_sec

	end
	
	print("\n\nAll calculation times:\n")
	calctimes_df = DataFrame(numareas=numareas_vec, numstates=numstates_vec, time=calctimes)
	print(calctimes_df)



df = ModelLikes.workprecision(prob, setups, abstols, reltols, solnames=solnames, save_everystep=false,numruns=1)


end # End of module Tst2









