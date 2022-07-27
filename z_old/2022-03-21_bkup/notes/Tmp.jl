module Tmp
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
using DataFrames          # for DataFrame()

# (1) List all function names here:
export say_hello, CparamsStructure, default_Cparams, sumy, sums, sumv, sumj, setup_DEC_Cmat, relative_probabilities_of_subsets, relative_probabilities_of_vicariants, discrete_maxent_distrib_of_smaller_daughter_ranges, array_in_array, setup_DEC_Cmat, update_Cijk_vals


#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst.jl")
"""
#######################################################


#######################################################
# (2) write the functions here
#######################################################

say_hello() = println("Hello dude!")



#######################################################
# Example use of maximum entropy on discrete case
# https://github.com/JuliaOpt/Convex.jl/issues/64
#######################################################
using Distributions  # for quantile
using Convex				 # for Convex.entropy(), maximize()
using SCS						 # for SCSSolve, solve (maximize(entropy()))





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
######################################################
Weight the possible range-sizes of the smaller daughter range

E.g., if you have subset sympatry from ancestor ABC,
the smaller daugher could be size 1, 2, or 3.

DEC says the weights of these would be 1 0 0
BayArea says the weights would be 0 0 1

With this function, input 
maxent_constraint_01=0.0001 = 1 0 0 
maxent_constraint_01=0.5    = 1 1 1 
maxent_constraint_01=0.9999 = 0 0 1 
######################################################

max_numareas=6
maxent_constraint_01 = 0.0001    # ranges from 0.0001 (all weight on ranges of size 1)
						                     # to 0.9999 (all weight on ranges of max size)
NA_val=NaN
relative_probabilities_of_subsets(max_numareas, maxent_constraint_01, NA_val)
"""

function relative_probabilities_of_subsets(max_numareas=6, maxent_constraint_01=0.5, NA_val=NaN)
	# Set up a matrix to hold the maxent distributions of relative prob of 
	# smaller daughter range sizes
	relprob_subsets_matrix = reshape(repeat([NA_val], max_numareas^2), (max_numareas,max_numareas))
	for i in 1:max_numareas
		ancestor_range_size = i
		maxent_result = discrete_maxent_distrib_of_smaller_daughter_ranges(ancestor_range_size, maxent_constraint_01)

		if (i == 1)
			relprob_subsets_matrix[i,1:i] .= maxent_result
		else
			relprob_subsets_matrix[i,1:i] = maxent_result
		end
	end
	
	return relprob_subsets_matrix
end



function relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01=0.5, NA_val=NaN)
	# Set up a matrix to hold the maxent distributions of relative prob of 
	# smaller daughter range sizes
	relprob_subsets_matrix = reshape(repeat([NA_val], max_numareas^2), (max_numareas,max_numareas))
	relprob_subsets_matrix[1,:] .= NA_val
	
	for i in 2:max_numareas
		ancestor_range_size = i
		tmpstates = collect(1:i)# .+ 0.0
		tmpstates_Floats = collect(1:i) .+ 0.0
		max_smaller_rangesize = median(tmpstates_Floats)
		possible_vicariance_smaller_rangesizes = tmpstates[tmpstates_Floats .< max_smaller_rangesize]
		vic_max_numareas = maximum(possible_vicariance_smaller_rangesizes)

		maxent_result = flat2(discrete_maxent_distrib_of_smaller_daughter_ranges(vic_max_numareas, maxent_constraint_01))

		if (i <= 3)
# 			print("\n")
# 			print(i)
# 
# 			print("\n")
# 			print(vic_max_numareas)
# 			
# 			print("\n")
# 			print(relprob_subsets_matrix[i,1:vic_max_numareas])
# 
# 			print("\n")
# 			print(maxent_result)
# 
# 			print("\n")
			relprob_subsets_matrix[i,1:vic_max_numareas] .= maxent_result
		else
			relprob_subsets_matrix[i,1:vic_max_numareas] = maxent_result
		end
	end
	
	return relprob_subsets_matrix
end




"""
numareas = 6
max_numareas=6
maxent_constraint_01=0.0
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

max_numareas=6
maxent_constraint_01=0.5
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

max_numareas=6
maxent_constraint_01=1.0
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

numareas = 2
max_numareas=6
maxent_constraint_01=0.0
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

max_numareas=6
maxent_constraint_01=0.5
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

max_numareas=6
maxent_constraint_01=1.0
discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

"""
function discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas=6, maxent_constraint_01=0.5)
	n = max_numareas
	x = 1:n

	#discrete_values_padded = cat(0, collect(x)[:], n+1; dims=1)
	discrete_values_padded = collect(x)
	maxent_constraint = quantile(discrete_values_padded, maxent_constraint_01)
	#maxent_constraint = 1

	# Vector of probabilities that must sum to 1
	p = Variable(length(x))
	probability_contraints = [0 <= p, p <= 1, sum(p) == 1];
	feature_constraints = sum(p'*x) == maxent_constraint

	# base or prior (e.g. Uniform)
	#h = pdf.(Uniform(1, n), x)
	
	# This solution updates p (look in p.values)
	problem = maximize(Convex.entropy(p), 0 <= p, p <= 1, sum(p) == 1, feature_constraints)
	sol = Convex.solve!(problem, SCS.SCSSolver(verbose=0))
	maxent_result = abs.(round.(p.value; digits=3))
	return maxent_result
end


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
max_numareas = length(areas_list)
maxent_constraint_01 = 0.0
maxent01symp = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
maxent01sub = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
maxent01jump = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
maxent_constraint_01 = 0.5
maxent01vic = relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01)
maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
predeclare_array_length=10000000
Carray = setup_DEC_Cmat(areas_list, states_list, Cparams)

"""

function setup_DEC_Cmat(areas_list, states_list, maxent01, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))); predeclare_array_length=Integer(min(length(states_list)*length(states_list)*round((length(states_list)/2)), 10000000)), min_precision=1e-9)
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
					if (ancstate == lstate == rstate)
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
						end # end if tmp_weightval > 0.0
					end # end if (ancstate == lstate == rstate)
				end # end if (y_wt > min_precision)
				
				# If one of the descendants is identical to the ancestor, 
				# (after we've excluded sympatry)
				# we can have jump dispersal or subset speciation
				if ( (ancstate == rstate) )
					# Subset sympatry
					if (s_wt > min_precision)
						# Check for subset sympatry: lstate smaller than rstate, lstate inside rstate
						if ((array_in_array(lstate, rstate) == true) && (lsize < rsize))
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
					
							end # if (tmp_weightval > 0.0)
							# end of jump dispersal
						end # end if ( (lsize == 1) && (array_in_array(lstate, rstate) == false) )
					end # end if (j_wt > min_precision)
				end # end if ( (ancstate == rstate) )
			
				# Vicariance
				if (v_wt > min_precision)
					# Check if the combined vector equals the ancestor vector
					combined_vector = sort(cat(lstate,rstate; dims=1))
					if ( combined_vector == sort(ancstate) )
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



#######################################################
# Update the Cijk_vals
#######################################################
function update_Cijk_vals(Carray, areas_list, states_list, maxent01, Cparams=default_Cparams(), dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list))) )

	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Carray_event_types = Carray.Carray_event_types;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals

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
	end # End update of s event weights
	
	# Update the "j" events (jump dispersal)
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
	end # End update of s event weights

	df1 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	
	row_weightvals_df = by(df1, :i, :weight => sum)
	row_weightvals = row_weightvals_df[!,2]
	
	# If i=1 is missing from row_weightvals_df, add it to row_weightvals
	if (in(1, row_weightvals_df[!,1]) == false)
		row_weightvals = cat([1], row_weightvals_df[!,2]; dims=1)
	end
	row_weightvals
	
	# Convert the weights to conditional event probabilities
	num_clado_events = length(Cijk_weights)
	Cijk_vals = collect(repeat([0.0], num_clado_events))
	for i in 1:length(states_list)
		TF = Carray_ivals .== i
		Cijk_vals[TF] = Cijk_weights[TF] ./ row_weightvals[i]
	end

	df2 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	row_weightvals_df = by(df2, :i, :weight => sum)
	
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






end # end module