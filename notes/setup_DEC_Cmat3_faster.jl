"""
include("/GitHub/PhyBEARS.jl/notes/setup_DEC_Cmat3_faster.jl")
"""


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

Cparams=CparamsStructure(0.9, 0.9, 0.9, 0.3)
Carray1b = setup_DEC_Cmat(areas_list, states_list, NaN, Cparams);
Carray3b = setup_DEC_Cmat3(areas_list, states_list, NaN, Cparams);
trt1b = prtC(Carray1b); sort!(trt1, :i)
trt3b = prtC(Carray3b); sort!(trt3, :i)
trt1b
trt3b
sum(prtC(Carray1b)[:,:val])
sum(prtC(Carray3b)[:,:val])

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
					rsize = length(rstate)
					lstate = ancstate
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
							# Record the range-copying event
							numC += 1
							Carray_event_types[numC] = "j"
							Carray_ivals[numC] = i
							Carray_jvals[numC] = i	# because one daughter is same as anc in jumps
							Carray_kvals[numC] = state_index	# rstate for jump
							Carray_pair[numC] = 2
							Cijk_weights[numC] = tmp_weightval * 2.0 
							row_weightvals[i] += tmp_weightval * 2.0
							continue
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

