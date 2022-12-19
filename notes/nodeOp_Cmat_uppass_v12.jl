# Uppass version of nodeOp_Cmat

# 1. For a given node, get the uppass probs at the bottom of the branch
#    (unless it's the root)
# 2. Run the Ds up from that
# 3. Run nodeOp_Cmat_get_condprobs to get uppass probs for either 
#    L or R descendant

function nodeOp_Cmat_uppass_v12!(res, current_nodeIndex, trdf, p_Ds_v12, solver_options)
	n = numstates = length(res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex])
	# Is this a root node?
	if (current_nodeIndex == res.root_nodeIndex)
		uppass_probs_just_below_node = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		res.anc_estimates_at_each_nodeIndex_branchTop .= normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0

	# Tip or internal nodes require passing probs up from branch bottom
	else
		# The uppass ancestral state probs will have been previously 
		# calculated at the branch bottom
		# BGB's "relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS"
		u0 = probs_at_branch_bottom = res.uppass_probs_at_each_nodeIndex_branchBot[current_nodeIndex] .+ 0.0
		time_start = trdf.node_age[trdf.ancNodeIndex[current_nodeIndex]]
		time_end = trdf.node_age[current_nodeIndex]
		tspan = [time_start, time_end]
		
		# Seems to work 
		(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)= branchOp_ClaSSE_Ds_v12(current_nodeIndex, res; u0, tspan, p_Ds_v12, solver_options=solver_options);
		# These are really conditional probabilities upwards, they don't 
		# have to add up to 1.0, unless normalized
		uppass_probs_just_below_node = sol_Ds.u[length(sol_Ds.u)]
		uppass_probs_just_below_node .= uppass_probs_just_below_node ./ sum(uppass_probs_just_below_node)
	end # END if (current_nodeIndex == res.root_nodeIndex)
			# END uppass from branch below
	
	# Combine info through node; Store results
	# Check if its a tip node
	if (trdf.nodeType[current_nodeIndex] == "tip")
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .* res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
		res.anc_estimates_at_each_nodeIndex_branchTop = res.anc_estimates_at_each_nodeIndex_branchTop ./ sum(res.anc_estimates_at_each_nodeIndex_branchTop)
		res.anc_estimates_at_each_nodeIndex_branchTop .= res.anc_estimates_at_each_nodeIndex_branchTop .+ 0.0
	
	# Internal nodes (not root)
	elseif (trdf.nodeType[current_nodeIndex] == "intern")
		# (Ignore direct ancestors for now)
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex][1] .= uppass_probs_just_below_node .+ 0.0
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .* res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] = res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] ./ sum(res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex])
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
		
		# Get uppass probs for Left and Right daughter branches
		node_above_Left_corner = trdf.leftNodeIndex[current_nodeIndex]
		node_above_Right_corner = trdf.rightNodeIndex[current_nodeIndex]
		# LEFT
		
		# Calculate the post-cladogenesis uppass probabilities for the Left node
		Ldownpass_likes = collect(repeat([1.0], n))
		Rdownpass_likes = res.normlikes_at_each_nodeIndex_branchTop[node_above_Right_corner]
		relprob_each_split_scenario = nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v12; use_Cijk_rates_t=true)

		for j in 1:n
			jTF = p_Ds_v12.p_indices.Carray_jvals .== j
			res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner][j] = sum(relprob_each_split_scenario[jTF])
		end
		res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner] .= res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner] ./ sum(res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner])
		
		# RIGHT
		# Calculate the post-cladogenesis uppass probabilities for the Left node
		Rdownpass_likes = collect(repeat([1.0], n))
		Ldownpass_likes = res.normlikes_at_each_nodeIndex_branchTop[node_above_Left_corner]
		relprob_each_split_scenario = nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v12; use_Cijk_rates_t=true)
		
		for k in 1:n
			kTF = p_Ds_v12.p_indices.Carray_jvals .== k
			res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner][k] = sum(relprob_each_split_scenario[kTF])
		end
		res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner] .= res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner] ./ sum(res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner])
		
		# Get ancestral range estimates for Left and Right daughter branches
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] .= Ldownpass_likes .* res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner]
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] .= res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] ./ sum(res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner])

		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] .= Rdownpass_likes .* res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner]
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] .= res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] ./ sum(res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner])

	end # END elseif internal nodes
end # END nodeOp_Cmat_uppass_v12!(res, current_nodeIndex, trdf, p_Ds_v12, solver_options)

# Get the conditional probabilities of all cladogenesis
# scenarios, conditional on ancestral range probs, Lprobs, Rprobs
#
# On an uppass, either Lprobs or Rprobs will be a vector of 1s
function nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v12; use_Cijk_rates_t=true)
	p = p_Ds_v12;
	
	# Figure out if you are solving for the left or right descendant
	# (if both are all 1.0s, assumes left)
	# (also assumes left if both are non-1st; but throws warning)
	left_or_right = ""
	if (all(Rdownpass_likes .== 1.0) == true)
		left_or_right = "right"
	elseif (all(Ldownpass_likes .== 1.0) == true)
		left_or_right = "left"
	else
		left_or_right = "left"
		txt = "WARNING in nodeOp_Cmat_get_condprobs(): Either Lprobs or Rprobs should be a vector of all 1.0. But this is not the case. Double-check your inputs."
		@warn txt
	end	
	
	
	# Go through each Ci (ancestral state index)
	# Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  if (use_Cijk_rates_t == true)
	  Cijk_vals = p.params.Cijk_rates_t # DIFFERENCE FOR V12
	else
		Cijk_vals = p.params.Cijk_vals
	end

	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	Carray_pair = p.p_indices.Carray_pair
	Carray_pair_floats = Float64.(Carray_pair)
	Carray_pair_ones = Carray_pair_floats .- 1.0
	
	# Recalculate probabilities of each cladogenesis event:
	anc_ivals = sort(unique(Carray_ivals))
	sumrates_by_i = collect(repeat([0.0], length(anc_ivals)))
	for i in (1:length(anc_ivals))
		ival = anc_ivals[i]
		iTF = Carray_ivals .== ival
		sumrates_by_i[i] = sum(Cijk_vals[iTF])
		p.params.Cijk_probs[iTF] .= Cijk_vals[iTF] ./ sumrates_by_i[i]
	end
	
	
	
	
	# See: calc_uppass_probs_v1.R
	# calc_uppass_scenario_probs_new2
	calc_uppass_scenario_probs_new2_CODE_IS_EG="""
	Lprobs = left_branch_downpass_likes[ COO_weights_columnar[[2]] + COOmat_0based_to_Qmat_1based]
	Rprobs = right_branch_downpass_likes[ COO_weights_columnar[[3]] + COOmat_0based_to_Qmat_1based]
	# Weights divided by the sum of the weights for that row
	# COO_weights_columnar[[1]] is the ancestral node, 0-based, cladogenesis-only
	# So, you will always add 1 to the index to get the rowsum conditional
	# on a particular ancestor (3 non-null states means 3 Rsp_rowsums)
	scenario_condprob = COO_weights_columnar[[4]] / Rsp_rowsums[ COO_weights_columnar[[1]] + 1 ]

	input_probs_each_split_scenario = cbind(ancprobs, Lprobs, Rprobs, scenario_condprob)
	# Multiply the probabilities through
	relprob_each_split_scenario = apply(X=input_probs_each_split_scenario, MARGIN=1, FUN=prod)
	sum_relprob_each_split_scenario = sum(relprob_each_split_scenario)
	prob_each_split_scenario = prob_each_split_scenario / sum(prob_each_split_scenario)
	prob_each_split_scenario = relprob_each_split_scenario
	"""
	
	# Calculate the relative probabilities, accounting for pairs
	ancprobs_by_scenario = uppass_probs_just_below_node[Carray_ivals]
	# Divide by 2.0 for the paired scenarios
	Lprobs_by_scenario_pt1 = Ldownpass_likes[Carray_jvals] ./ Carray_pair_floats
	# Multiply by 0.0 for the sympatry events already counted
	Lprobs_by_scenario_pt2 = Ldownpass_likes[Carray_kvals] .* Carray_pair_ones
	Lprobs = Lprobs_by_scenario_pt1 .+ Lprobs_by_scenario_pt2
	
	# Divide by 2.0 for the paired scenarios
	Rprobs_by_scenario_pt1 = Rdownpass_likes[Carray_kvals] ./ Carray_pair_floats
	# Multiply by 0.0 for the sympatry events already counted
	Rprobs_by_scenario_pt2 = Rdownpass_likes[Carray_jvals] .* Carray_pair_ones
	Rprobs = Rprobs_by_scenario_pt1 .+ Rprobs_by_scenario_pt2
	
	descprobs_by_scenario = p.params.Cijk_probs .* ancprobs_by_scenario .* Lprobs .* Rprobs
			
	# Because you are moving:
	# - FROM a vector of ancestral ranges
	# - TO a vector of cladogenesis scenarios
	# ...it *IS* legitimate to sum the conditional probabilities of
	#    all scenarios, then divide by the sum
	relprob_each_split_scenario = descprobs_by_scenario
	relprob_each_split_scenario .= descprobs_by_scenario ./ sum(descprobs_by_scenario)
	
  return(relprob_each_split_scenario)
end # END nodeOp_Cmat_get_condprobs = (tmpDs; tmp1, tmp2, p_Ds_v12) -> begin

