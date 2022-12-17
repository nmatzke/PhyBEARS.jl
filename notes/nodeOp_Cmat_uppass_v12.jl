# Uppass version of nodeOp_Cmat
nodeOp_Cmat_uppass_v12 = (uppass_probs_just_below_node, Lprobs, Rprobs, p_Ds_v12) -> begin
	p = p_Ds_v12;

	# Go through each Ci (ancestral state index)
	# Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  #Cijk_vals = p.params.Cijk_vals
  Cijk_vals = p.params.Cijk_rates_t # DIFFERENCE FOR V12
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	Carray_pair = p.p_indices.Carray_pair
	
	# Recalculate probabilities of each cladogenesis event:
	anc_ivals = sort(unique(Carray_ivals))
	sumrates_by_i = collect(repeat([0.0], length(anc_ivals)))
	for i in (1:length(anc_ivals))
		ival = anc_ivals[i]
		TF = Carray_ivals .== ival
		sumrates_by_i[i] = sum(Cijk_vals[TF])
		p.params.Cijk_probs[TF] .= Cijk_vals[TF] ./ sumrates_by_i[i]
	end
	
	# See: calc_uppass_probs_v1.R
	# calc_uppass_scenario_probs_new2
	# 
	
	
	# Error check
	txt="""
	input_probs_each_split_scenario = cbind(ancprobs, Lprobs, Rprobs, scenario_condprob)
	# Multiply the probabilities through
	relprob_each_split_scenario = apply(X=input_probs_each_split_scenario, MARGIN=1, FUN=prod)
	sum_relprob_each_split_scenario = sum(relprob_each_split_scenario)
	prob_each_split_scenario = prob_each_split_scenario / sum(prob_each_split_scenario)
	prob_each_split_scenario = relprob_each_split_scenario
	"""
	
	
	# Calculate likelihoods of states just before speciation
  @inbounds for i in 1:n
  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]

		yTFs = Carray_pair[Ci_eq_i] .== 1
		ysums = sum(Cijk_vals[Ci_eq_i][yTFs] .* tmp1[Cj_sub_i][yTFs] .* tmp2[Ck_sub_i][yTFs]) 
		
		tmpDs[i] = sum(Cijk_vals[Ci_eq_i]./Carray_pair[Ci_eq_i] .* tmp1[Cj_sub_i] .* tmp2[Ck_sub_i]) + sum((Carray_pair[Ci_eq_i] .== 2) .* Cijk_vals[Ci_eq_i]./Carray_pair[Ci_eq_i] .* tmp1[Ck_sub_i] .* tmp2[Cj_sub_i])

  end
  return(tmpDs)
end # END nodeOp_Cmat_v12 = (tmpDs; tmp1, tmp2, p_Ds_v12) -> begin

