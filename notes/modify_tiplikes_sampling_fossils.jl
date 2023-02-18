#######################################################
# Modify tiplikes by:
#######################################################
# psi, sampling rate for serially-sampled fossils
# rho (sampling_f)
# rho (tipsamp_f)

current_nodeIndex

res.sampling_f		# rho sampling probabilities for each state
res.tipsamp_f    	# sampling probabilities for each tip (a priori) -- the problem with this
									# is that it modifies Ei at the tip, for all states. Different tip branches
									# thus can't share the same Es calculation. Better to use states with different
									# rho to control for this (e.g. just a 'poorly sampled' and 'well-sampled' 
									# e.g. forests vs. not). But, leaving it in for experimental purposes.
inputs.setup.fossil_TF # is each node a fossil? If so, D multiplied by psi
inputs.setup.direct_TF # is each node a direct ancestor? 
												# NOTE: fossil direct ancestors should be represented by hooks, NOT
												#       direct ancestors, which are for recording/hypothesizing ancestral
												#       states. (Hooks are ultrashort branches that are not treated
												#       as speciation events.)

inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0         # zero out
inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0     # zero out
inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
inputs.res.sumLikes_at_node_at_branchTop

# Time-constant version of psi
function modify_tiplikes_sampling_fossils_v7(inputs)
	# Error checks
	taxa = inputs.trdf.taxa
	tipnames = sort(taxa[inputs.trdf.nodeType .== "tip"])
	check_tr_geog_tip_labels(tipnames, geog_df)

	dfnames = names(geog_df)
	area_column_nums = 2:length(dfnames)
	areas_txt_list = dfnames[area_column_nums]
	numareas = length(inputs.setup.areas_list)

	# Check if the number of areas in the geography file matches the number in the geog_df
	if (inputs.setup.numareas != numareas)
		txt = paste0(["STOP ERROR in tipranges_to_tiplikes(): inputs.setup.numareas=", numareas, ", but the number of areas in geog_df is ", numareas, ". Please fix and re-run."])
		error(txt)
	end

	# 	statenums = collect(1:length(states_list))
	# 	observed_statenums = collect(repeat([0], nrow(geog_df)))
	trdf_nodenums = collect(1:nrow(inputs.trdf))

	# Go through the geography file tips, match each one to the correct node of trdf,
	# then update the tip likelihoods at that node.

	for i in 1:nrow(geog_df)
		spname = geog_df[i,:tipnames]
		TF = spname .== inputs.trdf[!,:nodeName]
		nodeNum = trdf_nodenums[TF][1]
		
		# Check if it's a fossil
		if inputs.setup.fossil_TF[nodeNum] .== true
			# An m-type fossil (not a hook)
			
			
			# A k-type fossil (a hook)
			
			
			# Multiply Ds by psi_i for fossil tips
			inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .* inputs.p_Ds_v5.params.psi_vals
		end
		
		# Input likelihoods of 1 for the observed state, 0 otherwise
		inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0         # zero out
		inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0     # zero out
		inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
		inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
		inputs.res.sumLikes_at_node_at_branchTop[nodeNum] = 1.0
	end

	
end