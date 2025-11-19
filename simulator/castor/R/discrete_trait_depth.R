# calculate the phylogenetic depth at which a discrete (not necessarily binary) trait is conserved
# This function is similar to the consenTRAIT metric, with the important difference that it treats all states equally (i.e., presences are not more important than absences) and that it works with non-binary traits
# So even for binary traits this function will return different results than the consenTRAIT metric.
# That said, internally this function applies the consenTRAIT metric to each state and then averages the results over all possible states
discrete_trait_depth = function(tree, 
								tip_states, 				# vector of length Ntips, listing discrete tip states. In principle this may list arbitrary data formats, e.g. characters, integers etc
								min_fraction			= 0.9, 
								count_singletons		= TRUE,
								singleton_resolution	= 0,
								weighted				= FALSE, 
								Npermutations			= 0){
	Ntips = length(tree$tip.label)

	# map state space to integers 1,..,Nstates
	unique_states 		= sort(unique(tip_states))
	Nstates 			= length(unique_states)
	integer_tip_states 	= integer(Ntips)
	for(i in seq_len(Nstates)){
		integer_tip_states[tip_states==unique_states[i]] = i;
	}

	results = get_discrete_trait_depth_CPP(	Ntips 				= length(tree$tip.label),
											Nnodes 				= tree$Nnode,
											Nedges 				= nrow(tree$edge),
											Nstates				= Nstates,
											tree_edge			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
											edge_length 		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
											state_per_tip 		= integer_tip_states-1,
											threshold_fraction 	= min_fraction,
											count_singletons 	= count_singletons,
											weighted			= weighted,
											singleton_threshold = singleton_resolution,
											Npermutations 		= Npermutations,
											verbose 			= FALSE,
											verbose_prefix 		= "");
											

	return(list(unique_states			= unique_states,
				mean_depth				= results$tauD,
				var_depth				= results$varD,
				min_depth				= results$minD,
				max_depth				= results$maxD,
				Nmax_uniform			= results$Nmax_uniform,
				mean_depth_per_state	= results$tauD_per_state,
				var_depth_per_state 	= results$varD_per_state, 
				min_depth_per_state		= results$minD_per_state, 
				max_depth_per_state 	= results$maxD_per_state,
				Nmax_uniform_per_state	= results$Nmax_uniform_per_state,
				P						= results$Pvalue, 
				mean_random_depth		= results$mean_random_tauD))
}