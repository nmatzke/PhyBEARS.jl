consentrait_depth = function(	tree, 
								tip_states, 
								min_fraction				= 0.9, 
								count_singletons			= TRUE,
								singleton_resolution		= 0,
								weighted					= FALSE, 
								Npermutations				= 0){
	results = get_trait_depth_consenTRAIT_CPP(	Ntips 				= length(tree$tip.label),
												Nnodes 				= tree$Nnode,
												Nedges 				= nrow(tree$edge),
												tree_edge			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												edge_length 		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
												state_per_tip 		= tip_states,
												threshold_fraction 	= min_fraction,
												count_singletons 	= count_singletons,
												weighted			= weighted,
												singleton_threshold = singleton_resolution,
												Npermutations 		= Npermutations,
												verbose 			= FALSE,
												verbose_prefix 		= "");

	return(list(mean_depth 				= results$tauD,
				var_depth 				= results$varD, 
				min_depth				= results$minD, 
				max_depth 				= results$maxD,
				Npositives				= results$Npositives, 
				P						= results$Pvalue, 
				mean_random_depth		= results$mean_random_tauD,
				positive_clades			= results$positive_clades+1,	# 1D integer vector, listing clades found to be positive in the trait and included in the consenTRAIT statistic
				positives_per_clade		= results$positives_per_clade,	# number of positive tips per clade
				mean_depth_per_clade	= results$mean_depth_per_clade));
}