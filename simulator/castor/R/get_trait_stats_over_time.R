# calculate statistics (mean & standard deviation) of a numeric trait on a tree over time
# if tree$edge.length is missing, edges are assumed to have length 1.
get_trait_stats_over_time = function(	tree, 
										states,							# 1D numeric vector of size Nclades, listing the trait's value for each tip & node
										Ntimes				= NULL, 	# number of equidistant time points for which to calculate trait stats
										times				= NULL, 	# 1D array of time points in increasing order, for which to calculate trait stats
										include_quantiles 	= TRUE,		# include quantile information (e.g., median, CI95 and CI50) over time, i.e. in addition to means & stds. This increases computation time and memory.
										check_input			= TRUE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	
	# basic error checking
	if((!is.null(Ntimes)) && (!is.null(times))) stop("ERROR: Either Ntimes or times must be non-NULL, but not both")
	if(is.null(Ntimes) && is.null(times)) stop("ERROR: Both Ntimes and times are NULL; please specify one of the two")
	if(length(states)!=(Ntips+Nnodes)) stop("ERROR: Number of provided states does not much number of tips+nodes in the tree")
	if(check_input){
		if((!is.null(names(states))) && any(names(states[1:Ntips])!=tree$tip.label)) stop("ERROR: Names in states and tip labels in tree don't match (must be in the same order).")
		if((!is.null(names(states))) && (!is.null(tree$node.label)) && any(names(states[(Ntips+1):(Ntips+Nnodes)])!=tree$node.label)) stop("ERROR: Names in states and node labels in tree don't match (must be in the same order).")
	}
	
	if(is.null(times)){
		tree_span = castor::get_tree_span(tree)
		times = seq(from=0, to=(1-1e-8)*tree_span$max_distance, length.out=Ntimes)
	}else{
		Ntimes = length(times)
	}
	trait_stats = get_trait_stats_at_times_CPP(	Ntips 			= Ntips,
												Nnodes 			= Nnodes,
												Nedges 			= nrow(tree$edge),
												tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
												times			= times,
												states			= states,
												return_states	= include_quantiles)
	if(include_quantiles){
		trait_stats$states_matrix = matrix(trait_stats$states_matrix,nrow=Ntimes,byrow=TRUE)
		CI50lower 	= numeric(Ntimes)
		CI50upper 	= numeric(Ntimes)
		CI95lower 	= numeric(Ntimes)
		CI95upper 	= numeric(Ntimes)
		medians 	= numeric(Ntimes)
		for(i in seq_len(Ntimes)){
			quantiles = quantile(trait_stats$states_matrix[i,], probs=c(0.25, 0.75, 0.025, 0.975, 0.5), na.rm=TRUE)
			CI50lower[i] = quantiles[1]
			CI50upper[i] = quantiles[2]
			CI95lower[i] = quantiles[3]
			CI95upper[i] = quantiles[4]
			medians[i] 	 = quantiles[5]	
		}
	}

	
	return(list(Ntimes			= length(times),
				times			= times,
				clade_counts	= trait_stats$clade_counts,
				means 			= trait_stats$means,	# arithmetic means of trait
				stds 			= trait_stats$stds,		# standard deviations of trait
				medians			= (if(include_quantiles) medians else NULL),
				CI50lower		= (if(include_quantiles) CI50lower else NULL),
				CI50upper		= (if(include_quantiles) CI50upper else NULL),
				CI95lower		= (if(include_quantiles) CI95lower else NULL),
				CI95upper		= (if(include_quantiles) CI95upper else NULL)))
}