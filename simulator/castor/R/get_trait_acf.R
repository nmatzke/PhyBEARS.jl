get_trait_acf = function(	tree, 
							tip_states,
							Npairs				= 10000, # integer or numeric, maximum number of tip pairs to consider per tree. Set to Inf to not restrict the number of pairs.
							Nbins				= NULL, # integer, number of phylodistance bins to consider. If Nbins=NULL and phylodistances=NULL, then bins are chosen automatically to strike a good balance between accuracy and resolution.
							min_phylodistance 	= 0,	# non-negative numeric, minimum phylodistance to consider. Only relevant if phylodistances==NULL
							max_phylodistance 	= NULL,	# numeric, optional maximum phylodistance to consider. If NULL or negative, this will automatically set to the maximum phylodistance between any two tips.
							uniform_grid		= FALSE,# logical, specifying whether each phylodistance bin should be of equal size
							phylodistance_grid	= NULL){ # optional numeric vector listing phylodistances (left bin boundaries) in ascending order, within which to evaluate the ACF. Can be used as an alternative to Nbins. The first bin extends from phylodistances[1] to phylodistances[2], while the last bin extends from phylodistances[end] to max_phylodistance (if given) or Inf.
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode
	Nedges = nrow(tree$edge)
	if(!is.numeric(tip_states)) return(list(success=FALSE, error="tip_states must be numeric"))
	Npairs = min(Npairs,(length(tree$tip.label))^2)
	
	# prepare phylodistance grid
	grid_is_uniform = FALSE
	if(is.null(phylodistance_grid)){
		if(is.null(max_phylodistance) || (max_phylodistance<0)) max_phylodistance = max(find_farthest_tip_pair(tree)$distance)
		if(is.null(Nbins)){
			# determine number of phylodistance bins based on data
			Nbins = max(2,Npairs/10)
		}
		if(uniform_grid){
			phylodistance_grid = seq(from=min_phylodistance, to=max_phylodistance, length.out=(Nbins+1))[seq_len(Nbins)]
			grid_is_uniform = TRUE
		}else{
			# choose non-regular phylodistances grid, such that grid density is roughly proportional to the density of pairwise phylodistances
			tip_phylodistances 	= sort(sample_pairwise_tip_distances(tree, Npairs=Npairs))
			tip_phylodistances	= tip_phylodistances[(tip_phylodistances>=min_phylodistance) & (tip_phylodistances<=max_phylodistance)]
			if(length(tip_phylodistances)==0) return(list(success=FALSE, error=sprintf("Unable to determine suitable phylodistance grid: No phylodistances found in target interval")))
			phylodistance_grid 	= sort(unique(get_inhomogeneous_grid_1D_from_samples(Xstart=min_phylodistance, Xend=max_phylodistance, Ngrid=Nbins+1, randomX=tip_phylodistances)$grid[seq_len(Nbins)]))
			Nbins 				= length(phylodistance_grid)
		}
	}else{
		Nbins = length(phylodistance_grid)
		if(any(diff(phylodistance_grid)<=0)) return(list(success=FALSE, error="Provided phylodistance_grid must be strictly increasing"))
		if(is.null(max_phylodistance)){
			max_phylodistance = Inf
		}else{
			max_phylodistance = max(max_phylodistance,max(phylodistance_grid))
		}
		# checkif provided phylodistance grid is uniform (for efficiency purposes)
		if(length(phylodistance_grid)>1){
			grid_steps = diff(phylodistance_grid)
			if(min(grid_steps)>0.9999999*max(grid_steps)) grid_is_uniform = TRUE
		}else{
			grid_is_uniform = TRUE
		}
	}
	
	results = ACF_continuous_trait_CPP(	Ntips 				= Ntips,
										Nnodes 				= Nnodes,
										Nedges 				= Nedges,
										tree_edge			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
										edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
										state_per_tip 		= tip_states,
										max_Npairs			= Npairs,
										phylodistance_grid	= phylodistance_grid,
										max_phylodistance	= max_phylodistance,
										grid_is_uniform		= grid_is_uniform,
										verbose 			= FALSE,
										verbose_prefix 		= "")
										
	# explicitly define left & right bounds, as well as centers, of phylodistance grid cells
	left_phylodistances		= phylodistance_grid
	right_phylodistances	= c((if(Nbins>1) phylodistance_grid[2:Nbins] else NULL), (if(max_phylodistance==Inf) results$max_encountered_phylodistance else max_phylodistance))
	centroid_phylodistances	= 0.5*(left_phylodistances+right_phylodistances)

	return(list(success					= TRUE,
				phylodistances			= centroid_phylodistances,
				left_phylodistances		= left_phylodistances,
				right_phylodistances	= right_phylodistances,
				autocorrelations		= results$autocorrelations,
				mean_abs_differences	= results$mean_abs_deviations,
				mean_rel_differences	= results$mean_rel_deviations,
				Npairs_per_distance		= results$N_per_grid_point))
}