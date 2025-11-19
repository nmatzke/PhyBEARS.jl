geographic_acf = function(	trees, 
							tip_latitudes,
							tip_longitudes,
							Npairs				= 10000, # integer or numeric, maximum number of tip pairs to consider per tree. Set to Inf to not restrict the number of pairs.
							Nbins				= NULL, # integer, number of phylodistance bins to consider. If Nbins=NULL and phylodistances=NULL, then bins are chosen automatically to strike a good balance between accuracy and resolution.
							min_phylodistance 	= 0,	# non-negative numeric, minimum phylodistance to consider. Only relevant if phylodistances==NULL
							max_phylodistance 	= NULL,	# numeric, optional maximum phylodistance to consider. If NULL or negative, this will automatically set to the maximum phylodistance between any two tips.
							uniform_grid		= FALSE,# logical, specifying whether each phylodistance bin should be of equal size
							phylodistance_grid	= NULL){ # optional numeric vector listing phylodistances (left bin boundaries) in ascending order, within which to evaluate the ACF. Can be used as an alternative to Nbins. The first bin extends from phylodistances[1] to phylodistances[2], while the last bin extends from phylodistances[end] to max_phylodistance (if given) or Inf.

	# basic input checking
	if("phylo" %in% class(trees)){
		# trees[] is actually a single tree
		trees = list(trees)
		Ntrees = 1
		if(!(("list" %in% class(tip_latitudes)) && (length(tip_latitudes)==1))){
			tip_latitudes = list(tip_latitudes)
		}
		if(!(("list" %in% class(tip_longitudes)) && (length(tip_longitudes)==1))){
			tip_longitudes = list(tip_longitudes)
		}
	}else if("list" %in% class(trees)){
		# trees[] is a list of trees
		Ntrees = length(trees)
		if("list" %in% class(tip_latitudes)){
			if(length(tip_latitudes)!=Ntrees) return(list(success=FALSE,error=sprintf("Input list of tip_latitudes has length %d, but should be of length %d (number of trees)",length(tip_latitudes),Ntrees)))
		}else if("numeric" %in% class(tip_latitudes)){
			if(Ntrees!=1) return(list(success=FALSE,error=sprintf("Input tip_latitudes was given as a single vector, but expected a list of %d vectors (number of trees)",Ntrees)))
			if(length(tip_latitudes)!=length(trees[[1]]$tip.label)) return(list(success=FALSE,error=sprintf("Input tip_latitudes was given as a single vector of length %d, but expected length %d (number of tips in the input tree)",length(tip_latitudes),length(trees[[1]]$tip.label))))
			tip_latitudes = list(tip_latitudes)
		}
		if("list" %in% class(tip_longitudes)){
			if(length(tip_longitudes)!=Ntrees) return(list(success=FALSE,error=sprintf("Input list of tip_longitudes has length %d, but should be of length %d (number of trees)",length(tip_longitudes),Ntrees)))
		}else if("numeric" %in% class(tip_longitudes)){
			if(Ntrees!=1) return(list(success=FALSE,error=sprintf("Input tip_longitudes was given as a single vector, but expected a list of %d vectors (number of trees)",Ntrees)))
			if(length(tip_longitudes)!=length(trees[[1]]$tip.label)) return(list(success=FALSE,error=sprintf("ERROR: Input tip_longitudes was given as a single vector of length %d, but expected length %d (number of tips in the input tree)",length(tip_longitudes),length(trees[[1]]$tip.label))))
			tip_longitudes = list(tip_longitudes)
		}
	}else{
		return(list(success=FALSE,error=sprintf("Unknown data format '%s' for input trees[]: Expected a list of phylo trees or a single phylo tree",class(trees)[1])))
	}
	Ntotal_tip_pairs = 0 # approximate number of tip pairs that will be available for calculations. Used to estimate the proper resolution of the ACF if needed.
	for(tr in seq_len(Ntrees)){
		if(length(tip_latitudes[[tr]])!=length(trees[[tr]]$tip.label)) return(list(success=FALSE, error=sprintf("tip_latitudes of tree %d has length %d, but expected exactly %d (=Ntips) entries",tr,length(tip_latitudes[[tr]]),length(trees[[tr]]$tip.label))))
		if(length(tip_longitudes[[tr]])!=length(trees[[tr]]$tip.label)) return(list(success=FALSE, error=sprintf("tip_longitudes of tree %d has length %d, but expected exactly %d (=Ntips) entries",tr,length(tip_longitudes[[tr]]),length(trees[[tr]]$tip.label))))
		Ntotal_tip_pairs = Ntotal_tip_pairs + min(Npairs,(length(trees[[tr]]$tip.label))^2)
	}
	
	# prepare phylodistance grid
	grid_is_uniform = FALSE
	if(is.null(phylodistance_grid)){
		if(is.null(max_phylodistance) || (max_phylodistance<0)) max_phylodistance = max(sapply(seq_len(Ntrees), FUN=function(k) find_farthest_tip_pair(trees[[k]])$distance))
		if(is.null(Nbins)){
			# determine number of phylodistance bins based on data
			Nbins = max(2,Ntotal_tip_pairs/10)
		}
		if(uniform_grid){
			phylodistance_grid = seq(from=min_phylodistance, to=max_phylodistance, length.out=(Nbins+1))[1:Nbins]
			grid_is_uniform = TRUE
		}else{
			# choose non-regular phylodistances grid, such that grid density is roughly proportional to the density of pairwise phylodistances
			tip_phylodistances = vector(mode="list", Ntrees)
			for(tr in seq_len(Ntrees)){
				tip_phylodistances[[tr]] = sample_pairwise_tip_distances(trees[[tr]], Npairs=min(Npairs,(length(trees[[tr]]$tip.label))^2/2))
			}
			tip_phylodistances 	= sort(unlist(tip_phylodistances))
			tip_phylodistances	= tip_phylodistances[(tip_phylodistances>=min_phylodistance) & (tip_phylodistances<=max_phylodistance)]
			if(length(tip_phylodistances)==0) return(list(success=FALSE, error=sprintf("Unable to determine suitable phylodistance grid: No phylodistances found in target interval")))
			phylodistance_grid 	= sort(unique(get_inhomogeneous_grid_1D_from_samples(Xstart=min_phylodistance, Xend=max_phylodistance, Ngrid=Nbins+1, randomX=tip_phylodistances)$grid[1:Nbins]))
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

	# calculate ACF for each tree on the same phylodistance tree, and average all ACFs
	mean_autocorrelations 	= numeric(Nbins)
	SS_autocorrelations 	= numeric(Nbins) # sum of squared autocorrelations per phylodistance bin
	mean_geodistances		= numeric(Nbins)
	SS_geodistances			= numeric(Nbins)
	Npairs_per_distance		= integer(Nbins)
	max_encountered_phylodistance = 0
	for(tr in seq_len(Ntrees)){	
		tree = trees[[tr]]
		results = ACF_spherical_CPP(	Ntips 				= length(tree$tip.label),
										Nnodes 				= tree$Nnode,
										Nedges 				= nrow(tree$edge),
										tree_edge			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
										edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
										tip_latitudes 		= tip_latitudes[[tr]],
										tip_longitudes 		= tip_longitudes[[tr]],
										max_Npairs			= Npairs,
										phylodistance_grid	= phylodistance_grid,
										max_phylodistance	= max_phylodistance,
										grid_is_uniform		= grid_is_uniform,
										verbose 			= FALSE,
										verbose_prefix 		= "")
		valid_bins = which(is.finite(results$mean_autocorrelations))
		max_encountered_phylodistance		= max(max_encountered_phylodistance, results$max_encountered_phylodistance)
		mean_autocorrelations[valid_bins] 	= mean_autocorrelations[valid_bins] + results$mean_autocorrelations[valid_bins] * results$N_per_grid_point[valid_bins]
		SS_autocorrelations[valid_bins] 	= SS_autocorrelations[valid_bins] + results$SS_autocorrelations[valid_bins]
		mean_geodistances[valid_bins] 		= mean_geodistances[valid_bins] + results$mean_geodistances[valid_bins] * results$N_per_grid_point[valid_bins]
		SS_geodistances[valid_bins] 		= SS_geodistances[valid_bins] + results$SS_geodistances[valid_bins]
		Npairs_per_distance[valid_bins] 	= Npairs_per_distance[valid_bins] + results$N_per_grid_point[valid_bins]
	}
	mean_autocorrelations 	= mean_autocorrelations/Npairs_per_distance
	std_autocorrelations	= sqrt(SS_autocorrelations/Npairs_per_distance - mean_autocorrelations^2)
	mean_geodistances 		= mean_geodistances/Npairs_per_distance
	std_geodistances		= sqrt(SS_geodistances/Npairs_per_distance - mean_geodistances^2)
	
	# explicitly define left & right bounds, as well as centers, of phylodistance grid cells
	left_phylodistances		= phylodistance_grid
	right_phylodistances	= c((if(Nbins>1) phylodistance_grid[2:Nbins] else NULL), (if(max_phylodistance==Inf) max_encountered_phylodistance else max_phylodistance))
	centroid_phylodistances	= 0.5*(left_phylodistances+right_phylodistances)
	
	return(list(success 			= TRUE,
				phylodistances		= centroid_phylodistances,
				left_phylodistances	= left_phylodistances,
				right_phylodistances= right_phylodistances,
				autocorrelations	= mean_autocorrelations,
				std_autocorrelations= std_autocorrelations,
				mean_geodistances	= mean_geodistances,
				std_geodistances	= std_geodistances,
				Npairs_per_distance	= Npairs_per_distance))
}