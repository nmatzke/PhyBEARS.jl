# Fit a spherical Brownian Motion model with constant diffusivity for lineage dispersal, i.e. estimate the isotropic diffusivity D
#	of geographic dispersal assuming that geographic location evolves as a Brownian motion on a sphere (independently in each lineage).
# The diffusivity D is defined as follows: 
#	Under the planar approximation, the vector-valued displacement X \in R^2 between the two lineages is normally distributed with isotropic covariance matrix [2D*t 0; 0 2D*t], where t is the phylogenetic distance (e.g. in time units or ANI)
# 	In other words the squared distance between the two lineages has expectation Ndim * 2 * D * t (where Ndim=2)
# The diffusivity is defined in squared distance units per phylogenetic divergence. 
# The distance units are the same as used for the radius of the sphere. E.g., if you set radius=1 and phylogenetic edge lengths are measured in Myr, then diffusivity will be defined in squared radii per Myr.
# The input tree may include monofurcations and multifurcations, but note that multifurcations are internally first split into bifurcations
fit_sbm_const_biased_sampling = function(	trees, 					# either a single tree in phylo format, or a list of trees
							tip_latitudes, 			# either a 1D vector of size Ntips (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), listing geographical latitudes of the tips (in decimal degrees) of each tree
							tip_longitudes, 		# either a 1D vector of size Ntips (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), listing geographical longitudes of the tips (in decimal degrees) of each tree
							radius,					# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
							planar_approximation=FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
							only_basal_tip_pairs=FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
							min_MRCA_time=0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
							min_diffusivity=NULL,	# numeric, specifying the lower bound of allowed diffusivities. If omitted, it will be automatically chosen. Only relevant if planar_approximation==FALSE.
							max_diffusivity=NULL,	# numeric, specifying the upper bound of allowed diffusivities. If omitted, it will be automatically chosen. Only relevant if planar_approximation==FALSE.
							sampling_rate=NULL,		# function returning the sampling rate at any given geographic location (latitude, longitude). Should take 2 arguments, latitude (in decimal degrees, from -90 to 90) and longitude (in decimal degrees, from -180 to 180), and should return a probability within [0,1].
							Nbootstraps=0,			# integer, number of boostraps to perform. If <=0, no boostrapping is performed.	
							focal_diffusivities=NULL){	# optional integer vector, listing diffusivities D for which to calculate and return the loglikelihoods, e.g. for diagnostic purposes. This does not influence the actual fitting.
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
	if(is.null(min_diffusivity) || is.na(min_diffusivity)) min_diffusivity = NaN;
	if(is.null(max_diffusivity) || is.na(max_diffusivity)) max_diffusivity = NaN;
	if(is.null(Nbootstraps) || is.na(Nbootstraps) || (Nbootstraps<0)) Nbootstraps = 0;
	if(min_MRCA_time<0) min_MRCA_time = 0;
	with_sampling_rate = (!is.null(sampling_rate))
	if(planar_approximation && with_sampling_rate) return(list(success=FALSE, error="Accounting for sampling_rate is not available for the planar appoximation"))

	# loop through the trees and extract independet contrasts between sister clades
	phylogenetic_distances 	= numeric()
	latitudes1				= numeric()
	longitudes1				= numeric()
	latitudes2				= numeric()
	longitudes2				= numeric()
	for(i in 1:Ntrees){
		tree = trees[[i]]
		tip_latitudes_this_tree  = tip_latitudes[[i]]
		tip_longitudes_this_tree = tip_longitudes[[i]]
		
		# make sure tree does not have multifurcations
		tree = multifurcations_to_bifurcations(tree)$tree
		
		# extract independent pairs of sister tips
		tip_pairs = extract_independent_sister_tips(tree)
		if(only_basal_tip_pairs){
			# calculate number of nodes between tip pairs
			edge_counts = get_pairwise_distances(tree, A=tip_pairs[,1], B=tip_pairs[,2], as_edge_counts=TRUE, check_input=FALSE)
			# only keep tip pairs with at most 2 edges connecting them
			keep_pairs 	= which(edge_counts<=2)
			tip_pairs 	= tip_pairs[keep_pairs,,drop=FALSE]
		}
		if(min_MRCA_time>0){
			MRCAs 			= get_pairwise_mrcas(tree, tip_pairs[,1], tip_pairs[,2], check_input=FALSE)
			clade_heights	= castor::get_all_distances_to_root(tree)
			keep_pairs		= which(clade_heights[MRCAs]>=min_MRCA_time)
			tip_pairs 		= tip_pairs[keep_pairs,,drop=FALSE]
		}
		if(nrow(tip_pairs)==0) next; # no valid tip pairs found in this tree
		
		# calculate phylogenetic divergence between sister tips
		phylogenetic_distances = c(phylogenetic_distances, get_pairwise_distances(tree, A=tip_pairs[,1], B=tip_pairs[,2], check_input=FALSE))
		
		latitudes1 	= c(latitudes1, tip_latitudes_this_tree[tip_pairs[,1]])
		longitudes1	= c(longitudes1, tip_longitudes_this_tree[tip_pairs[,1]])
		latitudes2 	= c(latitudes2, tip_latitudes_this_tree[tip_pairs[,2]])
		longitudes2	= c(longitudes2, tip_longitudes_this_tree[tip_pairs[,2]])
	}

	# omit tip pairs with zero phylogenetic distance, because in that case the likelihood density is pathological
	valid_pairs 			= which(phylogenetic_distances>0)
	phylogenetic_distances 	= phylogenetic_distances[valid_pairs]
	latitudes1	 			= latitudes1[valid_pairs]
	longitudes1	 			= longitudes1[valid_pairs]
	latitudes2	 			= latitudes2[valid_pairs]
	longitudes2	 			= longitudes2[valid_pairs]
	NC 						= length(phylogenetic_distances)
	if(NC==0) return(list(success=FALSE, error="No valid tip pairs left for extracting independent contrasts"))
	
	# calculate geodesic distance for each tip pair
	geodesic_distances = radius * sapply(1:NC, FUN=function(p) geodesic_angle(latitudes1[p],longitudes1[p],latitudes2[p],longitudes2[p]))
		
	# prepare sampling_rate grid
	if(with_sampling_rate){
		NSRlat = 300; NSRlon = 300;
		sampling_lat_grid  = seq(-90,90,length.out=NSRlat)
		sampling_lon_grid  = seq(-180,180,length.out=NSRlon)
		sampling_rate_grid = matrix(1, nrow=NSRlat, ncol=NSRlon)
		for(r in 1:NSRlat){
			sampling_rate_grid[r,] = sapply(1:NSRlon, FUN=function(k) sampling_rate(sampling_lat_grid[r], sampling_lon_grid[k]))
		}
	}

	bootstrap_diffusivities = NULL
	if(planar_approximation){
		# ignore spherical structure and simply translate geodesic_angles to geodesic distances, then fit Brownian motion on a plane with isotropic diffusivity matrix
		Ndim 			= 2 # dimensionality of the Brownian motion
		diffusivity		= 0.5 * (1/(Ndim*NC)) * sum(geodesic_distances^2 / phylogenetic_distances);
		loglikelihood 	= -0.5*NC*Ndim*log(2*pi) - 0.5*Ndim*NC*log(2*diffusivity) - 0.5*Ndim*sum(log(phylogenetic_distances)) - (1/(4*diffusivity)) * sum(geodesic_distances^2 / phylogenetic_distances);
		if(is.nan(diffusivity)) return(list(success=FALSE, error="Fitted diffusivity is NaN"))
		if(is.nan(loglikelihood)) return(list(success=FALSE, error="Loglikelihood of fitted diffusivity is NaN"))
		if(Nbootstraps>0){
			bootstrap_diffusivities  = vector(mode="numeric", Nbootstraps)
			bootstrap_loglikelihoods = vector(mode="numeric", Nbootstraps)
			for(b in 1:Nbootstraps){
				bootstrap_distances 		= radius * sapply(1:NC, FUN=function(p) draw_SBM_geodesic_angle_CPP(phylogenetic_distances[p]*diffusivity/(radius^2)))
				bootstrap_diffusivities[b]	= 0.5 * (1/(Ndim*NC)) * sum(bootstrap_distances^2 / phylogenetic_distances);
				bootstrap_loglikelihoods[b]	= -0.5*NC*Ndim*log(2*pi) - 0.5*Ndim*NC*log(2*bootstrap_diffusivities[b]) - 0.5*Ndim*sum(log(phylogenetic_distances)) - (1/(4*bootstrap_diffusivities[b])) * sum(bootstrap_distances^2 / phylogenetic_distances);
			}
		}
	}else{
		# maximum-likelihood estimation based on full spherical Brownian motion model
		if(with_sampling_rate){
			# account for location-dependent sampling rate
			fit = fit_SBM_from_sampled_transitions_CPP(	radius 				= radius,
														time_steps			= phylogenetic_distances,
														old_thetas			= pi*latitudes1/180,
														old_phis			= pi*longitudes1/180,
														new_thetas			= pi*latitudes2/180,
														new_phis			= pi*longitudes2/180,
														Nlat				= NSRlat,
														Nlon				= NSRlon,
														sampling_rates		= as.vector(t(sampling_rate_grid)),
														max_error			= 1e-7,
														max_Legendre_terms	= 200,
														opt_epsilon			= 1e-10,
														max_iterations 		= 10000,
														min_diffusivity		= min_diffusivity,
														max_diffusivity		= max_diffusivity,
														Nbootstraps			= Nbootstraps);
		}else{
			# assume equal sampling rate everywhere
			fit = fit_SBM_diffusivity_from_transitions_CPP(	radius 					= radius,
															time_steps				= phylogenetic_distances,
															distances 				= geodesic_distances,
															max_error				= 1e-7,
															max_Legendre_terms		= 200,
															opt_epsilon				= 1e-10,
															max_iterations 			= 10000,
															min_diffusivity			= min_diffusivity,
															max_diffusivity			= max_diffusivity,
															Nbootstraps				= Nbootstraps,
															SBM_PD_functor			= list());
		}
		if(!fit$success) return(list(success=FALSE, error=fit$error))
		diffusivity 				= fit$fit_diffusivity
		loglikelihood				= fit$fit_loglikelihood
		bootstrap_diffusivities		= fit$bootstrap_diffusivities
		bootstrap_loglikleihoods	= fit$bootstrap_loglikelihoods
	}
	
	if(!is.null(focal_diffusivities)){
		if(with_sampling_rate){
			# account for location-dependent sampling rate
			focal_LLs = SBM_LLs_of_sampled_transitions_CPP(	radius 				= radius,
															time_steps			= phylogenetic_distances,
															old_thetas			= pi*latitudes1/180,
															old_phis			= pi*longitudes1/180,
															new_thetas			= pi*latitudes2/180,
															new_phis			= pi*longitudes2/180,
															diffusivities		= focal_diffusivities,
															Nlat				= NSRlat,
															Nlon				= NSRlon,
															sampling_rates		= as.vector(t(sampling_rate_grid)),
															max_error			= 1e-7,
															max_Legendre_terms	= 200)$loglikelihoods
		}else{
			# assume equal sampling rate everywhere
			focal_LLs = SBM_LLs_of_transitions_CPP(	radius 				= radius,
													time_steps			= phylogenetic_distances,
													distances 			= geodesic_distances,
													diffusivities		= focal_diffusivities,
													max_error			= 1e-7,
													max_Legendre_terms	= 200)$loglikelihoods
		}
	}else{
		focal_LLs = numeric();
	}

	return(list(success					= TRUE,
				diffusivity 			= diffusivity,
				loglikelihood			= loglikelihood,
				bootstrap_std			= (if(Nbootstraps>0) sd(bootstrap_diffusivities, na.rm=TRUE) else NULL),
				bootstrap_min			= (if(Nbootstraps>0) min(bootstrap_diffusivities, na.rm=TRUE) else NULL),
				bootstrap_max			= (if(Nbootstraps>0) max(bootstrap_diffusivities, na.rm=TRUE) else NULL),
				bootstrap_CI95lower		= (if(Nbootstraps>0) unname(quantile(bootstrap_diffusivities, probs=0.025, type=8, na.rm=TRUE)) else NULL),
				bootstrap_CI95upper		= (if(Nbootstraps>0) unname(quantile(bootstrap_diffusivities, probs=0.975, type=8, na.rm=TRUE)) else NULL),
				Ncontrasts				= NC,
				phylogenetic_distances 	= phylogenetic_distances,
				geodesic_distances		= geodesic_distances,
				focal_loglikelihoods	= focal_LLs))
}




