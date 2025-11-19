# Fit a spherical Brownian Motion model with constant diffusivity for lineage dispersal, i.e. estimate the isotropic diffusivity D
#	of geographic dispersal assuming that geographic location evolves as a Brownian motion on a sphere (independently in each lineage).
# The diffusivity D is defined as follows: 
#	Under the planar approximation, the vector-valued displacement X \in R^2 between the two lineages is normally distributed with isotropic covariance matrix [2D*t 0; 0 2D*t], where t is the phylogenetic distance (e.g. in time units or ANI)
# 	In other words the squared distance between the two lineages has expectation Ndim * 2 * D * t (where Ndim=2)
# The diffusivity is defined in squared distance units per phylogenetic divergence. 
# The distance units are the same as used for the radius of the sphere. E.g., if you set radius=1 and phylogenetic edge lengths are measured in Myr, then diffusivity will be defined in squared radii per Myr.
# The input tree may include monofurcations and multifurcations, but note that multifurcations are internally first split into bifurcations
fit_sbm_const = function(	trees, 								# either a single tree in phylo format, or a list of trees
							tip_latitudes, 						# either a 1D vector of size Ntips (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), listing geographical latitudes of the tips (in decimal degrees) of each tree
							tip_longitudes, 					# either a 1D vector of size Ntips (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), listing geographical longitudes of the tips (in decimal degrees) of each tree
							radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
							phylodistance_matrixes	= NULL,		# optional list of length Ntrees, each element being a phylodistance matrix for the tips of a tree. Hence phylodistance_matrixes[[tr]][i,j] is the phylogenetic (patristic) distance between tips i & j in tree tr. Can be used to specify distances between tips regardless of the edge lengths in the trees (i.e., trees are only used for topology). Some phylodistances may be NA or NaN; the corresponding tip pairs will be omitted from the fitting.
							clade_states			= NULL,		# optional, either an integer vector of length Ntips+Nnodes (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), specifying the discrete "state" of each tip and node in each tree. This can be used to limit independent contrasts to tip pairs whose total number of state-transitions (along their shortest path) is zero.
							planar_approximation	= FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
							only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
							only_distant_tip_pairs	= FALSE,	# logical, whether to only consider tip pairs located at distinct geographic locations
							min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
							max_MRCA_age			= Inf,		# numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set max_MRCA_age=Inf to disable this filter.
							max_phylodistance		= Inf,		# numeric, maximum allowed geodistance for an independent contrast to be included in the SBM fitting
							no_state_transitions	= FALSE,	# if TRUE, only tip pairs without state transitions along their shortest paths are considered. In particular, only tips in the same state are considered. Requires that clade_states[] is provided.
							only_state				= NULL,		# optional integer, specifying the state in which tip pairs (and their connecting ancestors) must be in order to be considered. Requires that clade_states[] is provided.
							min_diffusivity			= NULL,		# numeric, specifying the lower bound of allowed diffusivities. If omitted, it will be automatically chosen. Only relevant if planar_approximation==FALSE.
							max_diffusivity			= NULL,		# numeric, specifying the upper bound of allowed diffusivities. If omitted, it will be automatically chosen. Only relevant if planar_approximation==FALSE.
							Nbootstraps				= 0,		# integer, number of boostraps to perform. If <=0, no boostrapping is performed.
							NQQ						= 0,		# integer, optional number of simulations to perform for creating Q-Q plots of the theoretically expected distribution of geodistances vs the empirical distribution of geodistances (across independent contrasts). The resolution of the returned QQ plot will be equal to the number of independent contrasts used for fitting.
							SBM_PD_functor			= NULL,		# optional object, internally used SBM probability density functor
							focal_diffusivities		= NULL){	# optional integer vector, listing diffusivities D for which to calculate and return the loglikelihoods, e.g. for diagnostic purposes. This does not influence the actual fitting.
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
		if((!is.null(phylodistance_matrixes)) && (!(("list" %in% class(phylodistance_matrixes)) && (length(phylodistance_matrixes)==1)))){
			phylodistance_matrixes = list(phylodistance_matrixes)
		}
		if((!is.null(clade_states)) && (!(("list" %in% class(clade_states)) && (length(clade_states)==1)))){
			clade_states = list(clade_states)
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
		if(!is.null(phylodistance_matrixes)){
			if("list" %in% class(phylodistance_matrixes)){
				if(length(phylodistance_matrixes)!=Ntrees) return(list(success=FALSE,error=sprintf("Input phylodistance_matrixes list has length %d, but should be of length %d (number of trees)",length(phylodistance_matrixes),Ntrees)))
			}else if("matrix" %in% class(phylodistance_matrixes)){
				if(Ntrees!=1) return(list(success=FALSE,error=sprintf("Input phylodistance_matrixes was given as a single matrix, but expected a list of %d matrixes (number of trees)",Ntrees)))
				phylodistance_matrixes = list(phylodistance_matrixes)
			}else{
				return(list(success=FALSE, error="Unknown data type for phylodistance_matrixes: Exected a list of matrixes"))
			}
		}
		if(!is.null(clade_states)){
			if("list" %in% class(clade_states)){
				if(length(clade_states)!=Ntrees) return(list(success=FALSE,error=sprintf("Input list of clade_states has length %d, but should be of length %d (number of trees)",length(clade_states),Ntrees)))
			}else if(("numeric" %in% class(clade_states)) || ("integer" %in% class(clade_states))){
				if(Ntrees!=1) return(list(success=FALSE,error=sprintf("Input clade_states was given as a single vector, but expected a list of %d vectors (number of trees)",Ntrees)))
				if(length(clade_states)!=length(trees[[1]]$tip.label)) return(list(success=FALSE,error=sprintf("ERROR: Input clade_states was given as a single vector of length %d, but expected length %d (number of tips in the input tree)",length(clade_states),length(trees[[1]]$tip.label))))
				clade_states = list(clade_states)
			}
		}
	}else{
		return(list(success=FALSE,error=sprintf("Unknown data format '%s' for input trees[]: Expected a list of phylo trees or a single phylo tree",class(trees)[1])))
	}
	if(is.null(min_diffusivity) || is.na(min_diffusivity)) min_diffusivity = NaN;
	if(is.null(max_diffusivity) || is.na(max_diffusivity)) max_diffusivity = NaN;
	if(is.null(Nbootstraps) || is.na(Nbootstraps) || (Nbootstraps<0)) Nbootstraps = 0;
	if((!is.null(only_state)) && is.null(clade_states)) return(list(success=FALSE, error="Missing clade_states[], needed when only_state is specified"))
	if(no_state_transitions && is.null(clade_states)) return(list(success=FALSE, error="Missing clade_states[], needed when no_state_transitions=TRUE"))
	if(min_MRCA_time<0) min_MRCA_time = 0

	# loop through the trees and extract independet contrasts between sister clades
	phylodistances 		= numeric()
	geodistances		= numeric()
	tip_pairs_per_tree	= vector(mode="list", Ntrees)
	for(i in 1:Ntrees){
		tree 						= trees[[i]]
		tip_latitudes_this_tree  	= tip_latitudes[[i]]
		tip_longitudes_this_tree 	= tip_longitudes[[i]]
		clade_states_this_tree 		= clade_states[[i]]
		
		# make sure tree does not have multifurcations
		tree = multifurcations_to_bifurcations(tree)$tree
		
		# extract independent pairs of sister tips (independent contrasts)
		tip_pairs = extract_independent_sister_tips(tree)
		if(only_basal_tip_pairs){
			# calculate number of nodes between tip pairs
			edge_counts = get_pairwise_distances(tree, A=tip_pairs[,1], B=tip_pairs[,2], as_edge_counts=TRUE, check_input=FALSE)
			# only keep tip pairs with at most 2 edges connecting them
			keep_pairs 	= which(edge_counts<=2)
			tip_pairs 	= tip_pairs[keep_pairs,,drop=FALSE]
		}
		if((min_MRCA_time>0) || (max_MRCA_age<Inf)){
			MRCAs 			= get_pairwise_mrcas(tree, tip_pairs[,1], tip_pairs[,2], check_input=FALSE)
			clade_heights	= castor::get_all_distances_to_root(tree)
			tree_span		= max(clade_heights)
			keep_pairs		= which((clade_heights[MRCAs]>=min_MRCA_time) & (tree_span-clade_heights[MRCAs]<=max_MRCA_age))
			tip_pairs 		= tip_pairs[keep_pairs,,drop=FALSE]
		}
		
		# calculate geodistances between sister tips
		geodistances_this_tree = radius * geodesic_angles(tip_latitudes_this_tree[tip_pairs[,1]],tip_longitudes_this_tree[tip_pairs[,1]],tip_latitudes_this_tree[tip_pairs[,2]],tip_longitudes_this_tree[tip_pairs[,2]])

		# determine phylodistances between sister tips
		if(is.null(phylodistance_matrixes)){
			# calculate phylodistances from tree
			phylodistances_this_tree = get_pairwise_distances(tree, A=tip_pairs[,1], B=tip_pairs[,2], check_input=FALSE)
		}else{
			# get phylodistances from provided phylodistance matrixes
			phylodistance_matrix_this_tree = phylodistance_matrixes[[i]]
			if(nrow(phylodistance_matrix_this_tree)!=length(tree$tip.label)) return(list(success=FALSE,error=sprintf("ERROR: Input phylodistance_matrix #%d has %d rows, but expected %d rows (number of tips in the input tree)",i,nrow(phylodistance_matrix_this_tree),length(tree$tip.label))))
			if(ncol(phylodistance_matrix_this_tree)!=length(tree$tip.label)) return(list(success=FALSE,error=sprintf("ERROR: Input phylodistance_matrix #%d has %d columns, but expected %d columns (number of tips in the input tree)",i,ncol(phylodistance_matrix_this_tree),length(tree$tip.label))))
			phylodistances_this_tree = phylodistance_matrix_this_tree[tip_pairs]
		}

		# filter tip pairs
		# in particular, omit tip pairs with zero phylogenetic distance, because in that case the likelihood density is pathological
		# also omit tip pairs located at the same geographic location, if requested
		# also omit tips whose shortest path includes state transitions, if requested
		keep_pair = (is.finite(phylodistances_this_tree) & (phylodistances_this_tree>0))
		if(only_distant_tip_pairs) keep_pair = keep_pair & (geodistances_this_tree>0)
		if(max_phylodistance<Inf) keep_pair  = keep_pair & (phylodistances_this_tree<=max_phylodistance)
		if(no_state_transitions || (!is.null(only_state))){
			Ntransitions = count_transitions_between_clades(tree=tree, A=tip_pairs[,1], B=tip_pairs[,2], states=clade_states_this_tree, check_input=TRUE)
			if(no_state_transitions) keep_pair = keep_pair & (Ntransitions==0)
			if(!is.null(only_state)) keep_pair = keep_pair & (Ntransitions==0) & (clade_states_this_tree[tip_pairs[,1]]==only_state) & (clade_states_this_tree[tip_pairs[,2]]==only_state)
		}
		tip_pairs					= tip_pairs[keep_pair,,drop=FALSE]
		phylodistances_this_tree 	= phylodistances_this_tree[keep_pair]
		geodistances_this_tree		= geodistances_this_tree[keep_pair]
		tip_pairs_per_tree[[i]]		= tip_pairs
		
		if(nrow(tip_pairs)==0) next; # no valid tip pairs found in this tree

		phylodistances 	= c(phylodistances, phylodistances_this_tree)
		geodistances	= c(geodistances, geodistances_this_tree)
	}
	NC = length(phylodistances)
	if(NC==0) return(list(success=FALSE, error="No valid tip pairs left for extracting independent contrasts"))
	if(all(geodistances==0)) return(list(success=FALSE, error="The geodesic distances of all independent contrasts are zero"))

	bootstrap_diffusivities = NULL
	if(planar_approximation){
		# ignore spherical structure and simply translate geodesic_angles to geodesic distances, then fit Brownian motion on a plane with isotropic diffusivity matrix
		Ndim 			= 2 # dimensionality of the Brownian motion
		diffusivity		= 0.5 * (1/(Ndim*NC)) * sum(geodistances^2 / phylodistances);
		loglikelihood 	= -0.5*NC*Ndim*log(2*pi) - 0.5*Ndim*NC*log(2*diffusivity) - 0.5*Ndim*sum(log(phylodistances)) - (1/(4*diffusivity)) * sum(geodistances^2 / phylodistances);
		if(is.nan(diffusivity)) return(list(success=FALSE, error="Fitted diffusivity is NaN", Ncontrasts = NC))
		if(is.nan(loglikelihood)) return(list(success=FALSE, error="Loglikelihood of fitted diffusivity is NaN", Ncontrasts = NC))
		if(Nbootstraps>0){
			bootstrap_fit_diffusivities  = vector(mode="numeric", Nbootstraps)
			bootstrap_fit_loglikelihoods = vector(mode="numeric", Nbootstraps)
			bootstrap_loglikelihoods 	 = vector(mode="numeric", Nbootstraps)
			for(b in 1:Nbootstraps){
				bootstrap_distances 			= sapply(1:NC, FUN=function(p) sqrt(sum(stats::rnorm(n=2,mean=0,sd=sqrt(2*diffusivity*phylodistances[p]))^2))) # generate random distances according to 2-dimensional planar BM with the fitted diffusivity
				bootstrap_fit_diffusivities[b]	= 0.5 * (1/(Ndim*NC)) * sum(bootstrap_distances^2 / phylodistances);
				bootstrap_fit_loglikelihoods[b]	= -0.5*NC*Ndim*log(2*pi) - 0.5*Ndim*NC*log(2*bootstrap_fit_diffusivities[b]) - 0.5*Ndim*sum(log(phylodistances)) - (1/(4*bootstrap_fit_diffusivities[b])) * sum(bootstrap_distances^2 / phylodistances);
				bootstrap_loglikelihoods[b]		= -0.5*NC*Ndim*log(2*pi) - 0.5*Ndim*NC*log(2*diffusivity) - 0.5*Ndim*sum(log(phylodistances)) - (1/(4*diffusivity)) * sum(bootstrap_distances^2 / phylodistances);
			}
		}
	}else{
		# pre-calculate SBM probability density functor for efficiency
		if(is.null(SBM_PD_functor)) SBM_PD_functor = SBM_get_SBM_PD_functor_CPP(max_error = 1e-8, max_Legendre_terms = 200)
	
		# maximum-likelihood estimation based on full spherical Brownian motion model
		fit = fit_SBM_diffusivity_from_transitions_CPP(	radius 					= radius,
														time_steps				= phylodistances,
														distances 				= geodistances,
														max_error				= 1e-8,
														max_Legendre_terms		= 200,
														opt_epsilon				= 1e-10,
														max_iterations 			= 10000,
														min_diffusivity			= min_diffusivity,
														max_diffusivity			= max_diffusivity,
														Nbootstraps				= Nbootstraps,
														SBM_PD_functor			= SBM_PD_functor);
		if(!fit$success) return(list(success=FALSE, error=fit$error, Ncontrasts = NC))
		diffusivity 					= fit$fit_diffusivity
		loglikelihood					= fit$fit_loglikelihood
		bootstrap_fit_diffusivities		= fit$bootstrap_fit_diffusivities
		bootstrap_fit_loglikelihoods	= fit$bootstrap_fit_loglikelihoods
		bootstrap_loglikelihoods		= fit$bootstrap_loglikelihoods		
	}
	
	if(!is.null(focal_diffusivities)){
		focal_LLs = SBM_LLs_of_transitions_CPP(	radius 				= radius,
												time_steps			= phylodistances,
												distances 			= geodistances,
												diffusivities		= focal_diffusivities,
												max_error			= 1e-7,
												max_Legendre_terms	= 200)$loglikelihoods
	}else{
		focal_LLs = numeric();
	}
	
	# calculate data consistency with the fitted model, defined as the probability that the loglikelihoods of data generated by the model would deviate similarly or moe from the expected loglikelihood than the observed data does.
	if(Nbootstraps>0){
		mean_BLL = mean(bootstrap_loglikelihoods,na.rm=TRUE)
		consistency = sum(abs(bootstrap_loglikelihoods-mean_BLL)>=abs(loglikelihood-mean_BLL),na.rm=TRUE)/sum(!is.nan(bootstrap_loglikelihoods))
	}else{
		consistency = NA
	}
	
		
	#####################################
	# Calculate QQ-plot using simulations
	
	if(NQQ>0){
		sim_geodistances = numeric(NQQ * NC)
		next_g = 1 # index in sim_geodistances[] for placing the next simulated geodistance
		for(tr in 1:Ntrees){
			tip_pairs = tip_pairs_per_tree[[tr]]
			if(length(tip_pairs_per_tree[[tr]])>0){
				for(q in 1:NQQ){
					sim = castor::simulate_sbm(tree = trees[[tr]], radius = radius, diffusivity = diffusivity, root_latitude = NULL, root_longitude = NULL)
					if(!sim$success) return(list(success=FALSE, error=sprintf("Calculation of QQ failed at simulation %d for tree %d: Could not simulate SBM for the fitted model: %s",q,tr,sim$error), diffusivity=diffusivity, loglikelihood=loglikelihood, Ncontrasts=NC));
					sim_geodistances[next_g + c(1:nrow(tip_pairs))] = radius * geodesic_angles(sim$tip_latitudes[tip_pairs[,1]],sim$tip_longitudes[tip_pairs[,1]],sim$tip_latitudes[tip_pairs[,2]],sim$tip_longitudes[tip_pairs[,2]])
					next_g = next_g + nrow(tip_pairs)
				}
			}
		}
		probs  = c(1:NC)/NC
		QQplot = cbind(quantile(geodistances, probs=probs, na.rm=TRUE, type=8), quantile(sim_geodistances, probs=probs, na.rm=TRUE, type=8))
	}
	
	
	#####################################

	return(list(success					= TRUE,
				diffusivity 			= diffusivity,
				loglikelihood			= loglikelihood,
				Ncontrasts				= NC,
				phylodistances 			= phylodistances,
				geodistances			= geodistances,
				tip_pairs_per_tree		= tip_pairs_per_tree,
				focal_loglikelihoods	= focal_LLs,
				standard_error			= (if(Nbootstraps>0) sd(bootstrap_fit_diffusivities, na.rm=TRUE) else NULL),
				CI50lower				= (if(Nbootstraps>0) unname(quantile(bootstrap_fit_diffusivities, probs=0.25, type=8, na.rm=TRUE)) else NULL),
				CI50upper				= (if(Nbootstraps>0) unname(quantile(bootstrap_fit_diffusivities, probs=0.75, type=8, na.rm=TRUE)) else NULL),
				CI95lower				= (if(Nbootstraps>0) unname(quantile(bootstrap_fit_diffusivities, probs=0.025, type=8, na.rm=TRUE)) else NULL),
				CI95upper				= (if(Nbootstraps>0) unname(quantile(bootstrap_fit_diffusivities, probs=0.975, type=8, na.rm=TRUE)) else NULL),
				consistency				= consistency,
				QQplot					= (if(NQQ>0) QQplot else NULL),
				SBM_PD_functor			= SBM_PD_functor))
}




