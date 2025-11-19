# Fit a Spherical Brownian Motion model of diffusive geographic dispersal to two different sets of trees, and compare the fitted models
# This function can calculate the two-sided statistical significance of the diffusivity ratio (D1/D2), or equivalently, the statistical significance of the log_difference |log(D1)-log(D2)|
# The main inputs are:
# 	Two lists of trees
#	The corresponding geographic coordinates of the tips
fit_and_compare_sbm_const = function(	trees1, 				# either a single tree in phylo format, or a list of trees
										tip_latitudes1, 		# either a 1D vector of size Ntips (if trees1[] is a single tree) or a list of 1D vectors (if trees1[] is a list of trees), listing geographical latitudes of the tips (in decimal degrees) of each tree in trees1[]
										tip_longitudes1, 		# either a 1D vector of size Ntips (if trees1[] is a single tree) or a list of 1D vectors (if trees1[] is a list of trees), listing geographical longitudes of the tips (in decimal degrees) of each tree in trees1[]
										trees2, 				# either a single tree in phylo format, or a list of trees
										tip_latitudes2, 		# either a 1D vector of size Ntips (if trees2[] is a single tree) or a list of 1D vectors (if trees2[] is a list of trees), listing geographical latitudes of the tips (in decimal degrees) of each tree in trees2[]
										tip_longitudes2, 		# either a 1D vector of size Ntips (if trees2[] is a single tree) or a list of 1D vectors (if trees2[] is a list of trees), listing geographical longitudes of the tips (in decimal degrees) of each tree in trees2[]
										radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
										planar_approximation	= FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
										only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
										only_distant_tip_pairs	= FALSE,	# logical, whether to only consider tip pairs located at distinct geographic locations
										min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
										max_MRCA_age			= Inf,		# numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set max_MRCA_age=Inf to disable this filter.
										max_phylodistance		= Inf,		# numeric, maximum allowed geodistance for an independent contrast to be included in the SBM fitting
										min_diffusivity			= NULL,		# numeric, specifying the lower bound of allowed diffusivities. If omitted, it will be automatically chosen. Only relevant if planar_approximation==FALSE.
										max_diffusivity			= NULL,		# numeric, specifying the upper bound of allowed diffusivities. If omitted, it will be automatically chosen. Only relevant if planar_approximation==FALSE.
										Nbootstraps				= 0,		# integer, number of bootstraps to perform for estimating confidence intervals. Set to <=0 for no bootstrapping.
										Nsignificance			= 0,		# integer, number of random simulations to perform for testing the statistical significance of the log_difference of diffusivities (|log(D1)-log(D2)|). If <=0, the statistical significance is not calculated
										SBM_PD_functor			= NULL,		# optional object, internally used SBM probability density functor
										verbose					= FALSE,
										verbose_prefix 			= ""){	
	# fit SBM to each set of trees
	if(verbose) cat(sprintf("%sFitting SBM to first tree set..\n",verbose_prefix))
	fit1 = fit_sbm_const(	trees					= trees1, 
							tip_latitudes			= tip_latitudes1,
							tip_longitudes			= tip_longitudes1,
							radius					= radius,
							planar_approximation	= planar_approximation,
							only_basal_tip_pairs	= only_basal_tip_pairs,
							only_distant_tip_pairs	= only_distant_tip_pairs,
							min_MRCA_time			= min_MRCA_time,
							max_MRCA_age			= max_MRCA_age,
							max_phylodistance		= max_phylodistance,
							min_diffusivity			= min_diffusivity,
							max_diffusivity			= max_diffusivity,
							Nbootstraps				= Nbootstraps,
							SBM_PD_functor			= SBM_PD_functor)
	if(!fit1$success) return(list(success=FALSE, error=sprintf("Failed to fit SBM to tree set 1: %s",fit1$error)))
	if(is.null(SBM_PD_functor)) SBM_PD_functor = fit1$SBM_PD_functor

	if(verbose) cat(sprintf("%sFitting SBM to second tree set..\n",verbose_prefix))
	fit2 = fit_sbm_const(	trees					= trees2, 
							tip_latitudes			= tip_latitudes2,
							tip_longitudes			= tip_longitudes2,
							radius					= radius,
							planar_approximation	= planar_approximation,
							only_basal_tip_pairs	= only_basal_tip_pairs,
							only_distant_tip_pairs	= only_distant_tip_pairs,
							min_MRCA_time			= min_MRCA_time,
							max_MRCA_age			= max_MRCA_age,
							max_phylodistance		= max_phylodistance,
							min_diffusivity			= min_diffusivity,
							max_diffusivity			= max_diffusivity,
							Nbootstraps				= Nbootstraps,
							SBM_PD_functor			= SBM_PD_functor)
	if(!fit2$success) return(list(success=FALSE, error=sprintf("Failed to fit SBM to tree set 2: %s",fit2$error)))
	
	lin_difference	= abs(fit1$diffusivity - fit2$diffusivity)
	log_difference 	= abs(log(fit1$diffusivity) - log(fit2$diffusivity))
	if(Nsignificance>0){
		# calculate statistical significance of log_difference
		if(verbose) cat(sprintf("%sCalculating statistical significances..\n",verbose_prefix))
		# bring all data into a standard format
		if("phylo" %in% class(trees1)) trees1 = list(trees1)
		if("phylo" %in% class(trees2)) trees2 = list(trees2)
		Ntrees1 = length(trees1)
		Ntrees2 = length(trees2)
		if(!("list" %in% class(tip_latitudes1))) 	tip_latitudes1 	= list(tip_latitudes1)
		if(!("list" %in% class(tip_longitudes1))) 	tip_longitudes1 = list(tip_longitudes1)
		if(!("list" %in% class(tip_latitudes2))) 	tip_latitudes2 	= list(tip_latitudes2)
		if(!("list" %in% class(tip_longitudes2))) 	tip_longitudes2 = list(tip_longitudes2)
		if(verbose) cat(sprintf("%s  Fitting common SBM model to both tree sets..\n",verbose_prefix))
		fit_common = fit_sbm_const(	trees					= c(trees1,trees2), 
									tip_latitudes			= c(tip_latitudes1,tip_latitudes2), 
									tip_longitudes			= c(tip_longitudes1,tip_longitudes2), 
									radius					= radius,
									planar_approximation	= planar_approximation,
									only_basal_tip_pairs	= only_basal_tip_pairs,
									only_distant_tip_pairs	= only_distant_tip_pairs,
									min_MRCA_time			= min_MRCA_time,
									max_MRCA_age			= max_MRCA_age,
									max_phylodistance		= max_phylodistance,
									min_diffusivity			= min_diffusivity,
									max_diffusivity			= max_diffusivity,
									Nbootstraps				= 0,
									SBM_PD_functor			= SBM_PD_functor)
		if(verbose) cat(sprintf("%s  Assessing significance over %d SBM simulations..\n",verbose_prefix,Nsignificance))
		random_tip_latitudes1 	= vector(mode="list", Ntrees1)
		random_tip_longitudes1 	= vector(mode="list", Ntrees1)
		random_tip_latitudes2 	= vector(mode="list", Ntrees2)
		random_tip_longitudes2 	= vector(mode="list", Ntrees2)
		Ngreater_lin = 0
		Ngreater_log = 0
		Nsuccess 	 = 0
		for(r in 1:Nsignificance){
			# simulate the same SBM model on each of the input trees
			simulation_failed = FALSE
			for(tr in 1:Ntrees1){
				sim = simulate_sbm(tree = trees1[[tr]], radius = radius, diffusivity = fit_common$diffusivity)
				if(!sim$success){
					if(verbose) cat(sprintf("%s  WARNING: Simulation #%d of common SBM model failed for tree1 #%d\n",verbose_prefix,r,tr))
					simulation_failed = TRUE
					break
				}
				random_tip_latitudes1[[tr]]  = sim$tip_latitudes
				random_tip_longitudes1[[tr]] = sim$tip_longitudes
			}
			for(tr in 1:Ntrees2){
				sim = simulate_sbm(tree = trees2[[tr]], radius = radius, diffusivity = fit_common$diffusivity)
				if(!sim$success){
					if(verbose) cat(sprintf("%s  WARNING: Simulation #%d of common SBM model failed for tree2 #%d\n",verbose_prefix,r,tr))
					simulation_failed = TRUE
					break
				}
				random_tip_latitudes2[[tr]]  = sim$tip_latitudes
				random_tip_longitudes2[[tr]] = sim$tip_longitudes
			}
			if(simulation_failed) next; 
			
			# fit SBM separately to each tree set
			random_fit1 = fit_sbm_const(trees					= trees1, 
										tip_latitudes			= random_tip_latitudes1,
										tip_longitudes			= random_tip_longitudes1,
										radius					= radius,
										planar_approximation	= planar_approximation,
										only_basal_tip_pairs	= only_basal_tip_pairs,
										only_distant_tip_pairs	= only_distant_tip_pairs,
										min_MRCA_time			= min_MRCA_time,
										max_MRCA_age			= max_MRCA_age,
										max_phylodistance		= max_phylodistance,
										min_diffusivity			= min_diffusivity,
										max_diffusivity			= max_diffusivity,
										Nbootstraps				= 0,
										SBM_PD_functor			= SBM_PD_functor)
			if(!random_fit1$success){
				# fitting failed
				if(verbose) cat(sprintf("%s  WARNING: SBM fitting failed for random simulation #%d, for tree set 1\n",verbose_prefix,r))
				next;
			}
			random_fit2 = fit_sbm_const(trees					= trees2, 
										tip_latitudes			= random_tip_latitudes2,
										tip_longitudes			= random_tip_longitudes2,
										radius					= radius,
										planar_approximation	= planar_approximation,
										only_basal_tip_pairs	= only_basal_tip_pairs,
										only_distant_tip_pairs	= only_distant_tip_pairs,
										min_MRCA_time			= min_MRCA_time,
										max_MRCA_age			= max_MRCA_age,
										max_phylodistance		= max_phylodistance,
										min_diffusivity			= min_diffusivity,
										max_diffusivity			= max_diffusivity,
										Nbootstraps				= 0,
										SBM_PD_functor			= SBM_PD_functor)
			if(!random_fit2$success){
				# fitting failed
				if(verbose) cat(sprintf("%s  WARNING: SBM fitting failed for random simulation #%d, for tree set 2\n",verbose_prefix,r))
				next;
			}
			Nsuccess = Nsuccess + 1
			# compare the obtained lin_difference & log_difference to the original ones
			random_lin_difference = abs(random_fit1$diffusivity - random_fit2$diffusivity)
			random_log_difference = abs(log(random_fit1$diffusivity) - log(random_fit2$diffusivity))
			Ngreater_lin = Ngreater_lin + (random_lin_difference>=lin_difference)
			Ngreater_log = Ngreater_log + (random_log_difference>=log_difference)
		}
		lin_significance = Ngreater_lin / Nsuccess
		log_significance = Ngreater_log / Nsuccess
	}
		
	return(list(success				= TRUE,
				fit1				= fit1,
				fit2				= fit2,
				lin_difference		= lin_difference,
				log_difference 		= log_difference,
				lin_significance	= (if(Nsignificance>0) lin_significance else NULL),
				log_significance	= (if(Nsignificance>0) log_significance else NULL),
				fit_common			= (if(Nsignificance>0) fit_common else NULL)))

}