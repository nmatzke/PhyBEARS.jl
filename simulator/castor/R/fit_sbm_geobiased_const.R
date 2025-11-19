# fit a constant-diffusivity SBM model to a tree with given tip coordinates, while iteratively correcting for geographic sampling biases
fit_sbm_geobiased_const = function(	trees, 								# either a single tree in phylo format, or a list of trees
									tip_latitudes, 						# either a 1D vector of size Ntips (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), listing geographical latitudes of the tips (in decimal degrees) of each tree
									tip_longitudes, 					# either a 1D vector of size Ntips (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), listing geographical longitudes of the tips (in decimal degrees) of each tree
									radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
									reference_latitudes		= NULL,		# optional numeric vector of length NR, listing latitudes of reference coordinates based on which to calculate the geographic sampling density. If NULL, the geographic sampling density is estimated based on tip_latitudes[] and tip_longitudes[]
									reference_longitudes	= NULL,		# optional numeric vector of length NR, listing longitudes of reference coordinates based on which to calculate the geographic sampling density. If NULL, the geographic sampling density is estimated based on tip_latitudes[] and tip_longitudes[]
									only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
									only_distant_tip_pairs	= FALSE,	# logical, whether to only consider tip pairs located at distinct geographic locations
									min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
									max_MRCA_age			= Inf,		# numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set max_MRCA_age=Inf to disable this filter.
									max_phylodistance		= Inf,		# numeric, maximum allowed geodistance for an independent contrast to be included in the SBM fitting
									min_diffusivity			= NULL,		# numeric, specifying the lower bound of allowed diffusivities, for the first iteration. If omitted, it will be automatically chosen.
									max_diffusivity			= NULL,		# numeric, specifying the upper bound of allowed diffusivities, for the first iteration. If omitted, it will be automatically chosen.
									rarefaction				= 0.1,		# numeric, by how much to rarefy the simulated trees when geographically subsampling. This should be less than 1, so that geographic bias actually plays in. Note that the simulated trees will have the same size as the original tree.
									Nsims					= 100,		# integer, number of SBM simulatons to perform per iteration for assessing the effects of geographic bias. Smaller trees require larger Nsims (due to higher stochasticity). This must be at least 2, although values of 10 or even 100 are better.
									max_iterations			= 100,		# integer, maximum number of iterations (correction steps) to perform. Note that due to the stochastic nature of the SBM, exact convergence is not possible.
									Nbootstraps				= 0,		# integer, number of boostraps to perform, in the final iteration. If <=0, no boostrapping is performed.
									NQQ						= 0,		# integer, optional number of simulations to perform for creating Q-Q plots of the theoretically expected distribution of geodistances vs the empirical distribution of geodistances (across independent contrasts). The resolution of the returned QQ plot will be equal to the number of independent contrasts used for fitting.
									Nthreads				= 1,		# integer, number of parallel thread to use whenever possible
									include_simulations		= FALSE,	# logical, include the final simulated trees and coordinates in the returned results. This may be useful e.g. for checking the adequacy of the fitted model.
									SBM_PD_functor			= NULL,		# optional object, internally used SBM probability density functor
									verbose					= FALSE,
									verbose_prefix			= ""){
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
	Nsims 			= max(2,Nsims)
	max_iterations 	= max(1,max_iterations)
	rarefaction 	= max(0,min(1,rarefaction))
	
	# get flattened set of reference coordinates (across all tips)
	if(is.null(reference_latitudes) || is.null(reference_longitudes)){
		reference_latitudes  = unlist(tip_latitudes)
		reference_longitudes = unlist(tip_longitudes)
	}
	
	# fit a simple birth-death model to each tree, for simulating "similar" trees later on
	if(verbose) cat(sprintf("%sFitting basic birth-death models to %d trees..\n",verbose_prefix,Ntrees))
	BD_lambdas 	= rep(NA,Ntrees)
	BD_mus		= rep(NA,Ntrees)
	BD_rhos		= rep(NA,Ntrees)
	for(tr in seq_len(Ntrees)){
		if(length(trees[[tr]]$tip.label)<=4) next # this tree is too small
		BDfit = castor::fit_hbd_model_on_grid(	tree			= trees[[tr]],
												const_lambda	= TRUE,
												const_mu		= TRUE,
												guess_rho0		= rarefaction,
												min_rho0		= 0.0001*rarefaction,
												condition		= "auto",
												Ntrials			= 100,
												Nthreads		= Nthreads)
		if(!BDfit$success) next
		# figure out congruent BD model with the appropriate rho
		# the reason we don't fit with rho fixed at rarefaction, is that for some trees there is an upper bound on the biologically plausible values of rho (with higher values leading to mu<=0).
		BD_rhos[tr] 	= min(rarefaction,(if(BDfit$fitted_lambda<=BDfit$fitted_mu) 1 else BDfit$fitted_rho/(1-BDfit$fitted_mu/BDfit$fitted_lambda))) # set rho to the provided rarefaction, but avoid getting above a certain threshold where the congruent mu would be negative
		BD_lambdas[tr] 	= max(0,BDfit$fitted_lambda * BDfit$fitted_rho / BD_rhos[tr])
		BD_mus[tr] 		= max(0,BD_lambdas[tr] - (BDfit$fitted_lambda - BDfit$fitted_mu))
	}
	successfull_BD_fits = which(is.finite(BD_lambdas))
	if(length(successfull_BD_fits)==0) return(list(success=FALSE, error=sprintf("BD model fitting failed for all %d trees",Ntrees)))
	if(length(successfull_BD_fits)<Ntrees){
		# filter out failed trees
		if(verbose) cat(sprintf("%s  WARNING: Ignoring %d out of %d trees for which BD model fitting failed\n",verbose_prefix,Ntrees-length(successfull_BD_fits),Ntrees))
		BD_lambdas 		= BD_lambdas[successfull_BD_fits]
		BD_mus			= BD_mus[successfull_BD_fits]
		BD_rhos			= BD_rhos[successfull_BD_fits]
		trees			= trees[successfull_BD_fits]
		tip_latitudes	= tip_latitudes[successfull_BD_fits]
		tip_longitudes	= tip_longitudes[successfull_BD_fits]
		Ntrees			= length(successfull_BD_fits)
	}
	if(verbose) cat(sprintf("%s  Note: Congruents of fitted models have mean lambda = %g, mu = %g, r = %g, rho = %g\n",verbose_prefix,mean(BD_lambdas),mean(BD_mus),mean(BD_lambdas-BD_mus),mean(BD_rhos)))
	
	# fit SBM-const (without any correction)
	if(verbose) cat(sprintf("%sFitting diffusivity without any correction..\n",verbose_prefix))
	fit0 = fit_sbm_const(	trees					= trees,
							tip_latitudes			= tip_latitudes,
							tip_longitudes			= tip_longitudes,
							radius					= radius,
							only_basal_tip_pairs	= only_basal_tip_pairs,
							only_distant_tip_pairs	= only_distant_tip_pairs,
							min_MRCA_time			= min_MRCA_time,
							max_MRCA_age			= max_MRCA_age,
							max_phylodistance		= max_phylodistance,
							min_diffusivity			= min_diffusivity,
							max_diffusivity			= max_diffusivity,
							Nbootstraps				= Nbootstraps,
							NQQ						= 0,
							SBM_PD_functor			= SBM_PD_functor)
	if(!fit0$success) return(list(success=FALSE, error=sprintf("Could not fit SBM model in first round: %s",fit0$error)))
	if(verbose) cat(sprintf("%s  Fitted diffusivity: %g\n",verbose_prefix,fit0$diffusivity))

	
	root_ages 					= sapply(seq_len(Ntrees), FUN=function(tr) get_tree_span(trees[[tr]])$max_distance)
	correction_factor 			= 1 # multiplicative factor to multiply the last fit diffusivity with, in order to get the true diffusivity
	all_correction_factors 		= rep(NA,min(1000,max_iterations))
	all_diffusivity_estimates 	= rep(NA,min(1000,max_iterations)) # all_diffusivity_estimates[i] will be the geometric-mean diffusivity estimated for trees simulated assuming a correction factor all_correction_factors[i]
	Niterations			 		= 0
	stopping_criterion 			= NULL
	sims 						= vector(mode="list", Ntrees)
	Nsims_per_tree 				= numeric(Ntrees)
	tile_counts 				= NULL
	tile_latitudes 				= NULL
	tile_longitudes 			= NULL
	converged					= FALSE
	for(i in seq_len(max_iterations)){
		if(verbose) cat(sprintf("%sIteration %d\n",verbose_prefix,i))
		if(verbose) cat(sprintf("%s  Simulating %d replicate SBM models for each input tree (with D=%g)..\n",verbose_prefix,Nsims,fit0$diffusivity * correction_factor))
		for(tr in seq_len(Ntrees)){
			sims[[tr]] = simulate_geobiased_sbm(Nsims					= Nsims,
												Ntips					= length(trees[[tr]]$tip.label),
												radius					= radius,
												diffusivity				= fit0$diffusivity * correction_factor,
												lambda					= BD_lambdas[tr],
												mu						= BD_mus[tr],
												rarefaction 			= BD_rhos[tr],
												crown_age				= root_ages[tr],
												Nthreads				= Nthreads,
												omit_failed_sims		= TRUE,
												reference_latitudes		= reference_latitudes,
												reference_longitudes	= reference_longitudes,
												tile_counts				= tile_counts,
												tile_latitudes			= tile_latitudes,
												tile_longitudes			= tile_longitudes)
			if(!sims[[tr]]$success) return(list(success=FALSE, error=sprintf("Iteration %d failed for tree #%d: Simulations failed: %s",i,tr,sims[[tr]]$error)))
			Nsims_per_tree[tr] = length(sims[[tr]]$sims)
			# keep record of geographic sampling density (spherical tiles) to avoid costly re-computing
			if(is.null(tile_counts)){
				tile_counts 	= sims[[tr]]$tile_counts
				tile_latitudes 	= sims[[tr]]$tile_latitudes
				tile_longitudes = sims[[tr]]$tile_longitudes
			}
			if(length(sims[[tr]]$sims)==0){
				sims[[tr]]$success = FALSE
				next
			}
		}
		trees_with_valid_sims = which(sapply(seq_len(Ntrees), FUN=function(tr) sims[[tr]]$success))
		if(length(trees_with_valid_sims)==0) return(list(success=FALSE, error=sprintf("Iteration %d failed: Simulations failed completely for all trees",i)))
		
		aux_fit_SBM_to_simulation = function(r){
			# extract r-th simulation for each tree (recycle simulations cyclically if fewer than r are available e.g. due to some simulation failures)
			sim_trees 		= lapply(trees_with_valid_sims, FUN=function(tr) sims[[tr]]$sims[[1+(r-1) %% Nsims_per_tree[tr]]]$tree)
			sim_latitudes 	= lapply(trees_with_valid_sims, FUN=function(tr) sims[[tr]]$sims[[1+(r-1) %% Nsims_per_tree[tr]]]$latitudes)
			sim_longitudes 	= lapply(trees_with_valid_sims, FUN=function(tr) sims[[tr]]$sims[[1+(r-1) %% Nsims_per_tree[tr]]]$longitudes)
			sim_fit = fit_sbm_const(trees					= sim_trees, 
									tip_latitudes			= sim_latitudes,
									tip_longitudes			= sim_longitudes,
									radius 					= radius,
									only_basal_tip_pairs	= only_basal_tip_pairs,
									only_distant_tip_pairs	= only_distant_tip_pairs,
									min_MRCA_time			= min_MRCA_time,
									max_MRCA_age			= max_MRCA_age,
									max_phylodistance		= max_phylodistance,
									SBM_PD_functor			= SBM_PD_functor)
			return(if(sim_fit$success) sim_fit$diffusivity else NA)
		}
		if((Nsims>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
			if(verbose) cat(sprintf("%s  Fitting SBM diffusivity to %d simulations (parallelized)..\n",verbose_prefix,Nsims))
			sim_fit_diffusivities = unlist(parallel::mclapply(	seq_len(Nsims), 
																FUN = function(r){ aux_fit_SBM_to_simulation(r) }, 
																mc.cores = min(Nthreads, Nsims), 
																mc.preschedule = TRUE, 
																mc.cleanup = TRUE))
		}else{
			if(verbose) cat(sprintf("%s  Fitting SBM diffusivity to %d simulations (sequentially)..\n",verbose_prefix,Nsims))
			sim_fit_diffusivities = rep(NA, Nsims)
			for(r in seq_len(Nsims)){
				sim_fit_diffusivities[r] = aux_fit_SBM_to_simulation(r)
			}
		}

		valid_fits = which(is.finite(sim_fit_diffusivities))
		if(length(valid_fits)==0) return(list(success=FALSE, error=sprintf("Iteration %d failed: Could not fit SBM to any of the simulated datasets",i)))
		if(length(valid_fits)==1) return(list(success=FALSE, error=sprintf("Iteration %d failed: Could only fit SBM to one of the simulated datasets, but need at least 2 successful sims",i)))
		sim_fit_diffusivities = sim_fit_diffusivities[valid_fits]
		if(length(sim_fit_diffusivities)>=5) sim_fit_diffusivities = remove_outliers(X=sim_fit_diffusivities, outlier_prob=0.1)
		if(length(sim_fit_diffusivities)<2) return(list(success=FALSE, error=sprintf("Iteration %d failed: Nearly all sim-fitted diffusivities were filtered out as 'outliers'.",i)))
		mean_sim_fit_diffusivity 		= exp(mean(log(sim_fit_diffusivities))) # geometric mean
		se_sim_fit_log_diffusivity	 	= sd(log(sim_fit_diffusivities))/sqrt(length(sim_fit_diffusivities)) # standard error in log space
		all_correction_factors[i]		= correction_factor
		all_diffusivity_estimates[i] 	= mean_sim_fit_diffusivity
		Niterations 					= Niterations + 1
		if(verbose) cat(sprintf("%s  Geometric-mean diffusivity fitted to sims: %g (standard error of mean log-D %g)\n",verbose_prefix,mean_sim_fit_diffusivity,se_sim_fit_log_diffusivity))

		# check if approximate "convergence" was achieved
		if((abs(log(fit0$diffusivity)-log(mean_sim_fit_diffusivity))<min(se_sim_fit_log_diffusivity,0.1)) && (i>1) && (abs(log(correction_factor/all_correction_factors[i-1]))<0.1)){
			# mean_sim_fit_diffusivity is pretty close (by one standard deviation) from the fit0, so consider this a successful convergence
			if(verbose) cat(sprintf("%s  Achieved approximate convergence at iteration %d: true estimated diffusivity = %g\n",verbose_prefix,i,fit0$diffusivity*correction_factor))
			stopping_criterion = sprintf("Achieved approximate convergence at iteration %d (relative difference to uncorrected fit: %g)",i,abs(fit0$diffusivity-mean_sim_fit_diffusivity)/fit0$diffusivity)
			converged = TRUE
			break
		}
		
		# improve correction factor for next iteration
		if(i<=2){
			correction_factor = correction_factor * fit0$diffusivity/mean_sim_fit_diffusivity
		}else{
			# use linear interpolation (on a log scale) between the closest two iterations (closest to the target from below and above) if possible, to estimate improved correction factor
			aboves = which(all_diffusivity_estimates-fit0$diffusivity>=0)
			belows = which(all_diffusivity_estimates-fit0$diffusivity<=0)
			if((length(aboves)>0) && (length(belows)>0)){
				# found closest from below (left) and above (right)
				left  = belows[which.min(fit0$diffusivity-all_diffusivity_estimates[belows])]
				right = aboves[which.min(all_diffusivity_estimates[aboves]-fit0$diffusivity)]
				if((i!=left) && (i!=right)){
					# both points seem to be old, hence there's no point in trying them out again. Replace one of the two with the current iteration i.
					if((all_diffusivity_estimates[i]-fit0$diffusivity) * (all_diffusivity_estimates[left]-fit0$diffusivity)<0){
						# i and left landed on opposite sides of the target, so use those
						right = i
					}else{
						left = i
					}
				}
			}else{
				# can't bracket the target yet, so just use last 2 iterations
				left 	= i-1
				right	= i
			}
			correction_factor = exp(log(all_correction_factors[left]) + (log(fit0$diffusivity)-log(all_diffusivity_estimates[left])) * (log(all_correction_factors[right])-log(all_correction_factors[left]))/(log(all_diffusivity_estimates[right])-log(all_diffusivity_estimates[left])))
		}
		if(verbose) cat(sprintf("%s  Estimated correction factor: %g\n",verbose_prefix,correction_factor))
		if(!is.finite(correction_factor)) return(list(	success						= FALSE, 
														error						= sprintf("Sequence diverged (correction factor non-finite)"),
														Niterations					= Niterations, 
														uncorrected_fit_diffusivity = fit0$diffusivity, 
														all_correction_factors 		= all_correction_factors[seq_len(Niterations)],
														all_diffusivity_estimates 	= all_diffusivity_estimates[seq_len(Niterations)]))		
	}
	if(!converged) return(list(success=FALSE, error="Failed to converge", Nlat = sims[[1]]$Nlat, Nlon = sims[[1]]$Nlon, uncorrected_fit_diffusivity	= fit0$diffusivity, Ncontrasts = fit0$Ncontrasts, Ntrees = Ntrees))
	true_diffusivity = fit0$diffusivity * correction_factor
	
	# Calculate QQ-plot using simulations, based on estimated true diffusivity	
	if(NQQ>0){
		sim_geodistances = numeric(NQQ * fit0$Ncontrasts)
		next_g = 1 # index in sim_geodistances[] for placing the next simulated geodistance
		for(tr in 1:Ntrees){
			tip_pairs = fit0$tip_pairs_per_tree[[tr]]
			if(length(tip_pairs)>0){
				for(q in 1:NQQ){
					sim = castor::simulate_sbm(tree = trees[[tr]], radius = radius, diffusivity = true_diffusivity, root_latitude = NULL, root_longitude = NULL)
					if(!sim$success) return(list(success=FALSE, error=sprintf("Calculation of QQ failed at simulation %d for tree %d: Could not simulate SBM for the fitted model: %s",q,tr,sim$error), diffusivity=true_diffusivity));
					sim_geodistances[next_g + c(1:nrow(tip_pairs))] = radius * geodesic_angles(sim$tip_latitudes[tip_pairs[,1]],sim$tip_longitudes[tip_pairs[,1]],sim$tip_latitudes[tip_pairs[,2]],sim$tip_longitudes[tip_pairs[,2]])
					next_g = next_g + nrow(tip_pairs)
				}
			}
		}
		probs  = seq_len(fit0$Ncontrasts)/fit0$Ncontrasts
		QQplot = cbind(quantile(fit0$geodistances, probs=probs, na.rm=TRUE, type=8), quantile(sim_geodistances, probs=probs, na.rm=TRUE, type=8))
	}

	return(list(success						= TRUE,
				Nlat						= sims[[1]]$Nlat,
				Nlon						= sims[[1]]$Nlon,
				diffusivity					= true_diffusivity,
				correction_factor			= correction_factor,
				Niterations					= Niterations,
				stopping_criterion			= stopping_criterion,
				uncorrected_fit_diffusivity	= fit0$diffusivity,
				last_sim_fit_diffusivity	= mean_sim_fit_diffusivity,
				all_correction_factors		= all_correction_factors[seq_len(Niterations)], # correction factor at each iteration
				all_diffusivity_estimates	= all_diffusivity_estimates[seq_len(Niterations)], # mean diffusivity estimates for trees simulated at each iteration
				Ntrees						= Ntrees,
				lambda						= BD_lambdas,
				mu							= BD_mus,
				rarefaction					= rarefaction,
				Ncontrasts					= fit0$Ncontrasts,
				standard_error				= fit0$standard_error * correction_factor,
				CI50lower					= fit0$CI50lower * correction_factor,
				CI50upper					= fit0$CI50upper * correction_factor,
				CI95lower					= fit0$CI95lower * correction_factor,
				CI95upper					= fit0$CI95upper * correction_factor,
				QQplot						= (if(NQQ>0) QQplot else NULL),
				simulations					= (if(include_simulations) sims else NULL),
				SBM_PD_functor				= fit0$SBM_PD_functor))
}