# fit a homogenous birth-death pulled diversification rate (PDR) and rho*lambda0 on a grid to a given extant timetree, choosing the "best" grid size according to AIC or BIC
fit_hbd_pdr_on_best_grid_size = function(	tree, 
											oldest_age			= NULL,		# numeric, the oldest age to consider in the evaluation of the likelihood as well as for defining the age grid. Typically set to the stem age or root age. Can be NULL (equivalent to the root age).
											age0				= 0,		# non-negative numeric, youngest age (time before present) to consider when fitting and with respect to which rho is defined (rho(age0) is the fraction of lineages extant at age0 that are included in the tree)
											grid_sizes			= c(1,10),	# integer vector, listing the grid sizes to consider
											uniform_grid		= FALSE,	# logical, specifying whether the age grid should be uniform (equidistant age intervals). If FALSE, then the grid point density is chosen proportional to the square root of the LTT, hence resulting in higher resolution grid near the present.
											criterion			= "AIC",	# character, how to choose the optimal grid point. Options are "AIC" or "BIC".
											exhaustive			= TRUE,		# logical, whether to try all grid sizes for choosing the "best" one. If FALSE, the grid size is gradually increased until the selectin criterio (e.g., AIC) starts becoming worse, at which point the search is halted.
											min_PDR				= -Inf,		# numeric, lower bound for the fitted PDRs (applying to all grid points).
											max_PDR				= +Inf,		# numeric, upper bound for the fitted PDRs (applying to all grid points).
											min_rholambda0		= 1e-10,	# numeric, lower bound for the fitted rholambda0.
											max_rholambda0		= +Inf,		# numeric, upper bound for the fitted rholambda0.
											guess_PDR			= NULL,		# initial guess for the PDR. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant PDR at all ages), or a function handle (for generating guesses at each grid point; this function may also return NA at some time points).
											guess_rholambda0	= NULL,		# initial guess for rholambda0. Either NULL (an initial guess will be computed automatically) or a single strictly-positive numeric.
											fixed_PDR			= NULL,		# optional fixed PDR value. Either NULL (none of the PDRs are fixed), or a single scalar (all PDRs are fixed to this value) or a function handle specifying PDR for any arbitrary age (PDRs will be fixed at any age for which this function returns a finite number). The function PDR() need not return finite values for all times.
											fixed_rholambda0	= NULL,		# optional fixed value for rholambda0. If non-NULL and non-NA, then rholambda0 is not fitted. 
											splines_degree		= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by PDR and mu on the age grid.
											condition			= "auto",	# one of "crown" or "stem" or "none" or "auto" or "stem2" (or "stem3" etc) or "crown3" (or "crown4" etc), specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age!=root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
											relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
											Ntrials				= 1,
											Nbootstraps			= 0,		# (integer) optional number of parametric-bootstrap samples (random trees generated using the fitted PDR and rholambda0) for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100. Bootstrapping is only performed for the best grid size.
											Ntrials_per_bootstrap = NULL,	# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
											Nthreads			= 1,
											max_model_runtime	= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
											fit_control			= list(),
											verbose				= FALSE,
											verbose_prefix		= ""){
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input parameters..\n",verbose_prefix))
	root_age = get_tree_span(tree)$max_distance
	if(is.null(oldest_age)) oldest_age = root_age
	if(!is.null(guess_PDR)){
		if(class(guess_PDR) != "function"){
			if(length(guess_PDR)!=1){
				return(list(success=FALSE, error="Expecting either exactly one guess_PDR, or NULL, or a function handle"))
			}else{
				# convert guess_PDR to a function handle
				guess_PDR_value = guess_PDR
				guess_PDR = function(ages){ rep(guess_PDR_value, length(ages)) }
			}
		}
	}else{
		guess_PDR = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(fixed_PDR)){
		if(class(fixed_PDR) != "function"){
			if(length(fixed_PDR)!=1){
				return(list(success=FALSE, error="Expecting either exactly one fixed_PDR, or NULL, or a function handle"))
			}else{
				# convert fixed_PDR to a function handle
				fixed_PDR_value = fixed_PDR
				fixed_PDR = function(ages){ rep(fixed_PDR_value, length(ages)) }
			}
		}
	}else{
		fixed_PDR = function(ages){ rep(NA, length(ages)) }
	}
	if(length(min_PDR)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one min_PDR; instead, received %d",length(min_PDR))))
	if(length(max_PDR)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one max_PDR; instead, received %d",length(max_PDR))))
	if(!(criterion %in% c("AIC", "BIC"))) return(list(success=FALSE, error=sprintf("Invalid model selection criterion '%s'. Expected 'AIC' or 'BIC'",criterion)))
	Nmodels  = length(grid_sizes)
	
	# calculate tree LTT if needed
	if(!uniform_grid){
		LTT = count_lineages_through_time(tree=tree, Ntimes = max(100,10*max(grid_sizes)), regular_grid = TRUE, ultrametric=TRUE)
		LTT$ages = root_age - LTT$times
	}
	
	# determine order in which to examine models
	if(exhaustive){
		model_order = seq_len(Nmodels)
	}else{
		# examine models in the order of increasing grid sizes
		model_order = order(grid_sizes)
	}
	
	# fit HBD model on various grid sizes, keeping track of the "best" Ngrid
	if(verbose) cat(sprintf("%sFitting models with %s%d different grid sizes..\n",verbose_prefix,(if(exhaustive) "" else "up to "),Nmodels))
	AICs 		= rep(NA, times=Nmodels)
	BICs 		= rep(NA, times=Nmodels)
	best_fit	= NULL
	for(m in model_order){
		Ngrid = grid_sizes[m]
		if(uniform_grid || (Ngrid==1)){
			age_grid = seq(from=age0, to=oldest_age, length.out=Ngrid)
		}else{
			age_grid = get_inhomogeneous_grid_1D(Xstart = age0, Xend = oldest_age, Ngrid = Ngrid, densityX = rev(LTT$ages), densityY=sqrt(rev(LTT$lineages)), extrapolate=TRUE)
		}
		if(verbose) cat(sprintf("%s  Fitting model with grid size %d..\n",verbose_prefix,Ngrid))
		fit = fit_hbd_pdr_on_grid(	tree					= tree, 
									oldest_age				= oldest_age,
									age0					= age0,
									age_grid				= age_grid,
									min_PDR					= min_PDR,
									max_PDR					= max_PDR,
									min_rholambda0			= min_rholambda0,
									max_rholambda0			= max_rholambda0,
									guess_PDR				= guess_PDR(age_grid),
									guess_rholambda0		= guess_rholambda0,
									fixed_PDR				= fixed_PDR(age_grid),
									fixed_rholambda0		= fixed_rholambda0,
									splines_degree			= splines_degree,
									condition				= condition,
									relative_dt				= relative_dt,
									Ntrials					= Ntrials,
									Nbootstraps				= 0,
									Nthreads				= Nthreads,
									max_model_runtime		= max_model_runtime,
									fit_control				= fit_control,
									verbose					= FALSE,
									verbose_prefix			= paste0(verbose_prefix,"    "))
		if(!fit$success) return(list(success=FALSE, error=sprintf("Fitting model with grid size %d failed: %s",Ngrid,fit$error)))
		criterion_value = fit[[criterion]]
		if(is.null(best_fit)){
			best_fit = fit
			worsened = FALSE
		}else if(criterion_value<best_fit[[criterion]]){
			best_fit = fit
			worsened = FALSE
		}else{
			worsened = TRUE
		}
		AICs[m] = fit$AIC
		BICs[m] = fit$BIC
		if(verbose) cat(sprintf("%s  --> %s=%.10g. Best grid size so far: %d\n",verbose_prefix,criterion,criterion_value,length(best_fit$age_grid)))
		if((!exhaustive) && worsened) break; # model selection criterion became worse compared to the previous grid size, so stop search and keep best model found so far
	}

	if((Nbootstraps>0) && (!is.null(best_fit))){
		if(verbose) cat(sprintf("%s  Performing boostraps for best model, with grid size %d..\n",verbose_prefix,length(best_fit$age_grid)))
		best_fit = fit_hbd_pdr_on_grid(	tree					= tree, 
										oldest_age				= oldest_age,
										age0					= age0,
										age_grid				= best_fit$age_grid,
										min_PDR					= min_PDR,
										max_PDR					= max_PDR,
										min_rholambda0			= min_rholambda0,
										max_rholambda0			= max_rholambda0,
										guess_PDR				= guess_PDR(best_fit$age_grid),
										guess_rholambda0		= guess_rholambda0,
										fixed_PDR				= fixed_PDR(best_fit$age_grid),
										fixed_rholambda0		= fixed_rholambda0,
										splines_degree			= splines_degree,
										condition				= condition,
										relative_dt				= relative_dt,
										Ntrials					= Ntrials,
										Nbootstraps				= Nbootstraps,
										Ntrials_per_bootstrap	= Ntrials_per_bootstrap,
										Nthreads				= Nthreads,
										max_model_runtime		= max_model_runtime,
										fit_control				= fit_control,
										verbose					= FALSE,
										verbose_prefix			= paste0(verbose_prefix,"    "))		
	}
	
	return(list(success	 	= TRUE,
				best_fit 	= best_fit,
				grid_sizes	= grid_sizes,
				AICs		= AICs,
				BICs		= BICs))
}


