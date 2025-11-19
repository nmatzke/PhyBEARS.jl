# Fit a Spherical Brownian Motion (SBM) model with time-dependent diffusivity, by fitting the parameters of a functional form for the diffusivity
#
fit_sbm_parametric = function(	tree, 
								tip_latitudes, 						# numeric vector of size Ntips, listing geographical latitudes of the tips (in decimal degrees)
								tip_longitudes, 					# numeric vector of size Ntips, listing geographical longitudes of the tips (in decimal degrees)
								radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
								param_values,						# numeric vector of size NP, specifying fixed values for a some or all parameters. For fitted (i.e. non-fixed) parameters, use NaN or NA.
								param_guess,						# numeric vector of size NP, listing an initial guess for each parameter. For fixed parameters, guess values are ignored.
								diffusivity,						# function handle, mapping time & model_parameters to the current diffusivity, (time,param_values) --> diffusivity. Must be defined for all times in [0:root_age] and for all parameters within the imposed bounds. Must be vectorized in the time argument, i.e. return a vector the same size as time[].
								time_grid				= NULL,		# numeric vector of size NG>=1, listing times in ascending order, on which the diffusivity functional should be evaluated. This time grid must be fine enough to capture the possible variation in diffusivity() over time. If NULL or of length 1, then the diffusivity is assumed to be time-independent. This grid should cover the interval [0,root_age]; otherwise the diffusivity will be extrapolated as a constant where needed.
								clade_states			= NULL,		# optional, either an integer vector of length Ntips+Nnodes (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), specifying the discrete "state" of each tip and node in each tree. This can be used to limit independent contrasts to tip pairs whose total number of state-transitions (along their shortest path) is zero.
								planar_approximation	= FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
								only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
								only_distant_tip_pairs	= FALSE,	# logical, whether to only consider tip pairs located at distinct geographic locations
								min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
								max_MRCA_age			= Inf,		# numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set max_MRCA_age=Inf to disable this filter.
								max_phylodistance		= Inf,		# numeric, maximum allowed geodistance for an independent contrast to be included in the SBM fitting
								no_state_transitions	= FALSE,	# if TRUE, only tip pairs without state transitions along their shortest paths are considered. In particular, only tips in the same state are considered. Requires that clade_states[] is provided.
								only_state				= NULL,		# optional integer, specifying the state in which tip pairs (and their connecting ancestors) must be in order to be considered. Requires that clade_states[] is provided.
								param_min				= -Inf,		# numeric vector of size NP, specifying lower bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same lower bound is assumed for all params.
								param_max				= +Inf,		# numeric vector of size NP, specifying upper bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same upper bound is assumed for all params.
								param_scale				= NULL,		# numeric vector of size NP, specifying typical scales for the model parameters. For fixed parameters, scales are ignored. If NULL, scales are automatically estimated from other information (such as provided guess and bounds). May also be a single scalar, in which case the same scale is assumed for all params.
								Ntrials					= 1,		# number of fitting trials to perform, each time starting with random parameter values
								max_start_attempts		= 1,		# number of times to attempt finding a valid start point (per trial) before giving up. Randomly choosen start parameters may result in Inf/undefined objective, so this option allows the algorithm to keep looking for valid starting points.
								Nthreads				= 1,
								Nbootstraps				= 0,		# (integer) optional number of parametric-bootstrap samples for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
								Ntrials_per_bootstrap	= NULL,		# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
								NQQ						= 0,		# (integer) optional number of simulations to perform for creating Q-Q plots of the theoretically expected distribution of geodistances vs the empirical distribution of geodistances (across independent contrasts). The resolution of the returned QQ plot will be equal to the number of independent contrasts used for fitting.
								fit_control				= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
								SBM_PD_functor			= NULL,		# internally used SBM probability density functor
								focal_param_values		= NULL,		# optional 2D numeric matrix with NP columns and an arbitrary number of rows, specifying parameter combinations of particular interest for which the loglikelihood should be calculated for. Can be used e.g. to explore the shape of the loglikelihood function.
								verbose					= FALSE,	# boolean, specifying whether to print informative messages
								verbose_prefix			= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	Ntips = length(tree$tip.label)
	# basic input error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	if(tree$Nnode<2) return(list(success = FALSE, error="Input tree is too small"));
	root_age 			= get_tree_span(tree)$max_distance
	max_start_attempts 	= max(1,max_start_attempts)
	Ntrials 			= max(1,Ntrials)
	Nthreads 			= max(1,Nthreads)
	if((!is.null(time_grid)) && (time_grid[1]>tail(time_grid,1))) time_grid = rev(time_grid); # avoid common errors where time_grid is in reverse order
	Ntrials = pmax(1,Ntrials)
	if(is.null(time_grid)) time_grid = 0;
	if(("list" %in% class(tip_latitudes)) && (length(tip_latitudes)==Ntips)){
		tip_latitudes = unlist(tip_latitudes)
	}
	if(("list" %in% class(tip_latitudes)) && (length(tip_longitudes)==Ntips)){
		tip_longitudes = unlist(tip_longitudes)
	}
	if((!is.null(clade_states)) && ("list" %in% class(clade_states)) && (length(clade_states)==Ntips+tree$Nnode)){
		clade_states = unlist(clade_states)
	}
	if((!is.null(only_state)) && is.null(clade_states)) return(list(success=FALSE, error="Missing clade_states[], needed when only_state is specified"))
	if(no_state_transitions && is.null(clade_states)) return(list(success=FALSE, error="Missing clade_states[], needed when no_state_transitions=TRUE"))
	if(is.null(Nbootstraps) || is.na(Nbootstraps) || (Nbootstraps<0)) Nbootstraps = 0;
	
	# sanitize input params (reformat if needed)
	sanitized_params = sanitize_parameters_for_fitting(param_values, param_guess = param_guess, param_min = param_min, param_max = param_max, param_scale = param_scale)
	if(!sanitized_params$success) return(list(success=FALSE, error=sanitized_params$error))
	param_names		= sanitized_params$param_names
	param_values 	= sanitized_params$param_values
	param_guess 	= sanitized_params$param_guess
	param_min 		= sanitized_params$param_min
	param_max 		= sanitized_params$param_max
	param_scale 	= sanitized_params$param_scale
	fitted_params	= sanitized_params$fitted_params
	fixed_params	= sanitized_params$fixed_params
	NFP				= sanitized_params$NFP
	NP				= sanitized_params$NP
	
	# check if diffusivity functional is valid at least on the initial guess
	diffusivity_guess = diffusivity(time_grid,param_guess)
	if(!all(is.finite(diffusivity_guess))) return(list(success=FALSE, error=sprintf("Diffusivity is not a valid number for guessed parameters, at some times")));
	if(any(diffusivity_guess<0)) return(list(success=FALSE, error=sprintf("Diffusivity is negative for guessed parameters, at some times")));
	if(length(diffusivity_guess)!=length(time_grid)) return(list(success=FALSE, error=sprintf("The diffusivity function must return vectors of the same length as the input times")));


	###########################################
	# EXTRACT INDEPENDENT CONTRASTS FOR FITTING
	
	if(verbose) cat(sprintf("%sExtracting independent contrasts from tree..\n",verbose_prefix))
	ICs = get_SBM_independent_contrasts(tree					= tree,
										tip_latitudes			= tip_latitudes,
										tip_longitudes			= tip_longitudes,
										radius					= radius,
										clade_states			= clade_states,
										planar_approximation	= planar_approximation,
										only_basal_tip_pairs	= only_basal_tip_pairs,
										only_distant_tip_pairs	= only_distant_tip_pairs,
										min_MRCA_time			= min_MRCA_time,
										max_MRCA_age			= max_MRCA_age,
										max_phylodistance		= max_phylodistance,
										no_state_transitions	= no_state_transitions,
										only_state				= only_state)
	if(!ICs$success) return(list(success=FALSE, error=ICs$error))
	NC = ICs$NC

	################################
	# FITTING
	
	# pre-calculate SBM probability density functor for efficiency
	if(is.null(SBM_PD_functor)){
		if(verbose) cat(sprintf("%sPre-computing SBM probability density functor..\n",verbose_prefix))
		SBM_PD_functor = SBM_get_SBM_PD_functor_CPP(max_error = 1e-8, max_Legendre_terms = 200)
	}
	
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scale
	objective_function = function(fparam_values){
		params = param_values; params[fitted_params] = fparam_values * param_scale[fitted_params];
		if(any(is.nan(params)) || any(is.infinite(params))) return(Inf); # catch weird cases where params become NaN
		if(!is.null(param_names)) names(params) = param_names;
		diffusivities = diffusivity(time_grid,params)
		if(!(all(is.finite(diffusivities)))) return(Inf); # catch weird cases where diffusivity/rholambda0 become NaN
		if(length(time_grid)==1){
			# while age-grid has only one point (i.e., diffusivity is constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,root_age]
			input_time_grid 	= c(0,root_age);
			input_diffusivities	= c(diffusivities, diffusivities);
		}else{
			input_time_grid 	= time_grid;
			input_diffusivities	= diffusivities
		}
		results = TSBM_LL_of_transitions_CPP(	radius			= radius,
												MRCA_times		= ICs$MRCA_times,
												child_times1	= ICs$child_times1,
												child_times2	= ICs$child_times2,
												distances		= ICs$geodistances,
												time_grid		= input_time_grid,
												diffusivities	= input_diffusivities,
												splines_degree	= 1,
												SBM_PD_functor	= SBM_PD_functor);
		if(!results$success) return(Inf);
		LL = results$loglikelihood;
		if(is.na(LL) || is.nan(LL) || is.infinite(LL)) return(Inf);
		return(-LL);
	}
	
	# calculate loglikelihood for initial guess
	guess_loglikelihood = -objective_function(param_guess[fitted_params]/param_scale[fitted_params])
	
	# calculate loglikelihood for focal param values
	if((!is.null(focal_param_values)) && (nrow(focal_param_values)>0)){
		if(ncol(focal_param_values)!=NP) return(list(success=FALSE, error=sprintf("focal_param_values has %d columns, but expected exactly %d columns (=number of parameters)",ncol(focal_param_values),NP)))
		if(verbose) cat(sprintf("%sComputing loglikelihoods for focal param values..\n",verbose_prefix))
		focal_loglikelihoods = sapply(1:nrow(focal_param_values), FUN=function(r) -objective_function(focal_param_values[r,fitted_params]/param_scale[fitted_params]))
	}else{
		focal_loglikelihoods = NULL
	}

	# fit with various starting points
	fit_single_trial = function(trial){
		scales		 = param_scale[fitted_params]
		lower_bounds = param_min[fitted_params]
		upper_bounds = param_max[fitted_params]
		# randomly choose start values for fitted params (keep trying up to max_start_attempts times)
		Nstart_attempts = 0
		while(Nstart_attempts<max_start_attempts){
			# randomly choose start values for fitted params
			if(trial==1){
				start_values = param_guess[fitted_params]		
			}else{
				start_values = get_random_params(defaults=param_guess[fitted_params], lower_bounds=lower_bounds, upper_bounds=upper_bounds, scales=scales, orders_of_magnitude=4)
			}
			# check if start values yield NaN
			start_objective = objective_function(start_values/scales);
			Nstart_attempts = Nstart_attempts + 1
			if(is.finite(start_objective)) break;
		}
		# run fit
		if(is.finite(start_objective)){
			fit = stats::nlminb(start_values/scales, 
								objective	= objective_function, 
								lower		= lower_bounds/scales, 
								upper		= upper_bounds/scales, 
								control		= fit_control)
			return(list(objective_value=fit$objective, fparam_values = fit$par*scales, converged=(fit$convergence==0), Niterations=fit$iterations, Nevaluations=fit$evaluations[[1]], Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective));
		}else{
			return(list(objective_value=NA, fparam_values = NA, converged=FALSE, Niterations=0, Nevaluations=0, Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective));
		}
	}
	
	################################

	# run one or more independent fitting trials
    if((Ntrials>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
		# run trials in parallel using multiple forks
		# Note: Forks (and hence shared memory) are not available on Windows
		if(verbose) cat(sprintf("%sFitting %d model parameters (%d trials, parallelized)..\n",verbose_prefix,NFP,Ntrials))
		fits = parallel::mclapply(	1:Ntrials, 
									FUN = function(trial) fit_single_trial(trial), 
									mc.cores = min(Nthreads, Ntrials), 
									mc.preschedule = FALSE, 
									mc.cleanup = TRUE);
	}else{
		# run in serial mode
		if(verbose) cat(sprintf("%sFitting %d model parameters (%s)..\n",verbose_prefix,NFP,(if(Ntrials==1) "1 trial" else sprintf("%d trials",Ntrials))))
		fits = sapply(1:Ntrials,function(x) NULL)
		for(trial in 1:Ntrials){
			fits[[trial]] = fit_single_trial(trial)
		}
	}
	

	# extract information from best fit (note that some fits may have LL=NaN or NA)
	objective_values	= sapply(1:Ntrials, function(trial) fits[[trial]]$objective_value);
	valids				= which((!is.na(objective_values)) & (!is.nan(objective_values)) & (!is.null(objective_values)) & (!is.infinite(objective_values)));
	if(length(valids)==0) return(list(success=FALSE, error=sprintf("Fitting failed for all trials")));
	best 				= valids[which.min(sapply(valids, function(i) objective_values[i]))]
	objective_value		= -fits[[best]]$objective_value;
	loglikelihood		= objective_value;
	fitted_param_values = param_values; fitted_param_values[fitted_params] = fits[[best]]$fparam_values;
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
	if(!is.null(param_names)) names(fitted_param_values) = param_names;
	
	#######################################################################
	# estimate confidence intervals if needed, via parametric bootstrapping
	
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sEstimating confidence intervals using %d parametric bootstraps..\n",verbose_prefix,Nbootstraps))
		if(is.null(Ntrials_per_bootstrap)) Ntrials_per_bootstrap = max(1,Ntrials)
		bootstrap_params = matrix(NA,nrow=Nbootstraps,ncol=NP)
		bootstrap_LLs	 = rep(NA,times=Nbootstraps)
		for(b in 1:Nbootstraps){
			# simulate model with fitted parameters
			if(verbose) cat(sprintf("%s  Bootstrap #%d..\n",verbose_prefix,b))
			bootstrap = castor::simulate_sbm(	tree			= tree, 
												radius 			= radius,
												diffusivity 	= diffusivity(time_grid, fitted_param_values),
												time_grid		= time_grid,
												splines_degree	= 1,
												root_latitude	= NULL,
												root_longitude	= NULL)
			if(!bootstrap$success) return(list(success=FALSE, error=sprintf("Bootstrapping failed at bootstrap %d: Could not simulate SBM for the fitted model: %s",b,bootstrap$error), param_fitted=fitted_param_values, loglikelihood=loglikelihood, NFP=NFP, Ncontrasts=NC));

			# fit diffusivity using bootstrap-simulation
			fit = fit_sbm_parametric(	tree, 
										tip_latitudes			= bootstrap$tip_latitudes,
										tip_longitudes			= bootstrap$tip_longitudes,
										radius					= radius,
										param_values			= param_values,
										diffusivity				= diffusivity,
										time_grid				= time_grid,
										clade_states			= clade_states,
										planar_approximation	= planar_approximation,
										only_basal_tip_pairs	= only_basal_tip_pairs,
										only_distant_tip_pairs	= only_distant_tip_pairs,
										min_MRCA_time			= min_MRCA_time,
										max_MRCA_age			= max_MRCA_age,
										max_phylodistance		= max_phylodistance,
										no_state_transitions	= no_state_transitions,
										param_guess				= param_guess,
										param_min				= param_min,
										param_max				= param_max,
										param_scale				= param_scale,
										Ntrials					= Ntrials_per_bootstrap,
										max_start_attempts		= max_start_attempts,
										Nthreads				= Nthreads,
										Nbootstraps				= 0,
										fit_control				= fit_control,
										SBM_PD_functor			= SBM_PD_functor,
										focal_param_values		= matrix(fitted_param_values, ncol=NP, byrow=TRUE),
										verbose					= verbose,
										verbose_prefix			= paste0(verbose_prefix,"    "))
			if(!fit$success){
				if(verbose) cat(sprintf("%s  WARNING: Fitting failed for this bootstrap: %s\n",verbose_prefix,fit$error))
			}else{
				bootstrap_params[b,] = fit$param_fitted
				bootstrap_LLs[b]	 = fit$focal_loglikelihoods
			}
		}
		# calculate standard errors and confidence intervals from distribution of bootstrapped parameters
		standard_errors = sqrt(pmax(0, colMeans(bootstrap_params^2, na.rm=TRUE) - colMeans(bootstrap_params, na.rm=TRUE)^2))
		quantiles 	= sapply(1:ncol(bootstrap_params), FUN=function(p) quantile(bootstrap_params[,p], probs=c(0.25, 0.75, 0.025, 0.975, 0.5), na.rm=TRUE, type=8))
		CI50lower 	= quantiles[1,]
		CI50upper 	= quantiles[2,]
		CI95lower 	= quantiles[3,]
		CI95upper 	= quantiles[4,]
		medians		= quantiles[5,]
		mean_BLL	= mean(bootstrap_LLs,na.rm=TRUE)
		consistency = sum(abs(bootstrap_LLs-mean_BLL)>=abs(loglikelihood-mean_BLL),na.rm=TRUE)/sum(!is.nan(bootstrap_LLs))
	}
	
	#####################################
	# Calculate QQ-plot using simulations
	
	if(NQQ>0){
		if(verbose) cat(sprintf("%sCalculating QQ-plot using %d simulations..\n",verbose_prefix,NQQ))
		sim_geodistances = numeric(NQQ * NC)
		for(q in 1:NQQ){
			sim = castor::simulate_sbm(	tree = tree, radius = radius, diffusivity = diffusivity(time_grid, fitted_param_values), time_grid = time_grid, splines_degree = 1, root_latitude = NULL, root_longitude = NULL)
			if(!sim$success) return(list(success=FALSE, error=sprintf("Calculation of QQ failed at simulation %d: Could not simulate SBM for the fitted model: %s",q,sim$error), param_fitted=fitted_param_values, loglikelihood=loglikelihood, NFP=NFP, Ncontrasts=NC));
			sim_geodistances[(q-1)*NC + c(1:NC)] = radius * geodesic_angles(sim$tip_latitudes[ICs$tip_pairs[,1]],sim$tip_longitudes[ICs$tip_pairs[,1]],sim$tip_latitudes[ICs$tip_pairs[,2]],sim$tip_longitudes[ICs$tip_pairs[,2]])
		}
		probs  = c(1:NC)/NC
		QQplot = cbind(quantile(ICs$geodistances, probs=probs, na.rm=TRUE, type=8), quantile(sim_geodistances, probs=probs, na.rm=TRUE, type=8))
	}
	
	####################################
		
	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				param_fitted			= fitted_param_values,
				loglikelihood			= loglikelihood,
				NFP						= NFP,
				Ncontrasts				= NC,
				phylodistances			= ICs$phylodistances,
				geodistances			= ICs$geodistances,
				child_times1			= ICs$child_times1,
				child_times2			= ICs$child_times2,
				MRCA_times				= ICs$MRCA_times,
				AIC						= 2*NFP - 2*loglikelihood,
				BIC						= log(NC)*NFP - 2*loglikelihood,
				converged				= fits[[best]]$converged,
				Niterations				= fits[[best]]$Niterations,
				Nevaluations			= fits[[best]]$Nevaluations,
				guess_loglikelihood		= guess_loglikelihood,
				focal_loglikelihoods	= focal_loglikelihoods,
				trial_start_objectives	= -sapply(1:Ntrials, function(trial) fits[[trial]]$start_objective),
				trial_objective_values	= -objective_values,
				trial_Nstart_attempts	= sapply(1:Ntrials, function(trial) fits[[trial]]$Nstart_attempts),
				trial_Niterations		= sapply(1:Ntrials, function(trial) fits[[trial]]$Niterations),
				trial_Nevaluations		= sapply(1:Ntrials, function(trial) fits[[trial]]$Nevaluations),
				standard_errors			= (if(Nbootstraps>0) setNames(standard_errors, param_names) else NULL),
				medians					= (if(Nbootstraps>0) setNames(medians, param_names) else NULL),
				CI50lower				= (if(Nbootstraps>0) setNames(CI50lower, param_names) else NULL),
				CI50upper				= (if(Nbootstraps>0) setNames(CI50upper, param_names) else NULL),
				CI95lower				= (if(Nbootstraps>0) setNames(CI95lower, param_names) else NULL),
				CI95upper				= (if(Nbootstraps>0) setNames(CI95upper, param_names) else NULL),
				consistency				= (if(Nbootstraps>0) consistency else NULL),
				QQplot					= (if(NQQ>0) QQplot else NULL),
				SBM_PD_functor			= SBM_PD_functor))

}



