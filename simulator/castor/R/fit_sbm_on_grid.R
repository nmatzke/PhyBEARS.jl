# Fit a Spherical Brownian Motion (SBM) model with a diffusivity that varies polynomially between multiple grid points.
#
fit_sbm_on_grid = function(	tree, 
							tip_latitudes, 						# numeric vector of size Ntips, listing geographical latitudes of the tips (in decimal degrees)
							tip_longitudes, 					# numeric vector of size Ntips, listing geographical longitudes of the tips (in decimal degrees)
							radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
							clade_states			= NULL,		# optional, either an integer vector of length Ntips+Nnodes (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), specifying the discrete "state" of each tip and node in each tree. This can be used to limit independent contrasts to tip pairs whose total number of state-transitions (along their shortest path) is zero.
							planar_approximation	= FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
							only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
							only_distant_tip_pairs	= FALSE,	# logical, whether to only consider tip pairs located at distinct geographic locations
							min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
							max_MRCA_age			= Inf,		# numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set max_MRCA_age=Inf to disable this filter.
							max_phylodistance		= Inf,		# numeric, maximum allowed geodistance for an independent contrast to be included in the SBM fitting
							no_state_transitions	= FALSE,	# if TRUE, only tip pairs without state transitions along their shortest paths are considered. In particular, only tips in the same state are considered. Requires that clade_states[] is provided.
							only_state				= NULL,		# optional integer, specifying the state in which tip pairs (and their connecting ancestors) must be in order to be considered. Requires that clade_states[] is provided.
							time_grid				= 0,		# numeric vector of length NG, specifying the time points on which to fit the diffusivity. Between time points the diffusivity is assumed to vary polynomially (e.g., linearly or quadratically, see option splines_degree). If NULL or a single number, then the diffusivity is assumed to be time-independent. If the grid does not cover the entire interval [0,root_age], constant extrapolation is used beyond the time_grid's limits.
							guess_diffusivity		= NULL,		# optional first guess for the diffusivity. Either NULL (choose automatically), or a single numeric (same guess for all grid points) or a numeric vector of length NG. May also contain NAs (first guess will be automatically chosen at those grid points).
							min_diffusivity			= NULL,		# optional lower bound for the fitted diffusivity. Either NULL (choose automatically), or a single numeric (same lower bound at all grid points) or a numeric vector of length NG. Note that negative values are replaced by zeros. May also contain NAs (lower bound will be automatically chosen at those grid points).
							max_diffusivity			= Inf,		# optional upper bound for the fitted diffusivity. Either NULL (no upper bound), or a single numeric (same upper bound at all grid points) or a numeric vector of length NG. To disable a bound set it to Inf.
							Ntrials					= 1,		# number of fitting trials to perform, each time starting with random parameter values
							Nthreads				= 1,
							Nbootstraps				= 0,		# (integer) optional number of parametric-bootstrap samples for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
							Ntrials_per_bootstrap	= NULL,		# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
							NQQ						= 0,		# (integer) optional number of simulations to perform for creating Q-Q plots of the theoretically expected distribution of geodistances vs the empirical distribution of geodistances (across independent contrasts). The resolution of the returned QQ plot will be equal to the number of independent contrasts used for fitting.
							fit_control				= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
							SBM_PD_functor			= NULL,		# internally used SBM probability density functor
							verbose					= FALSE,	# boolean, specifying whether to print informative messages
							verbose_prefix			= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	Ntips = length(tree$tip.label)
	# basic input error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	if(tree$Nnode<2) return(list(success = FALSE, error="Input tree is too small"));
	root_age	= get_tree_span(tree)$max_distance
	Ntrials		= max(1,Ntrials)
	Nthreads	= max(1,Nthreads)
	if(min_MRCA_time<0) min_MRCA_time = Inf
	Ntrials = pmax(1,Ntrials)
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

	if(is.null(time_grid)){
		if((!is.null(min_diffusivity)) && (length(min_diffusivity)>1)) return(list(success = FALSE, error = sprintf("Invalid length of min_diffusivity[]; since no time grid was provided, you must provide a single (constant) min_diffusivity or none at all")));
		if((!is.null(max_diffusivity)) && (length(max_diffusivity)>1)) return(list(success = FALSE, error = sprintf("Invalid length of max_diffusivity[]; since no time grid was provided, you must provide a single (constant) max_diffusivity or none at all")));
		if((!is.null(guess_diffusivity)) && (length(guess_diffusivity)>1)) return(list(success = FALSE, error = sprintf("Invalid length of guess_diffusivity[]; since no time grid was provided, you must provide a single (constant) guess_diffusivity or none at all")));
		time_grid = 0 # single-point grid, means that diffusivity is assumed to be time-independent
		NG = 1
	}else{
		NG = length(time_grid)
		if(any(diff(time_grid)<=0)) return(list(success=FALSE, error="Time grid must be in strictly ascending order"))
		if((!is.null(guess_diffusivity)) && (length(guess_diffusivity)!=1) && (length(guess_diffusivity)!=NG)) return(list(success = FALSE, error = sprintf("Invalid length of guess_diffusivity[] (%d); since a time grid of size %d was provided, you must either provide one or %d guess_diffusivities",length(guess_diffusivity),NG,NG)));
		if((!is.null(min_diffusivity)) && (length(min_diffusivity)!=1) && (length(min_diffusivity)!=NG)) return(list(success = FALSE, error = sprintf("Invalid length of min_diffusivity[] (%d); since a time grid of size %d was provided, you must either provide one or %d min_diffusivities",length(min_diffusivity),NG,NG)));
		if((!is.null(max_diffusivity)) && (length(max_diffusivity)!=1) && (length(max_diffusivity)!=NG)) return(list(success = FALSE, error = sprintf("Invalid length of min_diffusivity[] (%d); since a time grid of size %d was provided, you must either provide one or %d min_diffusivities",length(min_diffusivity),NG,NG)));
	}
	
	# pre-calculate SBM probability density functor for efficiency
	if(is.null(SBM_PD_functor)){
		if(verbose) cat(sprintf("%sPre-computing SBM probability density functor..\n",verbose_prefix))
		SBM_PD_functor = SBM_get_SBM_PD_functor_CPP(max_error = 1e-8, max_Legendre_terms = 200)
	}

					
	##############################################
	# determine first guess and lower/upper bounds

	if(is.null(guess_diffusivity) || any(!is.finite(guess_diffusivity))){
		# fit time-independent SBM to get a rough estimate
		if(verbose) cat(sprintf("%sFitting const-diffusivity model for initual guess..\n",verbose_prefix))
		fit_const = fit_sbm_const(	tree,
									tip_latitudes		= tip_latitudes,
									tip_longitudes		= tip_longitudes,
									radius				= radius,
									clade_states		= clade_states,
									only_basal_tip_pairs= only_basal_tip_pairs,
									min_MRCA_time		= min_MRCA_time,
									max_MRCA_age		= max_MRCA_age,
									max_phylodistance	= max_phylodistance,
									no_state_transitions= no_state_transitions,
									only_state			= only_state,
									min_diffusivity		= (if(is.null(min_diffusivity) || all(!is.finite(min_diffusivity))) NULL else min(min_diffusivity,na.rm=TRUE)),
									max_diffusivity		= (if(is.null(max_diffusivity) || all(!is.finite(max_diffusivity))) NULL else max(max_diffusivity,na.rm=TRUE)),
									Nbootstraps			= 0,
									SBM_PD_functor		= SBM_PD_functor)
		if(!fit_const$success) return(list(success=FALSE, error=sprintf("Failed to fit constant-diffusivity model: %s",fit_const$error)))
		if(is.null(guess_diffusivity)){
			guess_diffusivity = rep(fit_const$diffusivity, NG)
		}else{
			if(length(guess_diffusivity)==1) guess_diffusivity = rep(guess_diffusivity, NG)
			guess_diffusivity[!is.finite(guess_diffusivity)] = fit_const$success
		}
	}else if(length(guess_diffusivity)==1){
		guess_diffusivity = rep(guess_diffusivity, NG)
	}
	if(is.null(min_diffusivity)){
		min_diffusivity = 0.0000001*guess_diffusivity
	}else{
		if(length(min_diffusivity)==1) min_diffusivity = rep(min_diffusivity, NG)
		min_diffusivity[!(is.finite(min_diffusivity))] = 0.0000001*guess_diffusivity[!(is.finite(min_diffusivity))]
		min_diffusivity[min_diffusivity<0] = 0
	}
	if(is.null(max_diffusivity)){
		max_diffusivity = rep(Inf,NG)
	}else{
		if(length(max_diffusivity)==1) max_diffusivity = rep(max_diffusivity, NG)
		max_diffusivity[!(is.finite(max_diffusivity))] = Inf
	}
	diffusivity_scale = abs(guess_diffusivity)
	diffusivity_scale[diffusivity_scale==0] = mean(diffusivity_scale)

	####################################
	# Fit diffusivity on the full grid
	
	
	if(verbose) cat(sprintf("%sFitting diffusivity at %d grid points..\n",verbose_prefix,NG))
	diffusivity_functor = function(times, params){
		if(NG==1){
			return(params[[1]])
		}else{
			return(approx(x=time_grid, y=unlist(params), xout=times, method="linear", rule=2)$y)
		}
	}
	fit = fit_sbm_parametric(	tree,
								tip_latitudes			= tip_latitudes,
								tip_longitudes			= tip_longitudes,
								radius 					= radius,
								param_values 			= rep(NA,NG),
								param_guess 			= guess_diffusivity,
								diffusivity 			= diffusivity_functor,
								time_grid 				= time_grid,
								clade_states			= clade_states,
								planar_approximation	= planar_approximation,
								only_basal_tip_pairs	= only_basal_tip_pairs,
								only_distant_tip_pairs	= only_distant_tip_pairs,
								min_MRCA_time			= min_MRCA_time,
								max_MRCA_age			= max_MRCA_age,
								max_phylodistance		= max_phylodistance,
								no_state_transitions	= no_state_transitions,
								only_state				= only_state,
								param_min				= min_diffusivity,
								param_max				= max_diffusivity,
								param_scale				= diffusivity_scale,
								Ntrials 				= Ntrials,
								Nthreads				= Nthreads,
								Nbootstraps				= Nbootstraps,
								Ntrials_per_bootstrap	= Ntrials_per_bootstrap,
								NQQ						= NQQ,
								fit_control				= fit_control,
								SBM_PD_functor			= SBM_PD_functor,
								verbose					= verbose,
								verbose_prefix			= sprintf("%s  ",verbose_prefix))
	if(!fit$success) return(list(success=FALSE, error=sprintf("Failed to fit SBM model: %s",fit$error)))
	
	
	####################################
		
	# return results
	return(list(success					= TRUE,
				objective_value			= fit$objective_value,
				objective_name			= fit$objective_name,
				time_grid				= time_grid,
				diffusivity				= fit$param_fitted,
				loglikelihood			= fit$loglikelihood,
				NFP						= fit$NFP,
				Ncontrasts				= fit$Ncontrasts,
				phylodistances			= fit$phylodistances,
				geodistances			= fit$geodistances,
				child_times1			= fit$child_times1,
				child_times2			= fit$child_times2,
				MRCA_times				= fit$MRCA_times,
				AIC						= fit$AIC,
				BIC						= fit$BIC,
				converged				= fit$converged,
				Niterations				= fit$Niterations,
				Nevaluations			= fit$Nevaluations,
				trial_start_objectives	= fit$trial_start_objectives,
				trial_objective_values	= fit$trial_objective_values,
				trial_Nstart_attempts	= fit$trial_Nstart_attempts,
				trial_Niterations		= fit$trial_Niterations,
				trial_Nevaluations		= fit$trial_Nevaluations,
				standard_errors			= fit$standard_errors,
				CI50lower				= fit$CI50lower,
				CI50upper				= fit$CI50upper,
				CI95lower				= fit$CI95lower,
				CI95upper				= fit$CI95upper,
				consistency				= fit$consistency,
				QQplot					= fit$QQplot,
				SBM_PD_functor			= SBM_PD_functor))
}



