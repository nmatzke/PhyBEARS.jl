# Fit a homogenous-birth-death cladogenic model-congruence-class to an ultrametric timetree, by estimating functional form parameters for the pulled speciation rate (PSR)
# An HBD congruence class is fully defined by a time-dependent pulled speciation rate (PSR)
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
fit_hbd_psr_parametric = function(	tree, 
									param_values,					# numeric vector of size NP, specifying fixed values for a some or all parameters. For fitted (i.e. non-fixed) parameters, use NaN or NA.
									param_guess			= NULL,		# numeric vector of size NP, listing an initial guess for each parameter. For fixed parameters, guess values are ignored.
									param_min			= -Inf,		# numeric vector of size NP, specifying lower bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same lower bound is assumed for all params.
									param_max			= +Inf,		# numeric vector of size NP, specifying upper bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same upper bound is assumed for all params.
									param_scale			= NULL,		# numeric vector of size NP, specifying typical scales for the model parameters. For fixed parameters, scales are ignored. If NULL, scales are automatically estimated from other information (such as provided guess and bounds). May also be a single scalar, in which case the same scale is assumed for all params.
									oldest_age			= NULL,		# either a numeric specifying the stem age or NULL (equivalent to the root age). This is similar to the "tot_time" option in the R function RPANDA::likelihood_bd
									age0				= 0,		# non-negative numeric, youngest age (time before present) to consider when fitting and with respect to which rholambda0 is defined (i.e. rholambda0 = rho(age0)*lambda(age0))
									PSR,							# function handle, mapping age & model_parameters to the current pulled speciation rate, (age,param_values) --> PSR. Must be defined for all ages in [age0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[].
									age_grid			= NULL,		# numeric vector of size NG>=1, listing ages in ascending order, on which the PSR functional should be evaluated. This age grid must be fine enough to capture the possible variation in PSR() over time. If NULL or of length 1, then the PSR is assumed to be time-independent.
									condition			= "auto",	# one of "crown" or "stem" or "stem2" (or "stem3", "stem4", .. etc) or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age<root_age, "stem2" when oldest_age>root_age, or "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
									relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
									Ntrials				= 1,		# number of fitting trials to perform, each time starting with random parameter values
									max_start_attempts	= 1,		# number of times to attempt finding a valid start point (per trial) before giving up. Randomly choosen start parameters may result in Inf/undefined objective, so this option allows the algorithm to keep looking for valid starting points.
									Nthreads			= 1,
									max_model_runtime	= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
									fit_control			= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
									verbose				= FALSE,	# boolean, specifying whether to print informative messages
									diagnostics			= FALSE,	# boolean, specifying whether to print detailed info (such as log-likelihood) at every iteration of the fitting. For debugging purposes mainly.
									verbose_prefix		= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	# basic input error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	if(tree$Nnode<1) return(list(success = FALSE, error="Input tree is too small"));
	if(age0<0) return(list(success = FALSE, error="age0 must be non-negative"));
	root_age = get_tree_span(tree)$max_distance
	if(is.null(oldest_age)) oldest_age = root_age;
	if(root_age<age0) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the root age (%g)",age0,root_age)));
	if(oldest_age<age0) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the oldest considered age (%g)",age0,oldest_age)));
	if((!(condition %in% c("crown","stem","auto"))) && (!startsWith(condition,"stem")) && (!startsWith(condition,"crown"))) return(list(success = FALSE, error = sprintf("Invalid condition '%s': Expected 'stem', 'stem2', 'stem<N>', 'crown', 'crown<N>', or 'auto'.",condition)));
	if(condition=="auto") condition = (if(abs(oldest_age-root_age)<=1e-10*root_age) "crown" else (if(oldest_age>root_age) "stem2" else "stem"))
	max_start_attempts 	= max(1,max_start_attempts)
	Ntrials 			= max(1,Ntrials)
	Nthreads 			= max(1,Nthreads)
	if(!is.null(age_grid)){
		if(any(diff(age_grid)<0)) age_grid = sort(age_grid) # avoid common errors where age_grid is in reverse order
		if((age_grid[1]>age0) || (tail(age_grid,1)<oldest_age)) return(list(success=FALSE, error=sprintf("Provided age_grid must cover the full considered age interval (age0=%.10g to oldest_age=%.10g)",age0,oldest_age)))
	}
	if(Ntrials<1) return(list(success = FALSE, error = sprintf("Ntrials must be at least 1")))
	if(is.null(age_grid)) age_grid = 0;
	if(is.null(max_model_runtime)) max_model_runtime = 0;
	
	# check if some of the functionals are actually fixed numbers
	model_fixed = TRUE
	if(is.numeric(PSR) && (length(PSR)==1)){
		# the provided PSR is actually a single number, so convert to a function
		PSR_value = PSR
		PSR = function(ages,params){ rep(PSR_value, times=length(ages)) }
	}else{
		model_fixed = FALSE
	}

	# trim tree at age0 if needed, while shifting time for the subsequent analyses (i.e. new ages will start counting at age0)
	if(age0>0){
		tree = trim_tree_at_height(tree,height=root_age-age0)$tree
		if(tree$Nnode<1) return(list(success = FALSE, error=sprintf("Tree is too small after trimming at age0 (%g)",age0)));
		if(!is.null(oldest_age)) oldest_age	= oldest_age - age0	
		if(!is.null(age_grid)) age_grid 	= age_grid - age0
		root_age = root_age - age0
	}
	
	# pre-compute some tree stats
	lineage_counter  = count_lineages_through_time(tree, Ntimes=max(3,log2(length(tree$tip.label))), include_slopes=TRUE, ultrametric=TRUE)
	sorted_node_ages = sort(get_all_branching_ages(tree));
	root_age 		 = tail(sorted_node_ages,1);
	age_epsilon		 = 1e-4*mean(tree$edge.length);
		
	# sanitize model parameters
	sanitized_params = sanitize_parameters_for_fitting(param_values, param_guess = param_guess, param_min = param_min, param_max = param_max, param_scale = param_scale)
	if(!sanitized_params$success) return(list(success=FALSE, error=sanitized_params$error))
	NP				= sanitized_params$NP
	NFP				= sanitized_params$NFP
	param_names		= sanitized_params$param_names
	param_values	= sanitized_params$param_values
	param_guess		= sanitized_params$param_guess
	param_min		= sanitized_params$param_min
	param_max		= sanitized_params$param_max
	param_scale		= sanitized_params$param_scale
	fitted_params	= sanitized_params$fitted_params
	fixed_params	= sanitized_params$fixed_params
	
	if((NFP>0) && model_fixed) return(list(success=FALSE, error="At least one model parameter is fitted, however all model aspects (PSR) are fixed"))
	
	# check if functionals are valid at least on the initial guess
	PSR_guess = PSR(age_grid+age0,param_guess)
	if(!all(is.finite(PSR_guess))) return(list(success=FALSE, error=sprintf("PSR is not a valid number for guessed parameters, at some ages")));
	if(length(PSR_guess)!=length(age_grid)) return(list(success=FALSE, error=sprintf("PSR function must return vectors of the same length as the input ages")));

	# set fit-control options, unless provided by the caller
	if(is.null(fit_control)) fit_control = list()
	if(is.null(fit_control$step.min)) fit_control$step.min 	= 0.001
	if(is.null(fit_control$x.tol)) fit_control$x.tol 		= 1e-10
	if(is.null(fit_control$rel.tol)) fit_control$rel.tol 	= 1e-10
	if(is.null(fit_control$iter.max)) fit_control$iter.max 	= 1000
	if(is.null(fit_control$eval.max)) fit_control$eval.max 	= 2 * fit_control$iter.max * NFP
		

	################################
	# FITTING
	
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scale
	objective_function = function(fparam_values,trial){
		params = param_values; params[fitted_params] = fparam_values * param_scale[fitted_params];
		if(any(is.nan(params)) || any(is.infinite(params))) return(Inf); # catch weird cases where params become NaN
		if(!is.null(param_names)) names(params) = param_names;
		PSRs = PSR(age_grid+age0,params)
		if(!all(is.finite(PSRs))) return(Inf); # catch weird cases where PSR becomes NaN
		if(length(age_grid)==1){
			# while age-grid has only one point (i.e., PSR is constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age);
			input_PSRs		= c(PSRs, PSRs);
		}else{
			input_age_grid 	= age_grid;
			input_PSRs		= PSRs
		}
		results = HBD_PSR_loglikelihood_CPP(branching_ages		= sorted_node_ages,
											oldest_age			= oldest_age,
											age_grid 			= input_age_grid,
											PSRs 				= input_PSRs,
											splines_degree		= 1,
											condition			= condition,
											relative_dt			= relative_dt,
											runtime_out_seconds	= max_model_runtime);
		loglikelihood = if((!results$success) || (!is.finite(results$loglikelihood))) -Inf else results$loglikelihood
		if(diagnostics){
			if(results$success){ cat(sprintf("%s  Trial %s: loglikelihood %.10g, model runtime %.5g sec\n",verbose_prefix,as.character(trial),loglikelihood,results$runtime)) }
			else{ cat(sprintf("%s  Trial %s: Model evaluation failed: %s\n",verbose_prefix,as.character(trial),results$error)) }
		}
		return(-loglikelihood)
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
			start_objective = objective_function(start_values/scales, trial)
			Nstart_attempts = Nstart_attempts + 1
			if(is.finite(start_objective)) break;
		}
		# run fit with chosen start_values
		if(is.finite(start_objective)){
			fit = stats::nlminb(start		= start_values/scales, 
								objective	= function(pars){ objective_function(pars, trial) },
								lower		= lower_bounds/scales, 
								upper		= upper_bounds/scales, 
								control		= fit_control)
			results = list(objective_value=fit$objective, fparam_values = fit$par*scales, converged=(fit$convergence==0), Niterations=fit$iterations, Nevaluations=fit$evaluations[[1]], Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective)
			if(diagnostics) cat(sprintf("%s  Trial %d: Final loglikelihood %.10g, Niterations %d, Nevaluations %d, converged = %d\n",verbose_prefix,trial,-results$objective_value,results$Niterations, results$Nevaluations, results$converged))
		}else{
			results = list(objective_value=NA, fparam_values = NA, converged=FALSE, Niterations=0, Nevaluations=0, Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective)
			if(diagnostics) cat(sprintf("%s  Trial %d: Start objective is non-finite. Skipping trial\n",verbose_prefix,trial))
		}
		return(results)
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
	objective_values	= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$objective_value))
	valids				= which((!is.na(objective_values)) & (!is.nan(objective_values)) & (!is.null(objective_values)) & (!is.infinite(objective_values)));
	if(length(valids)==0) return(list(success=FALSE, error=sprintf("Fitting failed for all trials")));
	best 				= valids[which.min(sapply(valids, function(i) objective_values[i]))]
	objective_value		= -fits[[best]]$objective_value;
	loglikelihood		= objective_value;
	fitted_param_values = param_values; fitted_param_values[fitted_params] = fits[[best]]$fparam_values;
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
	if(!is.null(param_names)) names(fitted_param_values) = param_names;
		
	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				param_fitted			= fitted_param_values,
				param_guess				= param_guess,
				NFP						= NFP,
				AIC						= 2*NFP - 2*loglikelihood,
				BIC						= log(sum((sorted_node_ages<=oldest_age) & (sorted_node_ages>=age0)))*NFP - 2*loglikelihood,
				converged				= fits[[best]]$converged,
				Niterations				= fits[[best]]$Niterations,
				Nevaluations			= fits[[best]]$Nevaluations,
				trial_start_objectives	= -unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$start_objective)),
				trial_objective_values	= -objective_values,
				trial_Nstart_attempts	= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Nstart_attempts)),
				trial_Niterations		= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Niterations)),
				trial_Nevaluations		= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Nevaluations))));
}


