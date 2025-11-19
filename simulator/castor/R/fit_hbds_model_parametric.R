# Fit a homogenous-birth-death-sampling cladogenic model to a timetree, by estimating parameters for functional forms of lambda & mu & psi & kappa
# An HBDS model is defined by a time-dependent speciation rate (lambda), a time-dependent extinction rate (mu), a time-dependent sampling rate (psi), and a time-dependent retention rate (kappa).
# Optionally, the model may include a finite set of "concentrated sampling attempts" that occurred at specific times t_1,..,t_NCSA with specific sampling probabilities psi_1,..,psi_NCSA and specific retenion probabilities kappa_1,..,kappa_NCSA
fit_hbds_model_parametric = function(tree, 
									param_values,						# numeric vector of size NP, specifying fixed values for a some or all parameters. For fitted (i.e. non-fixed) parameters, use NaN or NA.
									param_guess				= NULL,		# numeric vector of size NP, listing an initial guess for each parameter. For fixed parameters, guess values are ignored.
									param_min				= -Inf,		# numeric vector of size NP, specifying lower bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same lower bound is assumed for all params.
									param_max				= +Inf,		# numeric vector of size NP, specifying upper bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same upper bound is assumed for all params.
									param_scale				= NULL,		# numeric vector of size NP, specifying typical scales for the model parameters. For fixed parameters, scales are ignored. If NULL, scales are automatically estimated from other information (such as provided guess and bounds). May also be a single scalar, in which case the same scale is assumed for all params.
									root_age 				= NULL,		# optional numeric, the age of the root. Can be used to define a time offset, e.g. if the last tip was not actually sampled at the present. If NULL, this will be calculated from the tree and it will be assumed that the last tip was sampled at the present
									oldest_age				= NULL,		# either a numeric specifying the oldest age to consider in the likelihood. Can be NULL, in which case this is set to the root age.
									lambda					= 0,		# function handle, mapping age & model_parameters to the current speciation rate, (age,param_values) --> lambda. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[]. Can also be a single numeric.
									mu						= 0,		# function handle, mapping age & model_parameters to the current extinction rate, (age,param_values) --> mu. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[]. Can also be a single numeric.
									psi						= 0,		# function handle, mapping age & model_parameters to the current continuous (Poissonian) sampling rate, (age,param_values) --> psi. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[]. Can also be a single numeric.
									kappa					= 0,		# function handle, mapping age & model_parameters to the current retention probability (probability of a lineage remaining in the pool upon continuous sampling), (age,param_values) --> kappa. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[]. Can also be a single numeric.
									age_grid				= NULL,		# numeric vector of size NG>=1, listing ages in ascending order, on which the lambda, mu & kappa functionals should be evaluated. This age grid must be fine enough to capture the possible variation in lambda(), mu() and kappa() over time. If NULL or of length 1, then lambda & mu & kappa are assumed to be time-independent.
									CSA_ages				= NULL,		# optional numeric vector of length NCSA, listing the ages of concentrated sampling attempts in strictly ascending order. These ages are assumed to be known, i.e. they are not fitted
									CSA_probs				= NULL,		# function handle, mapping model_parameters to NCSA concentrated sampling probabilities, (param_values) --> (psi_1,..,psi_NCSA). Must be defined for all parameters within the imposed bounds. Must always return a numeric vector of length NCSA. Can also be a single numeric. Only relevant if NCSA>0.
									CSA_kappas				= 0,		# function handle, mapping model_parameters to NCSA concentrated retention probabilities, (param_values) --> (kappa_1,..,kappa_NCSA). Must be defined for all parameters within the imposed bounds. Must always return a numeric vector of length NCSA. Can also be a single numeric. Only relevant if NCSA>0.
									condition				= "auto",	# one of "crown" or "stem" or "none" or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age!=root_age, and "crown" when oldest_age==root_age. Note that "crown" really only makes sense when oldest_age==root_age.
									ODE_relative_dt			= 0.001,	# positive unitless number, relative integration time step for the ODE solvers. Relative to the typical time scales of the dynamics, as estimated from the theoretically maximum possible rate of change. Typical values are 0.001 - 0.1.
									ODE_relative_dy			= 1e-3,		# positive unitless mumber, unitless number, relative step for interpolating simulated values over time. So a ODE_relative_dy of 0.001 means that E is recorded and interpolated between points between which E differs by roughy 0.001. Typical values are 0.01-0.0001. A smaller E_value_step increases interpolation accuracy, but also increases memory requirements and adds runtime (scales with the tree's age span, not Ntips).
									CSA_age_epsilon			= NULL,		# non-negative numeric (in units of time), age radius around a concentrated sampling attempt, within which to assume that sampling events were due to the concentrated sampling attempt. If NULL, this is chosen automatically based on the anticipated scale of numerical rounding errors. Only relevant if NCSA>0.
									Ntrials					= 1,		# integer, number of fitting trials to perform, each time starting with random parameter values
									max_start_attempts		= 1,		# integer, number of times to attempt finding a valid start point (per trial) before giving up. Randomly choosen start parameters may result in Inf/undefined objective, so this option allows the algorithm to keep looking for valid starting points.
									Nthreads				= 1,		# integer, number of parallel threads to use for running multiple fitting trials
									max_model_runtime		= NULL,		# numeric, maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
									Nbootstraps				= 0,		# integer optional number of parametric-bootstrap samples for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
									Ntrials_per_bootstrap	= NULL,		# integer optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
									fit_control				= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
									focal_param_values		= NULL,		# optional 2D numeric matrix with NP columns and an arbitrary number of rows, specifying parameter combinations of particular interest for which the loglikelihood should be calculated for. Can be used e.g. to explore the shape of the loglikelihood function.
									verbose					= FALSE,	# boolean, specifying whether to print informative messages
									diagnostics				= FALSE,	# boolean, specifying whether to print detailed info (such as log-likelihood) at every iteration of the fitting. For debugging purposes mainly.
									verbose_prefix			= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	# basic input error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	if(tree$Nnode<2) return(list(success = FALSE, error="Tree is too small"));
	if((!is.null(oldest_age)) && (!is.null(age_grid)) && (tail(age_grid,1)<oldest_age)) return(list(success=FALSE, error=sprintf("Provided age grid must cover oldest_age (%g)",oldest_age)))
	if(is.null(CSA_ages)){
		NCSA = 0
		CSA_ages = numeric(0)
	}else{
		NCSA = length(CSA_ages)
		if(any(diff(CSA_ages)<=0)) return(list(success=FALSE, error=sprintf("CSA_ages must be in strictly ascending order")))
		if(any(CSA_ages<0)) return(list(success=FALSE, error=sprintf("CSA_ages must be non-negative")))
	}
	tree_span = get_tree_span(tree)$max_distance
	if(is.null(root_age)) root_age = tree_span
	if(is.null(oldest_age)) oldest_age = root_age;
	if(!(condition %in% c("crown","stem","auto","none"))) return(list(success = FALSE, error = sprintf("Invalid condition '%s': Extected 'stem', 'crown', 'none' or 'auto'.",condition)));
	if(condition=="auto") condition = (if(abs(oldest_age-root_age)<=1e-10*root_age) "crown" else "stem")
	if(is.null(age_grid)){
		age_grid = 0
	}else{
		if(any(diff(age_grid)<0)) age_grid = sort(age_grid) # avoid common errors where age_grid is in reverse order
	}
	if(is.null(Nbootstraps) || is.na(Nbootstraps) || (Nbootstraps<0)) Nbootstraps = 0;
	if(is.null(Ntrials_per_bootstrap)) Ntrials_per_bootstrap = max(1,Ntrials)
	if(is.null(max_model_runtime)) max_model_runtime = 0
	max_start_attempts 	= (if(is.null(max_start_attempts)) 1 else max(1,max_start_attempts))
	Ntrials  = (if(is.null(Ntrials)) 1 else max(1,Ntrials))
	Nthreads = (if(is.null(Nthreads)) 1 else max(1,Nthreads))
	
	# check if some of the functionals are actually fixed numbers
	model_fixed = TRUE
	if(is.numeric(lambda) && (length(lambda)==1)){
		# the provided lambda is actually a single number, so convert to a function
		lambda_value = lambda
		lambda = function(ages,params){ rep(lambda_value, times=length(ages)) }
	}else{
		model_fixed = FALSE
	}
	if(is.numeric(mu) && (length(mu)==1)){
		# the provided mu is actually a single number, so convert to a function
		mu_value = mu
		mu = function(ages,params){ rep(mu_value, times=length(ages)) }
	}else{
		model_fixed = FALSE
	}
	if(is.numeric(kappa) && (length(kappa)==1)){
		# the provided kappa is actually a single number, so convert to a function
		kappa_value = kappa
		kappa = function(ages,params){ rep(kappa_value, times=length(ages)) }
	}else{
		model_fixed = FALSE
	}
	if(is.numeric(psi) && (length(psi)==1)){
		# the provided psi is actually a single number, so convert to a function
		psi_value = psi
		psi = function(ages,params){ rep(psi_value, times=length(ages)) }
	}else{
		model_fixed = FALSE
	}
	if((NCSA>0) && is.numeric(CSA_probs) && (length(CSA_probs)==1)){
		CSA_prob_value = CSA_probs
		CSA_probs = function(params){ rep(CSA_prob_value, times=NCSA) }
	}else{
		model_fixed = FALSE
	}
	if((NCSA>0) && is.numeric(CSA_kappas) && (length(CSA_kappas)==1)){
		CSA_kappa_value = CSA_kappas
		CSA_kappas = function(params){ rep(CSA_kappa_value, times=NCSA) }
	}else{
		model_fixed = FALSE
	}

	# pre-compute some tree stats
	if((NCSA>0) && is.null(CSA_age_epsilon)){
		CSA_age_epsilon = (if(NCSA>1) min(tree_span,min(diff(CSA_ages))) else tree_span)/1000
	}
	tree_events = extract_HBDS_events_from_tree(tree, root_age = root_age, CSA_ages = CSA_ages, age_epsilon = CSA_age_epsilon)

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
	
	if((NFP>0) && model_fixed) return(list(success=FALSE, error="At least one model parameter is fitted, however all model aspects (lambda, mu, psi, kappa etc) are fixed"))
	
	# check if functionals are valid at least on the initial guess
	lambda_guess 	= lambda(age_grid,param_guess)
	mu_guess 		= mu(age_grid,param_guess)
	psi_guess 		= psi(age_grid,param_guess)
	kappa_guess 	= kappa(age_grid,param_guess)
	if(!all(is.finite(lambda_guess))) return(list(success=FALSE, error=sprintf("lambda is not a valid number for guessed parameters, at some ages")));
	if(!all(is.finite(mu_guess))) return(list(success=FALSE, error=sprintf("mu is not a valid number for guessed parameters, at some ages")));
	if(!all(is.finite(psi_guess))) return(list(success=FALSE, error=sprintf("psi is not a valid number for guessed parameters, at some ages")));
	if(!all(is.finite(kappa_guess))) return(list(success=FALSE, error=sprintf("kappa is not a valid number for guessed parameters, at some ages")));
	if(any(kappa_guess>1)) return(list(success=FALSE, error=sprintf("kappa is greater than 1 for guessed parameters, at some ages; expected probabilities between 0 and 1")));
	if(length(lambda_guess)!=length(age_grid)) return(list(success=FALSE, error=sprintf("lambda function must return vectors of the same length as the input ages")));
	if(length(mu_guess)!=length(age_grid)) return(list(success=FALSE, error=sprintf("mu function must return vectors of the same length as the input ages")));	
	if(length(psi_guess)!=length(age_grid)) return(list(success=FALSE, error=sprintf("psi function must return vectors of the same length as the input ages")));	
	if(length(kappa_guess)!=length(age_grid)) return(list(success=FALSE, error=sprintf("kappa function must return vectors of the same length as the input ages")));
	if(NCSA>0){
		CSA_probs_guess	= CSA_probs(param_guess)
		CSA_kappas_guess	= CSA_kappas(param_guess)
		if(any(!is.finite(CSA_probs_guess))) return(list(success=FALSE, error=sprintf("CSA_probs is not a valid number for guessed parameters, at some ages")));
		if(any(!is.finite(CSA_kappas_guess))) return(list(success=FALSE, error=sprintf("CSA_kappas is not a valid number for guessed parameters, at some ages")));
		if(length(CSA_probs_guess)!=NCSA) return(list(success=FALSE, error=sprintf("CSA_probs function must return vectors of the same length as the number of concentrated sampling attempts")));
		if(length(CSA_kappas_guess)!=NCSA) return(list(success=FALSE, error=sprintf("CSA_kappas function must return vectors of the same length as the number of concentrated sampling attempts")));
		if(any(CSA_probs_guess>1)) return(list(success=FALSE, error=sprintf("CSA_probs is greater than 1 for guessed parameters, at some ages; expected probabilities between 0 and 1")));
		if(any(CSA_kappas_guess>1)) return(list(success=FALSE, error=sprintf("CSA_kappas is greater than 1 for guessed parameters, at some ages; expected probabilities between 0 and 1")));
	}
	
	# set fit-control options, unless provided by the caller
	if(is.null(fit_control)) fit_control = list()
	if(is.null(fit_control$step.min)) fit_control$step.min = 0.001
	if(is.null(fit_control$x.tol)) fit_control$x.tol = 1e-8
	if(is.null(fit_control$iter.max)) fit_control$iter.max = 1000
	if(is.null(fit_control$eval.max)) fit_control$eval.max = 2 * fit_control$iter.max * NFP


	################################
	# FITTING
			
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scale
	objective_function = function(fparam_values,trial){
		params = param_values; params[fitted_params] = fparam_values * param_scale[fitted_params];
		if(any(is.nan(params)) || any(is.infinite(params))) return(Inf); # catch weird cases where params become NaN
		if(!is.null(param_names)) names(params) = param_names;
		lambdas 	= lambda(age_grid,params)
		mus 		= mu(age_grid,params)
		psis 		= psi(age_grid,params)
		kappas 		= kappa(age_grid,params)
		if(!(all(is.finite(lambdas)) && all(is.finite(mus)) && all(is.finite(psis)) && all(is.finite(kappas)))) return(Inf); # catch weird cases where lambda/mu/psi/kappa become NaN
		if(length(age_grid)==1){
			# while age-grid has only one point (i.e., lambda & mu & psi & kappa are constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age)
			input_lambdas	= c(lambdas, lambdas)
			input_mus		= c(mus, mus)
			input_psis		= c(psis, psis)
			input_kappas	= c(kappas, kappas)
		}else{
			input_age_grid 	= age_grid
			input_lambdas	= lambdas
			input_mus 		= mus
			input_psis 		= psis
			input_kappas	= kappas
		}
		results = get_HBDS_model_loglikelihood_CPP(	branching_ages					= tree_events$branching_ages,
													Ptip_ages						= tree_events$Ptip_ages,
													Pnode_ages						= tree_events$Pnode_ages,
													CSA_ages						= (if(NCSA==0) numeric(0) else CSA_ages),
													CSA_probs						= (if(NCSA==0) numeric(0) else CSA_probs(params)),
													CSA_kappas						= (if(NCSA==0) numeric(0) else CSA_kappas(params)),
													concentrated_tip_counts			= tree_events$concentrated_tip_counts,
													concentrated_node_counts		= tree_events$concentrated_node_counts,
													oldest_age						= oldest_age,
													age_grid						= input_age_grid,
													lambdas							= input_lambdas,
													mus								= input_mus,
													psis							= input_psis,
													kappas							= input_kappas,
													splines_degree					= 1,
													condition						= condition,
													relative_ODE_step				= ODE_relative_dt,
													E_value_step					= ODE_relative_dy,
													runtime_out_seconds				= max_model_runtime)
		loglikelihood = if((!results$success) || (!is.finite(results$loglikelihood))) -Inf else results$loglikelihood
		if(diagnostics){
			if(results$success){ cat(sprintf("%s  Trial %s: loglikelihood %.10g, model runtime %.5g sec\n",verbose_prefix,as.character(trial),loglikelihood,results$runtime)) }
			else{ cat(sprintf("%s  Trial %s: Model evaluation failed: %s\n",verbose_prefix,as.character(trial),results$error)) }
		}
		return(-loglikelihood);
	}
	
	# calculate loglikelihood for initial guess
	guess_loglikelihood = -objective_function(param_guess[fitted_params]/param_scale[fitted_params], trial=0)
	
	# calculate loglikelihood for focal param values
	if((!is.null(focal_param_values)) && (nrow(focal_param_values)>0)){
		if(ncol(focal_param_values)!=NP) return(list(success=FALSE, error=sprintf("focal_param_values has %d columns, but expected exactly %d columns (=number of parameters)",ncol(focal_param_values),NP)))
		if(verbose) cat(sprintf("%sComputing loglikelihoods for focal param values..\n",verbose_prefix))
		focal_loglikelihoods = sapply(1:nrow(focal_param_values), FUN=function(r) -objective_function(focal_param_values[r,fitted_params]/param_scale[fitted_params], trial="focal"))
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
	
	# determine number of data points used
	Ndata = sum(tree_events$branching_ages<=oldest_age) + sum(tree_events$Ptip_ages<=oldest_age) + sum(tree_events$Pnode_ages<=oldest_age) + sum(tree_events$concentrated_tip_counts[CSA_ages<=oldest_age]) + sum(tree_events$concentrated_node_counts[CSA_ages<=oldest_age])



	#######################################################################
	# estimate confidence intervals if needed, via parametric bootstrapping
	
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sEstimating confidence intervals using %d parametric bootstraps..\n",verbose_prefix,Nbootstraps))
		bootstrap_params = matrix(NA,nrow=Nbootstraps,ncol=NP)
		bootstrap_LLs	 = rep(NA,times=Nbootstraps)
		for(b in 1:Nbootstraps){
			# simulate model with fitted parameters
			if(verbose) cat(sprintf("%s  Bootstrap #%d..\n",verbose_prefix,b))
			bootstrap = generate_tree_hbds(	max_time			= root_age,
											include_extant		= FALSE,
											include_extinct		= FALSE,
											time_grid			= rev(root_age-age_grid),
											lambda				= rev(lambda(root_age-age_grid,fitted_param_values)),
											mu					= rev(mu(root_age-age_grid,fitted_param_values)),
											psi					= rev(psi(root_age-age_grid,fitted_param_values)),
											kappa				= rev(kappa(root_age-age_grid,fitted_param_values)),
											splines_degree		= 1,
											CSA_times			= rev(root_age-CSA_ages),
											CSA_probs			= rev(CSA_probs(fitted_param_values)),
											CSA_kappas			= rev(CSA_kappas(fitted_param_values)))
			if(!bootstrap$success){
				if(verbose) cat(sprintf("%s    WARNING: Bootstrap %d failed: Could not simulate HBDS for the fitted params: %s\n",verbose_prefix,b,bootstrap$error))
				next
				# return(list(success=FALSE, error=sprintf("Bootstrapping failed at bootstrap %d: Could not simulate HBDS for the fitted model: %s",b,bootstrap$error), param_fitted=fitted_param_values, loglikelihood=loglikelihood, NFP=NFP));
			}else{
				if(verbose) cat(sprintf("%s    Note: Bootstrap tree has %d tips and %d nodes\n",verbose_prefix,length(bootstrap$tree$tip.label),bootstrap$tree$Nnode))
			}
			
			# fit HBDS model using bootstrap-simulation
			fit = fit_hbds_model_parametric(tree					= bootstrap$tree, 
											param_values			= param_values,
											param_guess				= param_guess,
											param_min				= param_min,
											param_max				= param_max,
											param_scale				= param_scale,
											root_age 				= bootstrap$final_time-bootstrap$root_time,
											oldest_age				= oldest_age,
											lambda					= lambda,
											mu						= mu,
											psi						= psi,
											kappa					= kappa,
											age_grid				= age_grid,
											CSA_ages				= CSA_ages,
											CSA_probs				= CSA_probs,
											CSA_kappas				= CSA_kappas,
											condition				= condition,
											ODE_relative_dt			= ODE_relative_dt,
											ODE_relative_dy			= ODE_relative_dy,
											CSA_age_epsilon			= CSA_age_epsilon,
											Ntrials					= Ntrials,
											max_start_attempts		= max_start_attempts,
											Nthreads				= Nthreads,
											max_model_runtime		= max_model_runtime,
											Nbootstraps				= 0,
											fit_control				= fit_control,
											focal_param_values		= matrix(fitted_param_values, ncol=NP, byrow=TRUE),
											verbose					= verbose,
											diagnostics				= diagnostics,
											verbose_prefix			= paste0(verbose_prefix,"    "))
			if(!fit$success){
				if(verbose) cat(sprintf("%s  WARNING: Fitting failed for this bootstrap: %s\n",verbose_prefix,fit$error))
				next
			}
			bootstrap_params[b,] = fit$param_fitted
			bootstrap_LLs[b]	 = fit$focal_loglikelihoods
		}
		# calculate standard errors and confidence intervals from distribution of bootstrapped parameters
		standard_errors = sqrt(pmax(0, colMeans(bootstrap_params^2, na.rm=TRUE) - colMeans(bootstrap_params, na.rm=TRUE)^2))
		quantiles 	= sapply(1:ncol(bootstrap_params), FUN=function(p) quantile(bootstrap_params[,p], probs=c(0.25, 0.75, 0.025, 0.975, 0.5), na.rm=TRUE, type=8))
		CI50lower 	= quantiles[1,]
		CI50upper 	= quantiles[2,]
		CI95lower 	= quantiles[3,]
		CI95upper 	= quantiles[4,]
		medians 	= quantiles[5,]
		mean_BLL	= mean(bootstrap_LLs,na.rm=TRUE)
		consistency = sum(abs(bootstrap_LLs-mean_BLL)>=abs(loglikelihood-mean_BLL),na.rm=TRUE)/sum(!is.nan(bootstrap_LLs))
	}


	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				param_fitted			= fitted_param_values,
				param_guess				= param_guess,
				guess_loglikelihood		= guess_loglikelihood,
				focal_loglikelihoods	= focal_loglikelihoods,
				NFP						= NFP,
				Ndata					= Ndata,
				AIC						= 2*NFP - 2*loglikelihood,
				BIC						= log(Ndata)*NFP - 2*loglikelihood,
				condition				= condition,
				converged				= fits[[best]]$converged,
				Niterations				= fits[[best]]$Niterations,
				Nevaluations			= fits[[best]]$Nevaluations,
				trial_start_objectives	= -unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$start_objective)),
				trial_objective_values	= -objective_values,
				trial_Nstart_attempts	= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Nstart_attempts)),
				trial_Niterations		= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Niterations)),
				trial_Nevaluations		= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Nevaluations)),
				standard_errors			= (if(Nbootstraps>0) setNames(standard_errors, param_names) else NULL),
				medians					= (if(Nbootstraps>0) setNames(medians, param_names) else NULL),
				CI50lower				= (if(Nbootstraps>0) setNames(CI50lower, param_names) else NULL),
				CI50upper				= (if(Nbootstraps>0) setNames(CI50upper, param_names) else NULL),
				CI95lower				= (if(Nbootstraps>0) setNames(CI95lower, param_names) else NULL),
				CI95upper				= (if(Nbootstraps>0) setNames(CI95upper, param_names) else NULL),
				consistency				= (if(Nbootstraps>0) consistency else NULL)))
}


