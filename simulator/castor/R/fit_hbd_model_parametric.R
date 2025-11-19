# Fit a homogenous-birth-death cladogenic model to an ultrametric timetree, by estimating parameters for functional forms of lambda & mu
# An HBD model is defined by a time-dependent speciation rate (lambda), a time-dependent extinction rate (mu) and a rarefaction (rho0, sampling fraction)
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
fit_hbd_model_parametric = function(tree, 
									param_values,					# numeric vector of size NP, specifying fixed values for a some or all parameters. For fitted (i.e. non-fixed) parameters, use NaN or NA.
									param_guess				= NULL,		# numeric vector of size NP, listing an initial guess for each parameter. For fixed parameters, guess values are ignored.
									param_min				= -Inf,		# numeric vector of size NP, specifying lower bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same lower bound is assumed for all params.
									param_max				= +Inf,		# numeric vector of size NP, specifying upper bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same upper bound is assumed for all params.
									param_scale				= NULL,		# numeric vector of size NP, specifying typical scales for the model parameters. For fixed parameters, scales are ignored. If NULL, scales are automatically estimated from other information (such as provided guess and bounds). May also be a single scalar, in which case the same scale is assumed for all params.
									oldest_age				= NULL,		# either a numeric specifying the stem age or NULL (equivalent to the root age). This is similar to the "tot_time" option in the R function RPANDA::likelihood_bd
									age0					= 0,		# non-negative numeric, youngest age (time before present) to consider when fitting and with respect to which rho is defined (rho0:=rho(age0) is the fraction of lineages extant at age0 that are included in the tree)
									lambda,								# function handle, mapping age & model_parameters to the current speciation rate, (age,param_values) --> lambda. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[].
									mu						= 0,		# function handle, mapping age & model_parameters to the current extinction rate, (age,param_values) --> mu. Must be defined for all ages in [0:oldest_age] and for all parameters within the imposed bounds. Must be vectorized in the age argument, i.e. return a vector the same size as age[]. Can also be a single numeric.
									rho0					= 1,		# function handle, mapping model_parameters to the sampling fraction at age0 (aka. rarefaction), (param_values) --> rho. Must be defined for all parameters within the imposed bounds. Can also be a single numeric.
									age_grid				= NULL,		# numeric vector of size NG>=1, listing ages in ascending order, on which the lambda and mu functionals should be evaluated. This age grid must be fine enough to capture the possible variation in lambda() and mu() over time. If NULL or of length 1, then lambda & mu are assumed to be time-independent.
									condition				= "auto",	# one of "crown" or "stem" or "none" or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age>root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
									relative_dt				= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
									Ntrials					= 1,		# number of fitting trials to perform, each time starting with random parameter values
									max_start_attempts		= 1,		# number of times to attempt finding a valid start point (per trial) before giving up. Randomly choosen start parameters may result in Inf/undefined objective, so this option allows the algorithm to to keep looking for valid starting points.
									Nthreads				= 1,
									max_model_runtime		= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
									Nbootstraps				= 0,		# integer optional number of parametric-bootstrap samples for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
									Ntrials_per_bootstrap	= NULL,		# integer optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
									fit_algorithm 			= "nlminb",	# either "nlminb" or "subplex". What algorithm to use for fitting.
									fit_control				= list(),	# a named list containing options for the nlminb or subplex fitting routine (e.g. iter.max and rel.tol)
									focal_param_values		= NULL,		# optional 2D numeric matrix with NP columns and an arbitrary number of rows, specifying parameter combinations of particular interest for which the loglikelihood should be calculated for. Can be used e.g. to explore the shape of the loglikelihood function.
									verbose					= FALSE,	# boolean, specifying whether to print informative messages
									diagnostics				= FALSE,	# boolean, specifying whether to print detailed info (such as log-likelihood) at every iteration of the fitting. For debugging purposes mainly.
									verbose_prefix			= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	# basic input error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	if(!(fit_algorithm %in% c("subplex", "nlminb"))) return(list(success=FALSE, error=sprintf("ERROR: Invalid optimization algorithm '%s'",fit_algorithm)))
	if(tree$Nnode<2) return(list(success = FALSE, error="Tree is too small"));
	if(age0<0) return(list(success = FALSE, error="age0 must be non-negative"));
	if(is.null(Ntrials_per_bootstrap)) Ntrials_per_bootstrap = max(1,Ntrials)
	if((!is.null(oldest_age)) && (!is.null(age_grid)) && (tail(age_grid,1)<oldest_age)) return(list(success=FALSE, error=sprintf("Provided age grid must cover oldest_age (%g)",oldest_age)))
	if(is.null(max_model_runtime)) max_model_runtime = 0;
	if(is.null(age_grid)){
		age_grid = 0
	}else{
		if(any(diff(age_grid)<0)) age_grid = sort(age_grid) # avoid common errors where age_grid is in reverse order
		if(age_grid[1]>age0) return(list(success=FALSE, error=sprintf("Provided age grid must cover age0 (%g)",age0)))
	}
	max_start_attempts 	= max(1,max_start_attempts)
	Ntrials 			= max(1,Ntrials)
	Nthreads 			= max(1,Nthreads)
	original_Ntips 		= length(tree$tip.label)

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
	if(is.numeric(rho0) && (length(rho0)==1)){
		# the provided rho0 is actually a single number, so convert to a function
		rho0_value = rho0
		rho0 = function(params){ rho0_value }
	}else{
		model_fixed = FALSE
	}

	# trim tree at age0 if needed, while shifting time for the subsequent analyses (i.e. new ages will start counting at age0)
	if(age0>0){
		root_age = get_tree_span(tree)$max_distance
		if(root_age<age0) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the root age (%g)",age0,root_age)));
		if((!is.null(oldest_age)) && (oldest_age<age0)) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the oldest considered age (%g)",age0,oldest_age)));
		tree = trim_tree_at_height(tree,height=root_age-age0)$tree
		if(tree$Nnode<2) return(list(success = FALSE, error=sprintf("Tree is too small after trimming at age0 (%g)",age0)));
		if(!is.null(oldest_age)) oldest_age	= oldest_age - age0	
		if(!is.null(age_grid)) age_grid 	= age_grid - age0
		root_age = root_age - age0
	}

	# pre-compute some tree stats
	sorted_node_ages = sort(get_all_branching_ages(tree));
	root_age 		 = tail(sorted_node_ages,1);
	age_epsilon		 = 1e-4*mean(tree$edge.length);

	# more input error checking
	if(is.null(oldest_age)) oldest_age = root_age;
	if((!(condition %in% c("crown","stem","auto"))) && (!startsWith(condition,"stem")) && (!startsWith(condition,"crown"))) return(list(success = FALSE, error = sprintf("Invalid condition '%s': Expected 'stem', 'stem2', 'stem<N>', 'crown', 'crown<N>', or 'auto'.",condition)));
	if(condition=="auto") condition = (if(abs(oldest_age-root_age)<=1e-10*root_age) "crown" else "stem")
	if((!is.null(oldest_age)) && (!is.null(age_grid)) && (tail(age_grid,1)<oldest_age)) return(list(success=FALSE, error=sprintf("Provided age grid must cover oldest_age (%g)",oldest_age)))
	
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
	
	if((NFP>0) && model_fixed) return(list(success=FALSE, error="At least one model parameter is fitted, however all model aspects (lambda, mu, rho, kappa etc) are fixed"))
	
	# check if functionals are valid at least on the initial guess
	lambda_guess 	= lambda(age_grid+age0,param_guess)
	mu_guess 		= mu(age_grid+age0,param_guess)
	rho0_guess 		= rho0(param_guess)
	if(!all(is.finite(lambda_guess))) return(list(success=FALSE, error=sprintf("lambda is not a valid number for guessed parameters, at some ages")));
	if(!all(is.finite(mu_guess))) return(list(success=FALSE, error=sprintf("mu is not a valid number for guessed parameters, at some ages")));
	if(!is.finite(rho0_guess)) return(list(success=FALSE, error=sprintf("rho0 is not a valid number for guessed parameters")));
	if(length(lambda_guess)!=length(age_grid)) return(list(success=FALSE, error=sprintf("lambda function must return vectors of the same length as the input ages")));
	if(length(mu_guess)!=length(age_grid)) return(list(success=FALSE, error=sprintf("mu function must return vectors of the same length as the input ages")));
	
	# set fit-control options, unless provided by the caller
	if(is.null(fit_control)) fit_control = list()
	if(fit_algorithm=="nlminb"){
		if(is.null(fit_control$step.min)) fit_control$step.min 	= 0.001
		if(is.null(fit_control$x.tol)) fit_control$x.tol 		= (if(is.null(fit_control$xtol_abs)) 1e-10 else fit_control$xtol_abs)
		if(is.null(fit_control$rel.tol)) fit_control$rel.tol 	= (if(is.null(fit_control$ftol_rel)) 1e-10 else fit_control$ftol_rel)
		if(is.null(fit_control$iter.max)) fit_control$iter.max 	= 1000
		if(is.null(fit_control$eval.max)) fit_control$eval.max 	= (if(is.null(fit_control$maxeval)) 2 * fit_control$iter.max * NFP else fit_control$maxeval)
		# remove recycled subplex-specific options, because otherwise nlminb complains
		fit_control$maxeval 	= NULL
		fit_control$ftol_rel 	= NULL
		fit_control$xtol_abs	= NULL
	}else if(fit_algorithm=="subplex"){
		if(is.null(fit_control$maxeval)) fit_control$maxeval 	= (if(is.null(fit_control$eval.max)) 2 * fit_control$iter.max * NFP else fit_control$eval.max)
		if(is.null(fit_control$ftol_rel)) fit_control$ftol_rel 	= (if(is.null(fit_control$rel.tol)) 1e-10 else fit_control$rel.tol)
		if(is.null(fit_control$ftol_abs)) fit_control$ftol_abs 	= (if(is.null(fit_control$abs.tol)) 0 else fit_control$abs.tol)
		if(is.null(fit_control$xtol_abs)) fit_control$xtol_abs 	= (if(is.null(fit_control$x.tol)) 0 else fit_control$x.tol)
		if(is.null(fit_control$xtol_rel)) fit_control$xtol_rel 	= 1e-10
		if(is.null(fit_control$algorithm)) fit_control$algorithm= "NLOPT_LN_SBPLX"
	}
	
	sprint_param_values = function(params){
		if(is.null(names(params))){
			pnames = str(seq_len(length(params)))
		}else{
			pnames = names(params)
		}
		return(paste(sapply(seq_len(length(params)), FUN=function(p) sprintf("%s: %.6g",pnames[p],params[p])), collapse=", "))
	}
	
	################################
	# FITTING
	
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scale
	objective_function = function(fparam_values,trial){
		params = param_values; params[fitted_params] = fparam_values * param_scale[fitted_params];
		if(any(is.nan(params)) || any(is.infinite(params))) return(if(fit_algorithm == "optim") 1e100 else Inf); # catch weird cases where params become NaN
		if(!is.null(param_names)) names(params) = param_names;
		lambdas 	= lambda(age_grid+age0,params)
		mus 		= mu(age_grid+age0,params)
		input_rho0 	= rho0(params)
		if(!(all(is.finite(lambdas)) && all(is.finite(mus)) && is.finite(input_rho0))) return(Inf); # catch weird cases where lambda/mu/rho become NaN
		if(length(age_grid)==1){
			# while age-grid has only one point (i.e., lambda & mu are constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age);
			input_lambdas	= c(lambdas, lambdas);
			input_mus		= c(mus, mus);
		}else{
			input_age_grid 	= age_grid;
			input_lambdas	= lambdas
			input_mus 		= mus
		}
		results = HBD_model_loglikelihood_CPP(	branching_ages		= sorted_node_ages,
												oldest_age			= oldest_age,
												rarefaction			= input_rho0,
												age_grid 			= input_age_grid,
												lambdas 			= input_lambdas,
												mus 				= input_mus,
												splines_degree		= 1,
												condition			= condition,
												relative_dt			= relative_dt,
												runtime_out_seconds	= max_model_runtime);
		loglikelihood = if((!results$success) || (!is.finite(results$loglikelihood))) (if(fit_algorithm == "optim") -1e100 else -Inf) else results$loglikelihood
		if(diagnostics){
			if(results$success){ cat(sprintf("%s  Trial %s: loglikelihood %.10g, model runtime %.5g sec\n%s    Parameters: %s\n",verbose_prefix,as.character(trial),loglikelihood,results$runtime,verbose_prefix,sprint_param_values(params))) }
			else{ cat(sprintf("%s  Trial %s: Model evaluation failed: %s\n",verbose_prefix,as.character(trial),results$error)) }
		}
		return(-loglikelihood)
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
			start_values = param_guess[fitted_params]
			if(trial>1){
				boxed_left	= which((!is.infinite(lower_bounds)) & is.infinite(upper_bounds))
				boxed_right	= which((!is.infinite(upper_bounds)) & is.infinite(lower_bounds))
				boxed_dual  = which(!(is.infinite(lower_bounds) | is.infinite(upper_bounds))); # determine fitted params that are boxed, i.e. constrained to within finite lower & upper bounds
				unboxed 	= which(is.infinite(lower_bounds) & is.infinite(upper_bounds))
				if(length(boxed_dual)>0) 	start_values[boxed_dual] = lower_bounds[boxed_dual] + (upper_bounds[boxed_dual]-lower_bounds[boxed_dual]) * runif(n=length(boxed_dual),min=0,max=1)
				if(length(unboxed)>0) 	 	start_values[unboxed]	 = 10**runif(n=length(unboxed), min=-2, max=2) * start_values[unboxed]
				if(length(boxed_left)>0) 	start_values[boxed_left] = sapply(boxed_left, FUN=function(fp) random_semiboxed_left(lower_bound=lower_bounds[fp], default=start_values[fp], typical_scale=scales[fp], orders_of_magnitude=4))
				if(length(boxed_right)>0) 	start_values[boxed_right]= sapply(boxed_right, FUN=function(fp) -random_semiboxed_left(lower_bound=-upper_bounds[fp], default=-start_values[fp], typical_scale=scales[fp], orders_of_magnitude=4))
			}
			# make sure start fparams are within bounds
			start_values 	= pmax(lower_bounds,pmin(upper_bounds,start_values))
			start_objective = objective_function(start_values/scales, trial)
			Nstart_attempts = Nstart_attempts + 1
			if(is.finite(start_objective)) break;
		}
		# run fit
		if(is.finite(start_objective)){
			if(fit_algorithm == "nlminb"){
				fit = tryCatch({ stats::nlminb(start_values/scales, 
									objective	= function(pars){ objective_function(pars, trial) },
									lower		= lower_bounds/scales, 
									upper		= upper_bounds/scales, 
									control		= fit_control)
								}, error = function(e){ list(objective=NaN, par=NA, convergence=1, evaluations=NA, iterations=NA) })
				LL 				= -fit$objective;
				Nevaluations 	= fit$evaluations[[1]]
				Niterations		= fit$iterations
				converged		= (fit$convergence==0)
			}else if(fit_algorithm == "subplex"){
				fit = tryCatch({ nloptr::nloptr(x0 		= start_values/scales, 
												eval_f 	= function(pars){ objective_function(pars, trial) },
												lb 		= lower_bounds/scales, 
												ub 		= upper_bounds/scales, 
												opts 	= fit_control)
								}, error = function(e){ list(objective=NaN, solution=NA, status=1, iterations=NA) })
				LL 				= -fit$objective
				Nevaluations 	= NA
				Niterations 	= fit$iterations
				converged 		= (fit$status==0)
				fit$par 		= fit$solution
			}
			if(is.null(LL)){ LL = NaN; converged = FALSE; }
			results = list(objective_value=-LL, fparam_values = fit$par*scales, converged=converged, Niterations=Niterations, Nevaluations=Nevaluations, Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective)
			if(diagnostics) cat(sprintf("%s  Trial %d: Final loglikelihood %.10g, Niterations %d, Nevaluations %d, converged = %d\n",verbose_prefix,trial,LL,Niterations,Nevaluations,converged))
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
		
	# reverse any time shift due to earlier tree trimming
	age_grid 	= age_grid + age0
	root_age 	= root_age + age0
	oldest_age	= oldest_age + age0
	
	
	#######################################################################
	# estimate confidence intervals if needed, via parametric bootstrapping
	
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sEstimating confidence intervals using %d parametric bootstraps..\n",verbose_prefix,Nbootstraps))
		#   include a dummy age grid point at the end of age_grid if needed to cover the root age
		#   also include a dummy age grid point at the beginning if necessary, to cover present-day (age 0)
		sim_age_grid = age_grid
		if(tail(sim_age_grid,1)<root_age) sim_age_grid = c(sim_age_grid, root_age*1.01)
		if(sim_age_grid[1]>0) sim_age_grid = c(0,sim_age_grid)
		# calculate PSR on sim_age_grid
		sim_PSR = get_PSR_of_HBD_model(	oldest_age 		= tail(sim_age_grid,1),
										age_grid		= sim_age_grid,
										lambda			= lambda(sim_age_grid,fitted_param_values),
										mu				= mu(sim_age_grid,fitted_param_values),
										age0			= age0,
										rho0			= rho0(fitted_param_values),
										splines_degree	= 1,
										relative_dt		= relative_dt)
		if(!sim_PSR$success) return(list(success=FALSE, error=sprintf("Bootstrapping failed: Could not calculate PSR of fitted model: %s",sim_PSR$error)))
		bootstrap_params = matrix(NA,nrow=Nbootstraps,ncol=NP)
		bootstrap_LLs	 = rep(NA,times=Nbootstraps)
		for(b in 1:Nbootstraps){
			# simulate HBD model with fitted parameters
			if(verbose) cat(sprintf("%s  Bootstrap #%d..\n",verbose_prefix,b))
			bootstrap = castor::generate_tree_hbd_reverse(	Ntips			= original_Ntips,
															crown_age		= root_age,
															age_grid		= sim_PSR$ages,
															PSR				= sim_PSR$PSR,
															splines_degree	= 1,
															relative_dt		= relative_dt)
			if(!bootstrap$success){
				if(verbose) cat(sprintf("%s    WARNING: Bootstrap %d failed: Could not simulate HBD for the fitted params: %s\n",verbose_prefix,b,bootstrap$error))
				next
			}
			bootstrap_tree = bootstrap$trees[[1]]
			if(verbose) cat(sprintf("%s    Note: Bootstrap tree has %d tips and %d nodes\n",verbose_prefix,length(bootstrap_tree$tip.label),bootstrap_tree$Nnode))
			# fit HBD model using bootstrap-simulation
			fit = fit_hbd_model_parametric(	tree					= bootstrap_tree, 
											param_values			= param_values,
											param_guess				= param_guess,
											param_min				= param_min,
											param_max				= param_max,
											param_scale				= param_scale,
											oldest_age				= oldest_age,
											age0					= age0,
											lambda					= lambda,
											mu						= mu,
											rho0					= rho0,
											age_grid				= age_grid,
											condition				= condition,
											relative_dt				= relative_dt,
											Ntrials					= Ntrials_per_bootstrap,
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
				AIC						= 2*NFP - 2*loglikelihood,
				BIC						= log(sum((sorted_node_ages<=oldest_age) & (sorted_node_ages>=age0)))*NFP - 2*loglikelihood,
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


