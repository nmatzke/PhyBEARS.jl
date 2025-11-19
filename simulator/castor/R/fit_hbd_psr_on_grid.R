# Fit a homogenous-birth-death cladogenic model-congruence-class to an ultrametric timetree, by fitting the pulled speciation rate (PSR)
# An HBD model is defined by a time-dependent speciation rate (lambda), a time-dependent extinction rate (mu) and a rarefaction (rho, subsampling fraction)
# However, for each specific model and a given timetree there exists a continuum of alternative models that would all generate the same deterministic lineages-through-time (LTT) curve (when calculated backward in time), and all of these models actually have the same likelihood.
# Hence, each model is part of an "equivalence class" of models, and likelihood-based approaches can only discern between model classes, but not between the individual model members in a class
# It turns out that each HBD model-class is uniquely defined by its "pulled speciation rate" (PSR), which is also equal to the relative slope of the deterministic LTT
# This function thus fits model-classes, rather than models, by fitting the PSR on a discrete time grid.
fit_hbd_psr_on_grid = function(	tree, 
								oldest_age			= NULL,		# either a numeric specifying the stem age or NULL (equivalent to the root age). This is similar to the "tot_time" option in the R function RPANDA::likelihood_bd
								age0				= 0,		# non-negative numeric, youngest age (time before present) to consider when fitting
								age_grid			= NULL,		# either NULL, or a numeric vector of size NG, listing ages in ascending order, on which the PSR is defined as a piecewise linear curve. If NULL, the PSR is assumed to be time-independent.
								min_PSR				= 0,		# optional lower bound for the fitted PSRs. Either a single non-negative numeric (applying to all age-grid-points) or a non-negative numeric vector of size NG, specifying the lower bound at each age-grid point.
								max_PSR				= +Inf,		# optional upper bound for the fitted PSRs. Either a single non-negative numeric (applying to all age-grid-points) or a non-negative numeric vector of size NG, specifying the upper bound at each age-grid point.
								guess_PSR			= NULL,		# initial guess for the PSR. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant PSR at all ages) or a numeric vector of size NG specifying an initial guess for the PSR at each age-grid point (can include NAs)
								fixed_PSR			= NULL,		# optional fixed PSR values, on one or more of the age grid points. Either NULL (none of the PSRs are fixed), or a single scalar (all PSRs are fixed) or a numeric vector of size NG (some or all PSRs are fixed, can include NAs).
								splines_degree		= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by the PSR on the age grid.
								condition			= "auto",	# one of "crown" or "crown3" (or "crown4", .. etc) or "stem" or "stem2" (or "stem3", .. etc) or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age<root_age, "stem2" when oldest_age>root_age, or "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
								relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
								Ntrials					= 1,
								Nbootstraps				= 0,		# (integer) optional number of parametric-bootstrap samples (random trees generated using the fitted PSR) for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
								Ntrials_per_bootstrap	= NULL,		# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
								Nthreads				= 1,
								max_model_runtime		= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
								fit_control				= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
								verbose					= FALSE,		# boolean, specifying whether to print informative messages
								diagnostics				= FALSE,		# boolean, specifying whether to print detailed info (such as log-likelihood) at every iteration of the fitting. For debugging purposes mainly.
								verbose_prefix			= ""){			# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	original_Ntips = length(tree$tip.label)
	if((Nbootstraps>0) && (!is.null(Ntrials_per_bootstrap)) && (Ntrials_per_bootstrap<=0)) stop(sprintf("ERROR: Ntrials_per_bootstrap must be strictly positive, if bootstrapping is requested"))
	if(tree$Nnode<2) return(list(success = FALSE, error="Input tree is too small"));
	if(age0<0) return(list(success = FALSE, error="age0 must be non-negative"));
	root_age = get_tree_span(tree)$max_distance
	if(is.null(oldest_age)) oldest_age = root_age;
	if(root_age<age0) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the root age (%g)",age0,root_age)));
	if(oldest_age<age0) return(list(success=FALSE, error=sprintf("age0 (%g) is older than the oldest considered age (%g)",age0,oldest_age)));
	if((!is.null(age_grid)) && (length(age_grid)>1) && ((age_grid[1]>age0) || (tail(age_grid,1)<oldest_age))) return(list(success = FALSE, error=sprintf("Provided age-grid range (%g - %g) does not cover entire required age range (%g - %g)",age_grid[1],tail(age_grid,1),age0,oldest_age)));
	if((!(condition %in% c("crown","stem","auto"))) && (!startsWith(condition,"stem")) && (!startsWith(condition,"crown"))) return(list(success = FALSE, error = sprintf("Invalid condition '%s': Expected 'stem', 'stem2', 'stem<N>', 'crown', 'crown<N>', or 'auto'.",condition)));
	if(condition=="auto") condition = (if(abs(oldest_age-root_age)<=1e-10*root_age) "crown" else (if(oldest_age>root_age) "stem2" else "stem"))

	# trim tree at age0 if needed, while shifting time for the subsequent analyses (i.e. new ages will start counting at age0)
	if(age0>0){
		if(verbose) cat(sprintf("%sTrimming tree at age0=%g..\n",verbose_prefix,age0))
		tree = trim_tree_at_height(tree,height=root_age-age0)$tree
		if(tree$Nnode<2) return(list(success = FALSE, error=sprintf("Tree is too small after trimming at age0 (%g)",age0)));
		if(!is.null(oldest_age)) oldest_age	= oldest_age - age0	
		if(!is.null(age_grid)) age_grid 	= age_grid - age0
		root_age = root_age - age0
	}
								
	# pre-compute some tree stats
	if(verbose) cat(sprintf("%sPrecomputing some stats about the tree..\n",verbose_prefix))
	LTT0				= length(tree$tip.label);
	lineage_counter 	= count_lineages_through_time(tree, Ntimes=log2(LTT0), include_slopes=TRUE, ultrametric=TRUE)
	sorted_node_ages	= sort(get_all_branching_ages(tree));
	root_age 		 	= tail(sorted_node_ages,1);
	age_epsilon		 	= 1e-4*mean(tree$edge.length);

	# more error checking
	if(Ntrials<1) return(list(success = FALSE, error = sprintf("Ntrials must be at least 1")))
	if(is.null(oldest_age)) oldest_age = root_age;
	if(is.null(age_grid)){
		if((!is.null(guess_PSR)) && (length(guess_PSR)>1)) return(list(success = FALSE, error = sprintf("Invalid number of guessed PSRs; since no age grid was provided, you must provide a single (constant) guess_PSR or none at all")));
		age_grid = 0 # single-point grid, means that PSRs are assumed time-independent
		NG = 1
	}else{
		NG = length(age_grid)
		if((!is.null(guess_PSR)) && (length(guess_PSR)!=1) && (length(guess_PSR)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of guessed PSRs (%d); since an age grid of size %d was provided, you must either provide one or %d PSRs",length(guess_PSR),NG,NG)));
		if((length(age_grid)>1) && (age_grid[NG]>oldest_age-1e-5*(age_grid[NG]-age_grid[NG-1]))) age_grid[NG] = max(age_grid[NG],oldest_age); # if age_grid "almost" covers oldest_age (i.e. up to rounding errors), then fix the remaining difference
		if((length(age_grid)>1) && (age_grid[1]<1e-5*(age_grid[2]-age_grid[1]))) age_grid[1] = min(age_grid[1],0); # if age_grid "almost" covers present-day (i.e. up to rounding errors), then fix the remaining difference
	}
	if(is.null(max_model_runtime)) max_model_runtime = 0;
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree: Extected one of 0,1,2,3.")));
	if(NG==1) splines_degree = 1; # no point in using splines since PSR is assumed to be time-independent
	
	# reformat shape of input params to an internally standardized format
	if(verbose) cat(sprintf("%sPreparing for fitting..\n",verbose_prefix))
	if(length(min_PSR)==1) min_PSR = rep(min_PSR,times=NG);
	if(length(max_PSR)==1) max_PSR = rep(max_PSR,times=NG);
	min_PSR = pmax(0,min_PSR);
	max_PSR = pmax(min_PSR,max_PSR);
	if(is.null(guess_PSR)){
		guess_PSR = rep(NA,times=NG);
	}else if(length(guess_PSR)==1){
		guess_PSR = rep(guess_PSR,times=NG);
	}
	if(is.null(fixed_PSR)){
		fixed_PSR = rep(NA,times=NG);
	}else if(length(fixed_PSR)==1){
		fixed_PSR = rep(fixed_PSR,times=NG);
	}

	# verify that fixed params are within the imposed bounds
	if(any(fixed_PSR[!is.na(fixed_PSR)]<min_PSR[!is.na(fixed_PSR)]) || any(fixed_PSR[!is.na(fixed_PSR)]>max_PSR[!is.na(fixed_PSR)])){
		return(list(success = FALSE, error=sprintf("Some fixed PSRs are outside of their fitting bounds")));
	}
						
	#################################
	# PREPARE PARAMETERS TO BE FITTED
	
	# guess reasonable start params, if not provided
	default_guess_PSR = mean(lineage_counter$relative_slopes[lineage_counter$relative_slopes>=0]); # a reasonable guesstimate for the PSR is the relative LTT-slope. Omitting negative slopes likely due to numerical inaccuracies
	if((!is.finite(default_guess_PSR)) || (default_guess_PSR==0)) default_guess_PSR=log(LTT0)/root_age
	guess_PSR[!is.finite(guess_PSR)] = default_guess_PSR;
	guess_PSR = pmin(max_PSR, pmax(min_PSR, guess_PSR)); # make sure initial guess is within the imposed bounds
	
	# determine which parameters are to be fitted
	fixed_param_values 	= c(fixed_PSR); # may contain NAs, corresponding to non-fixed parameters
	fitted_params		= which(is.na(fixed_param_values))
	fixed_params		= which(!is.na(fixed_param_values))
	guess_param_values 	= c(guess_PSR); # should contain a valid numeric for each parameter, even if the parameter is fixed
	guess_param_values[fixed_params] = fixed_param_values[fixed_params] # make sure guessed param values are consistent with fixed param values
	min_param_values	= c(min_PSR);
	max_param_values	= c(max_PSR);
	NFP					= length(fitted_params)
		
	# determine typical parameter scales
	scale_PSR = abs(guess_PSR); scale_PSR[scale_PSR==0] = mean(scale_PSR); scale_PSR[scale_PSR==0]=default_guess_PSR;
	param_scales = c(scale_PSR);

	# set fit-control options, unless provided by the caller
	if(is.null(fit_control)) fit_control = list()
	if(is.null(fit_control$step.min)) fit_control$step.min 	= 0.001
	if(is.null(fit_control$rel.tol)) fit_control$rel.tol 	= 1e-10
	if(is.null(fit_control$x.tol)) fit_control$x.tol 		= 1e-10
	if(is.null(fit_control$iter.max)) fit_control$iter.max 	= 1000
	if(is.null(fit_control$eval.max)) fit_control$eval.max 	= 2 * fit_control$iter.max * NFP

	################################
	# FITTING
	
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scales
	objective_function = function(fparam_values, trial){
		param_values = fixed_param_values; param_values[fitted_params] = fparam_values * param_scales[fitted_params];
		if(any(is.nan(param_values)) || any(is.infinite(param_values))) return(Inf); # catch weird cases where params become NaN
		PSRs = param_values[1:NG]; 
		if(length(age_grid)==1){
			# while age-grid has only one point (i.e., PSRs are constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age);
			input_PSRs 		= c(PSRs, PSRs);
		}else{
			input_age_grid 	= age_grid;
			input_PSRs 		= PSRs
		}
		results = HBD_PSR_loglikelihood_CPP(branching_ages		= sorted_node_ages,
											oldest_age			= oldest_age,
											age_grid 			= input_age_grid,
											PSRs 				= input_PSRs,
											splines_degree		= splines_degree,
											condition			= condition,
											relative_dt			= relative_dt,
											runtime_out_seconds	= max_model_runtime);
		if(!results$success){
			loglikelihood = -Inf
		}else if(is.na(results$loglikelihood) || is.nan(results$loglikelihood)){
			loglikelihood = -Inf
		}else{		
			loglikelihood = results$loglikelihood
		}
		if(diagnostics){
			if(results$success){ cat(sprintf("%s  Trial %d: loglikelihood %.10g, model runtime %.5g sec\n",verbose_prefix,trial,loglikelihood,results$runtime)) }
			else{ cat(sprintf("%s  Trial %d: Model evaluation failed: %s\n",verbose_prefix,trial,results$error)) }
		}
		return(-loglikelihood)
	}
	

	# fit with various starting points
	fit_single_trial = function(trial){
		scales		 = param_scales[fitted_params]
		lower_bounds = min_param_values[fitted_params]
		upper_bounds = max_param_values[fitted_params]
		# randomly choose start values for fitted params
		start_values = guess_param_values[fitted_params]
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
		start_values = pmax(lower_bounds,pmin(upper_bounds,start_values))
		# run fit
		fit = stats::nlminb(start_values/scales, 
							objective	= function(fpars){ objective_function(fpars, trial) }, 
							lower		= lower_bounds/scales, 
							upper		= upper_bounds/scales, 
							control		= fit_control)
		LL 			= -fit$objective
		converged	= (fit$convergence==0)
		if(diagnostics) cat(sprintf("%s  Trial %d: Final loglikelihood %.10g, converged = %d\n",verbose_prefix,trial,LL,converged))
		return(list(objective_value=fit$objective, fparam_values = fit$par*scales, converged=converged, Niterations=fit$iterations, Nevaluations=fit$evaluations[1]));
	}
	
	################################

	# run one or more independent fitting trials
    if((Ntrials>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
		if(verbose) cat(sprintf("%sFitting %d model parameters (%d trials, parallelized)..\n",verbose_prefix,NFP,Ntrials))
		# run trials in parallel using multiple forks
		# Note: Forks (and hence shared memory) are not available on Windows
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
	loglikelihood		= objective_value
	fitted_param_values = fixed_param_values; fitted_param_values[fitted_params] = fits[[best]]$fparam_values;
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
	fitted_PSR			= fitted_param_values[1:NG]

	# reverse any time shift due to earlier tree trimming
	age_grid 	= age_grid + age0
	oldest_age 	= oldest_age + age0
	root_age 	= root_age + age0
	
	# calculate deterministic LTT of fitted congruence class on the age grid
	fitted_LTT = LTT0 * exp(-get_antiderivative_of_splines_function(age_grid, age0, fitted_PSR, splines_degree, age_grid));
	
	
	#######################################################################
	# estimate confidence intervals if needed, via parametric bootstrapping
	
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sEstimating confidence intervals using %d parametric bootstraps..\n",verbose_prefix,Nbootstraps))
		if(is.null(Ntrials_per_bootstrap)) Ntrials_per_bootstrap = max(1,Ntrials)
		#   include a dummy age grid point at the end of age_grid if needed (extrapolating the fitted PSR as a constant), so as to cover the root age
		#   also include a dummy age grid point at the beginning if necessary (extrapolating PSR as a constant), to cover present-day (age 0)
		sim_age_grid = age_grid
		sim_PSR		 = fitted_PSR
		if(tail(sim_age_grid,1)<root_age){
			# extrapolate sim_PSR (as a constant) all the way to the root (a bit more to avoid rounding errors)
			sim_age_grid = c(sim_age_grid, root_age*1.01)
			sim_PSR 	 = c(sim_PSR, tail(sim_PSR,1));
		}
		if(sim_age_grid[1]>0){
			# extrapolate sim_PSR (as a constant) all the way to the present-day (age 0)
			sim_age_grid = c(0,sim_age_grid)
			sim_PSR		 = c(sim_PSR[1],sim_PSR)
		}
		bootstrap_PSR = matrix(NA,nrow=Nbootstraps,ncol=NG)
		NBsucceeded	  = 0
		for(b in 1:Nbootstraps){
			# simulate model with fitted parameters
			if(verbose) cat(sprintf("%s  Bootstrap #%d..\n",verbose_prefix,b))
			bootstrap = castor::generate_tree_hbd_reverse(	Ntips			= original_Ntips,
															crown_age		= root_age,
															age_grid		= sim_age_grid, 
															PSR				= sim_PSR,
															splines_degree	= splines_degree,
															relative_dt		= relative_dt)
			if(!bootstrap$success) return(list(success=FALSE, error=sprintf("Bootstrapping failed: Could not generate tree for the fitted PSR: %s",bootstrap$error), age_grid=age_grid, fitted_PSR=fitted_PSR, loglikelihood=loglikelihood));
			bootstrap_tree = bootstrap$trees[[1]]

			# fit PSR to simulated tree
			fit = fit_hbd_psr_on_grid(	bootstrap_tree, 
										oldest_age			= oldest_age,
										age0				= age0,
										age_grid			= age_grid,
										min_PSR				= min_PSR,
										max_PSR				= max_PSR,
										guess_PSR			= guess_PSR,
										fixed_PSR			= fixed_PSR,
										splines_degree		= splines_degree,
										condition			= condition,
										relative_dt			= relative_dt,
										Ntrials				= Ntrials_per_bootstrap,
										Nbootstraps			= 0,
										Nthreads			= Nthreads,
										max_model_runtime	= max_model_runtime,
										fit_control			= fit_control,
										verbose				= verbose,
										diagnostics			= diagnostics,
										verbose_prefix		= paste0(verbose_prefix,"    "))
			if(!fit$success){
				if(verbose) cat(sprintf("%s  WARNING: Fitting failed for this bootstrap: %s\n",verbose_prefix,fit$error))
			}else{
				bootstrap_PSR[b,] = fit$fitted_PSR
				NBsucceeded = NBsucceeded + 1
			}
		}
		# calculate standard errors and confidence intervals from distribution of bootstrapped parameters
		standard_errors = sqrt(pmax(0, colMeans(bootstrap_PSR^2, na.rm=TRUE) - colMeans(bootstrap_PSR, na.rm=TRUE)^2))
		quantiles = sapply(1:ncol(bootstrap_PSR), FUN=function(p) quantile(bootstrap_PSR[,p], probs=c(0.25, 0.75, 0.025, 0.975), na.rm=TRUE, type=8))
		CI50lower = quantiles[1,]
		CI50upper = quantiles[2,]
		CI95lower = quantiles[3,]
		CI95upper = quantiles[4,]
	}
		
	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				fitted_PSR				= fitted_PSR,
				guess_PSR				= guess_param_values[1:NG],
				age_grid				= age_grid,
				fitted_LTT				= fitted_LTT,
				NFP						= NFP,
				AIC						= 2*NFP - 2*loglikelihood,
				BIC						= log(sum((sorted_node_ages<=oldest_age) & (sorted_node_ages>=age0)))*NFP - 2*loglikelihood,
				converged				= fits[[best]]$converged,
				Niterations				= fits[[best]]$Niterations,
				Nevaluations			= fits[[best]]$Nevaluations,
				bootstrap_estimates		= (if(Nbootstraps>0) bootstrap_PSR else NULL),
				standard_errors			= (if(Nbootstraps>0) standard_errors else NULL),
				CI50lower				= (if(Nbootstraps>0) CI50lower else NULL),
				CI50upper				= (if(Nbootstraps>0) CI50upper else NULL),
				CI95lower				= (if(Nbootstraps>0) CI95lower else NULL),
				CI95upper				= (if(Nbootstraps>0) CI95upper else NULL)))

}



completement = function(N, indices){
	pool = rep(TRUE,N);
	pool[indices] = FALSE;
	return(which(pool));
}


