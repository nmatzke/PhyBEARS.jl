# Fit a homogenous-birth-death cladogenic model-congruence-class to an ultrametric timetree, by fitting the pulled diversification rate (PDR)
# An HBD model is defined by a time-dependent speciation rate (lambda), a time-dependent extinction rate (mu) and a rarefaction (rho, subsampling fraction)
# However, for each specific model and a given timetree there exists a continuum of alternative models that would all generate the same deterministic lineages-through-time (LTT) curve (when calculated backward in time), and all of these models actually have the same likelihood.
# Hence, each model is part of an "equivalence class" of models, and likelihood-based approaches can only discern between model classes, but not between the individual model members in a class
# It turns out that each HBD model-class is uniquely defined by its "pulled diversification rate" (PDR) and the product rho(age0)*lambda(age0)=:rholambda0, where rho(age0) is the fraction of lineages extant at age0 that are represented in the tree.
# This function thus fits model-classes, rather than models, by fitting the PDR and the parameter rholambda0.
#
# This function can optionally fit the congruence class by considering only the age interval [age0, oldest_age].
# oldest_age is enforced by (mathematically) cutting the tree into multiple subtrees stemming at age oldest_age (omitting anything older) and treating each subtree as an independent realization of the same stochastic BD process
# age0 is enforced by (actually) trimming the tree at age0 (omitting anything younger than age0) and fitting the HBD class to the new (shorter) timetree. 
# In that case, the fitted rholambda0 is actually the product Phi(age0)*lambda(age0), where Phi(age0) is the probability that a lineage extant at age tau would be eventually included in the original timetree; Phi(age0) is essentially the sampling fraction at age0.
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
fit_hbd_pdr_on_grid = function(	tree, 
								oldest_age				= NULL,		# either a numeric specifying the stem age or NULL (equivalent to the root age). This is similar to the "tot_time" option in the R function RPANDA::likelihood_bd
								age0					= 0,		# non-negative numeric, youngest age (time before present) to consider when fitting and with respect to which rholambda0 is defined (i.e. rholambda0 = rho(age0)*lambda(age0))
								age_grid				= NULL,		# either NULL, or a numeric vector of size NG, listing ages in ascending order, on which the PDR is defined as a piecewise linear curve. If NULL, the PDR is assumed to be time-independent.
								min_PDR					= -Inf,		# optional lower bound for the fitted PDRs. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the lower bound at each age-grid point.
								max_PDR					= +Inf,		# optional upper bound for the fitted PDRs. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the upper bound at each age-grid point.
								min_rholambda0			= 1e-10,	# optional lower bound for the fitted rholambda0. Note that rholambda0 is always non-negative.
								max_rholambda0			= +Inf,		# optional upper bound for the fitted rholambda0
								guess_PDR				= NULL,		# initial guess for the PDR. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant PDR at all ages) or a numeric vector of size NG specifying an initial guess for the PDR at each age-grid point (can include NAs)
								guess_rholambda0		= NULL,		# initial guess for the product rho*lambda(0). Either NULL (an initial guess will be computed automatically) or a single strictly-positive numeric.
								fixed_PDR				= NULL,		# optional fixed PDR values, on one or more of the age grid points. Either NULL (none of the PDRs are fixed), or a single scalar (all PDRs are fixed) or a numeric vector of size NG (some or all PDRs are fixed, can include NAs).
								fixed_rholambda0		= NULL,		# optional fixed value for rholambda0. If non-NULL and non-NA, then rholambda0 is not fitted. 
								splines_degree			= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by the PDR on the age grid.
								condition				= "auto",	# one of "crown" or "stem" or "stem2" (or "stem3", "stem4", .. etc) or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age<root_age, "stem2" when oldest_age>root_age, or "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
								relative_dt				= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
								Ntrials					= 1,
								Nbootstraps				= 0,		# (integer) optional number of parametric-bootstrap samples (random trees generated using the fitted PDR & rholambda0) for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
								Ntrials_per_bootstrap	= NULL,		# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
								Nthreads				= 1,
								max_model_runtime		= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
								fit_control				= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
								verbose					= FALSE,	# boolean, specifying whether to print informative messages
								verbose_prefix			= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	original_Ntips = length(tree$tip.label)
	if(tree$Nnode<1) return(list(success = FALSE, error="Input tree is too small"));
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
		if(tree$Nnode<1) return(list(success = FALSE, error=sprintf("Tree is too small after trimming at age0 (%g)",age0)));
		if(!is.null(oldest_age)) oldest_age	= oldest_age - age0	
		if(!is.null(age_grid)) age_grid 	= age_grid - age0
		root_age = root_age - age0
	}
								
	# pre-compute some tree stats
	if(verbose) cat(sprintf("%sPrecomputing some stats about the tree..\n",verbose_prefix))
	LTT0				= length(tree$tip.label);
	lineage_counter 	= count_lineages_through_time(tree, Ntimes=max(3,log2(LTT0)), include_slopes=TRUE, ultrametric=TRUE)
	sorted_node_ages	= sort(get_all_branching_ages(tree));
	root_age 		 	= tail(sorted_node_ages,1);
	age_epsilon		 	= 1e-4*mean(tree$edge.length);

	# more error checking
	if(Ntrials<1) return(list(success = FALSE, error = sprintf("Ntrials must be at least 1")))
	if(is.null(age_grid)){
		if((!is.null(guess_PDR)) && (length(guess_PDR)>1)) return(list(success = FALSE, error = sprintf("Invalid number of guessed PDRs; since no age grid was provided, you must provide a single (constant) guess_PDR or none at all")));
		age_grid = 0 # single-point grid, means that PDRs are assumed time-independent
		NG = 1
	}else{
		NG = length(age_grid)
		if((!is.null(guess_PDR)) && (length(guess_PDR)!=1) && (length(guess_PDR)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of guessed PDRs (%d); since an age grid of size %d was provided, you must either provide one or %d PDRs",length(guess_PDR),NG,NG)));
		if((length(age_grid)>1) && (age_grid[NG]>oldest_age-1e-5*(age_grid[NG]-age_grid[NG-1]))) age_grid[NG] = max(age_grid[NG],oldest_age); # if age_grid "almost" covers oldest_age (i.e. up to rounding errors), then fix the remaining difference
		if((length(age_grid)>1) && (age_grid[1]<1e-5*(age_grid[2]-age_grid[1]))) age_grid[1] = min(age_grid[1],0); # if age_grid "almost" covers present-day (i.e. up to rounding errors), then fix the remaining difference
	}
	if(is.null(max_model_runtime)) max_model_runtime = 0;
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree: Extected one of 0,1,2,3.")));
	if(NG==1) splines_degree = 1; # no point in using splines since PDR is assumed to be time-independent
	
	# reformat shape of input params to an internally standardized format
	if(verbose) cat(sprintf("%sPreparing for fitting..\n",verbose_prefix))
	min_rholambda0 = max(0,min_rholambda0);
	max_rholambda0 = max(0,max_rholambda0);
	if(length(min_PDR)==1) min_PDR = rep(min_PDR,times=NG);
	if(length(max_PDR)==1) max_PDR = rep(max_PDR,times=NG);
	if(is.null(guess_rholambda0)) guess_rholambda0 = NA;
	if(is.null(fixed_rholambda0)) fixed_rholambda0 = NA;
	if(is.null(guess_PDR)){
		guess_PDR = rep(NA,times=NG);
	}else if(length(guess_PDR)==1){
		guess_PDR = rep(guess_PDR,times=NG);
	}
	if(is.null(fixed_PDR)){
		fixed_PDR = rep(NA,times=NG);
	}else if(length(fixed_PDR)==1){
		fixed_PDR = rep(fixed_PDR,times=NG);
	}

	# verify that fixed params are within the imposed bounds
	if((!is.na(fixed_rholambda0)) && ((fixed_rholambda0<min_rholambda0) || (fixed_rholambda0>max_rholambda0))){
		return(list(success = FALSE, error=sprintf("Fixed rholambda0 (%g) is outside of the requested bounds (%g - %g)",fixed_rholambda0,min_rholambda0,max_rholambda0)));
	}
	if(any(fixed_PDR[!is.na(fixed_PDR)]<min_PDR[!is.na(fixed_PDR)]) || any(fixed_PDR[!is.na(fixed_PDR)]>max_PDR[!is.na(fixed_PDR)])){
		return(list(success = FALSE, error=sprintf("Some fixed PDRs are outside of the requested bounds")));
	}
						
	#################################
	# PREPARE PARAMETERS TO BE FITTED
	
	# guess reasonable start params, if not provided
	default_guess_PDR = mean(lineage_counter$relative_slopes); # a reasonable guesstimate for the average PDR is the average of the relative LTT-slope
	guess_PDR[is.na(guess_PDR)] = default_guess_PDR;
	if(is.na(guess_rholambda0)) guess_rholambda0 = tail(lineage_counter$relative_slopes[lineage_counter$relative_slopes>0],1)
	if(is.null(guess_rholambda0) || (length(guess_rholambda0)==0) || (!is.finite(guess_rholambda0)) || (guess_rholambda0==0)) guess_rholambda0 = log(LTT0)/root_age
	
	# make sure initial guess is within the imposed bounds
	guess_PDR = pmin(max_PDR, pmax(min_PDR, guess_PDR));
	guess_rholambda0 = min(max_rholambda0, max(min_rholambda0, guess_rholambda0))
	
	# determine which parameters are to be fitted
	# convention: parameters are indexed as follows: [PDR[], rholambda0]
	fixed_param_values 	= c(fixed_PDR, fixed_rholambda0); # may contain NAs, corresponding to non-fixed parameters
	fitted_params		= which(is.na(fixed_param_values))
	fixed_params		= which(!is.na(fixed_param_values))
	guess_param_values 	= c(guess_PDR, guess_rholambda0); # should contain a valid numeric for each parameter, even if the parameter is fixed
	guess_param_values[fixed_params] = fixed_param_values[fixed_params] # make sure guessed param values are consistent with fixed param values
	min_param_values	= c(min_PDR,min_rholambda0);
	max_param_values	= c(max_PDR,max_rholambda0);
	NP					= length(fixed_param_values)
	NFP					= length(fitted_params);
	
	# determine typical parameter scales
	scale_PDR = abs(guess_PDR); scale_PDR[scale_PDR==0] = mean(scale_PDR);
	scale_rholambda0 = abs(guess_rholambda0);
	if(scale_rholambda0==0) scale_rholambda0 = log2(LTT0)/root_age;
	param_scales = c(scale_PDR,scale_rholambda0);

	# set fit-control options, unless provided by the caller
	if(is.null(fit_control)) fit_control = list()
	if(is.null(fit_control$step.min)) fit_control$step.min = 0.001
	if(is.null(fit_control$x.tol)) fit_control$x.tol = 1e-8
	if(is.null(fit_control$iter.max)) fit_control$iter.max = 1000
	if(is.null(fit_control$eval.max)) fit_control$eval.max = 2 * fit_control$iter.max * NFP
	

	################################
	# FITTING
	
	
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scales
	objective_function = function(fparam_values){
		param_values = fixed_param_values; param_values[fitted_params] = fparam_values * param_scales[fitted_params];
		if(any(is.nan(param_values)) || any(is.infinite(param_values))) return(Inf); # catch weird cases where params become NaN
		PDRs = param_values[1:NG]; 
		rholambda0 = param_values[NG+1];
		if(length(age_grid)==1){
			# while age-grid has only one point (i.e., PDRs are constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age);
			input_PDRs 		= c(PDRs, PDRs);
		}else{
			input_age_grid 	= age_grid;
			input_PDRs 		= PDRs
		}
		results = HBD_PDR_loglikelihood_CPP(branching_ages		= sorted_node_ages,
											oldest_age			= oldest_age,
											rholambda0 			= rholambda0,
											age_grid 			= input_age_grid,
											PDRs 				= input_PDRs,
											splines_degree		= splines_degree,
											condition			= condition,
											relative_dt			= relative_dt,
											runtime_out_seconds	= max_model_runtime,
											diff_PDR			= numeric(),
											diff_PDR_degree		= 0)
		if(!results$success) return(Inf)
		LL = results$loglikelihood
		if(is.na(LL) || is.nan(LL)) return(Inf)
		return(-LL)
	}
	

	fitted_grid_params = fitted_params[fitted_params!=(NG+1)]
	gradient_function = function(fparam_values){
		if(splines_degree!=1) return(NaN); # only implemented for splines_degree=1
		param_values = fixed_param_values; param_values[fitted_params] = fparam_values * param_scales[fitted_params];
		if(any(is.nan(param_values)) || any(is.infinite(param_values))) return(Inf); # catch weird cases where params become NaN
		PDRs = param_values[1:NG]; 
		rholambda0 = param_values[NG+1];
		if(NG==1){
			# while age-grid has only one point (i.e., PDRs are constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age);
			input_PDRs 		= c(PDRs, PDRs);
		}else{
			input_age_grid 	= age_grid;
			input_PDRs 		= PDRs
		}
		# calculate differentials of PDR in the directions of the various fitted parameters (fitted grid ages and PDRs)
		diff_PDR_degree = 1
		if(NG==1){
			diff_PDR = (if(is.na(fixed_PDR)) c(1,0,1,0) else numeric())
		}else{
			diff_PDR_all = derivatives_of_grid_curve_CPP(Xgrid=age_grid, Ygrid=PDRs) # yields a 3D array of size (2*NG)*NG*2, flattened in layer-row-major format, representing the derivatives of the PDR w.r.t. the grid ages and the PDR values on the grid points
			diff_PDR = unlist(lapply(fitted_grid_params, FUN=function(p) diff_PDR_all[((NG+p-1)*NG*(diff_PDR_degree+1) + 1):((NG+p-1)*NG*(diff_PDR_degree+1) + NG*(diff_PDR_degree+1))])) # extract only differentials along fitted parameters. Note that the first NG differentials are always omitted, because they correspond to the grid ages themselves, which are held constant in this case.
		}
		results = HBD_PDR_loglikelihood_CPP(branching_ages		= sorted_node_ages,
											oldest_age			= oldest_age,
											rholambda0 			= rholambda0,
											age_grid 			= input_age_grid,
											PDRs 				= input_PDRs,
											splines_degree		= splines_degree,
											condition			= condition,
											relative_dt			= relative_dt,
											runtime_out_seconds	= max_model_runtime,
											diff_PDR			= diff_PDR,
											diff_PDR_degree		= diff_PDR_degree);
		if(!results$success) return(rep(Inf,times=NFP));
		gradient_full = rep(NA, times=NP)
		gradient_full[NG+1] = results$dLL_drholambda0
		gradient_full[fitted_grid_params] = results$dLL_dPDR
		gradient = gradient_full[fitted_params] * param_scales[fitted_params]
		return(-gradient);
	}
		

	# fit with various starting points
	fit_single_trial = function(trial){
		scales		 = param_scales[fitted_params]
		lower_bounds = min_param_values[fitted_params]
		upper_bounds = max_param_values[fitted_params]
		# randomly choose start values for fitted params
		if(trial==1){
			start_values = guess_param_values[fitted_params]		
		}else{
			start_values = get_random_params(defaults=guess_param_values[fitted_params], lower_bounds=lower_bounds, upper_bounds=upper_bounds, scales=scales, orders_of_magnitude=4)
		}

		# check if start values yield NaN
		start_LL = objective_function(start_values/scales)
		if(!is.finite(start_LL)) return(list(objective_value=NA, fparam_values = rep(NA,times=NFP), converged=FALSE, Niterations=0, Nevaluations=1));
												
		# fine-tune some fitting controls based on the initial model evaluation
		if(is.null(fit_control$rel.tol)) fit_control$rel.tol = max(1e-30,min(1e-5,0.0001/abs(start_LL)))		
				
		# run fit
		fit = stats::nlminb(start_values/scales, 
							objective	= objective_function, 
							gradient	= (if(splines_degree==1) gradient_function else NULL),
							lower		= lower_bounds/scales, 
							upper		= upper_bounds/scales, 
							control		= fit_control)
		return(list(objective_value=fit$objective, fparam_values = fit$par*scales, converged=(fit$convergence==0), Niterations=fit$iterations, Nevaluations=fit$evaluations[1]));
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
	loglikelihood		= objective_value
	fitted_param_values = fixed_param_values; fitted_param_values[fitted_params] = fits[[best]]$fparam_values;
	fitted_PDR			= fitted_param_values[1:NG]
	fitted_rholambda0	= fitted_param_values[NG+1] 
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
	
	# reverse any time shift due to earlier tree trimming
	age_grid 	= age_grid + age0
	oldest_age 	= oldest_age + age0
	root_age 	= root_age + age0
	
	
	#######################################################################
	# estimate confidence intervals if needed, via parametric bootstrapping
	
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sEstimating confidence intervals using %d parametric bootstraps..\n",verbose_prefix,Nbootstraps))
		if(verbose) cat(sprintf("%s  Calculating pulled speciation rate from PDR, for simulating trees..\n",verbose_prefix))
		# first calculate the PSR from the PDR & rholambda0
		#   include a dummy age grid point at the end of age_grid if needed (extrapolating the fitted PDR as a constant), so as to cover the root age
		#   also include a dummy age grid point at the beginning if necessary (extrapolating PDR as a constant), to cover present-day (age 0)
		sim_age_grid = age_grid
		sim_PDR		 = fitted_PDR
		if(tail(sim_age_grid,1)<root_age){
			# extrapolate sim_PDR (as a constant) all the way to the root (a bit more to avoid rounding errors)
			sim_age_grid = c(sim_age_grid, root_age*1.01)
			sim_PDR 	 = c(sim_PDR, tail(sim_PDR,1));
		}
		if(sim_age_grid[1]>0){
			# extrapolate sim_PDR (as a constant) all the way to the present-day (age 0)
			sim_age_grid = c(0,sim_age_grid)
			sim_PDR		 = c(sim_PDR[1],sim_PDR)
		}
		sim = get_PSR_from_PDR_HBD(	age0			= age0,
									oldest_age 		= tail(sim_age_grid,1),
									age_grid		= sim_age_grid,
									PDR				= sim_PDR,
									rholambda0		= fitted_rholambda0,
									splines_degree	= splines_degree,
									relative_dt		= relative_dt,
									include_nLTT0	= TRUE);
		if(!sim$success) return(list(success=FALSE, error=sprintf("Bootstrapping failed: Could not calculate PSR corresponding to fitted PDR: %s",sim$error), age_grid=age_grid, fitted_PDR=fitted_PDR, fitted_rholambda0=fitted_rholambda0, loglikelihood=loglikelihood));
		if(is.null(Ntrials_per_bootstrap)) Ntrials_per_bootstrap = max(1,Ntrials)
		bootstrap_params = matrix(NA,nrow=Nbootstraps,ncol=NG+1)
		NBsucceeded		 = 0
		for(b in 1:Nbootstraps){
			# simulate model with fitted parameters
			if(verbose) cat(sprintf("%s  Bootstrap #%d..\n",verbose_prefix,b))
			bootstrap = castor::generate_tree_hbd_reverse(	Ntips			= LTT0/sim$nLTT0,
															crown_age		= root_age,
															age_grid		= sim$ages, 
															PSR				= sim$PSR,
															splines_degree	= 1,
															relative_dt		= relative_dt)
			if(!bootstrap$success) return(list(success=FALSE, error=sprintf("Bootstrapping failed: Could not generate tree for the fitted PDR: %s",bootstrap$error), age_grid=age_grid, fitted_PDR=fitted_PDR, fitted_rholambda0=fitted_rholambda0, loglikelihood=loglikelihood));
			bootstrap_tree = bootstrap$trees[[1]]

			# fit PSR to simulated tree
			fit = fit_hbd_pdr_on_grid(	tree				= bootstrap_tree, 
										oldest_age			= oldest_age,
										age0				= age0,
										age_grid			= age_grid,
										min_PDR				= min_PDR,
										max_PDR				= max_PDR,
										min_rholambda0		= min_rholambda0,
										max_rholambda0		= max_rholambda0,
										guess_PDR			= guess_PDR,
										guess_rholambda0	= guess_rholambda0,
										fixed_PDR			= fixed_PDR,
										fixed_rholambda0	= fixed_rholambda0,
										splines_degree		= 1,
										condition			= condition,
										relative_dt			= relative_dt,
										Ntrials				= Ntrials_per_bootstrap,
										Nbootstraps			= 0,
										Nthreads			= Nthreads,
										max_model_runtime	= max_model_runtime,
										fit_control			= fit_control,
										verbose				= verbose,
										verbose_prefix		= paste0(verbose_prefix,"    "))
			if(!fit$success){
				if(verbose) cat(sprintf("%s  WARNING: Fitting failed for this bootstrap: %s\n",verbose_prefix,fit$error))
			}else{
				bootstrap_params[b,] = c(fit$fitted_PDR, fit$fitted_rholambda0)
				NBsucceeded = NBsucceeded + 1
			}
		}
		# calculate standard errors and confidence intervals from distribution of bootstrapped parameters
		standard_errors_flat = sqrt(pmax(0, colMeans(bootstrap_params^2, na.rm=TRUE) - colMeans(bootstrap_params, na.rm=TRUE)^2))
		standard_errors = list(PDR=standard_errors_flat[1:NG], rholambda0=standard_errors_flat[NG+1])
		quantiles = sapply(1:ncol(bootstrap_params), FUN=function(p) quantile(bootstrap_params[,p], probs=c(0.25, 0.75, 0.025, 0.975, 0.5), na.rm=TRUE, type=8))
		CI50lower = list(PDR=quantiles[1,1:NG], rholambda0=quantiles[1,NG+1])
		CI50upper = list(PDR=quantiles[2,1:NG], rholambda0=quantiles[2,NG+1])
		CI95lower = list(PDR=quantiles[3,1:NG], rholambda0=quantiles[3,NG+1])
		CI95upper = list(PDR=quantiles[4,1:NG], rholambda0=quantiles[4,NG+1])
		medians   = list(PDR=quantiles[5,1:NG], rholambda0=quantiles[5,NG+1])
		bootstrap_estimates = list(PDR=bootstrap_params[,1:NG], rholambda0=bootstrap_params[,NG+1])
	}

	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				fitted_PDR				= fitted_PDR,
				fitted_rholambda0		= fitted_rholambda0,
				guess_PDR				= guess_param_values[1:NG],
				guess_rholambda0		= guess_param_values[NG+1],
				age_grid				= age_grid,
				NFP						= NFP,
				AIC						= 2*NFP - 2*loglikelihood,
				BIC						= log(sum((sorted_node_ages<=oldest_age) & (sorted_node_ages>=age0)))*NFP - 2*loglikelihood,
				converged				= fits[[best]]$converged,
				Niterations				= fits[[best]]$Niterations,
				Nevaluations			= fits[[best]]$Nevaluations,
				bootstrap_estimates		= (if(Nbootstraps>0) bootstrap_estimates else NULL),
				standard_errors			= (if(Nbootstraps>0) standard_errors else NULL),
				medians					= (if(Nbootstraps>0) medians else NULL),
				CI50lower				= (if(Nbootstraps>0) CI50lower else NULL),
				CI50upper				= (if(Nbootstraps>0) CI50upper else NULL),
				CI95lower				= (if(Nbootstraps>0) CI95lower else NULL),
				CI95upper				= (if(Nbootstraps>0) CI95upper else NULL)))
}




