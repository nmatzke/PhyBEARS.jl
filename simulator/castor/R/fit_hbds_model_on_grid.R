# Fit a homogenous-birth-death-sampling cladogenic model to a timetree, by estimating lambda & mu & psi & kappa on a temporal grid
# An HBDS model is defined by a time-dependent speciation rate (lambda), a time-dependent extinction rate (mu), a time-dependent sampling rate (psi), and a time-dependent retention rate (kappa).
# Optionally, the model may include a finite set of "concentrated sampling attempts" that occurred at specific times t_1,..,t_NCSA with specific sampling probabilities psi_1,..,psi_NCSA and specific retenion probabilities kappa_1,..,kappa_NCSA
fit_hbds_model_on_grid = function(	tree, 
									root_age 				= NULL,		# optional numeric, the age of the root. Can be used to define a time offset, e.g. if the last tip was not actually sampled at the present. If NULL, this will be calculated from the tree and it will be assumed that the last tip was sampled at the present (age 0)
									oldest_age				= NULL,		# numeric, specifying the oldest age to consider. Can also be NULL, in which case this is set to the root age
									age_grid				= NULL,		# either NULL, or a numeric vector of size NG, listing ages in ascending order, on which lambda, mu, psi, kappa are defined as piecewise linear curves. If NULL, all model variables are assumed to be time-independent.
									CSA_ages				= NULL,		# optional numeric vector of length NCSA, listing the ages of concentrated sampling attempts in strictly ascending order. These ages are assumed to be known, i.e. they are not fitted
									min_lambda				= 0,		# lower bound for the fitted lambdas. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the lower bound at each age-grid point.
									max_lambda				= +Inf,		# upper bound for the fitted lambdas. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the upper bound at each age-grid point.
									min_mu					= 0,		# lower bound for the fitted mus. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the lower bound at each age-grid point.
									max_mu					= +Inf,		# upper bound for the fitted mus. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the upper bound at each age-grid point.
									min_psi					= 0,		# lower bound for the fitted psis. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the lower bound at each age-grid point.
									max_psi					= +Inf,		# upper bound for the fitted psis. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the upper bound at each age-grid point.
									min_kappa				= 0,		# lower bound for the fitted kappas. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the lower bound at each age-grid point.
									max_kappa				= 1,		# upper bound for the fitted kappas. Either a single numeric (applying to all age-grid-points) or a numeric vector of size NG, specifying the upper bound at each age-grid point. Note that kappa is a probability, and hence this upper bound should not be greater than 1.
									min_CSA_probs			= 0,		# lower bound for the fitted CSA_probs (sampling probabilities during CSAs). Either a single numeric (applying to all CSAs) or a numeric vector of size NCSA, specifying the lower bound at each CSA.
									max_CSA_probs			= 1,		# upper bound for the fitted CSA_probs. Either a single numeric (applying to all CSAs) or a numeric vector of size NCSA, specifying the upper bound at each CSA. Note that CSA_probs are probabilities, and hence this upper bound should not be greater than 1.
									min_CSA_kappas			= 0,		# lower bound for the fitted CSA_kappas (retention probabilities during CSAs). Either a single numeric (applying to all CSAs) or a numeric vector of size NCSA, specifying the lower bound at each CSA.
									max_CSA_kappas			= 1,		# upper bound for the fitted CSA_kappas. Either a single numeric (applying to all CSAs) or a numeric vector of size NCSA, specifying the upper bound at each CSA. Note that CSA_kappas are probabilities, and hence this upper bound should not be greater than 1.
									guess_lambda			= NULL,		# initial guess for the lambda. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant lambda at all ages) or a numeric vector of size NG specifying an initial guess for the lambda at each age-grid point (can include NAs)
									guess_mu				= NULL,		# initial guess for the mu. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant mu at all ages) or a numeric vector of size NG specifying an initial guess for the mu at each age-grid point (can include NAs)
									guess_psi				= NULL,		# initial guess for the psi. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant psi at all ages) or a numeric vector of size NG specifying an initial guess for the psi at each age-grid point (can include NAs)
									guess_kappa				= NULL,		# initial guess for the kappa. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant kappa at all ages) or a numeric vector of size NG specifying an initial guess for the kappa at each age-grid point (can include NAs)
									guess_CSA_probs			= NULL,		# initial guess for the CSA_probs. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant sampling probability at all CSAs) or a numeric vector of size NCSA specifying an initial guess for the sampling probability at each CSA (can include NAs)
									guess_CSA_kappas		= NULL,		# initial guess for the CSA_kappas. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant retention probability at all CSAs) or a numeric vector of size NCSA specifying an initial guess for the retention probability at each CSA (can include NAs)
									fixed_lambda			= NULL,		# optional fixed lambda values, on one or more of the age grid points. Either NULL (none of the lambdas are fixed), or a single scalar (all lambdas are fixed) or a numeric vector of size NG (some or all lambdas are fixed, can include NAs).
									fixed_mu				= NULL,		# optional fixed mu values, on one or more of the age grid points. Either NULL (none of the mus are fixed), or a single scalar (all mus are fixed) or a numeric vector of size NG (some or all mus are fixed, can include NAs).
									fixed_psi				= NULL,		# optional fixed psi values, on one or more of the age grid points. Either NULL (none of the psis are fixed), or a single scalar (all psis are fixed) or a numeric vector of size NG (some or all psis are fixed, can include NAs).
									fixed_kappa				= NULL,		# optional fixed kappa values, on one or more of the age grid points. Either NULL (none of the kappas are fixed), or a single scalar (all kappas are fixed) or a numeric vector of size NG (some or all kappas are fixed, can include NAs).
									fixed_CSA_probs			= NULL,		# optional fixed CSA_prob values, on one or more of the CSAs. Either NULL (none of the CSA_probs are fixed), or a single scalar (all CSA_probs are fixed) or a numeric vector of size NCSA (some or all CSA_probs are fixed, can include NAs).
									fixed_CSA_kappas		= NULL,		# optional fixed CSA_kappa values, on one or more of the CSAs. Either NULL (none of the CSA_kappas are fixed), or a single scalar (all CSA_kappas are fixed) or a numeric vector of size NCSA (some or all CSA_kappas are fixed, can include NAs).
									fixed_age_grid			= NULL,		# optional age grid on which fixed_lambda, fixed_mu, fixed_psi and fixed_kappa are defined (rather than on the age_grid). If provided, then any fixed_lambda, fixed_mu, fixed_psi and fixed_kappa must be defined on the entire fixed_age_grid. This may be useful if you want to fit some parameters on a coarse grid, but want to specify (fix) some other parameters on a much finer grid. Entries in fixed_age_grid[] must be in ascending order.
									const_lambda			= FALSE,	# logical, whether to enforce a constant (time-independent) fitted speciation rate. Only relevant for those lambdas that are fitted (i.e. fixed lambda values are kept as is).
									const_mu				= FALSE,	# logical, whether to enforce a constant (time-independent) fitted extinction rate. Only relevant for those lambdas that are fitted (i.e. fixed lambda values are kept as is).
									const_psi				= FALSE,	# logical, whether to enforce a constant (time-independent) fitted sampling rate psi. Only relevant for those psis that are fitted (i.e. fixed psi values are kept as is).
									const_kappa				= FALSE,	# logical, whether to enforce a constant (time-independent) fitted retention probability kappa. Only relevant for those kappas that are fitted (i.e. fixed kappas values are kept as is).
									const_CSA_probs			= FALSE,	# logical, whether to enforce a constant (CSA-independent) fitted sampling probability during CSAs. Only relevant for those CSA_probs that are fitted (i.e. fixed CSA_probs values are kept as is).
									const_CSA_kappas		= FALSE,	# logical, whether to enforce a constant (CSA-independent) fitted retention probability during CSAs. Only relevant for those CSA_kappas that are fitted (i.e. fixed CSA_kappas values are kept as is).
									splines_degree			= 1,		# integer, either 0, 1 or 2 or 3, specifying the degree for the splines defined by lambda and mu on the age grid.
									condition				= "auto",	# one of "crown" or "stem" or "none" or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age!=root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
									ODE_relative_dt			= 0.001,	# positive unitless number, relative integration time step for the ODE solvers. Relative to the typical time scales of the dynamics, as estimated from the theoretically maximum possible rate of change. Typical values are 0.001 - 0.1.
									ODE_relative_dy			= 1e-3,		# positive unitless mumber, unitless number, relative step for interpolating simulated values over time. So a ODE_relative_dy of 0.001 means that E is recorded and interpolated between points between which E differs by roughy 0.001. Typical values are 0.01-0.0001. A smaller E_value_step increases interpolation accuracy, but also increases memory requirements and adds runtime (scales with the tree's age span, not Ntips).
									CSA_age_epsilon			= NULL,		# non-negative numeric (in units of time), age radius around a concentrated sampling attempt, within which to assume that sampling events were due to the concentrated sampling attempt. If NULL, this is chosen automatically based on the anticipated scale of numerical rounding errors. Only relevant if NCSA>0.
									Ntrials					= 1,		# integer, number of fitting trials to run
									max_start_attempts		= 1,		# integer, number of times to attempt finding a valid start point (per trial) before giving up. Randomly choosen start parameters may result in Inf/undefined objective, so this option allows the algorithm to keep looking for valid starting points.
									Nthreads				= 1,		# integer, number of parallel threads to use when performing multiple fitting trials
									max_model_runtime		= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
									Nbootstraps				= 0,		# integer optional number of parametric-bootstrap samples for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
									Ntrials_per_bootstrap	= NULL,		# integer optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
									fit_control				= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
									focal_param_values		= NULL,		# optional list, each element of which is a named list containing lambda, mu, psi, kappa, CSA_probs and CSA_kappas of particular interest and for which to calculate the loglikelihood. Can be used e.g. to explore the shape of the loglikelihood function.
									verbose					= FALSE,	# boolean, specifying whether to print informative messages
									diagnostics				= FALSE,	# boolean, specifying whether to print detailed info (such as log-likelihood) at every iteration of the fitting. For debugging purposes mainly.
									verbose_prefix			= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	if(tree$Nnode<2) return(list(success = FALSE, error="Input tree is too small"));
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree: Expected one of 0,1,2,3.")));
	if(is.null(CSA_ages)){
		NCSA = 0
		CSA_ages = numeric(0)
		if((!is.null(guess_CSA_probs)) && (length(guess_CSA_probs)>0)) return(list(success=FALSE, error=sprintf("No guessed CSA_probs must be provided, since no CSA ages were provided")))
		if((!is.null(guess_CSA_kappas)) && (length(guess_CSA_kappas)>0)) return(list(success=FALSE, error=sprintf("No guessed CSA_kappas must be provided, since no CSA ages were provided")))
		if((!is.null(fixed_CSA_probs)) && (length(fixed_CSA_probs)>0)) return(list(success=FALSE, error=sprintf("No fixed CSA_probs must be provided, since no CSA ages were provided")))
		if((!is.null(fixed_CSA_kappas)) && (length(fixed_CSA_kappas)>0)) return(list(success=FALSE, error=sprintf("No fixed CSA_kappas must be provided, since no CSA ages were provided")))
	}else{
		NCSA = length(CSA_ages)
		if(any(diff(CSA_ages)<=0)) return(list(success=FALSE, error=sprintf("CSA_ages must be in strictly ascending order")))
		if(any(CSA_ages<0)) return(list(success=FALSE, error=sprintf("CSA_ages must be non-negative")))
		if((!is.null(guess_CSA_probs)) && (length(guess_CSA_probs)!=1) && (length(guess_CSA_probs)!=NCSA)) return(list(success = FALSE, error = sprintf("Invalid number of guessed CSA_probs (%d); expected either a single value or %d values",length(guess_CSA_probs),NCSA)));
		if((!is.null(guess_CSA_kappas)) && (length(guess_CSA_kappas)!=1) && (length(guess_CSA_kappas)!=NCSA)) return(list(success = FALSE, error = sprintf("Invalid number of guessed CSA_kappas (%d); expected either a single value or %d values",length(guess_CSA_kappas),NCSA)));
		if((!is.null(fixed_CSA_probs)) && (length(fixed_CSA_probs)!=1) && (length(fixed_CSA_probs)!=NCSA)) return(list(success = FALSE, error = sprintf("Invalid number of fixed CSA_probs (%d); expected either a single value or %d values",length(fixed_CSA_probs),NCSA)));
		if((!is.null(fixed_CSA_kappas)) && (length(fixed_CSA_kappas)!=1) && (length(fixed_CSA_kappas)!=NCSA)) return(list(success = FALSE, error = sprintf("Invalid number of fixed CSA_kappas (%d); expected either a single value or %d values",length(fixed_CSA_kappas),NCSA)));

	}
	tree_span = get_tree_span(tree)$max_distance
	if(is.null(root_age)) root_age = tree_span
	if(is.null(oldest_age)) oldest_age = root_age;
	if(!(condition %in% c("crown","stem","auto","none"))) return(list(success = FALSE, error = sprintf("Invalid condition '%s': Extected 'stem', 'crown', 'none' or 'auto'.",condition)));
	if(condition=="auto") condition = (if(abs(oldest_age-root_age)<=1e-10*root_age) "crown" else "stem")
	if(is.null(max_model_runtime)) max_model_runtime = 0
	if(is.null(Ntrials_per_bootstrap)) Ntrials_per_bootstrap = max(1,Ntrials)
	max_start_attempts 	= (if(is.null(max_start_attempts)) 1 else max(1,max_start_attempts))
	Ntrials  = (if(is.null(Ntrials)) 1 else max(1,Ntrials))
	Nthreads = (if(is.null(Nthreads)) 1 else max(1,Nthreads))
	original_fixed_lambda 	= fixed_lambda
	original_fixed_mu 		= fixed_mu
	original_fixed_psi	 	= fixed_psi
	original_fixed_kappa 	= fixed_kappa
	has_fixed_lambda 		= (!is.null(fixed_lambda)) && any(is.finite(fixed_lambda))
	has_fixed_mu 			= (!is.null(fixed_mu)) && any(is.finite(fixed_mu))
	has_fixed_psi 			= (!is.null(fixed_psi)) && any(is.finite(fixed_psi))
	has_fixed_kappa 		= (!is.null(fixed_kappa)) && any(is.finite(fixed_kappa))

	if(is.null(age_grid) || (length(age_grid)<=1)){
		if((!is.null(guess_lambda)) && (length(guess_lambda)>1)) return(list(success = FALSE, error = sprintf("Invalid number of guessed lambdas; since no age grid was provided, you must provide a single (constant) guess_lambda or none at all")));
		if((!is.null(guess_mu)) && (length(guess_mu)>1)) return(list(success = FALSE, error = sprintf("Invalid number of guessed mus; since no age grid was provided, you must provide a single (constant) guess_mu or none at all")));
		if((!is.null(guess_psi)) && (length(guess_psi)>1)) return(list(success = FALSE, error = sprintf("Invalid number of guessed psis; since no age grid was provided, you must provide a single (constant) guess_psi or none at all")));
		if((!is.null(guess_kappa)) && (length(guess_kappa)>1)) return(list(success = FALSE, error = sprintf("Invalid number of guessed kappas; since no age grid was provided, you must provide a single (constant) guess_kappa or none at all")));
		age_grid = 0 # single-point grid, means that lambda, mu, psi, kappa are assumed time-independent
		NG = 1
	}else{
		NG = length(age_grid)
		if((splines_degree>0) && (!is.null(oldest_age)) && (tail(age_grid,1)<oldest_age)) return(list(success=FALSE, error=sprintf("Provided age grid must cover oldest_age (%g) when splines_degree>0",oldest_age)))
		if(any(diff(age_grid)<=0)) return(list(success=FALSE, error=sprintf("age_grid must list ages in strictly ascending order")))
		if((!is.null(guess_lambda)) && (length(guess_lambda)!=1) && (length(guess_lambda)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of guessed lambdas (%d); since an age grid of size %d was provided, you must either provide one or %d lambdas",length(guess_lambda),NG)));
		if((!is.null(guess_mu)) && (length(guess_mu)!=1) && (length(guess_mu)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of guessed mus (%d); since an age grid of size %d was provided, you must either provide one or %d mus",length(guess_mu),NG)));
		if((!is.null(guess_psi)) && (length(guess_psi)!=1) && (length(guess_psi)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of guessed psis (%d); since an age grid of size %d was provided, you must either provide one or %d psis",length(guess_psi),NG)));
		if((!is.null(guess_kappa)) && (length(guess_kappa)!=1) && (length(guess_mu)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of guessed kappas (%d); since an age grid of size %d was provided, you must either provide one or %d kappas",length(guess_kappa),NG)));
		if((length(age_grid)>1) && (age_grid[NG]>oldest_age-1e-5*(age_grid[NG]-age_grid[NG-1]))) age_grid[NG] = max(age_grid[NG],oldest_age); # if age_grid "almost" covers oldest_age (i.e. up to rounding errors), then fix the remaining difference
		if((length(age_grid)>1) && (age_grid[1]<1e-5*(age_grid[2]-age_grid[1]))) age_grid[1] = min(age_grid[1],0); # if age_grid "almost" covers present-day (i.e. up to rounding errors), then fix the remaining difference
	}
	if(NG==1) splines_degree = 1; # no point in using splines since lambda & mu & psi & kappa are assumed to be time-independent
	if((NCSA>0) && is.null(CSA_age_epsilon)) CSA_age_epsilon = (if(NCSA>1) min(tree_span,min(diff(CSA_ages))) else tree_span)/1000
	if(!is.null(fixed_age_grid)){
		if(fixed_age_grid[1]>tail(fixed_age_grid,1)) return(list(success=FALSE, error="fixed_age_grid must be strictly increasing"))
		if(tail(fixed_age_grid,1)<oldest_age) return(list(success=FALSE, error=sprintf("fixed_age_grid must cover oldest_age (%g)",oldest_age)))
		if(fixed_age_grid[1]>0) return(list(success=FALSE, error=sprintf("fixed_age_grid must cover present-day (age 0)")))
		if(verbose && (length(fixed_age_grid)<NG)) cat(sprintf("%sWARNING: fixed_age_grid should generally be finer than age_grid",verbose_prefix))
		NFG = length(fixed_age_grid)
	}
	
	# reformat shape of input params to an internally standardized format
	if(length(min_lambda)==1) min_lambda 	= rep(min_lambda,times=NG);
	if(length(max_lambda)==1) max_lambda 	= rep(max_lambda,times=NG);
	if(length(min_mu)==1) min_mu 			= rep(min_mu,times=NG);
	if(length(max_mu)==1) max_mu 			= rep(max_mu,times=NG);
	if(length(min_psi)==1) min_psi 			= rep(min_psi,times=NG);
	if(length(max_psi)==1) max_psi 			= rep(max_psi,times=NG);
	if(length(min_kappa)==1) min_kappa 		= rep(min_kappa,times=NG);
	if(length(max_kappa)==1) max_kappa 		= rep(max_kappa,times=NG);
	if(length(min_CSA_kappas)==1) min_CSA_kappas 	= rep(min_CSA_kappas,times=NCSA)
	if(length(max_CSA_kappas)==1) max_CSA_kappas 	= rep(max_CSA_kappas,times=NCSA)
	if(length(min_CSA_probs)==1) min_CSA_probs 		= rep(min_CSA_probs,times=NCSA)
	if(length(max_CSA_probs)==1) max_CSA_probs 		= rep(max_CSA_probs,times=NCSA)
	if(is.null(guess_lambda)){
		guess_lambda = rep(NA,times=NG);
	}else if(length(guess_lambda)==1){
		guess_lambda = rep(guess_lambda,times=NG);
	}
	if(is.null(guess_mu)){
		guess_mu = rep(NA,times=NG);
	}else if(length(guess_mu)==1){
		guess_mu = rep(guess_mu,times=NG);
	}
	if(is.null(guess_psi)){
		guess_psi = rep(NA,times=NG);
	}else if(length(guess_psi)==1){
		guess_psi = rep(guess_psi,times=NG);
	}
	if(is.null(guess_kappa)){
		guess_kappa = rep(NA,times=NG);
	}else if(length(guess_kappa)==1){
		guess_kappa = rep(guess_kappa,times=NG);
	}
	if(is.null(guess_CSA_probs)){
		guess_CSA_probs = rep(NA,times=NCSA);
	}else if(length(guess_CSA_probs)==1){
		guess_CSA_probs = rep(guess_CSA_probs,times=NCSA);
	}
	if(is.null(guess_CSA_kappas)){
		guess_CSA_kappas = rep(NA,times=NCSA);
	}else if(length(guess_CSA_kappas)==1){
		guess_CSA_kappas = rep(guess_CSA_kappas,times=NCSA);
	}
	if(!has_fixed_lambda){
		fixed_lambda = rep(NA,times=NG)
	}else{
		if(is.null(fixed_age_grid)){
			if(length(fixed_lambda)==1){
				fixed_lambda = rep(fixed_lambda,times=NG)
			}else if(length(fixed_lambda)!=NG){
				return(list(success=FALSE, error=sprintf("fixed_lambda has invalid length (%d); expected either length 1 or %d",length(fixed_lambda),NG)))
			}
		}else{
			if(length(fixed_lambda)==1){
				full_fixed_lambda = rep(fixed_lambda,times=NFG)
				fixed_lambda = rep(fixed_lambda,times=NG)
			}else if(length(fixed_lambda)!=NFG){
				return(list(success=FALSE, error=sprintf("fixed_lambda has invalid length (%d); expected either length 1 or %d",length(fixed_lambda),NFG)))
			}else if(any(!is.finite(fixed_lambda))){
				return(list(success=FALSE, error=sprintf("fixed_lambda must be defined and finite on the entire fixed_age_grid, since the latter was provided")))
			}else{
				full_fixed_lambda 	= fixed_lambda
				fixed_lambda 		= approx(x=fixed_age_grid, y=fixed_lambda, xout=age_grid, rule=2)$y
			}
		}
	}
	if(!has_fixed_mu){
		fixed_mu = rep(NA,times=NG)
	}else{
		if(is.null(fixed_age_grid)){
			if(length(fixed_mu)==1){
				fixed_mu = rep(fixed_mu,times=NG)
			}else if(length(fixed_mu)!=NG){
				return(list(success=FALSE, error=sprintf("fixed_mu has invalid length (%d); expected either length 1 or %d",length(fixed_mu),NG)))
			}
		}else{
			if(length(fixed_mu)==1){
				full_fixed_mu = rep(fixed_mu,times=NFG)
				fixed_mu = rep(fixed_mu,times=NG)
			}else if(length(fixed_mu)!=NFG){
				return(list(success=FALSE, error=sprintf("fixed_mu has invalid length (%d); expected either length 1 or %d",length(fixed_mu),NFG)))
			}else if(any(!is.finite(fixed_mu))){
				return(list(success=FALSE, error=sprintf("fixed_mu must be defined and finite on the entire fixed_age_grid, since the latter was provided")))
			}else{
				full_fixed_mu 	= fixed_mu
				fixed_mu 		= approx(x=fixed_age_grid, y=fixed_mu, xout=age_grid, rule=2)$y
			}
		}
	}
	if(!has_fixed_psi){
		fixed_psi = rep(NA,times=NG)
	}else{
		if(is.null(fixed_age_grid)){
			if(length(fixed_psi)==1){
				fixed_psi = rep(fixed_psi,times=NG)
			}else if(length(fixed_psi)!=NG){
				return(list(success=FALSE, error=sprintf("fixed_psi has invalid length (%d); expected either length 1 or %d",length(fixed_psi),NG)))
			}
		}else{
			if(length(fixed_psi)==1){
				full_fixed_psi = rep(fixed_psi,times=NFG)
				fixed_psi = rep(fixed_psi,times=NG)
			}else if(length(fixed_psi)!=NFG){
				return(list(success=FALSE, error=sprintf("fixed_psi has invalid length (%d); expected either length 1 or %d",length(fixed_psi),NFG)))
			}else if(any(!is.finite(fixed_psi))){
				return(list(success=FALSE, error=sprintf("fixed_psi must be defined and finite on the entire fixed_age_grid, since the latter was provided")))
			}else{
				full_fixed_psi 	= fixed_psi
				fixed_psi 		= approx(x=fixed_age_grid, y=fixed_psi, xout=age_grid, rule=2)$y
			}
		}
	}
	if(!has_fixed_kappa){
		fixed_kappa = rep(NA,times=NG)
	}else{
		if(is.null(fixed_age_grid)){
			if(length(fixed_kappa)==1){
				fixed_kappa = rep(fixed_kappa,times=NG)
			}else if(length(fixed_kappa)!=NG){
				return(list(success=FALSE, error=sprintf("fixed_kappa has invalid length (%d); expected either length 1 or %d",length(fixed_kappa),NG)))
			}
		}else{
			if(length(fixed_kappa)==1){
				full_fixed_kappa = rep(fixed_kappa,times=NFG)
				fixed_kappa = rep(fixed_kappa,times=NG)
			}else if(length(fixed_kappa)!=NFG){
				return(list(success=FALSE, error=sprintf("fixed_kappa has invalid length (%d); expected either length 1 or %d",length(fixed_kappa),NFG)))
			}else if(any(!is.finite(fixed_kappa))){
				return(list(success=FALSE, error=sprintf("fixed_kappa must be defined and finite on the entire fixed_age_grid, since the latter was provided")))
			}else{
				full_fixed_kappa 	= fixed_kappa
				fixed_kappa 		= approx(x=fixed_age_grid, y=fixed_kappa, xout=age_grid, rule=2)$y
			}
		}
	}
	if(is.null(fixed_CSA_probs)){
		fixed_CSA_probs = rep(NA,times=NCSA);
	}else if(length(fixed_CSA_probs)==1){
		fixed_CSA_probs = rep(fixed_CSA_probs,times=NCSA);
	}
	if(is.null(fixed_CSA_kappas)){
		fixed_CSA_kappas = rep(NA,times=NCSA);
	}else if(length(fixed_CSA_kappas)==1){
		fixed_CSA_kappas = rep(fixed_CSA_kappas,times=NCSA);
	}

	# pre-compute sampling & branching ages
	LTT			= count_lineages_through_time(tree, Ntimes=log2(length(tree$tip.label)), include_slopes=FALSE)
	LTT$ages	= root_age - LTT$times
	tree_events = extract_HBDS_events_from_tree(tree, root_age = root_age, CSA_ages = CSA_ages, age_epsilon = CSA_age_epsilon)
	if(any((tree_events$concentrated_node_counts+tree_events$concentrated_tip_counts>0) & (!is.na(fixed_CSA_probs)) & (fixed_CSA_probs==0))) return(list(success=FALSE, error=sprintf("Some CSAs include sampled tips and/or nodes, even though their sampling probability is fixed to 0")))
	if(any((tree_events$concentrated_node_counts>0) & (!is.na(fixed_CSA_kappas)) & (fixed_CSA_kappas==0))) return(list(success=FALSE, error=sprintf("Some CSAs include ancestral sampled nodes, even though their kappa is fixed to 0")))
									
	#################################
	# PREPARE PARAMETERS TO BE FITTED

	# guess reasonable start params, if not provided
	default_guess_lambda				= (1/root_age) * sum(1/stats::approx(x=LTT$ages, y=LTT$lineages, xout=tree_events$branching_ages, yleft=0, yright=1)$y, na.rm=TRUE) # a good approximation for the pulled speciation rate is 1/root_age * sum_b 1/LTT(i-th branching_age)
	default_guess_mu 					= default_guess_lambda/10
	default_guess_psi 					= (1/root_age) * sum(1/stats::approx(x=LTT$ages, y=LTT$lineages, xout=tree_events$Ptip_ages, yleft=0, yright=1)$y, na.rm=TRUE) # a good approximation for the pulled sampling rate is 1/root_age * sum_b 1/LTT(i-th Ptip_sampling_age)
	default_guess_kappa					= (if((length(tree_events$Ptip_ages) + length(tree_events$Pnode_ages))==0) 0 else length(tree_events$Pnode_ages)/(length(tree_events$Ptip_ages) + length(tree_events$Pnode_ages)))
	guess_lambda[is.na(guess_lambda)] 	= default_guess_lambda
	guess_mu[is.na(guess_mu)] 			= default_guess_mu
	guess_psi[is.na(guess_psi)] 		= default_guess_psi
	guess_kappa[is.na(guess_kappa)] 	= default_guess_kappa
	if(NCSA>0){
		guess_CSA_probs[is.na(guess_CSA_probs)]		= tree_events$concentrated_tip_counts[is.na(guess_CSA_probs)] / stats::approx(x=LTT$ages, y=LTT$lineages, xout=CSA_ages[is.na(guess_CSA_probs)], yleft=tail(LTT$lineages,1), yright=1)$y
		guess_CSA_kappas[is.na(guess_CSA_kappas)]	= tree_events$concentrated_node_counts[is.na(guess_CSA_kappas)]/(tree_events$concentrated_tip_counts[is.na(guess_CSA_kappas)] + tree_events$concentrated_tip_counts[is.na(guess_CSA_kappas)])
	}
	
	# make sure guessed param values are consistent with fixed param values
	guess_lambda[is.finite(fixed_lambda)] 		= fixed_lambda[is.finite(fixed_lambda)]
	default_guess_mu[is.finite(fixed_mu)] 		= fixed_mu[is.finite(fixed_mu)]
	default_guess_psi[is.finite(fixed_psi)] 	= fixed_psi[is.finite(fixed_psi)]
	default_guess_kappa[is.finite(fixed_kappa)] = fixed_kappa[is.finite(fixed_kappa)]
	if(NCSA>0){
		guess_CSA_probs[is.finite(fixed_CSA_probs)] 	= fixed_CSA_probs[is.finite(fixed_CSA_probs)]
		guess_CSA_kappas[is.finite(fixed_CSA_kappas)] 	= fixed_CSA_kappas[is.finite(fixed_CSA_kappas)]
	}	
		
	# make sure initial guess is within the imposed bounds
	guess_lambda 		= pmin(max_lambda, pmax(min_lambda, guess_lambda))
	guess_mu	 		= pmin(max_mu, pmax(min_mu, guess_mu))
	guess_psi	 		= pmin(max_psi, pmax(min_psi, guess_psi))
	guess_kappa	 		= pmin(max_kappa, pmax(min_kappa, guess_kappa))
	guess_CSA_probs		= pmin(max_CSA_probs, pmax(min_CSA_probs, guess_CSA_probs))
	guess_CSA_kappas	= pmin(max_CSA_kappas, pmax(min_CSA_kappas, guess_CSA_kappas))	
	
	# determine which parameters are to be fitted
	# convention: parameters are indexed as follows: [lambda[], mu[], psi[], kappa[], CSA_probs[], CSA_kappas[]]
	fixed_param_values_flat = c(fixed_lambda, fixed_mu, fixed_psi, fixed_kappa, fixed_CSA_probs, fixed_CSA_kappas) # may contain NAs, corresponding to non-fixed parameters
	fitted_params			= which(is.na(fixed_param_values_flat))
	fixed_params			= which(!is.na(fixed_param_values_flat))
	guess_param_values_flat = c(guess_lambda, guess_mu, guess_psi, guess_kappa, guess_CSA_probs, guess_CSA_kappas) # should contain a valid numeric for each parameter, even if the parameter is fixed
	guess_param_values_flat[fixed_params] = fixed_param_values_flat[fixed_params] # make sure guessed param values are consistent with fixed param values
	min_param_values_flat	= c(min_lambda,min_mu,min_psi,min_kappa,min_CSA_probs,min_CSA_kappas)
	max_param_values_flat	= c(max_lambda,max_mu,max_psi,min_kappa,max_CSA_probs,max_CSA_kappas)

	# determine free (i.e. independent) fitted parameters
	# for example, if lambda is enforced to be time-independent, this reduces the number of free parameters
	# free2fitted[frp] (where frp=1,..,Nfree) will be a list of fitted parameter indices represented by the frp-th free parameter
	# fitted2free[fp] will be the index of the free parameter representing the fp-th fitted parameter
	NP				= 4*NG + 2*NCSA				# total number of parameters
	NFlambda		= sum(is.na(fixed_lambda)) 	# number of non-fixed lambda
	NFmu			= sum(is.na(fixed_mu))		# number of non-fixed mu
	NFpsi			= sum(is.na(fixed_psi))
	NFkappa			= sum(is.na(fixed_kappa))
	NFCSA_probs		= sum(is.na(fixed_CSA_probs))
	NFCSA_kappas	= sum(is.na(fixed_CSA_kappas))
	fitted2free		= (if(NFlambda==0) c() else (if(const_lambda) rep(1,times=NFlambda) else c(1:NFlambda)))
	fitted2free		= c(fitted2free, (if(NFmu==0) c() else (if(length(fitted2free)==0) 0 else tail(fitted2free,1))+(if(const_mu) rep(1,times=NFmu) else c(1:NFmu))))
	fitted2free		= c(fitted2free, (if(NFpsi==0) c() else (if(length(fitted2free)==0) 0 else tail(fitted2free,1))+(if(const_psi) rep(1,times=NFpsi) else c(1:NFpsi))))
	fitted2free		= c(fitted2free, (if(NFkappa==0) c() else (if(length(fitted2free)==0) 0 else tail(fitted2free,1))+(if(const_kappa) rep(1,times=NFkappa) else c(1:NFkappa))))
	fitted2free		= c(fitted2free, (if(NFCSA_probs==0) c() else (if(length(fitted2free)==0) 0 else tail(fitted2free,1))+(if(const_CSA_probs) rep(1,times=NFCSA_probs) else c(1:NFCSA_probs))))
	fitted2free		= c(fitted2free, (if(NFCSA_kappas==0) c() else (if(length(fitted2free)==0) 0 else tail(fitted2free,1))+(if(const_CSA_kappas) rep(1,times=NFCSA_kappas) else c(1:NFCSA_kappas))))
	Nfree			= length(unique(fitted2free)); # number of free (i.e. independently) fitted parameters
	free2fitted		= lapply(1:Nfree, FUN=function(frp) which(fitted2free==frp))
		
	# determine typical parameter scales
	scale_lambda = abs(guess_lambda); scale_lambda[scale_lambda==0] = mean(scale_lambda);
	scale_mu 	 = abs(guess_mu);
	if(all(scale_mu==0)){ 
		scale_mu[] = scale_lambda; 
	}else{ 
		scale_mu[scale_mu==0] = mean(scale_mu); 
	}
	scale_psi = abs(guess_psi)
	if(all(scale_psi==0)){ 
		scale_psi[] = scale_lambda/10; 
	}else{ 
		scale_psi[scale_psi==0] = mean(scale_psi); 
	}
	scale_kappa	= abs(guess_kappa)
	if(all(scale_kappa==0)){ 
		scale_kappa[] = 1
	}else{ 
		scale_kappa[scale_kappa==0] = mean(scale_kappa); 
	}
	scale_CSA_probs = abs(guess_CSA_probs)
	if(all(scale_CSA_probs==0)){ 
		scale_CSA_probs[] = 1
	}else{ 
		scale_CSA_probs[scale_CSA_probs==0] = mean(scale_CSA_probs); 
	}
	scale_CSA_kappas = abs(guess_CSA_kappas)
	if(all(scale_CSA_kappas==0)){ 
		scale_CSA_kappas[] = 1
	}else{ 
		scale_CSA_kappas[scale_CSA_kappas==0] = mean(scale_CSA_kappas); 
	}
	param_scales = c(scale_lambda,scale_mu,scale_psi,scale_kappa,scale_CSA_probs,scale_CSA_kappas)
	
	# define auxiliary function for obtaining full parameter list from rescaled free fitted parameters
	# input: fparam_values[] is a 1D vector of length NFP, listing rescaled values for the free fitted parameters
	# output: param_values[] will be a 1D vector of length NP, listing all model parameter values
	fparam_scales = (if(Nfree==0) numeric(0) else sapply(1:Nfree, FUN = function(frp) mean(param_scales[fitted_params[free2fitted[[frp]]]])))
	expand_free_fitted_params = function(fparam_values){
		fparam_values = fparam_values * fparam_scales;
		param_values = fixed_param_values_flat; 
		param_values[fitted_params] = pmax(min_param_values_flat[fitted_params], pmin(max_param_values_flat[fitted_params], fparam_values[fitted2free]))
		return(param_values)
	}

	# set fit-control options, unless provided by the caller
	if(is.null(fit_control)) fit_control = list()
	if(is.null(fit_control$step.min)) fit_control$step.min = 0.001
	if(is.null(fit_control$x.tol)) fit_control$x.tol = 1e-8
	if(is.null(fit_control$iter.max)) fit_control$iter.max = 1000
	if(is.null(fit_control$eval.max)) fit_control$eval.max = 2 * fit_control$iter.max * Nfree


	################################
	# FITTING
	
	
	# objective function: negated log-likelihood
	# input argument is either the full set of parameters, or the subset of fitted parameters rescaled according to param_scales
	objective_function = function(fparam_values, param_values, trial){
		if(is.null(param_values)){
			param_values_flat = expand_free_fitted_params(fparam_values);
			if(any(is.nan(param_values_flat)) || any(is.infinite(param_values_flat))) return(Inf); # catch weird cases where params become NaN
			param_values = unflatten_hbds_params(param_values_flat, NG, NCSA)
		}
		CSA_probs	= param_values$CSA_probs
		CSA_kappas	= param_values$CSA_kappas
		if(length(age_grid)==1){
			# while age-grid has only one point (i.e., lambda & mu are constant over time), we need to provide a least 2 grid points to the loglikelihood calculator, spanning the interval [0,oldest_age]
			input_age_grid 	= c(0,oldest_age)
			input_lambdas	= c(param_values$lambda, param_values$lambda)
			input_mus		= c(param_values$mu, param_values$mu)
			input_psis		= c(param_values$psi, param_values$psi)
			input_kappas	= c(param_values$kappa, param_values$kappa)
		}else{
			input_age_grid 	= age_grid
			input_lambdas	= param_values$lambda
			input_mus 		= param_values$mu
			input_psis 		= param_values$psi
			input_kappas	= param_values$kappa
			if(tail(age_grid,1)<oldest_age){
				# age_grid does not cover the oldest_age, so extend all profiles as a constant
				input_age_grid 	= c(input_age_grid,oldest_age)
				input_lambdas	= c(input_lambdas, tail(input_lambdas,1))
				input_mus		= c(input_mus, tail(input_mus,1))
				input_psis		= c(input_psis, tail(input_psis,1))
				input_kappas	= c(input_kappas, tail(input_kappas,1))
			}
			if(age_grid[1]>0){
				# age_grid does not cover the present-day, so extend all profiles as a constant
				input_age_grid 	= c(0,input_age_grid)
				input_lambdas	= c(input_lambdas[1], input_lambdas)
				input_mus		= c(input_mus[1], input_mus)
				input_psis		= c(input_psis[1], input_psis)
				input_kappas	= c(input_kappas[1], input_kappas)
			}
		}
		if(!is.null(fixed_age_grid)){
			# use parameter profiles defined on fixed_age_grid (if available) or interpolate onto fixed_age_grid
			input_lambdas 	= (if(has_fixed_lambda) full_fixed_lambda else evaluate_spline(Xgrid=input_age_grid, Ygrid=input_lambdas, splines_degree=splines_degree, Xtarget=fixed_age_grid, extrapolate="const"))
			input_mus 		= (if(has_fixed_mu) full_fixed_mu else evaluate_spline(Xgrid=input_age_grid, Ygrid=input_mus, splines_degree=splines_degree, Xtarget=fixed_age_grid, extrapolate="const"))
			input_psis 		= (if(has_fixed_psi) full_fixed_psi else evaluate_spline(Xgrid=input_age_grid, Ygrid=input_psis, splines_degree=splines_degree, Xtarget=fixed_age_grid, extrapolate="const"))
			input_kappas 	= (if(has_fixed_kappa) full_fixed_kappa else evaluate_spline(Xgrid=input_age_grid, Ygrid=input_kappas, splines_degree=splines_degree, Xtarget=fixed_age_grid, extrapolate="const"))
			input_age_grid	= fixed_age_grid
		}
		results = get_HBDS_model_loglikelihood_CPP(	branching_ages					= tree_events$branching_ages,
													Ptip_ages						= tree_events$Ptip_ages,
													Pnode_ages						= tree_events$Pnode_ages,
													CSA_ages						= CSA_ages,
													CSA_probs						= CSA_probs,
													CSA_kappas						= CSA_kappas,
													concentrated_tip_counts			= tree_events$concentrated_tip_counts,
													concentrated_node_counts		= tree_events$concentrated_node_counts,
													oldest_age						= oldest_age,
													age_grid						= input_age_grid,
													lambdas							= input_lambdas,
													mus								= input_mus,
													psis							= input_psis,
													kappas							= input_kappas,
													splines_degree					= splines_degree,
													condition						= condition,
													relative_ODE_step				= ODE_relative_dt,
													E_value_step					= ODE_relative_dy,
													runtime_out_seconds				= max_model_runtime)
		loglikelihood = if((!results$success) || (!is.finite(results$loglikelihood))) -Inf else results$loglikelihood
		if(diagnostics && (trial>=0)){
			if(results$success){ cat(sprintf("%s  Trial %d: loglikelihood %.10g, model runtime %.5g sec\n",verbose_prefix,trial,loglikelihood,results$runtime)) }
			else{ cat(sprintf("%s  Trial %d: Model evaluation failed: %s\n",verbose_prefix,trial,results$error)) }
		}
		return(-loglikelihood)
	}
	
	# calculate loglikelihood for initial guess
	if(Nfree>0){
		guess_fparam_values = sapply(1:Nfree, FUN = function(frp) guess_param_values_flat[fitted_params[free2fitted[[frp]][1]]])
		guess_loglikelihood = -objective_function(fparam_values=guess_fparam_values/fparam_scales, param_values=NULL, trial=0)
	}else{
		guess_fparam_values = numeric(0)
		guess_loglikelihood = -objective_function(fparam_values=numeric(0), param_values=NULL, trial=0)
	}
	
	# calculate loglikelihood for focal param values
	if((!is.null(focal_param_values)) && (length(focal_param_values)>0)){
		if(verbose) cat(sprintf("%sComputing loglikelihoods for focal param values..\n",verbose_prefix))
		focal_loglikelihoods = numeric(length(focal_param_values))
		for(f in 1:length(focal_loglikelihoods)){
			if(is.null(focal_param_values[[f]]$lambda)) return(list(success=FALSE, error=sprintf("Badly formatted input focal_param_values[%d]: Expected a named list containing lambda, mu, psi, kappa, CSA_probs and CSA_kappas",f)))
			focal_loglikelihoods[f] = -objective_function(fparam_values=NULL, param_values=focal_param_values[[f]], trial=-1)
		}
	}else{
		focal_loglikelihoods = NULL
	}
	
	
	# fit with various starting points
	fit_single_trial = function(trial){
		# determine bounds & start params
		lower_bounds = sapply(1:Nfree, FUN = function(frp) max(min_param_values_flat[fitted_params[free2fitted[[frp]]]]))
		upper_bounds = sapply(1:Nfree, FUN = function(frp) min(max_param_values_flat[fitted_params[free2fitted[[frp]]]]))
		start_values = sapply(1:Nfree, FUN = function(frp) guess_param_values_flat[fitted_params[free2fitted[[frp]][1]]])
		Nstart_attempts = 0
		while(Nstart_attempts<max_start_attempts){
			# randomly choose start values for fitted params
			if(trial==1){
				start_values = guess_fparam_values
			}else{
				start_values = get_random_params(	defaults=guess_fparam_values,
													lower_bounds=lower_bounds, 
													upper_bounds=upper_bounds, 
													scales=fparam_scales, 
													orders_of_magnitude=4)
			}
			# check if start values yield NaN
			start_objective = objective_function(fparam_values=start_values/fparam_scales, param_values=NULL, trial=trial);
			Nstart_attempts = Nstart_attempts + 1
			if(is.finite(start_objective)) break;
		}
		# run fit
		if(is.finite(start_objective)){
			fit = tryCatch({ stats::nlminb(start_values/fparam_scales, 
								objective	= function(pars) objective_function(fparam_values=pars, param_values=NULL, trial=trial), 
								lower		= lower_bounds/fparam_scales, 
								upper		= upper_bounds/fparam_scales, 
								control		= fit_control)
							}, error = function(e){ list(objective=NA, par=NA, convergence=1, iterations=NA, evaluations=NA) })
			results = list(objective_value=fit$objective, fparam_values = fit$par, converged=(fit$convergence==0), Niterations=fit$iterations, Nevaluations=fit$evaluations[1], Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective)
			if(diagnostics){
				if(is.finite(fit$objective)){
					cat(sprintf("%s  Trial %d: Final loglikelihood %.10g, Niterations %d, Nevaluations %d, converged = %d\n",verbose_prefix,trial,-results$objective_value, results$Niterations, results$Nevaluations, results$converged))
				}else{
					cat(sprintf("%s  Trial %d: Fitting failed for unknown reason\n",verbose_prefix,trial))
				}
			}
		}else{
			results = list(objective_value=NA, fparam_values = NA, converged=FALSE, Niterations=0, Nevaluations=0, Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective)
			if(diagnostics) cat(sprintf("%s  Trial %d: Start objective is non-finite. Skipping trial\n",verbose_prefix,trial))
		}
		return(results)
	}
	
	################################

	# run one or more independent fitting trials
	if(Nfree>0){
		if((Ntrials>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
			# run trials in parallel using multiple forks
			# Note: Forks (and hence shared memory) are not available on Windows
			if(verbose) cat(sprintf("%sFitting %d free parameters (%d trials, parallelized)..\n",verbose_prefix,Nfree,Ntrials))
			fits = parallel::mclapply(	1:Ntrials, 
										FUN = function(trial) fit_single_trial(trial), 
										mc.cores = min(Nthreads, Ntrials), 
										mc.preschedule = FALSE, 
										mc.cleanup = TRUE);
		}else{
			# run in serial mode
			if(verbose) cat(sprintf("%sFitting %d free parameters (%s)..\n",verbose_prefix,Nfree,(if(Ntrials==1) "1 trial" else sprintf("%d trials",Ntrials))))
			fits = sapply(1:Ntrials,function(x) NULL)
			for(trial in 1:Ntrials){
				fits[[trial]] = fit_single_trial(trial)
			}
		}
	
		# extract information from best fit (note that some fits may have LL=NaN or NA)
		objective_values		 = unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$objective_value))
		valids					 = which((!is.na(objective_values)) & (!is.nan(objective_values)) & (!is.null(objective_values)) & (!is.infinite(objective_values)));
		if(length(valids)==0) return(list(success=FALSE, error=sprintf("Fitting failed for all trials")));
		best_fit				 = fits[[valids[which.min(sapply(valids, function(i) objective_values[i]))]]]
		objective_value			 = -best_fit$objective_value
		loglikelihood			 = objective_value
		fitted_param_values_flat = expand_free_fitted_params(best_fit$fparam_values)
		if(is.null(objective_value) || any(is.na(fitted_param_values_flat)) || any(is.nan(fitted_param_values_flat))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));
	}else{
		fitted_param_values_flat = fixed_param_values_flat
		objective_value 		 = guess_loglikelihood
		loglikelihood			 = guess_loglikelihood
		best_fit 				 = list(converged=TRUE, Niterations = 0, Nevaluations = 0, Nstart_attempts = 0, start_objective=loglikelihood)
	}
	fitted_param_values	= unflatten_hbds_params(fitted_param_values_flat, NG, NCSA)
	guess_param_values	= unflatten_hbds_params(guess_param_values_flat, NG, NCSA)

	# determine number of data points used
	Ndata = sum(tree_events$branching_ages<=oldest_age) + sum(tree_events$Ptip_ages<=oldest_age) + sum(tree_events$Pnode_ages<=oldest_age) + sum(tree_events$concentrated_tip_counts[CSA_ages<=oldest_age]) + sum(tree_events$concentrated_node_counts[CSA_ages<=oldest_age])
	
	#######################################################################
	# estimate confidence intervals if needed, via parametric bootstrapping
	
	if(Nfree==0) Nbootstraps = 0
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sEstimating confidence intervals using %d parametric bootstraps..\n",verbose_prefix,Nbootstraps))
		bootstrap_LLs	 		= rep(NA,times=Nbootstraps)
		bootstrap_params_flat 	= matrix(NA,nrow=Nbootstraps,ncol=NP)
		for(b in 1:Nbootstraps){
			# simulate model with fitted parameters
			if(verbose) cat(sprintf("%s  Bootstrap #%d..\n",verbose_prefix,b))
			bootstrap = castor::generate_tree_hbds(	max_time			= root_age,
													include_extant		= FALSE,
													include_extinct		= FALSE,
													time_grid			= rev(root_age - age_grid),
													lambda				= rev(fitted_param_values$lambda),
													mu					= rev(fitted_param_values$mu),
													psi					= rev(fitted_param_values$psi),
													kappa				= rev(fitted_param_values$kappa),
													splines_degree		= splines_degree,
													CSA_times			= rev(root_age - CSA_ages),
													CSA_probs			= rev(fitted_param_values$CSA_probs),
													CSA_kappas			= rev(fitted_param_values$CSA_kappas),
													no_full_extinction	= TRUE)
			
			if(!bootstrap$success){
				if(verbose) cat(sprintf("%s    WARNING: Bootstrap #%d failed: %s\n",verbose_prefix,b,bootstrap$error))
				next
			}
			# fit model to simulated tree
			fit = fit_hbds_model_on_grid(	tree					= bootstrap$tree, 
											root_age 				= bootstrap$root_age,
											oldest_age				= oldest_age,
											age_grid				= age_grid,
											CSA_ages				= CSA_ages,
											min_lambda				= min_lambda,
											max_lambda				= max_lambda,
											min_mu					= min_mu,
											max_mu					= max_mu,
											min_psi					= min_psi,
											max_psi					= max_psi,
											min_kappa				= min_kappa,
											max_kappa				= max_kappa,
											min_CSA_probs			= min_CSA_probs,
											max_CSA_probs			= max_CSA_probs,
											min_CSA_kappas			= min_CSA_kappas,
											max_CSA_kappas			= max_CSA_kappas,
											guess_lambda			= guess_lambda,
											guess_mu				= guess_mu,
											guess_psi				= guess_psi,
											guess_kappa				= guess_kappa,
											guess_CSA_probs			= guess_CSA_probs,
											guess_CSA_kappas		= guess_CSA_kappas,
											fixed_lambda			= original_fixed_lambda,
											fixed_mu				= original_fixed_mu,
											fixed_psi				= original_fixed_psi,
											fixed_kappa				= original_fixed_kappa,
											fixed_CSA_probs			= fixed_CSA_probs,
											fixed_CSA_kappas		= fixed_CSA_kappas,
											fixed_age_grid			= fixed_age_grid,
											const_lambda			= const_lambda,
											const_mu				= const_mu,
											const_psi				= const_psi,
											const_kappa				= const_kappa,
											const_CSA_probs			= const_CSA_probs,
											const_CSA_kappas		= const_CSA_kappas,
											splines_degree			= splines_degree,
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
											focal_param_values		= list(fitted_param_values),
											verbose					= verbose,
											diagnostics				= diagnostics,
											verbose_prefix			= paste0(verbose_prefix,"    "))
			if(!fit$success){
				if(verbose) cat(sprintf("%s  WARNING: Fitting failed for this bootstrap: %s\n",verbose_prefix,fit$error))
			}else{
				bootstrap_params_flat[b,] = c(fit$param_fitted$lambda, fit$param_fitted$mu, fit$param_fitted$psi, fit$param_fitted$kappa, fit$param_fitted$CSA_probs, fit$param_fitted$CSA_kappas)
				bootstrap_LLs[b]	 	  = fit$focal_loglikelihoods
			}
		}
		# calculate standard errors and confidence intervals from distribution of bootstrapped parameters
		standard_errors = unflatten_hbds_params(sqrt(pmax(0, colMeans(bootstrap_params_flat^2, na.rm=TRUE) - colMeans(bootstrap_params_flat, na.rm=TRUE)^2)), NG, NCSA)
		quantiles_flat 	= sapply(1:ncol(bootstrap_params_flat), FUN=function(p) quantile(bootstrap_params_flat[,p], probs=c(0.25, 0.75, 0.025, 0.975, 0.5), na.rm=TRUE, type=8))
		CI50lower 		= unflatten_hbds_params(quantiles_flat[1,], NG, NCSA)
		CI50upper 		= unflatten_hbds_params(quantiles_flat[2,], NG, NCSA)
		CI95lower 		= unflatten_hbds_params(quantiles_flat[3,], NG, NCSA)
		CI95upper 		= unflatten_hbds_params(quantiles_flat[4,], NG, NCSA)
		medians 		= unflatten_hbds_params(quantiles_flat[5,], NG, NCSA)
		mean_BLL		= mean(bootstrap_LLs,na.rm=TRUE)
		consistency 	= sum(abs(bootstrap_LLs-mean_BLL)>=abs(loglikelihood-mean_BLL), na.rm=TRUE)/sum(!is.nan(bootstrap_LLs))
	}
		
	# return results
	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				guess_loglikelihood		= guess_loglikelihood,
				param_guess				= guess_param_values,
				param_fitted			= fitted_param_values,
				age_grid				= age_grid,
				CSA_ages				= CSA_ages,
				NFP						= Nfree,
				Ndata					= Ndata,
				AIC						= 2*Nfree - 2*loglikelihood,
				BIC						= log(Ndata)*Nfree - 2*loglikelihood,
				condition				= condition,
				converged				= best_fit$converged,
				Niterations				= best_fit$Niterations,
				Nevaluations			= best_fit$Nevaluations,
				focal_loglikelihoods	= focal_loglikelihoods,
				standard_errors			= (if(Nbootstraps>0) standard_errors else NULL),
				CI50lower				= (if(Nbootstraps>0) CI50lower else NULL),
				CI50upper				= (if(Nbootstraps>0) CI50upper else NULL),
				CI95lower				= (if(Nbootstraps>0) CI95lower else NULL),
				CI95upper				= (if(Nbootstraps>0) CI95upper else NULL),
				consistency				= (if(Nbootstraps>0) consistency else NULL)))
}


unflatten_hbds_params = function(params, NG, NCSA){
	lambda		= params[1:NG]
	mu			= params[(NG+1):(2*NG)]
	psi			= params[(2*NG+1):(3*NG)]
	kappa		= params[(3*NG+1):(4*NG)]
	CSA_probs	= (if(NCSA==0) numeric(0) else params[(4*NG+1):(4*NG+NCSA)])
	CSA_kappas	= (if(NCSA==0) numeric(0) else params[(4*NG+NCSA+1):(4*NG+2*NCSA)])
	return(list(lambda=lambda, mu=mu, psi=psi, kappa=kappa, CSA_probs=CSA_probs, CSA_kappas=CSA_kappas))
}

