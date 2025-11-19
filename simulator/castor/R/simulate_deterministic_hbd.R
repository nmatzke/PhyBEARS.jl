# Predict various deterministic features of a homogenous birth-death cladogenic model, backward in time.
# The speciation rate lambda and extinction rate mu are specified on a discrete age-grid, and assumed to vary linearly (or polynomially, as splines) between grid points (see "degree" argument).
# This function calculates, among others, the following features over time:
#	Deterministic LTT curve
#	Deterministic total diversity
#	Deterministic shadow diversity
#   Pulled diversification rate (PDR)
# Alternatively, the lambda may be omitted and instead the PDR may be provided together with only the present birth_rate (lambda0).
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
#
simulate_deterministic_hbd = function(	LTT0, 						# number of extant species represented in the tree at age0, i.e. after rarefaction. This is equal to the value of the LTT at anchor_age.
										oldest_age,					# numeric, specifying how far back (time before present) to simulate the model
										age0			= 0,		# numeric, specifying the age (time before present) at which LTT0, lambda0 and rho are specified. Typically this is 0 (i.e., present-day), but may also be somewhere in the past (>0)
										rho0			= 1,		# numeric within (0,1], specifying the fraction of extant diversity represented in the tree at age0. Can also be NULL, which is equivalent to setting rarefaction=1.
										age_grid		= NULL,		# either NULL, or empty, or a numeric vector of size NG, listing ages in ascending order, on which birth/mu are specified. If NULL or empty, then lambda and mu mut be a single scalar. The returned time series will be defined on an age-grid that may be finer than this grid. If of size >=2, the age_grid must cover oldest_age and age0.
										lambda			= NULL,		# either NULL, or a single numeric (constant speciation rate over time), or a numeric vector of size NG (listing speciation rates at each age in grid_ages[]).
										mu				= NULL,		# either a single numeric (constant extinction rate over time), or a numeric vector of size NG (listing extinction rates at each age in grid_ages[]), or NULL (in which case mu_over_lambda must be provided).
										mu_over_lambda	= NULL,		# ratio mu/lambda. Either a single numeric (constant over time), or a numeric vector of size NG (listing mu/lambda at each age in grid_ages[]), or NULL (in which case mu must be provided).
										PDR				= NULL,		# either NULL, or a single numeric (constant PDR over time), or a numeric vector of size NG (listing PDR at each age in grid_ages[]). Only needed if lambda is NULL.
										lambda0			= NULL,		# either NULL, or a single numeric specifying the speciation rate at age0. Only needed if lambda is NULL.
										splines_degree	= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by lambda, mu and PDR on the age grid.
										relative_dt		= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.	
										allow_unreal	= FALSE){	# logical, specifying whether BD models with unrealistic parameters (e.g., negative mu or negative Pmissing) should be supported. This may be desired e.g. when examining model congruence classes with negative mu.
	# check validity of input variables
	if(is.null(rho0)) rho0 = 1;
	if(is.null(mu) && is.null(mu_over_lambda)) return(list(success = FALSE, error = sprintf("Missing either mu or mu_over_lambda")))
	if(!(is.null(mu) || is.null(mu_over_lambda))) return(list(success = FALSE, error = sprintf("Expected either mu or mu_over_lambda, but not both")))
	if(is.null(lambda)){
		if(is.null(PDR)) return(list(success = FALSE, error = sprintf("PDR must be provided when lambda are omitted")))
		if(is.null(lambda0)) return(list(success = FALSE, error = sprintf("lambda0 must be provided when lambda are omitted")))
	}else{
		if(!is.null(lambda0)) return(list(success = FALSE, error = sprintf("lambda0 must not be explicitly provided when lambda are provided (due to potential ambiguity)")))
		if(!is.null(PDR)) return(list(success = FALSE, error = sprintf("PDR must not be explicitly provided when lambda are provided (due to potential ambiguity)")))
	}
	if(is.null(age_grid) || (length(age_grid)<=1)){
		if((!is.null(lambda)) && (length(lambda)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of lambda values; since no age grid was provided, you must either provide a single (constant) lambda or none")))
		if((!is.null(mu)) && (length(mu)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of mu values; since no age grid was provided, you must provide a single (constant) mu")))
		if((!is.null(mu_over_lambda)) && (length(mu_over_lambda)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of mu_over_lambda; since no age grid was provided, you must provide a single (constant) mu_over_lambda")))
		# create dummy age grid
		NG 			= 2;
		age_grid	= seq(from=0,to=1.01*oldest_age,length.out=NG)
		if(!is.null(lambda)) lambda = rep(lambda,times=NG);
		if(!is.null(PDR)) PDR = rep(PDR,times=NG);
		if(!is.null(mu)) mu = rep(mu,times=NG);
		if(!is.null(mu_over_lambda)) mu_over_lambda = rep(mu_over_lambda,times=NG);
	}else{
		NG = length(age_grid);
		if(age_grid[1]>tail(age_grid,1))return(list(success = FALSE, error = sprintf("Values in age_grid must be strictly increasing")))
		if((age_grid[1]>oldest_age) || (age_grid[NG]<oldest_age)) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including oldest_age (%g)",oldest_age)))
		if((age_grid[1]>age0) || (age_grid[NG]<age0)) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including age0 (%g)",age0)))
		if((!is.null(lambda)) && (length(lambda)!=1) && (length(lambda)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of lambda; since an age grid of size %d was provided, you must either provide zero, one or %d lambda",NG,NG)))
		if((!is.null(mu)) && (length(mu)!=1) && (length(mu)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of mu; since an age grid of size %d was provided, you must either provide one or %d mu",NG,NG)))
		if((!is.null(mu_over_lambda)) && (length(mu_over_lambda)!=1) && (length(mu_over_lambda)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of mu_over_lambda; since an age grid of size %d was provided, you must either provide one or %d mu_over_lambda",NG,NG)))
		if((!is.null(lambda)) && (length(lambda)==1)) lambda = rep(lambda,times=NG);
		if((!is.null(PDR)) && (length(PDR)==1)) PDR = rep(PDR,times=NG);
		if((!is.null(mu)) && (length(mu)==1)) mu = rep(mu,times=NG);
		if((!is.null(mu_over_lambda)) && (length(mu_over_lambda)==1)) mu_over_lambda = rep(mu_over_lambda,times=NG);
	}
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
		
	# simulate model backward in time
	census_age = age_grid[1]
	simulation = simulate_deterministic_HBD_model_CPP(	census_age			= census_age,
														oldest_age			= oldest_age,
														age_grid 			= age_grid,
														lambdas 			= (if(is.null(lambda)) numeric() else lambda),
														mus 				= (if(is.null(mu)) numeric() else mu),
														mu_over_lambda		= (if(is.null(mu_over_lambda)) numeric() else mu_over_lambda),
														PDRs	 			= (if(is.null(PDR)) numeric() else PDR),
														anchor_age			= age0,
														anchor_rho	 		= rho0,
														anchor_lambda		= (if(is.null(lambda0)) NaN else lambda0),
														anchor_LTT 			= LTT0,
														splines_degree		= splines_degree,
														relative_dt			= relative_dt,
														allow_unreal		= allow_unreal)
	if(!simulation$success) return(list(success = FALSE, error = sprintf("Could not simulate model: %s",simulation$error)))
	rholambda0 = simulation$anchor_rho*simulation$anchor_lambda;
	census_rholambda = simulation$census_rho*simulation$census_lambda;

	return(list(success							= TRUE,
				ages							= simulation$refined_age_grid, # potentially refined ages grid, on which all returned variables are defined
				total_diversity					= simulation$total_diversity,
				shadow_diversity				= simulation$shadow_diversity, # shadow diversity, defined w.r.t. census_age
				Pmissing						= simulation$Pmissing,
				Pextinct						= simulation$Pextinct,
				LTT								= simulation$LTT,
				lambda							= simulation$lambda,
				mu								= simulation$mu,
				diversification_rate			= simulation$diversification_rate,
				PDR								= simulation$PDR,				# pulled diversification rate, defined w.r.t. census_age
				PND								= simulation$PND,				# pulled normalized diversity, defined w.r.t. census_age
				SER								= census_rholambda - simulation$PDR, 	# shadow extinction rate, defined w.r.t. census_age
				PER								= simulation$census_lambda - simulation$PDR, # pulled extinction rate, defined w.r.t. census_age
				PSR								= simulation$lambda * (1-simulation$Pmissing), # pulled speciation rate, defined w.r.t. census_age
				rholambda0						= rholambda0)); # product rho*lambda at age0
}

