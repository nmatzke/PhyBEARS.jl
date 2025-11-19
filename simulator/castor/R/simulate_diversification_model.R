# Simulate a deterministic speciation/extinction model
# The model is parameterized similarly to generate_random_tree(), but the returned value is a time series of diversity (number of clades) over time
simulate_diversification_model = function(	times,									# times at which to calculate diversities, in ascending order
											parameters				= list(), 		# named list of model parameters. For entries and default values see the main function body below
											added_rates_times		= NULL,			# numeric vector of size NAR, or empty or NULL
											added_birth_rates_pc	= NULL,			# numeric vector of size NAR, or empty or NULL
											added_death_rates_pc 	= NULL,			# numeric vector of size NAR, or empty or NULL
											added_periodic			= FALSE,		# (logical) if TRUE, added pc birth & death rates are extended periodically if needed. If FALSE, they are extended with zeros.
											start_time				= NULL,			# Beginning (earliest) time of the simulation (<=times[1]). If NULL, this is set to times[1]. If reverse==FALSE, the tree is simulated forwards from this time point.
											final_time				= NULL,			# Final (latest) time of the simulation (>=times[NT]). If NULL, this is set to times[NT]. If reverse==TRUE, the tree is simulated backwards from this time point.
											start_diversity			= 1,			# diversity of extant clades at start_time. Only relevant if reverse==FALSE.
											final_diversity			= NULL,			# diversity of extant clades at final_time (prior to any rarefaction or collapsing at some resolution). Only relevant if reverse==TRUE.
											reverse					= FALSE,		# (boolean) if true, then the tree model is integrated in backward time direction. In that case, start_diversity is interpreted as the true diversity at times.back()
											include_coalescent		= FALSE,		# (boolean) if true, the coalescent diversities are also calculated (in addition to the total diversities)
											include_event_rates		= FALSE,		# (boolean) include birth & death rates in returned values
											include_Nevents			= FALSE,		# (boolean) include an estimate of the total birth (speciation) and death (extinction) events during the simulation
											max_runtime				= NULL){		# (numeric) max allowed runtime in seconds. If NULL or <=0, this option is ignored
	NT = length(times)
	
	# basic error checking
	if((!is.null(start_time)) && (!reverse) && (start_time>times[1])) stop(sprintf("start_time must be equal to or smaller than the first requested time point (got start_time=%g, first requested time point = %g)",start_time,times[1]))
	if((!is.null(final_time)) && reverse && (final_time<times[NT])) stop(sprintf("final_time must be equal to or larger than the last requested time point (got final_time=%.10g, last requested time point = %.10g)",final_time,times[NT]))
	if(is.null(start_time)) start_time = times[1];
	if(is.null(final_time)) final_time = times[NT];
	if(reverse && is.null(final_diversity)) stop(sprintf("Missing final_diversity, needed for simulating tree in reverse"))
	if((!reverse) && is.null(start_diversity)) stop(sprintf("Missing start_diversity, needed for simulating tree forwards"))
	if(is.null(start_diversity)) start_diversity = NaN; # not relevant, but assign a dummy value needed for the CPP function
	if(is.null(final_diversity)) final_diversity = NaN; # not relevant, but assign a dummy value needed for the CPP function
							
	# set default model parameters
	if(is.null(parameters$birth_rate_intercept)) 	parameters$birth_rate_intercept = 0;
	if(is.null(parameters$birth_rate_factor)) 		parameters$birth_rate_factor = 0;
	if(is.null(parameters$birth_rate_exponent)) 	parameters$birth_rate_exponent = 1;
	if(is.null(parameters$death_rate_intercept)) 	parameters$death_rate_intercept = 0;
	if(is.null(parameters$death_rate_factor))		parameters$death_rate_factor = 0;
	if(is.null(parameters$death_rate_exponent)) 	parameters$death_rate_exponent = 1;
	if(is.null(parameters$resolution)) 				parameters$resolution = 0;
	if(is.null(parameters$rarefaction)) 			parameters$rarefaction = 1;
						
	# run simulation
	simulation = simulate_deterministic_diversity_growth_CPP(	birth_rate_intercept 		= parameters$birth_rate_intercept,
																birth_rate_factor 			= parameters$birth_rate_factor,
																birth_rate_exponent 		= parameters$birth_rate_exponent,
																death_rate_intercept 		= parameters$death_rate_intercept,
																death_rate_factor 			= parameters$death_rate_factor,
																death_rate_exponent 		= parameters$death_rate_exponent,
																resolution					= parameters$resolution,
																rarefaction					= parameters$rarefaction,
																Nsplits 					= 2,
																additional_rates_times		= (if(is.null(added_rates_times)) numeric() else added_rates_times),
																additional_birth_rates_pc	= (if(is.null(added_birth_rates_pc)) numeric() else added_birth_rates_pc),
																additional_death_rates_pc	= (if(is.null(added_death_rates_pc)) numeric() else added_death_rates_pc),
																additional_periodic			= added_periodic,
																times 						= times,
																start_time					= start_time,
																final_time					= final_time,
																start_diversity				= start_diversity,
																final_diversity				= final_diversity,
																reverse						= reverse,
																include_coalescent			= include_coalescent,
																include_probabilities		= TRUE,
																include_birth_rates			= include_event_rates,
																include_death_rates			= include_event_rates,
																include_Nevents				= include_Nevents,
																runtime_out_seconds			= (if(is.null(max_runtime)) 0 else max_runtime));

	# get probabilities of survival/discovery/representation
	# note that the simulation may return empty vectors (e.g. if this calculation was not requested) or NULL (calculation not available)
	Psurvival 		= simulation$Psurvival	 # probability of survival to the present (regardless of discovery)
	Pdiscovery		= simulation$Pdiscovery
	Prepresentation	= simulation$Prepresentation # probability of survival to the present and discovery
	if((!is.null(Psurvival)) && (length(Psurvival)==0)) Psurvival = NULL
	if((!is.null(Pdiscovery)) && (length(Pdiscovery)==0)) Pdiscovery = NULL
	if((!is.null(Prepresentation)) && (length(Prepresentation)==0)) Prepresentation = NULL
	if(is.null(Psurvival) && (!is.null(Prepresentation)) && (!is.null(Pdiscovery))){
		Psurvival = Prepresentation/Pdiscovery
	}else if(is.null(Prepresentation) && (!is.null(Psurvival)) && (!is.null(Pdiscovery))){
		Prepresentation = Psurvival * Pdiscovery
	}else if(is.null(Pdiscovery) && (!is.null(Psurvival)) && (!is.null(Prepresentation))){
		Pdiscovery = Prepresentation/Psurvival
	}

	return(list(success					= simulation$success,
				coalescent_diversities	= (if(include_coalescent) simulation$coalescent_diversities else NULL),
				total_diversities		= simulation$total_diversities,
				Psurvival				= Psurvival,
				Pdiscovery				= Pdiscovery,
				Prepresentation			= Prepresentation,
				birth_rates				= (if(include_event_rates) simulation$birth_rates else NULL),
				death_rates				= (if(include_event_rates) simulation$death_rates else NULL),
				Nbirths					= simulation$Nbirths,
				Ndeaths					= simulation$Ndeaths,
				error					= (if(!simulation$success) simulation$error else NULL)));
}
