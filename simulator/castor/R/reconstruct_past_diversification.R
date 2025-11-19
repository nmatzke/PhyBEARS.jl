# Estimate past diversity, birth & death rates, based on a time series of diversities (coalescent or not)
# This is a non-parametric algorithm that assumes knowledge (or alternatively, constancy) of the birth rates at all time.
# Note: This function is currently only implemented for bifurcating trees. 
# Input:
#	A time series of diversities. If coalescent==TRUE, these should be the diversities visible in a coalescent phylogenetic tree (after rarefaction).
# 	Corresponding assumed per-capita birth rates. Alternatively, these can be assumed to be constant and be estimated directly from the time series.
# Output:
#   The estimated true past diversities (N(tau))
#   The estimated corresponding death (extinction) rates (delta(tau))
#   The probability of a size-1 clade surviving from each time point to the present (P(tau)), including after rarefaction (if applicable)
#   The estimated total number of speciation & extinction events, over the considered age interval
reconstruct_past_diversification = function(times, 						# numeric vector of size Ntimes, listing times in ascending order
											diversities,				# numeric vector of size Ntimes, listing diversities (coalescent or not) at each time point. If coalescent==TRUE, then these should correspond to the diversities visible after rarefaction.
											birth_rates_pc				= NULL,	# 1D numeric array of size Ntimes, listing known or assumed per-capita birth rates. Can also be of size 1, in which case the same per-capita birth rate is assumed throughout. Can also be empty, in which case a constant per-capita birth rate is assumed and estimated from the last slope of the coalescent_diversity curve.
											rarefaction					= NULL,	# (numeric) Optional rarefaction fraction. If coalescent==TRUE, this rarefaction is assumed to have been applied to the final coalescent diversities. Otherwise it is merely used to calculate survival_chances. Set to 1 to assume no rarefaction was performed.
											discovery_fractions			= NULL,	# 1D numeric array of size Ntimes, listing the fractions of discovered lineages over time. Can be used as an alternative to rarefaction, for example if discovery of extant species is non-random or biased
											discovery_fraction_slopes	= NULL,	# 1D numeric array of size Ntimes, listing the 1st derivative of discovery_fractions (w.r.t. time) over time. If NULL and needed, will be estimated from discovery_fractions via basic finite differences
											max_age						= NULL, # (numeric) Optional max distance from the end time to be considered. If NULL or <=0 or Inf, all provided times are considered.
											coalescent					= FALSE,
											smoothing_span				= 0,	# (integer) Optional sliding window size (number of time points) for smoothening the diversities time series via Savitzky-Golay-filter. If <=2, no smoothing is done. Smoothening the coalescent diversity can reduce the noise in the non-parametric reconstruction. 
											smoothing_order				= 1){	# (integer) Polynomial order of the Savitzky-Golay smoothing filter.
											
	Ntimes = length(times)
	max_time = times[Ntimes]

	# basic error checking
	if(Ntimes<2) return(list(success = FALSE, error="Time series is too small"))
	if(!(smoothing_order %in% c(1,2,3,4))) stop(sprintf("ERROR: smoothing_order must be an integer between 1 and 4 (got %d)",smoothing_order))
	if((!coalescent) && (is.null(birth_rates_pc) || length(birth_rates_pc)==0)) stop("ERROR: For non-coalescent diversity time series birth_rates_pc must be explicitly provided (either for all time points, or as a single constant value)")
	if(is.null(rarefaction) && is.null(discovery_fractions)) rarefaction = 1
	if((!is.null(rarefaction)) && (!is.null(discovery_fractions))) stop("ERROR: Either rarefaction or discovery_fractions must be NULL")
	if(is.null(max_age) || (max_age<0) || (max_age==Inf)) max_age = 0;
	
	if((!is.null(discovery_fractions)) && is.null(discovery_fraction_slopes)){
		# estimate slopes of discovery_fractions via finite differences
		discovery_fraction_slopes = diff(discovery_fractions)/diff(times);
		discovery_fraction_slopes = c(discovery_fraction_slopes[1],discovery_fraction_slopes) # add one more identical point to the end to keep the same size as the original time series
	}

	if(coalescent && is.null(discovery_fractions)){
		reconstruction = reconstruct_past_diversity_from_coalescent_CPP(times,
																		raw_coalescent_diversities	= diversities,
																		birth_rates_pc	= (if(is.null(birth_rates_pc)) numeric() else birth_rates_pc),
																		rarefaction		= rarefaction,
																		max_age			= max_age,
																		smoothing_span	= smoothing_span,
																		smoothing_order	= smoothing_order);
	}else if(coalescent && (!is.null(discovery_fractions))){
		reconstruction = reconstruct_past_diversity_from_biased_coalescent_CPP(	times,
																				raw_coalescent_diversities	= diversities,
																				birth_rates_pc				= (if(is.null(birth_rates_pc)) numeric() else birth_rates_pc),
																				discovery_fractions			= discovery_fractions,
																				discovery_fraction_slopes	= discovery_fraction_slopes,
																				max_age						= max_age,
																				smoothing_span				= smoothing_span,
																				smoothing_order				= smoothing_order);
	}else{
		reconstruction = reconstruct_past_diversifications_CPP(	times,
																raw_diversities = diversities,
																birth_rates_pc	= (if(is.null(birth_rates_pc)) numeric() else birth_rates_pc),
																rarefaction		= rarefaction, # only used to calculate Prepresentation
																Nsplits			= 2,
																max_age			= max_age,
																smoothing_span	= smoothing_span,
																smoothing_order	= smoothing_order)
	}
	if(!reconstruction$success) return(list(success=FALSE, error=reconstruction$error))
	
	# get probabilities of survival/discovery/representation
	# note that some C++ routines may return empty vectors (e.g. if this calculation was not requested) or NULL (calculation not available)
	Psurvival 		= reconstruction$Psurvival	 # probability of survival to the present (regardless of discovery)
	Pdiscovery		= (if(is.null(discovery_fractions)) reconstruction$Pdiscovery else discovery_fractions)
	Prepresentation	= reconstruction$Prepresentation # probability of survival to the present and discovery
	if((!is.null(Psurvival)) && (length(Psurvival)==0)) Psurvival = NULL
	if((!is.null(Pdiscovery)) && (length(Pdiscovery)==0)) Psurvival = NULL
	if((!is.null(Prepresentation)) && (length(Prepresentation)==0)) Psurvival = NULL
	if(is.null(Psurvival) && (!is.null(Prepresentation)) && (!is.null(Pdiscovery))){
		Psurvival = pmax(0,pmin(1,Prepresentation/Pdiscovery))
	}else if(is.null(Prepresentation) && (!is.null(Psurvival)) && (!is.null(Pdiscovery))){
		Prepresentation = pmax(0,pmin(1,Psurvival * Pdiscovery))
	}else if(is.null(Pdiscovery) && (!is.null(Psurvival)) && (!is.null(Prepresentation))){
		Pdiscovery = pmax(0,pmin(1,Prepresentation/Psurvival))
	}
	
	return(list(success							= TRUE,
				Ntimes							= Ntimes,
				total_diversities 				= (if(coalescent) reconstruction$total_diversities else diversities),
				coalescent_diversities 			= (if(coalescent) diversities else reconstruction$coalescent_diversities),
				birth_rates 					= reconstruction$birth_rates, # these are absolute (total) lineage creation rates, not per-capita rates. To obtain per-capita rates, divide these by total_diversities[]
				death_rates 					= reconstruction$death_rates, # these are absolute (total) lineage loss rates, not per-capita rates. To obtain per-capita rates, divide these by total_diversities[]
				Psurvival		 				= Psurvival,
				Pdiscovery						= Pdiscovery,
				Prepresentation	 				= Prepresentation,
				total_births 					= reconstruction$total_births,
				total_deaths 					= reconstruction$total_deaths,
				last_birth_rate_pc				= (if(coalescent) reconstruction$last_birth_rate_pc else reconstruction$birth_rates[Ntimes]/reconstruction$true_diversities[Ntimes]),
				last_death_rate_pc				= (if(coalescent) reconstruction$last_death_rate_pc else reconstruction$death_rates[Ntimes]/reconstruction$true_diversities[Ntimes]),
				pulled_diversification_rates	= (if(coalescent) reconstruction$pulled_diversification_rates else NULL),	# does not depend on the provided birth_rates_pc
				pulled_extinction_rates			= (if(coalescent) reconstruction$pulled_extinction_rates else NULL),
				pulled_total_diversities		= (if(coalescent) reconstruction$pulled_total_diversities else NULL),
				diversification_rates			= (if(coalescent) NULL else reconstruction$diversification_rates)))
}


sliding_window_average1D = function(signal,Naverage){
	if(Naverage<=1){
		return(signal)
	}else{
		Ntimes = length(signal);
		Nleft = floor(Naverage/2)
		Nright = max(0,Naverage - Nleft - 1)
		average = sapply(1:Ntimes, FUN = function(n) mean(signal[max(1,(n-Nleft)):min(Ntimes,(n+Nright))]))
		return(average)
	}
}