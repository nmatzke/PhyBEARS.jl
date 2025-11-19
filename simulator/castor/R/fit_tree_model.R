# Fit a model of random phylogenetic tree generation to an existing phylogenetic tree, based on the waiting times between speciation and extinction events
# In the model, tips (species) are born and killed at Poissonian rates, each of which is a power-law function of the number of extant tips
# The tip split at each speciation event, as well as the tip killed, are chosen randomly uniformly among all extant tips
# The model is parameterized via the various power-law parameters for the birth/death rates
# Any parameter set to NULL will be fitted (within reasonable ranges). 
#   To fix a parameter, set its value explicitly.
#   To fit a parameter within specific lower and upper bounds (constraints), set its value to a 2-tuple (lower,upper)
# The tree need not be ultrametric; any tips not extending all the way to the crown (tips with maximum distance from root), will be interpreted as extinct
fit_tree_model = function(	tree, 
							parameters				= list(),	# NULL, or a named list of parameter values to be fixed or constrained within an interval. Parameters not included in this list are assumed "free" (non-fixed).
							first_guess				= list(),	# NULL, or a named list of initial guesses for a subset of parameters
							min_age					= 0,		# min distance from the tree tips, for a node/tip to be considered in the fitting. Set to NULL or <=0 for no constraint. Must be <=max_age.
							max_age	 				= 0, 		# max distance from the tree tips, for a node/tip to be considered in the fitting. Set to NULL or <=0 or Inf for no constraint.
							age_centile				= NULL,		# fraction of youngest nodes/tips to consider for the fitting. This can be used as an alternative to max_age. E.g. if set to 0.6, then the 60% youngest nodes/tips are considered.
							Ntrials					= 1,
							Nthreads				= 1,
							coalescent				= FALSE,
							discovery_fraction		= NULL,		# optional functional mapping age --> discovery fraction (=probability of a lineage at age tau, that has an extant descendant today, being discovered today). For example, discovery_fraction(0) equals the fraction of extant lineages represented in the tree. If this is provided, then parameters$rarefaction is fixed to 1, and the dscovery_fraction is applied after simulation. Only relevant for coalescent trees.
							fit_control				= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
							min_R2					= -Inf,
							min_wR2					= -Inf,
							grid_size				= 100,
							max_model_runtime		= NULL,		# maximum time (in seconds) to allocate for each evaluation of a model. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
							objective				= 'LL'){ 	# either 'LL' (log-likelihood of waiting times) or 'R2' (R2 of diversity over time) or 'wR2' (weighted R2), 'lR2' (R2 on a logarithmic scale), 'MRD' (mean relative deviation)
	Ntips  				= length(tree$tip.label)
	Nnodes 				= tree$Nnode
	Nedges 				= nrow(tree$edge)
	max_model_runtime 	= (if(is.null(max_model_runtime)) 0 else max_model_runtime);

	if((Ntips<2) || (Nnodes<2)){
		# tree is trivial (~empty)
		return(list(success = FALSE, Nspeciations = 0, Nextinctions = 0, error="Tree is too small"))
	}
	if(is.null(min_age) || (min_age<0)) min_age = 0;
	if(is.null(max_age) || (max_age<0) || (max_age==Inf)) max_age = NULL;
	
	# basic error checking
	if(!(objective %in% c("LL","R2","wR2","lR2",'MRD'))) stop(sprintf("ERROR: Unknown fitting objective '%s'",objective))
	if((!is.null(age_centile)) && (!is.null(max_age))) stop("ERROR: Either a max_age or an age_centile must be specified, but not both")
	if(is.null(age_centile) && is.null(max_age)) stop("ERROR: Both max_age and age_centile are unspecified; please specify one of the two")
	if((!is.null(age_centile)) && (age_centile<0 || age_centile>1)) stop(sprintf("ERROR: age_centile must be between 0 and 1; was set to %g instead",age_centile))
	if(is.null(min_R2)) min_R2   = -Inf
	if(is.null(min_wR2)) min_wR2 = -Inf
	if((!is.null(discovery_fraction)) && (!is.null(parameters$rarefaction)) && (!is.na(parameters$rarefaction)) && (parameters$rarefaction!=discovery_fraction(0))) stop(sprintf("Inconsistent provided rarefaction (%g) and discovery_fraction at age 0 (%g). You can just set parameters$rarefaction to NULL, to automatically adjust to the discovery_fraction at age 0",parameters$rarefaction, discovery_fraction(0)))
	if((!is.null(discovery_fraction)) && (!coalescent)) stop("Provided non-null discovery_fraction for non-coalescent tree")
	if(!is.null(discovery_fraction)){
		# do not fit or consider rarefaction during simulations, so formally set to 1
		# discovery_fraction to the entire time series will be applied after simulations
		parameters$rarefaction  = 1
		first_guess$rarefaction = 1
		todays_discovery_fraction = discovery_fraction(0)
	}else{
		todays_discovery_fraction = NULL
	}
	
	# determine age interval if needed
	if(is.null(max_age)){
		distances_to_root = castor::get_all_distances_to_root(tree)
		max_age = stats::quantile(max(distances_to_root)-distances_to_root, age_centile)
	}else{
		age_centile = NA;
	}
	if(is.null(max_age)) max_age = 0;

	# get waiting times between speciation and extinction events, as observed in the tree
	events = get_speciation_extinction_events_CPP(	Ntips,
													Nnodes,
													Nedges,
													tree_edge	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
													edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													min_age		= min_age,
													max_age 	= max_age,
													only_clades	= numeric(),
													omit_clades = numeric());
	speciation_time_range = if(events$Nspeciations==0) NA else (max(events$speciation_times)-min(events$speciation_times))
	extinction_time_range = if(events$Nextinctions==0) NA else (max(events$extinction_times)-min(events$extinction_times))
	if((events$Nspeciations==0) && (events$Nextinctions==0)){
		# no speciation/extinction events found
		return(list(success			= FALSE,
					error			= "Tree contains no speciation or extinction events",
					Nspeciations	= events$Nspeciations, 
					Nextinctions	= events$Nextinctions,
					objective		= objective))
	}
													
	# get diversity-vs-time curve on a time grid
	grid_times = seq(from=(if(max_age>0) max(0,events$max_distance_from_root-max_age) else 0.0), to=(if(min_age>0) max(0,events$max_distance_from_root-min_age) else events$max_distance_from_root*(1-1e-3/grid_size)), length.out=grid_size)
	tree_diversities_on_grid = count_clades_at_times_CPP(	Ntips		= Ntips,
															Nnodes		= Nnodes,
															Nedges		= Nedges,
															tree_edge	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
															edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
															times		= grid_times,
															degree		= 1);
	tree_diversity_slopes_on_grid = first_derivative(grid_times, tree_diversities_on_grid)
	
	# pre-calculate discovery fraction on time grid & events
	if(!is.null(discovery_fraction)){
		discovery_fraction_at_grid_times = discovery_fraction(max_age-grid_times)
		if(objective=='LL') discovery_fraction_at_speciation_times = discovery_fraction(max_age-events$speciation_times)
	}
															
															
	#################################
	# PREPARE PARAMETERS TO BE FITTED
	
	parameter_names = c("birth_rate_intercept", "birth_rate_factor", "birth_rate_exponent", "death_rate_intercept", "death_rate_factor", "death_rate_exponent", "resolution", "rarefaction", "extant_diversity")
	if((!is.null(parameters$birth_rate_factor)) && (parameters$birth_rate_factor==0) && (is.null(parameters$birth_rate_exponent) || (length(parameters$birth_rate_exponent)==2))) parameters$birth_rate_exponent = 0; # this parameter cannot possibly be determined
	if((!is.null(parameters$death_rate_factor)) && (parameters$death_rate_factor==0) && (is.null(parameters$death_rate_exponent) || (length(parameters$death_rate_exponent)==2))) parameters$death_rate_exponent = 0; # this parameter cannot possibly be determined
	min_params = list(); 
	max_params = list();
	param_fixed = list()
	param_constrained = list()
	for(param_name in parameter_names){
		value = parameters[[param_name]]
		if(is.null(value) || (length(value)==2)){
			param_fixed[[param_name]] = FALSE
		}else{
			param_fixed[[param_name]] = TRUE
		}
		if((!is.null(value)) && (length(value)==2)){
			param_constrained[[param_name]] = TRUE
			min_params[[param_name]] = value[1]
			max_params[[param_name]] = value[2]
		}else{
			param_constrained[[param_name]] = FALSE
			min_params[[param_name]] = 0
			max_params[[param_name]] = 1e+50
		}
	}
	max_params$rarefaction = max(0,min(1,max_params$rarefaction)) # special case, parameter must be within [0,1]
	max_params$resolution = max(0,min(max_age,max_params$resolution)) # special case, parameter should be within [0,max_age]
	if((param_fixed$resolution) && (parameters$resolution<=0) && (param_fixed$rarefaction) && (!param_fixed$extant_diversity)){
		# extant diversity can be immediately calculated based on known rarefaction and Ntips
		parameters$extant_diversity = (Ntips/parameters$rarefaction)/(if(is.null(todays_discovery_fraction) || (!coalescent)) 1 else todays_discovery_fraction)
		param_fixed$extant_diversity = TRUE
		param_constrained$extant_diversity = FALSE
	}
	fitted_parameter_names = names(param_fixed)[!unlist(param_fixed,use.names=FALSE)]
	if(is.null(fitted_parameter_names) || (length(fitted_parameter_names)==0)) stop("ERROR: All model parameters are fixed")
	NFP = length(fitted_parameter_names);
	include_speciations	= !((!is.null(parameters$birth_rate_intercept)) && (!is.null(parameters$birth_rate_factor)) && (parameters$birth_rate_intercept==0) && (parameters$birth_rate_factor==0))
	include_extinctions	= !((!is.null(parameters$death_rate_intercept)) && (!is.null(parameters$death_rate_factor)) && (parameters$death_rate_intercept==0) && (parameters$death_rate_factor==0))
	fit_birth_rates		= !(param_fixed$birth_rate_intercept && param_fixed$birth_rate_factor && param_fixed$birth_rate_exponent)
	fit_death_rates		= !(param_fixed$death_rate_intercept && param_fixed$death_rate_factor && param_fixed$death_rate_exponent)
	
	# prevent zero birth rates
	if(param_fixed$birth_rate_intercept && (parameters$birth_rate_intercept==0) && (events$Nspeciations>0)) min_params$birth_rate_factor = 1e-6*events$Nspeciations/(Ntips*speciation_time_range)
		
	# very rough first guess values
	guessed_birth_rate_exponent = (if(param_fixed$birth_rate_exponent) parameters$birth_rate_exponent else 1)
	guessed_rarefaction = (if(param_fixed$rarefaction) parameters$rarefaction else 1)
	guessed_extant_diversity = (Ntips/guessed_rarefaction)/(if(is.null(todays_discovery_fraction) || (!coalescent)) 1 else todays_discovery_fraction)
	start_params = list(	birth_rate_intercept 	= events$Nspeciations/speciation_time_range,
							birth_rate_factor 		= (tree_diversities_on_grid[grid_size]-tree_diversities_on_grid[grid_size-2])/(guessed_rarefaction*(tree_diversities_on_grid[grid_size]**guessed_birth_rate_exponent)*(grid_times[grid_size]-grid_times[grid_size-2])),
							birth_rate_exponent 	= guessed_birth_rate_exponent, 
							death_rate_intercept	= (if(coalescent) events$Nspeciations/speciation_time_range else events$Nextinctions/extinction_time_range),
							death_rate_factor		= (if(coalescent) log(events$Nspeciations)/speciation_time_range else log(events$Nextinctions)/extinction_time_range),
							death_rate_exponent		= 1,
							resolution				= max(min_age,grid_times[length(grid_times)]-grid_times[length(grid_times)-1]),
							rarefaction				= guessed_rarefaction,
							extant_diversity		= guessed_extant_diversity)
	# adopt provided first guesses if available
	if(!is.null(first_guess)){
		for(param_name in parameter_names){
			if(!is.null(first_guess[[param_name]])){
				start_params[[param_name]] = first_guess[[param_name]];
			}
		}
	}
	# constrain first guesses if needed
	for(param_name in parameter_names){
		if(param_fixed[[param_name]]){
			start_params[[param_name]] = parameters[[param_name]]
		}else if(is.na(start_params[[param_name]])){
			start_params[[param_name]] = 0.5*(min_params[[param_name]]+max_params[[param_name]])
		}else{
			start_params[[param_name]] = max(min_params[[param_name]],min(max_params[[param_name]],start_params[[param_name]]))
		}
	}


	################################
	# AUXILIARY FUNCTION DEFINITIONS
	
	# auxiliary function for mapping fitted --> all parameters
	expand_fitted_model_parameters = function(fitted_params){
		expanded_params = parameters;
		if(length(fitted_params)==0) return(expanded_params)
		for(param_name in fitted_parameter_names) expanded_params[[param_name]] = fitted_params[[param_name]];
		return(expanded_params)
	}
	
	# objective function: negated log-likelihood of Poisson process with variable borth/death probability rates
	objective_function = function(fitted_params){
		params = expand_fitted_model_parameters(fitted_params);
		if(any(is.nan(unlist(params)))) return(Inf); # catch weird cases where params become NaN
		if(objective %in% c('R2','wR2','lR2','MRD')){
			objective_value = get_objective_value_for_tree_model(params,objective)
		}else if(objective=='LL'){	
			objective_value = 0;
			if(coalescent){
				# for coalescent trees, deaths (extinctions) are not observed
				# however death rates influence the apparent birth rates (as observed in the tree)
				# so only the apparent waiting times between births are available for fitting
				model_predictions = simulate_deterministic_diversity_growth_CPP(birth_rate_intercept 		= params$birth_rate_intercept,
																				birth_rate_factor 			= params$birth_rate_factor,
																				birth_rate_exponent 		= params$birth_rate_exponent,
																				death_rate_intercept 		= params$death_rate_intercept,
																				death_rate_factor 			= params$death_rate_factor,
																				death_rate_exponent 		= params$death_rate_exponent,
																				resolution					= params$resolution,
																				rarefaction					= params$rarefaction,
																				Nsplits 					= 2,
																				additional_rates_times		= numeric(),
																				additional_birth_rates_pc	= numeric(),
																				additional_death_rates_pc	= numeric(),
																				additional_periodic			= FALSE,
																				times 						= events$speciation_times,
																				start_time					= events$speciation_times[1],
																				final_time					= events$max_distance_from_root,
																				start_diversity				= 1,
																				final_diversity				= params$extant_diversity,
																				reverse						= coalescent,
																				include_coalescent			= coalescent,
																				include_probabilities		= TRUE,
																				include_birth_rates			= TRUE,
																				include_death_rates			= FALSE,
																				include_Nevents				= FALSE,
																				runtime_out_seconds			= max_model_runtime)
				if(!model_predictions$success){
					objective_value = -Inf; # simulation of deterministic model failed
				}else{
					if(!is.null(discovery_fraction)) model_predictions$Prepresentation = model_predictions$Prepresentation * discovery_fraction_at_speciation_times
					apparent_birth_rates = model_predictions$birth_rates * pmax(1e-16,model_predictions$Prepresentation);
					valids				 = which(apparent_birth_rates>0)
					objective_value 	 = (sum(log(apparent_birth_rates[valids])) - sum(apparent_birth_rates[valids]*events$speciation_waiting_times[valids]))/length(valids);
				}
				
			}else{
				# tree is not coalescent, i.e. both births (speciations) & deaths (extinctions) are visible
				if(fit_birth_rates){
					birth_rates   	= rep(params$birth_rate_intercept, events$Nspeciations) + params$birth_rate_factor * (events$speciation_diversities ** params$birth_rate_exponent)
					objective_value = objective_value + sum(log(birth_rates)) - sum(birth_rates*events$speciation_waiting_times);
				}
				if(fit_death_rates){
					death_rates		= rep(params$death_rate_intercept, events$Nextinctions) + params$death_rate_factor * (events$extinction_diversities ** params$death_rate_exponent)			
					objective_value = objective_value + sum(log(death_rates)) - sum(death_rates*events$extinction_waiting_times);
				}
			}
		}
		if(is.na(objective_value) || is.nan(objective_value)) objective_value = -Inf
		return(-objective_value);
	}
		
	get_objective_value_for_tree_model = function(params,objective){
		simulation = simulate_deterministic_diversity_growth_CPP(	birth_rate_intercept 		= params$birth_rate_intercept,
																	birth_rate_factor 			= params$birth_rate_factor,
																	birth_rate_exponent 		= params$birth_rate_exponent,
																	death_rate_intercept 		= params$death_rate_intercept,
																	death_rate_factor 			= params$death_rate_factor,
																	death_rate_exponent 		= params$death_rate_exponent,
																	resolution					= params$resolution,
																	rarefaction					= params$rarefaction,
																	Nsplits 					= 2,
																	additional_rates_times		= numeric(),
																	additional_birth_rates_pc	= numeric(),
																	additional_death_rates_pc	= numeric(),
																	additional_periodic			= FALSE,
																	times 						= grid_times,
																	start_time					= grid_times[1],
																	final_time					= events$max_distance_from_root,
																	start_diversity				= tree_diversities_on_grid[1],
																	final_diversity				= params$extant_diversity,
																	reverse						= coalescent,
																	include_coalescent			= coalescent,
																	include_probabilities		= FALSE,
																	include_birth_rates			= FALSE,
																	include_death_rates			= FALSE,
																	include_Nevents				= FALSE,
																	runtime_out_seconds			= max_model_runtime);
		if(!simulation$success) return(NaN);
		predicted_diversities = (if(coalescent) simulation$coalescent_diversities else simulation$total_diversities);
		if(!is.null(discovery_fraction)) predicted_diversities = predicted_diversities * discovery_fraction_at_grid_times
		if(objective=='wR2'){
			objective_value = 1.0 - mean((predicted_diversities/tree_diversities_on_grid - 1)**2)	
		}else if(objective=='R2'){
			objective_value = 1.0 - mean((predicted_diversities-tree_diversities_on_grid)**2)/variance(tree_diversities_on_grid);
		}else if(objective=='lR2'){
			valids = which((predicted_diversities>0) & (tree_diversities_on_grid>0))
			objective_value = 1.0 - mean((log(predicted_diversities[valids])-log(tree_diversities_on_grid[valids]))**2)/variance(log(tree_diversities_on_grid[valids]));				
		}else if(objective=='MRD'){
			valids = which(tree_diversities_on_grid>0)
			objective_value = -mean(abs(predicted_diversities[valids]/tree_diversities_on_grid[valids] - 1))	
		}
		return(objective_value);
	}

	# fit with various starting points
	fit_single_trial = function(trial){
		lower_bounds = unlist(min_params[fitted_parameter_names])
		upper_bounds = unlist(max_params[fitted_parameter_names])
		# randomly choose start values
		initial_fitted_params = unlist(start_params[fitted_parameter_names])
		if(trial>1){
			constrained  	= which(unlist(param_constrained[fitted_parameter_names]))
			unconstrained  	= which(!unlist(param_constrained[fitted_parameter_names]))
			if(length(constrained)>0) initial_fitted_params[constrained] 		= lower_bounds[constrained] + (upper_bounds[constrained]-lower_bounds[constrained]) * runif(n=length(constrained),min=0,max=1)
			if(length(unconstrained)>0) initial_fitted_params[unconstrained]	= 10**runif(n=length(unconstrained), min=-2, max=2) * initial_fitted_params[unconstrained]
		}
		initial_fitted_params = pmax(lower_bounds,pmin(upper_bounds,initial_fitted_params))
		# run fit
		fit 	= stats::nlminb(initial_fitted_params, objective=objective_function, lower=lower_bounds, upper=upper_bounds, control=fit_control)
		params 	= expand_fitted_model_parameters(fit$par);
		if(objective=='R2'){
			R2  = -fit$objective;
		}else{
			R2 = get_objective_value_for_tree_model(params,'R2')
		}
		if(objective=='wR2'){
			wR2  = -fit$objective;
		}else{
			wR2 = get_objective_value_for_tree_model(params,'wR2')
		}
		if(objective=='lR2'){
			lR2  = -fit$objective;
		}else{
			lR2 = get_objective_value_for_tree_model(params,'lR2')
		}
		if(objective=='MRD'){
			MRD  = fit$objective;
		}else{
			MRD = -get_objective_value_for_tree_model(params,'MRD')
		}
		return(list(objective_value=-fit$objective, params = fit$par, R2=R2, wR2=wR2, lR2=lR2, MRD=MRD));
	}
	
	################################

	# run one or more independent fitting trials
    if((Ntrials>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
		# run trials in parallel using multiple forks
		# Note: Forks (and hence shared memory) are not available on Windows
		fits = parallel::mclapply(	1:Ntrials, 
									FUN = function(trial) fit_single_trial(trial), 
									mc.cores = min(Nthreads, Ntrials), 
									mc.preschedule = FALSE, 
									mc.cleanup = TRUE);
	}else{
		# run in serial mode
		fits = sapply(1:Ntrials,function(x) NULL)
		for(trial in 1:Ntrials){
			fits[[trial]] = fit_single_trial(trial)
		}
	}

	# extract information from best fit (note that some fits may have LL=NaN or NA)
	returned_value_on_error = list(	success					= FALSE, 
									error					= "Fitting failed for all trials", 
									start_parameters		= start_params,
									Nspeciations			= events$Nspeciations, 
									Nextinctions			= events$Nextinctions,
									grid_times				= grid_times,
									tree_diversities		= tree_diversities_on_grid,
									fitted_parameter_names	= fitted_parameter_names,
									objective				= objective)
	objective_values	= sapply(1:Ntrials, function(trial) fits[[trial]]$objective_value)
	R2s			 		= sapply(1:Ntrials, function(trial) fits[[trial]]$R2)
	wR2s			 	= sapply(1:Ntrials, function(trial) fits[[trial]]$wR2)
	valids				= which((!is.na(objective_values)) & (!is.nan(objective_values)) & (!is.null(objective_values)) & (!is.infinite(objective_values)) & (R2s>=min_R2) & (wR2s>=min_wR2))
	if(length(valids)==0) return(returned_value_on_error); # fitting failed for all trials
	best 				= valids[which.max(sapply(valids, function(i) objective_values[i]))]
	objective_value		= fits[[best]]$objective_value;
	fitted_params 		= fits[[best]]$params;
	parameters		 	= expand_fitted_model_parameters(fitted_params);
	if(is.null(objective_value) || any(is.na(fitted_params)) || any(is.nan(fitted_params))) return(returned_value_on_error); # fitting failed
	locally_fitted_parameters = setNames(lapply(1:NFP, FUN=function(fp) sapply(valids, FUN=function(v) fits[[v]]$params[[fp]])), fitted_parameter_names)
	
	# calculate deterministic diversities on grid_times, based on the fitted birth-death model
	simulation = simulate_deterministic_diversity_growth_CPP(	birth_rate_intercept 		= parameters$birth_rate_intercept,
																birth_rate_factor 			= parameters$birth_rate_factor,
																birth_rate_exponent 		= parameters$birth_rate_exponent,
																death_rate_intercept 		= parameters$death_rate_intercept,
																death_rate_factor 			= parameters$death_rate_factor,
																death_rate_exponent 		= parameters$death_rate_exponent,
																resolution					= parameters$resolution,
																rarefaction					= parameters$rarefaction,
																Nsplits 					= 2,
																additional_rates_times		= numeric(),
																additional_birth_rates_pc	= numeric(),
																additional_death_rates_pc	= numeric(),
																additional_periodic			= FALSE,
																times 						= grid_times,
																start_time					= grid_times[1],
																final_time					= events$max_distance_from_root,
																start_diversity				= tree_diversities_on_grid[1],
																final_diversity				= parameters$extant_diversity,
																reverse						= coalescent,
																include_coalescent			= coalescent,
																include_probabilities		= FALSE,
																include_birth_rates			= FALSE,
																include_death_rates			= FALSE,
																include_Nevents				= FALSE,
																runtime_out_seconds			= 0);
	model_diversities_on_grid = (if(coalescent) simulation$coalescent_diversities else simulation$total_diversities)
	if(!is.null(discovery_fraction)) model_diversities_on_grid = model_diversities_on_grid * discovery_fraction_at_grid_times
	
	# return results
	if(!is.null(discovery_fraction)) parameters$rarefaction	= todays_discovery_fraction # indicate effective rarefaction (=discovery_fraction at age 0)
	return(list(success						= TRUE,
				objective_value				= (if(objective=='MRD') -objective_value else objective_value), 
				parameters					= parameters, 
				start_parameters			= start_params,
				R2							= fits[[best]]$R2,
				wR2							= fits[[best]]$wR2,
				lR2							= fits[[best]]$lR2,
				MRD							= fits[[best]]$MRD,
				Nspeciations				= events$Nspeciations, 
				Nextinctions				= events$Nextinctions,
				grid_times					= grid_times,
				tree_diversities			= tree_diversities_on_grid,
				model_diversities			= model_diversities_on_grid,
				fitted_parameter_names 		= fitted_parameter_names,
				locally_fitted_parameters	= locally_fitted_parameters,
				objective					= objective,
				Ntips						= Ntips,
				Nnodes						= Nnodes,
				min_age						= min_age,
				max_age						= max_age,
				age_centile					= age_centile))
}


variance = function(X){
	return(sum((X-mean(X))**2)/length(X));
}

back = function(X){
	return(X[length(X)])
}


sprint_tree_parameters = function(parameters, separator){
	return(sprintf("birth_rate_intercept = %g%sbirth_rate_factor = %g%sbirth_rate_exponent = %g%sdeath_rate_intercept = %g%sdeath_rate_factor = %g%sdeath_rate_exponent = %g%sresolution = %g%srarefaction = %g",parameters$birth_rate_intercept,separator,parameters$birth_rate_factor,separator,parameters$birth_rate_exponent,separator,parameters$death_rate_intercept,separator,parameters$death_rate_factor,separator,parameters$death_rate_exponent,separator,parameters$rarefaction,separator,parameters$resolution))
}


first_derivative = function(times, values){
	slopes = diff(values)/diff(times);
	slopes = c(slopes[1],slopes) # add one more identical point to the end to keep the same size as the original time series
	return(slopes)
}
