# simulate a Multiple State Speciation and Extinction (MuSSE) model, whereby tree generation is coupled with the evolution of a discrete trait
# birth & death rates (=speciation & extintion rates) are evolving according to a discrete-state continuous-time Markov model
# Transitions between states occur exclusively along edges, i.e. anagenetically, according to rates specified by a transition_matrix.
# The simulation is halted as soon as Ntips>=max_tips (if max_tips>0) and/or time>=max_time (if max_time>0) and/or time>=max_time_eq+equilibrium_time (if max_time_eq>=0)
simulate_musse = function(	Nstates,							# number of discrete possible states for the trait
							NPstates				= NULL,		# optional number of proxy states, for hiding the original states (i.e. according to a Hidden State Speciation Extinction model)
							proxy_map				= NULL,		# optional 1D integer vector of size Nstates, mapping states to proxy-states, in the case of a HiSSE model. Hence, proxy_map[s] is an integer in 1:NPstates, specifying which proxy-state the state s belongs to. Only relevant if NPstates!=NULL and NPstates!=Nstates
							parameters				= list(), 	# named list of MuSSE model parameters. For names and default values see the main function body below.
							start_state				= NULL,		# integer between 1 and Nstates, specifying the state of the first lineage. If NULL, the start state is chosen randomly.
							max_tips				= NULL,
							max_extant_tips			= NULL,
							max_Psampled_tips		= NULL,		# integer, maximum number of Poissonian-sampled tips
							max_time				= NULL,
							max_time_eq				= NULL,
							max_events				= NULL,		# integer, specifying the max number of speciation/extinction/transition events prior to halting the simulation. Set to NULL to not impose any limit on the number of events.
							sampling_fractions		= NULL,		# numeric vector of size NPstates, listing present-day sampling fractions of extant tips depending on state. sampling_fractions[p] = probability of including an extant species in the tree, if its proxy state is p
							reveal_fractions		= NULL,		# numeric vector of size NPstates, listing reveal fractions depending on state. reveal_fractions[p] = probability of knowing a tip's proxy state, if its proxy state is p
							sampling_rates			= NULL,		# optional numeric vector of size NPstates, listing Poissonian sampling rates depending on state. sampling_fractions[p] = rate of sampling a lineage over time if its proxy state is p. Can also be a single number of NULL.
							coalescent 				= TRUE,
							as_generations			= FALSE,	# if FALSE, then edge lengths correspond to time. If TRUE, then edge lengths correspond to generations (hence if coalescent==false, all edges will have unit length).
							no_full_extinction		= TRUE,		# if true, then extinction of the entire tree is prevented. This is done by temporarily disabling extinctions when the number of extant tips is 1.
							tip_basename			= "",		# basename for tips (e.g. "tip."). 
							node_basename			= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
							include_event_times		= FALSE,
							include_rates			= FALSE,
							include_labels			= TRUE){	# whether to include tip-labels and node-labels as names in the returned state vectors (e.g. tip_states and node_states). Setting this to FALSE may slightly increase computational time/memory efficiency.
	# check that parameters are properly named
	if(is.null(parameters$transition_matrix) && (!is.null(parameters$transition_matrix_A))) parameters$transition_matrix = parameters$transition_matrix_A;
	parameter_names = c("transition_matrix", "birth_rates", "death_rates")
	invalids = setdiff(names(parameters),parameter_names)
	if(length(invalids)>0) stop(sprintf("ERROR: Unknown parameter '%s'",invalids[1]))

	# simulate as a special case of dSSE
	results = simulate_dsse(Nstates					= Nstates,
							NPstates				= NPstates,
							proxy_map				= proxy_map,
							parameters				= list(	transition_matrix_A = parameters$transition_matrix, 
															transition_matrix_C = parameters$transition_matrix_C, 
															birth_rates 		= parameters$birth_rates, 
															death_rates 		= parameters$death_rates),
							start_state				= start_state,
							max_tips				= max_tips,
							max_extant_tips			= max_extant_tips,
							max_Psampled_tips		= max_Psampled_tips,
							max_time				= max_time,
							max_time_eq				= max_time_eq,
							max_events				= max_events,
							sampling_fractions		= sampling_fractions,
							reveal_fractions		= reveal_fractions,
							sampling_rates			= sampling_rates,
							coalescent 				= coalescent,
							as_generations			= as_generations,
							no_full_extinction		= no_full_extinction,
							tip_basename			= tip_basename,
							node_basename			= node_basename,
							include_event_times		= include_event_times,
							include_rates			= include_rates,
							include_labels			= include_labels)	

	return(results);
	
}