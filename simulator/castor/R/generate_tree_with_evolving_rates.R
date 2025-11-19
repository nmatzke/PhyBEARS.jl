# generate a random phylogenetic tree, by randomly splitting tips at a certain rate and randomly killing tips at a certain rate
# birth & death rates (=speciation & extintion rates) are themselves evolving under one of the following models:
#	'BM': Constrained Brownian Motion model (constrained into an interval via reflection)
#	'OU': Constrained Ornstein-Uhlenbeck model (constrained into an interval by cutting)
#	'Mk': Discrete-state continuous-time Markov model, as defined by the transition rate matrix
# The simulation is halted as soon as Ntips>=max_tips (if max_tips>0) and/or time>=max_time (if max_time>0) and/or time>=max_time_eq+equilibrium_time (if max_time_eq>=0)
generate_tree_with_evolving_rates = function(parameters				= list(), 	# named list of model parameters. For entries and default values (different for each rate_model) see the main function body below
											 rate_model				= 'BM',		# how speciation/extinction rates are evolving along lineages over time. Options are 'BM' (Brownian motion) or 'OU' (Ornstein-Uhlenbeck) or 'Mk'.
											 max_tips				= NULL, 
											 max_time				= NULL,
											 max_time_eq			= NULL,
											 coalescent 			= TRUE,
											 as_generations			= FALSE,	# if FALSE, then edge lengths correspond to time. If TRUE, then edge lengths correspond to generations (hence if coalescent==false, all edges will have unit length).
											 tip_basename			= "",		# basename for tips (e.g. "tip."). 
											 node_basename			= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
											 include_event_times	= FALSE,
											 include_rates			= FALSE){
	if(is.null(max_tips) && is.null(max_time) && is.null(max_time_eq)) stop("ERROR: At least one of max_tips and/or max_time and/or max_time_eq must be non-NULL")
	discrete_state_model = (rate_model=="Mk")
	
	# set default model parameters
	if(rate_model=="BM"){
		parameter_names = c("birth_rate_diffusivity", "min_birth_rate_pc", "max_birth_rate_pc", "death_rate_diffusivity", "min_death_rate_pc", "max_death_rate_pc", "root_birth_rate_pc", "root_death_rate_pc", "rarefaction")
		if(is.null(parameters$birth_rate_diffusivity)) 	parameters$birth_rate_diffusivity = 1;
		if(is.null(parameters$min_birth_rate_pc)) 		parameters$min_birth_rate_pc = 0;
		if(is.null(parameters$max_birth_rate_pc)) 		parameters$max_birth_rate_pc = 1;
		if(is.null(parameters$death_rate_diffusivity)) 	parameters$death_rate_diffusivity = 1;
		if(is.null(parameters$min_death_rate_pc))		parameters$min_death_rate_pc = 0;
		if(is.null(parameters$max_death_rate_pc)) 		parameters$max_death_rate_pc = 1;
		if(is.null(parameters$root_birth_rate_pc))		parameters$root_birth_rate_pc = runif(1, parameters$min_birth_rate_pc, parameters$max_birth_rate_pc);
		if(is.null(parameters$root_death_rate_pc))		parameters$root_death_rate_pc = runif(1, parameters$min_death_rate_pc, parameters$max_death_rate_pc);
		parameters$root_birth_rate_pc = max(parameters$min_birth_rate_pc, min(parameters$max_birth_rate_pc, parameters$root_birth_rate_pc))
		parameters$root_death_rate_pc = max(parameters$min_death_rate_pc, min(parameters$max_death_rate_pc, parameters$root_death_rate_pc))
	}else if(rate_model=="OU"){
		parameter_names = c("birth_rate_pc_mean", "birth_rate_pc_decay_rate", "birth_rate_pc_std", "min_birth_rate_pc", "max_birth_rate_pc", "death_rate_pc_mean", "death_rate_pc_decay_rate", "death_rate_pc_std", "min_death_rate_pc", "max_death_rate_pc", "root_birth_rate_pc", "root_death_rate_pc", "rarefaction")
		if(is.null(parameters$birth_rate_pc_mean)) 		 parameters$birth_rate_pc_mean = 1;
		if(is.null(parameters$birth_rate_pc_decay_rate)) parameters$birth_rate_pc_decay_rate = 1;
		if(is.null(parameters$birth_rate_pc_std)) 		 parameters$birth_rate_pc_std = 1;
		if(is.null(parameters$min_birth_rate_pc)) 		 parameters$min_birth_rate_pc = 0;
		if(is.null(parameters$max_birth_rate_pc)) 		 parameters$max_birth_rate_pc = 1;
		if(is.null(parameters$death_rate_pc_mean)) 		 parameters$death_rate_pc_mean = 1;
		if(is.null(parameters$death_rate_pc_decay_rate)) parameters$death_rate_pc_decay_rate = 1;
		if(is.null(parameters$death_rate_pc_std)) 		 parameters$death_rate_pc_std = 1;
		if(is.null(parameters$min_death_rate_pc))		 parameters$min_death_rate_pc = 0;
		if(is.null(parameters$max_death_rate_pc)) 		 parameters$max_death_rate_pc = 1;
		if(is.null(parameters$root_birth_rate_pc))		 parameters$root_birth_rate_pc = rnorm(1, mean=parameters$birth_rate_pc_mean, sd=parameters$birth_rate_pc_std)
		if(is.null(parameters$root_death_rate_pc))		 parameters$root_death_rate_pc = rnorm(1, mean=parameters$death_rate_pc_mean, sd=parameters$death_rate_pc_std)
		parameters$root_birth_rate_pc = max(parameters$min_birth_rate_pc, min(parameters$max_birth_rate_pc, parameters$root_birth_rate_pc))
		parameters$root_death_rate_pc = max(parameters$min_death_rate_pc, min(parameters$max_death_rate_pc, parameters$root_death_rate_pc))
	}else if(rate_model=="Mk"){
		parameter_names = c("Nstates", "transition_matrix", "state_birth_rates", "state_death_rates", "start_state", "rarefaction")
		if(is.null(parameters$Nstates)) 			parameters$Nstates = 1;
		if(is.null(parameters$transition_matrix))	parameters$transition_matrix = matrix(0,nrow=parameters$Nstates, ncol=parameters$Nstates);
		# pc birth rates corresponding to each state
		if(is.null(parameters$state_birth_rates)){
			parameters$state_birth_rates = rep(1, times=parameters$Nstates);
		}else if(length(parameters$state_birth_rates)==1){
			parameters$state_birth_rates = rep(parameters$state_birth_rates, times=parameters$Nstates);
		}
		# pc death rates corresponding to each state
		if(is.null(parameters$state_death_rates)){
			parameters$state_death_rates = rep(1, times=parameters$Nstates);
		}else if(length(parameters$state_death_rates)==1){
			parameters$state_death_rates = rep(parameters$state_death_rates, times=parameters$Nstates);
		}
		if(is.null(parameters$start_state)) parameters$start_state = sample.int(n=parameters$Nstates,size=1)
	}else{
		stop(sprintf("ERROR: Invalid rate_model '%s'",rate_model))
	}
	if(is.null(parameters$rarefaction)) parameters$rarefaction = 1;
	
	# check if some passed parameters are not recognized
	invalids = setdiff(names(parameters),parameter_names)
	if(length(invalids)>0) stop(sprintf("ERROR: Unknown parameter '%s' provided for model '%s'",invalids[1],rate_model))
	
	if(parameters$rarefaction<=0 || parameters$rarefaction>1) stop("ERROR: rarefaction parameter must be between 0 (non-inclusive) and 1 (inclusive).")
	
	if(rate_model=="BM"){
		results = generate_random_tree_BM_rates_CPP(max_tips					= (if(is.null(max_tips)) -1 else max_tips),
													max_time					= (if(is.null(max_time)) -1 else max_time),
													max_time_since_equilibrium	= (if(is.null(max_time_eq)) -1 else max_time_eq),
													birth_rate_diffusivity 		= parameters$birth_rate_diffusivity, 
													min_birth_rate_pc 			= parameters$min_birth_rate_pc,
													max_birth_rate_pc 			= parameters$max_birth_rate_pc, 
													death_rate_diffusivity 		= parameters$death_rate_diffusivity,
													min_death_rate_pc			= parameters$min_death_rate_pc,
													max_death_rate_pc			= parameters$max_death_rate_pc,
													root_birth_rate_pc			= parameters$root_birth_rate_pc,
													root_death_rate_pc			= parameters$root_death_rate_pc,
													coalescent					= coalescent,
													Nsplits						= 2,
													as_generations				= as_generations,
													include_event_times			= include_event_times,
													include_rates				= include_rates);
	}else if(rate_model=="OU"){
		results = generate_random_tree_OU_rates_CPP(max_tips					= (if(is.null(max_tips)) -1 else max_tips),
													max_time					= (if(is.null(max_time)) -1 else max_time),
													max_time_since_equilibrium	= (if(is.null(max_time_eq)) -1 else max_time_eq),
													birth_rate_pc_mean 			= parameters$birth_rate_pc_mean, 
													birth_rate_pc_decay_rate	= parameters$birth_rate_pc_decay_rate, 
													birth_rate_pc_std			= parameters$birth_rate_pc_std, 
													min_birth_rate_pc 			= parameters$min_birth_rate_pc,
													max_birth_rate_pc 			= parameters$max_birth_rate_pc, 
													death_rate_pc_mean 			= parameters$death_rate_pc_mean, 
													death_rate_pc_decay_rate	= parameters$death_rate_pc_decay_rate, 
													death_rate_pc_std			= parameters$death_rate_pc_std, 
													min_death_rate_pc			= parameters$min_death_rate_pc,
													max_death_rate_pc			= parameters$max_death_rate_pc,
													root_birth_rate_pc			= parameters$root_birth_rate_pc,
													root_death_rate_pc			= parameters$root_death_rate_pc,
													coalescent					= coalescent,
													Nsplits						= 2,
													as_generations				= as_generations,
													include_event_times			= include_event_times,
													include_rates				= include_rates);
	}else{
		results = generate_random_tree_Mk_rates_CPP(max_tips					= (if(is.null(max_tips)) -1 else max_tips),
													max_extant_tips				= -1,
													max_sampled_tips			= -1,
													max_time					= (if(is.null(max_time)) -1 else max_time),
													max_time_since_equilibrium	= (if(is.null(max_time_eq)) -1 else max_time_eq),
													max_events					= -1,
													Nstates						= parameters$Nstates,
													start_state					= max(1,min(parameters$Nstates, parameters$start_state)) - 1,
													state_birth_rates			= parameters$state_birth_rates, 
													state_death_rates			= parameters$state_death_rates,
													state_sampling_rates		= rep(0,times=parameters$Nstates),
													transition_matrix_A			= as.vector(t(parameters$transition_matrix)), # flatten in row-major format
													transition_matrix_C			= numeric(), # no cladogenic transitions included in this model
													as_generations				= as_generations,
													no_full_extinction			= TRUE,
													include_extant				= TRUE,
													include_extinct				= (!coalescent),
													include_event_times			= include_event_times,
													include_rates				= include_rates);
	
	}
	if(!results$success) return(list(success=FALSE, error=results$error)); # something went wrong
	results	= flatten_list_first_level(results) # flatten 1st-level list structure
	Ntips	= results$Ntips
	Nnodes 	= results$Nnodes
	tree = list(Nnode 		= Nnodes,
				tip.label 	= paste(tip_basename, 1:Ntips, sep=""),
				node.label 	= (if(is.null(node_basename)) NULL else paste(node_basename, 1:Nnodes, sep="")),
				edge 		= matrix(results$tree_edge,ncol=2,byrow=TRUE) + 1L,
				edge.length = results$edge_length,
				root 		= results$root+1L)
	class(tree) = "phylo"
	attr(tree,"order") = NULL
	clade_states = results$clade_states
	root_time = results$root_time
	birth_rates_pc = results$birth_rates_pc
	death_rates_pc = results$death_rates_pc
	
	
	# rarefy if needed
	Nrarefied = 0;
	if(parameters$rarefaction<1){
		rarefaction_depth = parameters$rarefaction*Ntips;
		if(rarefaction_depth<2){
			return(list(success=FALSE, error=sprintf("Rarefaction (%g) is too low for the generated tree (%d tips)", parameters$rarefaction,Ntips)))
		}
		keep_tips 	= sample.int(n=Ntips, size=rarefaction_depth, replace=FALSE)
		rarefaction = castor::get_subtree_with_tips(tree, only_tips=keep_tips, omit_tips=FALSE, collapse_monofurcations=TRUE)
		tree 		= rarefaction$subtree
		Nrarefied 	= Ntips - length(tree$tip.label)
		root_time = root_time + rarefaction$root_shift; # update root time, in case root has changed
		if(include_rates){
			birth_rates_pc	= c(birth_rates_pc[rarefaction$new2old_tip],birth_rates_pc[Ntips+rarefaction$new2old_node])
			death_rates_pc	= c(death_rates_pc[rarefaction$new2old_tip],death_rates_pc[Ntips+rarefaction$new2old_node])
		}
		Ntips 	= length(tree$tip.label)
		Nnodes 	= tree$Nnode
		if(!is.null(clade_states)) clade_states = clade_states[rarefaction$new2old_tip]
	}
	
	return(list(success				= TRUE,
				tree				= tree,
				root_time			= root_time,
				final_time			= results$final_time,
				equilibrium_time	= results$equilibrium_time,
				Nbirths		 		= sum(results$Nbirths),
				Ndeaths				= sum(results$Ndeaths),
				Nrarefied			= Nrarefied, # number of tips removed via rarefaction at the end
				states				= (if(is.null(clade_states) || (!discrete_state_model)) NULL else clade_states+1L), # only relevant for discrete-state rate models
				start_state			= (if(discrete_state_model) parameters$start_state else NULL), # only relevant for discrete-state rate models
				root_birth_rate_pc	= (if(discrete_state_model) NULL else parameters$root_birth_rate_pc),
				root_death_rate_pc	= (if(discrete_state_model) NULL else parameters$root_death_rate_pc),
				birth_times			= (if(include_event_times) results$birth_times else NULL),
				death_times			= (if(include_event_times) results$death_times else NULL),
				birth_rates_pc		= (if(include_rates) birth_rates_pc else NULL),
				death_rates_pc		= (if(include_rates) death_rates_pc else NULL)));
	
}