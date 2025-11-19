# simulate a Discrete-State Speciation and Extinction (dSSE) model, whereby tree generation is coupled with the evolution of a discrete trait
# birth & death rates (=speciation & extintion rates) are evolving according to a discrete-state continuous-time Markov model.
# Transitions between states can occur either along an edge (anagenetically) and/or during speciation events (cladogenetically).
# Anagenetic transition rates are specified through a transition_matrix_A, while cladogenic transition proabilities are specified through a transition_matrix_C
# The simulation is halted as soon as Ntips>=max_tips (if max_tips>0) and/or time>=max_time (if max_time>0) and/or time>=max_time_eq+equilibrium_time (if max_time_eq>=0)
simulate_dsse = function(	Nstates,							# number of discrete possible states for the trait
							NPstates				= NULL,		# optional number of proxy states, for hiding the original states (i.e. according to a Hidden State Speciation Extinction model)
							proxy_map				= NULL,		# optional 1D integer vector of size Nstates, mapping states to proxy-states, in the case of a HiSSE model. Hence, proxy_map[s] is an integer in 1:NPstates, specifying which proxy-state the state s belongs to. Only relevant if NPstates!=NULL and NPstates!=Nstates
							parameters				= list(), 	# named list of dSSE model parameters. For names and default values see the main function body below.
							start_state				= NULL,		# integer between 1 and Nstates, specifying the state of the first lineage. If NULL, the root state is chosen randomly.
							max_tips				= NULL, 	# integer, specifying the max number of extant+Psampled (if coalescent=TRUE) or extant+extinct+Psampled (if coalescent=FALSE) tips in the simulated tree (prior to any present-day sampling)
							max_extant_tips			= NULL,		# integer, specifying the max number of extant (non-sampled) tips in the simulated tree (prior to any present-day subsampling)
							max_Psampled_tips		= NULL,		# integer, specifying the max number of Poissonian-sampled tips in the simulated tree
							max_time				= NULL,
							max_time_eq				= NULL,
							max_events				= NULL,		# integer, specifying the max number of speciation/extinction/transition events prior to halting the simulation. Set to NULL to not impose any limit on the number of events.
							sampling_fractions		= NULL,		# numeric vector of size NPstates, listing present-day sampling fractions depending on state. sampling_fractions[p] = probability of including an extant species in the tree, if its proxy state is p
							reveal_fractions		= NULL,		# numeric vector of size NPstates, listing reveal fractions depending on state. reveal_fractions[p] = probability of knowing a tip's proxy state, if its proxy state is p
							sampling_rates			= NULL,		# numeric vector of size NPstates, listing Poissonian sampling rates depending on state. sampling_rates[p] = rate of sampling a lineage if its proxy state is p. If NULL, Poissonian sampling is omitted. Can also be a single number.
							coalescent 				= TRUE,
							as_generations			= FALSE,	# if FALSE, then edge lengths correspond to time. If TRUE, then edge lengths correspond to generations (hence if coalescent==false, all edges will have unit length).
							no_full_extinction		= TRUE,		# if true, then extinction of the entire tree is prevented. This is done by temporarily disabling extinctions when the number of extant tips is 1.
							tip_basename			= "",		# basename for tips (e.g. "tip."). 
							node_basename			= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
							include_event_times		= FALSE,
							include_rates			= FALSE,
							include_labels			= TRUE){	# whether to include tip-labels and node-labels as names in the returned state vectors (e.g. tip_states and node_states). Setting this to FALSE may slightly increase computational time/memory efficiency.
	# basic input checking
	if(is.null(max_tips) && is.null(max_extant_tips) && is.null(max_Psampled_tips) && is.null(max_time) && is.null(max_time_eq)) stop("ERROR: At least one of max_tips and/or max_time and/or max_time_eq must be non-NULL")
	is_hisse_model = !(is.null(NPstates) || (NPstates==0) || (NPstates==Nstates))
	if(is_hisse_model && is.null(proxy_map)) stop("ERROR: Missing proxy_map, needed for HiSSE model")
	if(is_hisse_model && (length(proxy_map)!=Nstates)) stop("ERROR: proxy_map has length %d, but should have length %d (Nstates)",length(proxy_map),Nstates)
	if(is_hisse_model && (length(unique(proxy_map))!=NPstates)) stop("ERROR: Not all %d proxy states are represented in proxy_map",NPstates)
	if((!is_hisse_model) && (!is.null(proxy_map)) & ((length(proxy_map)!=Nstates) || (any(proxy_map!=(1:Nstates))))) stop("ERROR: Non-trivial proxy_map contradicts non-HiSSE model")
	if(!is_hisse_model) NPstates = Nstates
	
	# set default model parameters
	parameter_names = c("transition_matrix_A", "transition_matrix_C", "birth_rates", "death_rates")
	if(is.null(parameters$transition_matrix_A)) parameters$transition_matrix_A = matrix(0,nrow=Nstates, ncol=Nstates);
	if(is.null(parameters$transition_matrix_C)) parameters$transition_matrix_C = diag(Nstates);
	# pc birth rates corresponding to each state
	if(is.null(parameters$birth_rates)){
		parameters$birth_rates = rep(1, times=Nstates);
	}else if(length(parameters$birth_rates)==1){
		parameters$birth_rates = rep(parameters$birth_rates, times=Nstates);
	}else if(length(parameters$birth_rates)!=Nstates){
		stop(sprintf("ERROR: Invalid number of birth_rates; expected %d, but got %d",Nstates,length(parameters$birth_rates)))
	}
	# pc death rates corresponding to each state
	if(is.null(parameters$death_rates)){
		parameters$death_rates = rep(1, times=Nstates);
	}else if(length(parameters$death_rates)==1){
		parameters$death_rates = rep(parameters$death_rates, times=Nstates);
	}else if(length(parameters$death_rates)!=Nstates){
		stop(sprintf("ERROR: Invalid number of death_rates; expected %d, but got %d",Nstates,length(parameters$death_rates)))
	}
	# start state
	if(is.null(start_state)){
		start_state = sample.int(n=Nstates,size=1)
	}else if(!(start_state %in% (1:Nstates))){
		stop(sprintf("ERROR: Invalid start_state (%d): Must be an integer between 1 and %d",start_state,Nstates))
	}
	# prepare present-day sampling fractions per proxy state
	if(is.null(sampling_fractions) || (length(sampling_fractions)==0)){
		sampling_fractions = rep(1,NPstates);
	}else if(length(sampling_fractions)==1){
		sampling_fractions = rep(sampling_fractions,NPstates);
	}else if(length(sampling_fractions)!=NPstates){
		stop(sprintf("ERROR: Invalid number of sampling fractions (%d), expected either 0, 1 or %d (NPstates)",length(sampling_fractions),NPstates))
	}
	# prepare reveal fractions  = probability of knowing a tip's proxy state, depending on its proxy state
	if(is.null(reveal_fractions) || (length(reveal_fractions)==0)){
		reveal_fractions = rep(1,NPstates);
	}else if(length(reveal_fractions)==1){
		reveal_fractions = rep(reveal_fractions,NPstates);
	}else if(length(reveal_fractions)!=NPstates){
		stop(sprintf("ERROR: Invalid number of reveal fractions (%d), expected either 0, 1 or %d (NPstates)",length(reveal_fractions),NPstates))
	}
	# prepare Poissonian sampling rates corresponding to each proxy state
	if(is.null(sampling_rates)){
		sampling_rates = rep(0, times=NPstates);
	}else if(length(sampling_rates)==1){
		sampling_rates = rep(sampling_rates, times=NPstates);
	}else if(length(sampling_rates)!=NPstates){
		stop(sprintf("ERROR: Invalid number of sampling_rates; expected %d, but got %d",NPstates,length(sampling_rates)))
	}
	

	# check biological/physical validity of model parameters	
	if(any((sampling_fractions<0) | (sampling_fractions>1))) stop("ERROR: sampling_fractions must be between 0 and 1.")
	if(any(abs(rowSums(parameters$transition_matrix_A))>1e-6*max(abs(parameters$transition_matrix_A)))) stop("ERROR: Anagenetic transition rate matrix does not seem to be valid; some row sums are not zero.")
	if(any(parameters$transition_matrix_C<0)) stop("ERROR: Cladogenic transition probability matrix does not seem to be valid; some entries are negative.")
	if(any(abs(rowSums(parameters$transition_matrix_C)-1)>1e-6)) stop("ERROR: Cladogenic transition probability matrix does not seem to be valid; some row sums differ from 1.")

	# check if some passed parameters are not recognized
	invalids = setdiff(names(parameters),parameter_names)
	if(length(invalids)>0) stop(sprintf("ERROR: Unknown parameter '%s'",invalids[1]))
	
	include_extinct = (!coalescent)
	results = generate_random_tree_Mk_rates_CPP(max_tips					= (if(is.null(max_tips)) -1 else max_tips),
												max_extant_tips				= (if(is.null(max_extant_tips)) -1 else max_extant_tips),
												max_sampled_tips			= (if(is.null(max_Psampled_tips)) -1 else max_Psampled_tips),
												max_time					= (if(is.null(max_time)) -1 else max_time),
												max_time_since_equilibrium	= (if(is.null(max_time_eq)) -1 else max_time_eq),
												max_events					= (if(is.null(max_events)) -1 else max_events),
												Nstates						= Nstates,
												start_state					= max(1,min(Nstates, start_state)) - 1,
												state_birth_rates			= parameters$birth_rates, 
												state_death_rates			= parameters$death_rates,
												state_sampling_rates		= (if(is_hisse_model) sampling_rates[proxy_map] else sampling_rates),
												transition_matrix_A			= as.vector(t(parameters$transition_matrix_A)), # flatten in row-major format
												transition_matrix_C			= as.vector(t(parameters$transition_matrix_C)), # flatten in row-major format
												as_generations				= as_generations,
												no_full_extinction			= no_full_extinction,
												include_extant				= any(sampling_fractions>0),
												include_extinct				= include_extinct,
												include_event_times			= include_event_times,
												include_rates				= include_rates);

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
	class(tree) = "phylo";
	attr(tree,"order") = NULL
	root_time	 	= results$root_time
	clade_states 	= results$clade_states
	extinct_tips 	= results$extinct_tips+1
	extant_tips  	= results$extant_tips+1
	Psampled_tips 	= results$sampled_tips+1
	birth_rates	 	= results$birth_rates_pc
	death_rates	 	= results$death_rates_pc
	
	
	# sample extant tips if needed
	NnonsampledExtant = 0
	if(!(all(sampling_fractions==1) || all(sampling_fractions==0))){
		old2new_clade = rep(0,Ntips+Nnodes) # will be populated later
		if(length(unique(sampling_fractions))==1){
			remove_extant_tips = sample(extant_tips, size=(1-sampling_fractions[1])*length(extant_tips), replace=FALSE);
		}else{
			extant_tip_pstates = clade_states[extant_tips]+1L
			if(is_hisse_model) extant_tip_pstates = proxy_map[extant_tip_pstates]
			remove_extant_tip = logical(length(extant_tips))
			for(pstate in 1:NPstates){
				extant_tips_with_pstate = which(extant_tip_pstates==pstate)
				remove_extant_tip[extant_tips_with_pstate] = as.logical(rbinom(n=length(extant_tips_with_pstate), size=1, prob=1-sampling_fractions[pstate])) # probability of keeping a given extant tip in the tree is sampling_fractions[pstate]
			}
			remove_extant_tips = which(remove_extant_tip)
		}
		if(length(remove_extant_tips)>=Ntips-1) return(list(success=FALSE, error=sprintf("Sampling fractions are too low for the generated tree (%d tips)",Ntips)))
		rarefaction = castor::get_subtree_with_tips(tree, omit_tips=extant_tips[remove_extant_tips], collapse_monofurcations=TRUE, force_keep_root=FALSE)
		tree 		= rarefaction$subtree
		NnonsampledExtant 	= length(remove_extant_tips)
		root_time 	= root_time + rarefaction$root_shift # update root time, in case root has changed
		Ntips 		= length(tree$tip.label)
		Nnodes 		= tree$Nnode
		old2new_clade[rarefaction$new2old_clade] = c(1:(Ntips+Nnodes))
		if(!is.null(clade_states)) clade_states = clade_states[rarefaction$new2old_clade]
		if(!is.null(birth_rates)) birth_rates = birth_rates[rarefaction$new2old_clade]
		if(!is.null(death_rates)) death_rates = death_rates[rarefaction$new2old_clade]
		if(!is.null(extinct_tips)) extinct_tips = old2new_clade[extinct_tips]
		if(!is.null(extant_tips)) extant_tips = old2new_clade[extant_tips]
		if(!is.null(Psampled_tips)) Psampled_tips = old2new_clade[Psampled_tips]
	}

	tip_states  = (if(is.null(clade_states)) NULL else clade_states[1:Ntips]+1L)
	node_states = (if(is.null(clade_states)) NULL else clade_states[(Ntips+1):(Ntips+Nnodes)]+1L)
	tip_proxy_states = NULL; node_proxy_states = NULL;
	if(is_hisse_model){
		if(!is.null(tip_states)) tip_proxy_states = proxy_map[tip_states]
		if(!is.null(node_states)) node_proxy_states = proxy_map[node_states]
	}

	# extract tip/node states & proxy states (if is_hisse_model)
	tip_states  = (if(is.null(clade_states)) NULL else clade_states[1:Ntips]+1L)
	node_states = (if(is.null(clade_states)) NULL else clade_states[(Ntips+1):(Ntips+Nnodes)]+1L)
	tip_proxy_states = NULL; node_proxy_states = NULL;
	if(is_hisse_model){
		if(!is.null(tip_states)) tip_proxy_states = proxy_map[tip_states]
		if(!is.null(node_states)) node_proxy_states = proxy_map[node_states]
	}
	
	# make some tip states (or proxy states, if is_hisse_model) unknown (i.e. assign state=NA)
	if((!is.null(tip_states)) && any(reveal_fractions<1)){
		tip_known = logical(Ntips)
		for(state in 1:NPstates){
			if(is_hisse_model){
				tips_with_state = which(tip_proxy_states==state)
			}else{
				tips_with_state = which(tip_states==state)
			}
			tip_known[tips_with_state] = as.logical(rbinom(n=length(tips_with_state), size=1, prob=reveal_fractions[state]))
		}
		if(is_hisse_model){
			tip_proxy_states[!tip_known] = NA
		}else{
			tip_states[!tip_known] = NA
		}		
	}
	
	# add labels to tip & node states
	if(include_labels){
		if(!is.null(tip_states)) names(tip_states) = tree$tip.label
		if((!is.null(node_states)) && (!is.null(tree$node.label))) names(node_states) = tree$node.label
		if(!is.null(tip_proxy_states)) names(tip_proxy_states) = tree$tip.label
		if((!is.null(node_proxy_states)) && (!is.null(tree$node.label))) names(node_proxy_states) = tree$node.label
	}
	
	return(list(success				= TRUE,
				tree				= tree,
				root_time			= root_time,
				final_time			= results$final_time,
				equilibrium_time	= results$equilibrium_time,
				Nbirths		 		= results$Nbirths,  	# number of birth events in each state
				Ndeaths				= results$Ndeaths,  	# number of death events in each state
				NPsamplings			= results$Nsamplings,
				Ntransitions_A		= matrix(results$Ntransitions_A,ncol=Nstates,byrow=TRUE), # number of anagenetic transition events between each pair of states
				Ntransitions_C		= matrix(results$Ntransitions_C,ncol=Nstates,byrow=TRUE), # number of cladogenic transition events between each pair of states
				NnonsampledExtant	= NnonsampledExtant, # number of extant tips omited at the end
				tip_states			= tip_states,
				node_states			= node_states,
				tip_proxy_states	= tip_proxy_states,
				node_proxy_states	= node_proxy_states,
				start_state			= start_state,
				extant_tips			= extant_tips,
				extinct_tips		= extinct_tips,
				Psampled_tips		= Psampled_tips,
				birth_times			= (if(include_event_times) results$birth_times else NULL),
				death_times			= (if(include_event_times) results$death_times else NULL),
				sampling_times		= (if(include_event_times) results$sampling_times else NULL),
				clade_birth_rates	= (if(include_rates) birth_rates else NULL),
				clade_death_rates	= (if(include_rates) death_rates else NULL)));
	
}