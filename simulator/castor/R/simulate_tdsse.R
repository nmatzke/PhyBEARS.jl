# simulate a time-dependent Discrete-State Speciation and Extinction (tdSSE) model, whereby tree generation is coupled with the evolution of a discrete trait
# birth & death rates (=speciation & extintion rates) are evolving according to a discrete-state continuous-time Markov model with time-dependent transition rates
# Transitions between states can occur either along an edge (anagenetically) and/or during speciation events (cladogenetically).
# Anagenetic transition rates are specified through a time-dependent transition_matrix_A, while cladogenic transition proabilities are specified through a time-dependent transition_matrix_C
# The simulation is halted as soon as Ntips>=max_tips (if max_tips>0) and/or time>=max_time (if max_time>0)
simulate_tdsse = function(	Nstates,							# number of discrete possible states for the trait
							NPstates				= NULL,		# optional number of proxy states, for hiding the original states (i.e. according to a Hidden State Speciation Extinction model)
							proxy_map				= NULL,		# optional 1D integer vector of size Nstates, mapping states to proxy-states, in the case of a HiSSE model. Hence, proxy_map[s] is an integer in 1:NPstates, specifying which proxy-state the state s belongs to. Only relevant if NPstates!=NULL and NPstates!=Nstates
							time_grid				= NULL,		# numeric vector listing grid times in ascending order. The time grid should generally cover the maximum possible simulation time, otherwise it will be polynomially extrapolated (according to splines_degree).
							parameters				= list(), 	# named list of dSSE model parameters. For names and default values see the main function body below.
							splines_degree			= 1,		# polynomial degree of time-dependent model parameters (birth_rates, death_rates, transition_rates) between time-grid points
							start_state				= NULL,		# integer between 1 and Nstates, specifying the state of the first lineage. If NULL, the root state is chosen randomly.
							max_tips				= NULL, 	# integer, specifying the max number of tips in the simulated tree (prior to any subsampling)
							max_time				= NULL,
							max_events				= NULL,		# integer, specifying the max number of speciation/extinction/transition events prior to halting the simulation. Set to NULL to not impose any limit on the number of events.
							sampling_fractions		= NULL,		# numeric vector of size NPstates, listing sampling fractions depending on state. sampling_fractions[p] = probability of including a species in the tree, if its proxy state is p
							reveal_fractions		= NULL,		# numeric vector of size NPstates, listing reveal fractions depending on state. reveal_fractions[p] = probability of knowing a tip's proxy state, if its proxy state is p
							coalescent 				= TRUE,
							as_generations			= FALSE,	# if FALSE, then edge lengths correspond to time. If TRUE, then edge lengths correspond to generations (hence if coalescent==false, all edges will have unit length).
							no_full_extinction		= TRUE,		# if true, then extinction of the entire tree is prevented. This is done by temporarily disabling extinctions when the number of extant tips is 1.
							Nsplits					= 2,	 	# number of children generated at each diversification event. If set to 2, a bifurcating tree is generated. If >2, the tree will be multifurcating.
							tip_basename			= "",		# basename for tips (e.g. "tip."). 
							node_basename			= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
							include_birth_times		= FALSE,
							include_death_times		= FALSE,
							include_labels			= TRUE){	# whether to include tip-labels and node-labels as names in the returned state vectors (e.g. tip_states and node_states). Setting this to FALSE may slightly increase computational time/memory efficiency.
	# basic input checking
	if(is.null(max_tips) && is.null(max_time)) stop("ERROR: At least one of max_tips and/or max_time must be non-NULL")
	is_hisse_model = !(is.null(NPstates) || (NPstates==0) || (NPstates==Nstates))
	if(is_hisse_model && is.null(proxy_map)) stop("ERROR: Missing proxy_map, needed for HiSSE model")
	if(is_hisse_model && (length(proxy_map)!=Nstates)) stop("ERROR: proxy_map has length %d, but should have length %d (Nstates)",length(proxy_map),Nstates)
	if(is_hisse_model && (length(unique(proxy_map))!=NPstates)) stop("ERROR: Not all %d proxy states are represented in proxy_map",NPstates)
	if((!is_hisse_model) && (!is.null(proxy_map)) & ((length(proxy_map)!=Nstates) || (any(proxy_map!=(1:Nstates))))) stop("ERROR: Non-trivial proxy_map contradicts non-HiSSE model")
	if(!is_hisse_model) NPstates = Nstates
	if(is.null(time_grid)) stop("ERROR: Missing time_grid")
	if(length(time_grid)<splines_degree+1) stop(sprintf("time_grid length must be at least splines_degree+1"))
	if((length(time_grid)>1) && (tail(time_grid,1)<=time_grid[1])) stop("Values in time_grid must be in strictly ascending order")
	NG = length(time_grid)
	
	# check parameters and set to defaults if needed
	parameter_names = c("transition_matrix_A", "transition_matrix_C", "birth_rates", "death_rates")
	# pc birth rates corresponding to each state
	if(is.null(parameters$birth_rates)){
		parameters$birth_rates = matrix(1, nrow=Nstates, ncol=NG);
	}else if(length(parameters$birth_rates)==1){
		parameters$birth_rates = matrix(parameters$birth_rates, nrow=Nstates, ncol=NG);
	}else if((!(is.array(parameters$birth_rates))) && (!is.matrix(parameters$birth_rates))){
		stop("ERROR: birth_rates must either be empty, a single number or a 2D matrix")
	}else if(nrow(parameters$birth_rates)!=Nstates){
		stop(sprintf("ERROR: Invalid number of rows in birth_rates; expected %d (Nstates), but got %d",Nstates,nrow(parameters$birth_rates)))
	}else if(ncol(parameters$birth_rates)!=NG){
		stop(sprintf("ERROR: Invalid number of columns in birth_rates; expected %d (Ngrid), but got %d",NG,ncol(parameters$birth_rates)))
	}
	# pc death rates corresponding to each state
	if(is.null(parameters$death_rates)){
		parameters$death_rates = matrix(1, nrow=Nstates, ncol=NG);
	}else if(length(parameters$death_rates)==1){
		parameters$death_rates = matrix(parameters$death_rates, nrow=Nstates, ncol=NG);
	}else if((!(is.array(parameters$death_rates))) && (!is.matrix(parameters$death_rates))){
		stop("ERROR: death_rates must either be empty, a single number or a 2D matrix")
	}else if(nrow(parameters$death_rates)!=Nstates){
		stop(sprintf("ERROR: Invalid number of rows in death_rates; expected %d (Nstates), but got %d",Nstates,nrow(parameters$death_rates)))
	}else if(ncol(parameters$death_rates)!=NG){
		stop(sprintf("ERROR: Invalid number of columns in death_rates; expected %d (Ngrid), but got %d",NG,ncol(parameters$death_rates)))
	}
	# start state
	if(is.null(start_state)){
		start_state = sample.int(n=Nstates,size=1)
	}else if(!(start_state %in% (1:Nstates))){
		stop(sprintf("ERROR: Invalid start_state (%d): Must be an integer between 1 and %d",start_state,Nstates))
	}
	# sampling fractions
	if(is.null(sampling_fractions) || (length(sampling_fractions)==0)){
		sampling_fractions = rep(1,NPstates);
	}else if(length(sampling_fractions)==1){
		sampling_fractions = rep(sampling_fractions,NPstates);
	}else if(length(sampling_fractions)!=NPstates){
		stop(sprintf("ERROR: Invalid number of sampling fractions (%d), expected either 0, 1 or %d (NPstates)",length(sampling_fractions),NPstates))
	}
	# reveal fractions  = probability of knowing a tip's state, depending on its actual state
	if(is.null(reveal_fractions) || (length(reveal_fractions)==0)){
		reveal_fractions = rep(1,NPstates);
	}else if(length(reveal_fractions)==1){
		reveal_fractions = rep(reveal_fractions,NPstates);
	}else if(length(reveal_fractions)!=NPstates){
		stop(sprintf("ERROR: Invalid number of reveal fractions (%d), expected either 0, 1 or %d (NPstates)",length(reveal_fractions),NPstates))
	}
	# anagenetic transition rates
	if(!is.null(parameters$transition_matrix_A)){
		if((!is.array(parameters$transition_matrix_A)) && (!is.matrix(parameters$transition_matrix_A))) stop("ERROR: transition_matrix_A must be a 2D matrix or a 3D array")
		if(nrow(parameters$transition_matrix_A)!=Nstates) stop(sprintf("ERROR: transition_matrix_A must have %d (Nstates) rows, but instead found %d",Nstates,nrow(parameters$transition_matrix_A)))
		if(ncol(parameters$transition_matrix_A)!=Nstates) stop(sprintf("ERROR: transition_matrix_A must have %d (Nstates) columns, but instead found %d",Nstates,ncol(parameters$transition_matrix_A)))
		if(is.matrix(parameters$transition_matrix_A)){
			# time-independent transition rates
			parameters$transition_matrix_A = array(parameters$transition_matrix_A,dim=c(Nstates,Nstates,NG))
		}else if(dim(parameters$transition_matrix_A)[3]!=NG){
			stop(sprintf("ERROR: transition_matrix_A must have %d layers (Ngrid), but instead found %d",NG,dim(parameters$transition_matrix_A)[3]))
		}
		if(any(abs(rowSums(parameters$transition_matrix_A[,,1]))>1e-6*max(abs(parameters$transition_matrix_A[,,1])))) stop("ERROR: Anagenetic transition rate matrix does not seem to be valid; at time %g some row sums are not zero.",time_grid[1])
	}
	# cladogenic transition rates
	if(!is.null(parameters$transition_matrix_C)){
		if((!is.array(parameters$transition_matrix_C)) && (!is.matrix(parameters$transition_matrix_C))) stop("ERROR: transition_matrix_C must be a 2D matrix or a 3D array")
		if(nrow(parameters$transition_matrix_C)!=Nstates) stop(sprintf("ERROR: transition_matrix_C must have %d (Nstates) rows, but instead found %d",Nstates,nrow(parameters$transition_matrix_C)))
		if(ncol(parameters$transition_matrix_C)!=Nstates) stop(sprintf("ERROR: transition_matrix_C must have %d (Nstates) columns, but instead found %d",Nstates,ncol(parameters$transition_matrix_C)))
		if(is.matrix(parameters$transition_matrix_C)){
			# time-independent transition rates
			parameters$transition_matrix_C = array(parameters$transition_matrix_C,dim=c(Nstates,Nstates,NG))
		}else if(dim(parameters$transition_matrix_C)[3]!=NG){
			stop(sprintf("ERROR: transition_matrix_C must have %d layers (Ngrid), but instead found %d",NG,dim(parameters$transition_matrix_C)[3]))
		}
		if(any(parameters$transition_matrix_C<0)) stop("ERROR: Cladogenic transition probability matrix does not seem to be valid; some entries are negative.")
		if(any(abs(rowSums(parameters$transition_matrix_C[,,1])-1)>1e-6)) stop("ERROR: Cladogenic transition probability matrix does not seem to be valid; at time %g some row sums differ from 1.",time_grid[1])
	}
	# sampling fractions
	if(any(sampling_fractions<=0 | sampling_fractions>1)) stop("ERROR: sampling_fractions must be between 0 (non-inclusive) and 1 (inclusive).")

	# check if some passed parameters are not recognized
	invalids = setdiff(names(parameters),parameter_names)
	if(length(invalids)>0) stop(sprintf("ERROR: Unknown parameter '%s'",invalids[1]))

	results = generate_random_tree_tdSSE_CPP(	max_tips					= (if(is.null(max_tips)) -1 else max_tips),
												max_time					= (if(is.null(max_time)) -1 else max_time),
												max_events					= (if(is.null(max_events)) -1 else max_events),
												Nstates						= Nstates,
												start_state					= max(1,min(Nstates, start_state)) - 1,
												time_grid					= time_grid,
												state_birth_rates			= as.vector(t(parameters$birth_rates)), # flatten in row-major format
												state_death_rates			= as.vector(t(parameters$death_rates)),
												transition_matrix_A			= (if(is.null(parameters$transition_matrix_A)) numeric() else as.vector(aperm(parameters$transition_matrix_A,c(3,2,1)))), # flatten in layer-row-major format
												transition_matrix_C			= (if(is.null(parameters$transition_matrix_C)) numeric() else as.vector(aperm(parameters$transition_matrix_C,c(3,2,1)))), # flatten in layer-row-major format
												splines_degree				= splines_degree,
												coalescent					= coalescent,
												Nsplits						= Nsplits,
												as_generations				= as_generations,
												no_full_extinction			= no_full_extinction,
												include_birth_times			= include_birth_times,
												include_death_times			= include_death_times)

	if(!results$success) return(list(success=FALSE, error=results$error)); # something went wrong
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
	
	
	# sub-sample tips (rarefy) if needed
	Nrarefied = 0;
	if(any(sampling_fractions<1)){
		if(length(unique(sampling_fractions))==1){
			keep_tips = sample.int(n=Ntips, size=sampling_fractions[1]*Ntips, replace=FALSE);
		}else{
			tip_pstates = results$clade_states[1:Ntips]+1L
			if(is_hisse_model) tip_pstates = proxy_map[tip_pstates]
			keep_tip = logical(Ntips)
			for(pstate in 1:NPstates){
				tips_with_pstate = which(tip_pstates==pstate)
				keep_tip[tips_with_pstate] = as.logical(rbinom(n=length(tips_with_pstate), size=1, prob=sampling_fractions[pstate])) # probability of keeping a tip in the tree is sampling_fractions[pstate]
			}
			keep_tips = which(keep_tip)
		}
		if(length(keep_tips)<2) return(list(success=FALSE, error=sprintf("Sampling fractions are too low for the generated tree (%d tips)",Ntips)))
		rarefaction = castor::get_subtree_with_tips(tree, only_tips=keep_tips, collapse_monofurcations=TRUE, force_keep_root=FALSE)
		tree 		= rarefaction$subtree
		Nrarefied 	= Ntips - length(tree$tip.label)
		results$root_time = results$root_time + rarefaction$root_shift; # update root time, in case root has changed
		Ntips 	= length(tree$tip.label)
		Nnodes 	= tree$Nnode
		if(!is.null(results$clade_states)) results$clade_states = results$clade_states[rarefaction$new2old_clade]
	}

	tip_states  = (if(is.null(results$clade_states)) NULL else results$clade_states[1:Ntips]+1L)
	node_states = (if(is.null(results$clade_states)) NULL else results$clade_states[(Ntips+1):(Ntips+Nnodes)]+1L)
	tip_proxy_states = NULL; node_proxy_states = NULL;
	if(is_hisse_model){
		if(!is.null(tip_states)) tip_proxy_states = proxy_map[tip_states]
		if(!is.null(node_states)) node_proxy_states = proxy_map[node_states]
	}

	# extract tip/node states & proxy states (if is_hisse_model)
	tip_states  = (if(is.null(results$clade_states)) NULL else results$clade_states[1:Ntips]+1L)
	node_states = (if(is.null(results$clade_states)) NULL else results$clade_states[(Ntips+1):(Ntips+Nnodes)]+1L)
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
				root_time			= results$root_time,
				final_time			= results$final_time,
				Nbirths		 		= results$Nbirths,  	# number of birth events in each state
				Ndeaths				= results$Ndeaths,  	# number of death events in each state
				Ntransitions_A		= matrix(results$Ntransitions_A,ncol=Nstates,byrow=TRUE), # number of anagenetic transition events between each pair of states
				Ntransitions_C		= matrix(results$Ntransitions_C,ncol=Nstates,byrow=TRUE), # number of cladogenic transition events between each pair of states
				Nrarefied			= Nrarefied, # number of tips removed via rarefaction at the end
				tip_states			= tip_states,
				node_states			= node_states,
				tip_proxy_states	= tip_proxy_states,
				node_proxy_states	= node_proxy_states,
				start_state			= start_state,
				birth_times			= (if(include_birth_times) results$birth_times else NULL),
				death_times			= (if(include_death_times) results$death_times else NULL)));
	
}