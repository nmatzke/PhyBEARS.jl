# Perform random simulations of a fixed-rates continuous-time Markov model of discrete character evolution
# Starting with a specified vector of root_probabilities, and moving from root to tips, each node is assigned a random state according to its parent's state and according to the markov transition matrix.
# Optionally, multiple independent simulations can be performed using the same model (e.g. as part of some Monte Carlo integration)
simulate_mk_model = function(	tree, 
								Q, 
								root_probabilities	= "stationary", 
								include_tips		= TRUE, 
								include_nodes		= TRUE, 
								Nsimulations		= 1,
								drop_dims			= TRUE){
	if(ncol(Q)!=nrow(Q)) stop(sprintf("ERROR: Transition matrix is not quadratic (has %d rows and %d columns)",nrow(Q),ncol(Q)))
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	Nstates = ncol(Q);
	if(is.character(root_probabilities)){
		if(root_probabilities=="stationary"){
			root_probabilities = get_stationary_distribution(Q);
		}else if(root_probabilities=="flat"){
			root_probabilities = rep(1.0/Nstates, times=Nstates);
		}else{
			stop(sprintf("ERROR: Unknown root_probabilities '%s'",root_probabilities))
		}
	}else{
		if(length(root_probabilities)!=Nstates) stop(sprintf("ERROR: root_probabilities has length (%d) different from the number states (%d)",length(root_probabilities),Nstates))
	}
	results = simulate_fixed_rates_Markov_model_CPP(Ntips				= Ntips,
													Nnodes				= Nnodes,
													Nedges				= nrow(tree$edge),
													Nstates				= Nstates,
													tree_edge 			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
													edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													transition_matrix	= as.vector(t(Q)),				# flatten in row-major format
													root_probabilities	= root_probabilities,
													include_tips		= include_tips,
													include_nodes		= include_nodes,
													Nsimulations		= Nsimulations);

	tip_states  = NULL
	node_states = NULL
	if(include_tips) tip_states   = 1L + as.integer(if(drop_dims && Nsimulations==1) results$tip_states else matrix(results$tip_states, ncol=Ntips, byrow=TRUE));
	if(include_nodes) node_states = 1L + as.integer(if(drop_dims && Nsimulations==1) results$node_states else matrix(results$node_states, ncol=Nnodes, byrow=TRUE));
	return(list(tip_states=tip_states, node_states=node_states));
}

