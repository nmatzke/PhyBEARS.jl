# Perform random simulation of an Ornstein-Uhlenbeck model of continuous trait evolution along a tree, moving from root to tips
# The root's state is drawn randomly from the OU stationary distribution
# Optionally, multiple independent simulations can be performed using the same model (e.g. as part of some Monte Carlo integration)
simulate_ou_model = function(	tree, 
								stationary_mean, 
								spread, 
								decay_rate, 
								include_tips	= TRUE, 
								include_nodes	= TRUE, 
								Nsimulations	= 1,
								drop_dims		= TRUE){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	results = simulate_Ornstein_Uhlenbeck_on_tree_CPP(Ntips				= Ntips,
													Nnodes				= Nnodes,
													Nedges				= nrow(tree$edge),
													tree_edge 			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
													edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													stationary_mean		= stationary_mean,
													stationary_std		= spread,
													decay_rate			= decay_rate,
													include_tips		= include_tips,
													include_nodes		= include_nodes,
													Nsimulations		= Nsimulations);

	tip_states  = NULL
	node_states = NULL
	if(include_tips) tip_states = (if(drop_dims && Nsimulations==1) results$tip_states else matrix(results$tip_states, ncol=Ntips, byrow=TRUE));
	if(include_nodes) node_states = (if(drop_dims && Nsimulations==1) results$node_states else matrix(results$node_states, ncol=Nnodes, byrow=TRUE));
	return(list(tip_states=tip_states, node_states=node_states));
}