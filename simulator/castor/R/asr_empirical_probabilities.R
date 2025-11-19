# For each node in a tree, calculate the empirical frequencies or probabilities of states of a discrete trait, based on the states of descending tips
# This may be a very crude reconstruction of ancestral state probabilities
# Requirements:
#   The tree must be rooted; the root should be the unique node with no parent
#   The tree can include multifurcations as well as monofurcations
asr_empirical_probabilities = function(	tree, 
										tip_states,
										Nstates=NULL,
										probabilities=TRUE,
										check_input=TRUE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;

	# basic error checking
	if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
	if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be integers"))	
	if(is.null(Nstates)) Nstates = max(tip_states);
	if(check_input){
		min_tip_state = min(tip_states)
		max_tip_state = max(tip_states)
		if((min_tip_state<1) || (max_tip_state>Nstates)) stop(sprintf("ERROR: tip_states must be integers between 1 and %d, but found values between %d and %d",Nstates,min_tip_state,max_tip_state))
		if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
	}
			
	results = get_empirical_state_frequencies_per_node_CPP(	Ntips 		= Ntips,
															Nnodes		= Nnodes,
															Nedges		= nrow(tree$edge),
															Nstates		= Nstates,
															tree_edge 	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
															tip_states	= tip_states-1);
	ancestral_likelihoods = matrix(results$frequencies_per_node, ncol=Nstates, byrow=TRUE) # unflatten
	if(probabilities) ancestral_likelihoods = ancestral_likelihoods/rowSums(ancestral_likelihoods); 
	return(list(success=TRUE, ancestral_likelihoods=ancestral_likelihoods));
}