# Reconstruction of continuous ancestral states via subtree averaging
# Note that reconstructed states are local estimates, i.e. they only take into account the tips descending from the reconstructed node
# The function returns the estimated ancestral states (=averages) as well as the corresponding standard deviations
# Requirements:
# 	Tree can be multifurcating, and can also include nodes with a single child
#	Tree must be rooted.
asr_subtree_averaging = function(	tree, 
									tip_states, # numeric vector of size Ntips
									check_input	= TRUE){
	Ntips  = length(tree$tip.label)

	# basic error checking
	if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
	if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be numeric"))
	if(check_input){
		if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
	}

	results = get_mean_state_per_node_CPP(	Ntips		= Ntips,
											Nnodes		= tree$Nnode,
											Nedges		= nrow(tree$edge),
											tree_edge	= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
											edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
											tip_states	= tip_states)
	return(list(success				= TRUE,
				ancestral_states	= results$means,
				ancestral_stds		= results$stds,
				ancestral_counts	= results$counts)); # number of tips considered for (descending from) each node
}