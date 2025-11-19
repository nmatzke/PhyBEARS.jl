# calculate relative evolutionary divergence (RED) for each node
# The RED of a node measures the relative placement of the node between the root and the node's descending tips
# The root's RED is always 0, and the RED of tips is always 1.
# Parks, D.H., Chuvochina, M. et al (2018). DOI:10.1101/256800
get_reds = function(tree){
	REDs = get_relative_evolutionary_divergences_CPP(	Ntips			= length(tree$tip.label),
														Nnodes			= tree$Nnode,
														Nedges			= nrow(tree$edge),
														tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
														edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length))
	return(REDs);
}