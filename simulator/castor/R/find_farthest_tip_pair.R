# Given a phylogenetic tree, find the two most distant tips in the tree
# Note that while the tree must be rooted, the answer does not actually depend on the rooting
find_farthest_tip_pair = function(tree, as_edge_counts=FALSE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode

	results = get_farthest_tip_pair_CPP(Ntips					= Ntips,
										Nnodes					= Nnodes,
										Nedges					= nrow(tree$edge),
										tree_edge				= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
										edge_length				= (if(as_edge_counts || is.null(tree$edge.length)) numeric() else tree$edge.length))

	return(list(tip1		= (results$farthest_tip1 + 1),
				tip2		= (results$farthest_tip2 + 1),
				distance	= results$max_tip_distance))
}
