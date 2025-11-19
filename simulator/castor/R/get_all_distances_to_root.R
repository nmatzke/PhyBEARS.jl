# calculate distance from root, for each clade (tips+nodes)
# distance from root = cumulative branch length from root to the clade
get_all_distances_to_root = function(tree, as_edge_count=FALSE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	distances = get_distances_from_root_CPP(Ntips			= Ntips,
											Nnodes			= Nnodes,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
											edge_length		= (if(as_edge_count || is.null(tree$edge.length)) numeric() else tree$edge.length))
	return(distances)
}