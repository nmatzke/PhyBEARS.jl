# calculate minimum and maximum distance of any tip from the root
# distance from root = cumulative branch length from root to the clade
get_tree_span = function(tree, as_edge_count=FALSE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	if(Nnodes==0) return(list(min_distance = 0, max_distance = 0))
	
	results = get_min_max_tip_distance_from_root_CPP(	Ntips			= Ntips,
														Nnodes			= Nnodes,
														Nedges			= nrow(tree$edge),
														tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
														edge_length		= (if(as_edge_count || is.null(tree$edge.length)) numeric() else tree$edge.length));

	return(list(min_distance = results$min_distance,
				max_distance = results$max_distance));
}