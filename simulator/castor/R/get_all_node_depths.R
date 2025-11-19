# calculate mean phylogenetic depth of each node (=mean distance to each descending tip)
get_all_node_depths = function(tree, as_edge_count=FALSE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	depths = get_mean_depth_per_node_CPP(	Ntips			= Ntips,
											Nnodes			= Nnodes,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
											edge_length		= (if(as_edge_count || is.null(tree$edge.length)) numeric() else tree$edge.length));
	return(depths);
}