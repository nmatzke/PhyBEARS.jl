# extend the terminal edges (edges leading to tips) so that each tip has the same fixed distance (new_height) from the root
# if a tip already extends beyond the specified new_height, its incoming edge remains unchanged
# this is a quick-and-dirty way to make the tree ultrametric
# if new_height<0 or new_height==NULL, then it is set to the max_distance_to_root of the input tree
extend_tree_to_height = function(tree, new_height=NULL){
	results = extend_tree_to_height_CPP(	Ntips 			= length(tree$tip.label),
											Nnodes			= tree$Nnode,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge)) - 1,
											edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
											new_height 		= (if(is.null(new_height)) -1.0 else new_height));
	tree$edge.length = results$new_edge_length;
	return(list(tree=tree, max_extension=results$max_extension))
}

