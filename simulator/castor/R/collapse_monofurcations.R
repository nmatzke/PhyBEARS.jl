# remove monofurcations (nodes with only one child) by connecting their incoming and outgoing edges
# if force_keep_root==TRUE, then the root node is always kept even if it only has one child
collapse_monofurcations = function(tree, force_keep_root=TRUE, as_edge_counts=FALSE){
	Ntips 	= length(tree$tip.label)
	Nnodes	= tree$Nnode
	Nedges	= nrow(tree$edge)
	
	results = get_tree_with_collapsed_monofurcations_CPP(	Ntips 			= length(tree$tip.label),
															Nnodes			= tree$Nnode,
															Nedges			= nrow(tree$edge),
															tree_edge		= as.vector(t(tree$edge)) - 1,
															edge_length		= (if(is.null(tree$edge.length) || as_edge_counts) numeric() else tree$edge.length),
															force_keep_root = force_keep_root,
															force_keep_nodes= integer());
	# reformat results into a valid "phylo" object
	Nnodes_new	 	= results$Nnodes_new
	new2old_node	= results$new2old_node + 1; # switch to 1-based indices
	collapsed_tree = list(	Nnode 		= Nnodes_new,
							tip.label 	= tree$tip.label,
							node.label 	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_node]),
							edge 		= matrix(results$new_tree_edge,ncol=2,byrow=TRUE) + 1,
							edge.length = (if(is.null(tree$edge.length)) NULL else results$new_edge_length),
							root 		= results$new_root+1)
	class(collapsed_tree) = "phylo";
	attr(collapsed_tree,"order") = NULL
	
	return(list(tree			= collapsed_tree, 
				new2old_node	= new2old_node, 
				Nnodes_removed	= Nnodes-Nnodes_new));
}