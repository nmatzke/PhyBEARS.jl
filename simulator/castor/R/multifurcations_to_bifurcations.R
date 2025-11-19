# replace multifurcations with multiple bifurcations in a phylogenetic tree
# Note that monofurcations are kept.
# All tips and nodes in the input tree retain their original indices, however the returned tree may include additional nodes and edges. Edge indices may change.
multifurcations_to_bifurcations = function(	tree, 
											dummy_edge_length	= 0, 
											new_node_basename	= "node.", 
											new_node_start_index= NULL){ 	# if NULL, it will be set to Nnodes+1. Only relevant if the input tree already includes node labels.
	Ntips 	= length(tree$tip.label);
	Nnodes	= tree$Nnode;
	Nedges	= nrow(tree$edge);
	
	results = multifurcations_to_bifurcations_CPP(	Ntips,
													Nnodes,
													Nedges,
													tree_edge 			= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													edge_length 		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													dummy_edge_length	= dummy_edge_length);
	
	# reformat results into a valid "phylo" object
	Nnew_nodes = results$Nnew_nodes;
	new_tree = list(Nnode 		= Nnew_nodes,
					tip.label 	= tree$tip.label,
					edge 		= matrix(results$new_tree_edge,ncol=2,byrow=TRUE) + 1,
					edge.length = results$new_edge_length,
					root 		= tree$root)
	if(!is.null(tree$node.label)){
		if(is.null(new_node_start_index)) new_node_start_index = Nnodes+1;
		new_tree$node.label = c(tree$node.label, if(Nnew_nodes<=Nnodes) c() else paste(new_node_basename, new_node_start_index+(0:(Nnew_nodes-Nnodes-1)), sep=""))
	}
	class(new_tree) = "phylo";
	attr(new_tree,"order") = NULL

	return(list(tree 			= new_tree, 
				old2new_edge 	= results$old2new_edge+1,
				Nnodes_added 	= (new_tree$Nnode-Nnodes)));
}