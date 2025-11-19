# Merge specific nodes with their parent nodes, thus creating multifurcations
merge_nodes_to_multifurcations = function(tree, 
									nodes_to_merge,
									merge_with_parents=FALSE, # if FALSE, the specified nodes will be merged with their children (if these are nodes). If TRUE, the specified nodes will be merged with their parents.
									keep_ancestral_ages=FALSE){
	Ntips 	= length(tree$tip.label)
	Nnodes	= tree$Nnode
	Nedges	= nrow(tree$edge)
	nodes_to_merge = map_tip_or_node_names_to_indices(tree, A=nodes_to_merge, type='node', list_title='nodes_to_merge', check_input=TRUE)
	
	results = merge_nodes_to_multifurcations_CPP(	Ntips 				= length(tree$tip.label),
													Nnodes				= tree$Nnode,
													Nedges				= nrow(tree$edge),
													tree_edge			= as.vector(t(tree$edge)) - 1,
													edge_length			= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													nodes_to_merge 		= nodes_to_merge-1,
													merge_with_parents	= merge_with_parents,
													keep_ancestral_ages = keep_ancestral_ages)

	# reformat results into a valid "phylo" object
	Nnodes_new	 	= results$Nnodes_new
	Nedges_new	 	= results$Nedges_new
	old2new_node	= results$old2new_node + 1; # switch to 1-based indices
	new2old_node	= results$new2old_node + 1; # switch to 1-based indices
	new_tree = list(Nnode 		= Nnodes_new,
					tip.label 	= tree$tip.label,
					node.label 	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_node]),
					edge 		= matrix(results$new_tree_edge,ncol=2,byrow=TRUE) + 1,
					edge.length = results$new_edge_length,
					root 		= results$new_root+1)
	class(new_tree) = "phylo";
	attr(new_tree,"order") = NULL
	
	return(list(tree			= new_tree,
				new2old_node	= new2old_node,
				old2new_node	= old2new_node,
				Nnodes_removed	= Nnodes-Nnodes_new,
				Nedges_removed	= Nedges-Nedges_new));
}