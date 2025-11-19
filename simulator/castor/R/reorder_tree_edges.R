# reorder edges in a tree for postorder (tips-->root) or pre-order (root-->tips) traversal
# to achieve ape's default "cladewise" order use root_to_tips=TRUE and depth_first_search=TRUE
reorder_tree_edges = function(tree, root_to_tips=TRUE, depth_first_search=TRUE, index_only=FALSE){ 
	Ntips 	= length(tree$tip.label);
	Nnodes	= tree$Nnode;
	Nedges  = nrow(tree$edge);
	
	new2old_edge = 1 + sort_tree_edges_root_to_tips_CPP(Ntips				= Ntips,
														Nnodes				= Nnodes,
														Nedges				= Nedges,
														depth_first_search 	= depth_first_search,
														root_to_tips		= root_to_tips,
														tree_edge			= as.vector(t(tree$edge)) - 1); # flatten in row-major format and adjust clade indices to 0-based
	if(index_only){
		return(new2old_edge);
	}else{
		tree$edge = tree$edge[new2old_edge,]
		if(!is.null(tree$edge.label))  tree$edge.label  = tree$edge.label[new2old_edge];
		if(!is.null(tree$edge.length)) tree$edge.length = tree$edge.length[new2old_edge];
		if(root_to_tips && depth_first_search){
			attr(tree,"order") = "cladewise"
		}else{
			attr(tree,"order") = NULL
		}
		return(tree);
	}
}