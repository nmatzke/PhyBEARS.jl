# Date (make ultrametric) a rooted tree based on relative evolutionary divergences (RED)
date_tree_red = function(	tree, 
							anchor_node	= NULL,		# index of node to use as anchor (1,..,Nnodes). If NULL or negative, the root is taken as anchor.
							anchor_age	= 1){		# age of anchor node.
	Ntips  = length(tree$tip.label);
	
	results = date_tree_via_RED_CPP(Ntips,
									Nnodes		= tree$Nnode,
									Nedges		= nrow(tree$edge),
									tree_edge	= as.vector(t(tree$edge)) - 1,
									edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
									anchor_node	= (if(is.null(anchor_node)) -1 else (anchor_node-1)),
									anchor_age 	= anchor_age);
	if(!results$success) return(list(success = FALSE, error = results$error));
	tree$edge.length = results$edge_times;

	return(list(success	= TRUE,
				tree	= tree,
				REDs	= results$node_REDs))
}
