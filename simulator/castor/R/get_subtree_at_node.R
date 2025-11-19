# extract the subtree descending from a specific node, and place that node as the root of the extracted subtree
get_subtree_at_node = function(tree, node){ 
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	Nedges = nrow(tree$edge);

	# figure out node
	if(is.character(node)){
		if(is.null(tree$node.label)) stop("ERROR: Tree must have node labels when specifying node as character")
		node_i = match(node, tree$node.label)
		if(is.na(node_i)) stop(sprintf("ERROR: node '%s' not found in tree nodes",node))
		node = node_i
	}else if(is.numeric(node) && (as.integer(node)==node)){
		node = as.integer(node)
		if((node<1) || (node>Nnodes)) stop(sprintf("ERROR: node must be between 1 and %d (=Nnodes), but instead is %d",Nnodes,node));
	}else{
		stop("ERROR: node must be a character or integer")
	}
	
	# extract subtree
	results = get_subtree_at_node_CPP(	Ntips			= Ntips,
										Nnodes			= Nnodes,
										Nedges 			= Nedges,
										tree_edge 		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
										new_root_node	= node-1);
	
	# reformat results into a valid "phylo" object
	Ntips_kept  	= results$Ntips_kept
	Nnodes_kept 	= results$Nnodes_kept
	new2old_clade 	= results$new2old_clade + 1; # switch to 1-based indices
	new2old_edge	= results$new2old_edge + 1;
	subtree = list(	Nnode 		= Nnodes_kept,
					tip.label 	= (if(is.null(tree$tip.label)) NULL else tree$tip.label[new2old_clade[1:Ntips_kept]]),
					node.label 	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips]),
					edge 		= matrix(results$new_tree_edge,ncol=2,byrow=TRUE) + 1,
					edge.length = (if(is.null(tree$edge.length)) NULL else tree$edge.length[new2old_edge]),
					root 		= results$new_root+1)
	class(subtree) = "phylo";
	attr(subtree,"order") = NULL

	return(list(subtree			= subtree, 
				new2old_tip		= new2old_clade[1:Ntips_kept], 
				new2old_node	= new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips,
				new2old_edge	= new2old_edge));
}
