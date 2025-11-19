# extract the subtrees descending from specific nodes, placing those focal nodes as the roots of the extracted subtrees
# this function is optimized for extracting lots of subtrees from large trees
get_subtrees_at_nodes = function(tree, nodes){ 
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	Nedges = nrow(tree$edge);
	NS	   = length(nodes)

	# figure out nodes
	if(is.character(nodes)){
		if(is.null(tree$node.label)) stop("ERROR: Tree must have node labels when specifying nodes as characters")
		node_indices = match(nodes, tree$node.label)
		if(any(is.na(node_indices))){
			invalids = which(is.na(node_indices))
			stop(sprintf("ERROR: %d nodes (e.g., '%s') were not found in tree nodes",length(invalids),nodes[invalids[1]]))
		}
		nodes = node_indices
	}else if(is.numeric(nodes) && all(as.integer(nodes)==nodes)){
		nodes = as.integer(nodes)
		if(any((nodes<1) | (nodes>Nnodes))) stop(sprintf("ERROR: nodes must be between 1 and %d (=Nnodes)",Nnodes));
	}else{
		stop("ERROR: node must be a character or integer")
	}
	
	# extract subtrees
	results = get_subtrees_at_nodes_CPP(Ntips			= Ntips,
										Nnodes			= Nnodes,
										Nedges 			= Nedges,
										tree_edge 		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
										new_root_nodes	= nodes-1);
	
	# reformat results into a valid "phylo" objects
	subtrees 		= vector(mode = "list", length = NS)
	new2old_tips 	= vector(mode = "list", length = NS)
	new2old_nodes	= vector(mode = "list", length = NS)
	new2old_edges	= vector(mode = "list", length = NS)
	if(NS>0){
		for(n in 1:NS){
			Ntips_kept  		= results$Ntips_kept[[n]]
			Nnodes_kept 		= results$Nnodes_kept[[n]]
			new2old_clade 		= results$new2old_clade[[n]] + 1; # switch to 1-based indices
			new2old_edges[[n]]	= results$new2old_edge[[n]] + 1;
			new2old_tips[[n]]	= new2old_clade[1:Ntips_kept]
			new2old_nodes[[n]]	= new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips
			subtrees[[n]] = list(	Nnode 		= Nnodes_kept,
									tip.label 	= (if(is.null(tree$tip.label)) NULL else tree$tip.label[new2old_clade[1:Ntips_kept]]),
									node.label 	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips]),
									edge 		= matrix(results$new_tree_edge[[n]],ncol=2,byrow=TRUE) + 1,
									edge.length = (if(is.null(tree$edge.length)) NULL else tree$edge.length[new2old_edges[[n]]]),
									root 		= results$new_root[[n]]+1,
									root.edge	= (if(results$stem_edges[n]<0) NULL else tree$edge.length[1+results$stem_edges[n]]))
			class(subtrees[[n]]) = "phylo";
			attr(subtrees[[n]],"order") = NULL
		}
	}
	
	if(is.null(tree$edge.length)){
		stem_lengths = NULL
	}else{
		stem_lengths = c(0,tree$edge.length)[2+results$stem_edges]
	}

	return(list(subtrees		= subtrees,
				new2old_tips	= new2old_tips,
				new2old_nodes	= new2old_nodes,
				new2old_edges	= new2old_edges));
}
