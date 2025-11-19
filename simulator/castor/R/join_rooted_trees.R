join_rooted_trees = function(	tree1, 					# rooted phylogenetic tree of class "phylo"
								tree2,					# rooted phylogenetic tree of class "phylo"
								target_edge1,			# integer, edge index in tree1 onto which tree2 is to be joined. If <=0, then this refers to the edge leading into the root of tree1.
								target_edge_length1,	# numeric, length of the edge segment in tree1 from the joining-point to the next child node, i.e. how far from the child of target_edge1 should the joining occur.
								root_edge_length2){		# numeric, length of the edge leading into the root of tree2, i.e. the distance from the joining point to the root of tree2
	Ntips1  = length(tree1$tip.label);
	Nnodes1 = tree1$Nnode;
	Nedges1 = nrow(tree1$edge);
	Ntips2  = length(tree2$tip.label);
	Nnodes2 = tree2$Nnode;
	Nedges2 = nrow(tree2$edge);
	
	# join trees
	results = join_rooted_trees_CPP(Ntips1				= Ntips1,
									Nnodes1				= Nnodes1,
									Nedges1				= Nedges1,
									tree_edge1			= as.vector(t(tree1$edge))-1,	# flatten in row-major format and make indices 0-based,
									edge_length1		= (if(is.null(tree1$edge.length)) numeric() else tree1$edge.length),
									Ntips2				= Ntips2,
									Nnodes2				= Nnodes2,
									Nedges2				= Nedges2,
									tree_edge2			= as.vector(t(tree2$edge))-1,
									edge_length2		= (if(is.null(tree2$edge.length)) numeric() else tree2$edge.length),
									target_edge1		= target_edge1-1,
									target_edge_length1	= target_edge_length1,
									root_edge_length2	= root_edge_length2)
		
	# reformat results into a valid "phylo" object
	Ntips  			= results$Ntips
	Nnodes	 		= results$Nnodes
	Nclades			= Ntips+Nnodes
	clade1_to_clade	= results$clade1_to_clade + 1
	clade2_to_clade	= results$clade2_to_clade + 1
	tip_labels		= character(Ntips)
	tip_labels[clade1_to_clade[seq_len(Ntips1)]] = tree1$tip.label
	tip_labels[clade2_to_clade[seq_len(Ntips2)]] = tree2$tip.label
	if(is.null(tree1$node.label) && is.null(tree2$node.label)){
		node_labels = NULL
	}else{
		node_labels = character(Nnodes)
		if(!is.null(tree1$node.label)) node_labels[clade1_to_clade[Ntips1+seq_len(Nnodes1)] - Ntips] = tree1$node.label
		if(!is.null(tree2$node.label)) node_labels[clade2_to_clade[Ntips2+seq_len(Nnodes2)] - Ntips] = tree2$node.label
	}
	tree = list(Nnode 		= Nnodes,
				tip.label 	= tip_labels,
				node.label 	= node_labels,
				edge 		= matrix(results$tree_edge,ncol=2,byrow=TRUE) + 1,
				edge.length = results$edge_length,
				root 		= results$root+1)
	class(tree) 		= "phylo"
	attr(tree,"order") 	= NULL

	return(list(tree				= tree,
				clade1_to_clade		= clade1_to_clade, 
				clade2_to_clade		= clade2_to_clade));
}
