# split a rooted tree at a certain height (distance from the root), yielding multiple sub-trees each rooted at the splitting time
split_tree_at_height = function(tree, 
								height			= 0, 	# distance from root at which to split. If zero, the original tree is returned
								by_edge_count	= FALSE){
	Ntips  			= length(tree$tip.label);
	Nnodes		 	= tree$Nnode;
	Nedges 			= nrow(tree$edge);
	clade_labels	= c(tree$tip.label, tree$node.label)

	# split
	results = split_tree_at_height_CPP(	Ntips,
										Nnodes,
										Nedges,
										tree_edge 		= as.vector(t(tree$edge)) - 1,
										edge_length		= (if(by_edge_count || is.null(tree$edge.length)) numeric() else tree$edge.length),
										root_edge		= (if(is.null(tree$root.edge)) 0 else tree$root.edge),
										split_height	= height)
	if(!results$success){
		return(list(success=FALSE, error=results$error));
	}else if(results$Nsubtrees==0){
		return(list(Nsubtrees		= 0,
					subtrees		= list(),
					clade2subtree	= rep(0,Ntips+Nnodes)));
	}else{
		# reformat results into a list of valid "phylo" objects
		Nsubtrees 		= results$Nsubtrees
		NStips  		= results$NStips 	# 1D array listing Ntips for each subtree
		NSnodes	 		= results$NSnodes 	# 1D array listing Nnodes for each subtree
		NSclades		= NStips+NSnodes
		NSedges			= results$NSedges 	# 1D array listing Nedges for each subtree
		new2old_clade	= results$new2old_clade + 1; 	# 1D array mapping new clade indices (per tree) to old clade indices, stored serially for all subtrees
		new2old_edge	= results$new2old_edge + 1;		# 1D array mapping new edge indices (per tree) to old edge indices, stored serially for all subtrees
		clade2subtree	= results$clade2subtree + 1;	# 1D array mapping old clade indices to their new subtree. A negative value means that the clade was not included in any of the subtrees.
		subtrees 		= vector(mode="list", Nsubtrees)
		clade_offset	= 0;
		edge_offset		= 0;
		for(i in 1:Nsubtrees){
			subtree = list(	Nnode 		= NSnodes[i],
							tip.label 	= (if(NStips[i]==0) vector("character") else clade_labels[new2old_clade[(1+clade_offset):(clade_offset+NStips[i])]]),
							node.label 	= (if(is.null(tree$node.label) || (NSnodes[i]==0)) NULL else clade_labels[new2old_clade[(clade_offset+NStips[i]+1):(clade_offset+NSclades[i])]]),
							edge 		= (if(NSedges[i]==0) matrix(nrow=0,ncol=2) else matrix(results$subtree_edges[(1+2*edge_offset):(2*edge_offset+2*NSedges[i])],ncol=2,byrow=TRUE) + 1),
							edge.length = (if(is.null(tree$edge.length)) NULL else (if(NSedges[i]==0) vector("double") else tree$edge.length[new2old_edge[(1+edge_offset):(edge_offset+NSedges[i])]])),
							root 		= results$new_roots[i]+1,
							root.edge	= results$root_edges[i])
			class(subtree) 			= "phylo";
			attr(subtree,"order") 	= NULL
			subtrees[[i]] 	= list(	tree 			= subtree,
									new2old_clade 	= (if(NSclades[i]==0) vector("integer") else new2old_clade[(1+clade_offset):(clade_offset+NSclades[i])]),
									new2old_edge 	= (if(NSedges[i]==0) vector("integer") else new2old_edge[(1+edge_offset):(edge_offset+NSedges[i])]));
			clade_offset 	= clade_offset + NSclades[i];
			edge_offset	 	= edge_offset + NSedges[i];
		}
		return(list(Nsubtrees		= Nsubtrees,
					subtrees		= subtrees,
					clade2subtree	= clade2subtree));
	}
}
