# root (or re-root) a tree based on a specified outgroup tip
# Note that the number of tips & nodes remains the same
# If update_indices==FALSE, then tip & node indices also remain the same
root_via_outgroup = function(tree, outgroup, update_indices=TRUE){ 
	Ntips 	= length(tree$tip.label);
	Nnodes	= tree$Nnode;
	
	# figure out outgroup tip index
	if(is.character(outgroup)){
		if(is.null(tree$tip.label)) stop("ERROR: Tree must have tip labels when specifying outgroup as character")
		outgroup = match(outgroup, tree$tip.label)
		if(is.na(outgroup)) stop(sprintf("ERROR: outgroup '%s' not found in tree tips",outgroup))
	}else if(is.numeric(outgroup) && (as.integer(outgroup)==outgroup)){
		outgroup = as.integer(outgroup)
		if((outgroup<1) || (outgroup>Ntips)) stop(sprintf("ERROR: outgroup must be between 1 and %d (=Ntips), but instead is %d",Ntips,outgroup));
	}else{
		stop("ERROR: outgroup must be a character or integer")
	}
	
	# determine parent of outgroup (note that edge directions may not make sense if the tree is unrooted)
	# since the outgroup is a tip, there is exactly one edge connected to it, namely the edge connecting it to its parent
	new_root_node = match(outgroup, tree$edge[,2]);
	if(is.na(new_root_node)) match(outgroup, tree$edge[,1])
	if(is.na(new_root_node)) stop("ERROR: Could not determine parent of outgroup; maybe the tree is a forest?") # something went wrong
	
	
	new_edges = root_tree_at_node_CPP(	Ntips			= Ntips,
										Nnodes			= Nnodes,
										Nedges			= nrow(tree$edge),
										tree_edge		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
										new_root_node	= new_root_node-1);
	new_edges = matrix(new_edges+1, ncol=2, byrow=TRUE); # unflatten returned table and shift clade indices to 1-based
	tree$edge = new_edges;
	tree$root = Ntips + new_root_node;
	
	# update node indices if required
	correct_root_node = 1; # correct index that the root node should have
	if(update_indices && (new_root_node!=correct_root_node)){
		# swap indices with wrong node
		temp_root_index = 0;
		tree$edge[tree$edge==(Ntips+new_root_node)] 	= temp_root_index;
		tree$edge[tree$edge==(Ntips+correct_root_node)] = Ntips+new_root_node;
		tree$edge[tree$edge==temp_root_index] 			= Ntips+correct_root_node;
		tree$root 										= Ntips+correct_root_node;
		
		if(!is.null(tree$node.label)){
			root_label 							= tree$node.label[new_root_node];
			tree$node.label[new_root_node] 		= tree$node.label[correct_root_node];
			tree$node.label[correct_root_node] 	= root_label;
		}
	}
	attr(tree,"order") = NULL
	
	return(tree);
}