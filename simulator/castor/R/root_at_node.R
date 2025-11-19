# adjust edge directions such that the specified node becomes the new root
# Note that the number of tips & nodes remains the same
# If update_indices==FALSE, then tip & node indices also remain the same
root_at_node = function(tree, new_root_node, update_indices=TRUE){ 
	Ntips 	= length(tree$tip.label);
	Nnodes	= tree$Nnode;
	
	# figure out root node
	if(is.character(new_root_node)){
		if(is.null(tree$node.label)) stop("ERROR: Tree must have node labels when specifying new_root_node as character")
		new_root_node = match(new_root_node, tree$node.label)
		if(is.na(new_root_node)) stop(sprintf("ERROR: new_root_node '%s' not found in tree nodes",new_root_node))
	}else if(is.numeric(new_root_node) && (as.integer(new_root_node)==new_root_node)){
		new_root_node = as.integer(new_root_node)
		if((new_root_node<1) || (new_root_node>Nnodes)) stop(sprintf("ERROR: new_root_node must be between 1 and %d (=Nnodes), but instead is %d",Nnodes,new_root_node));
	}else{
		stop("ERROR: new_root_node must be a character or integer")
	}
	
	new_edges = root_tree_at_node_CPP(	Ntips			= Ntips,
										Nnodes			= Nnodes,
										Nedges			= nrow(tree$edge),
										tree_edge		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
										new_root_node	= new_root_node-1);
	new_edges = matrix(new_edges+1, ncol=2, byrow=TRUE); # unflatten returned table and shift clade indices to 1-based
	tree$edge = new_edges;
	tree$root = Ntips + new_root_node;
	tree$root.edge = 0
	
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