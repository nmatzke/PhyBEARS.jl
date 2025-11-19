# for each node, get a list of tips whose MRCA is that node
get_tips_for_mrcas = function(tree, mrca_nodes, check_input=TRUE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	Nmrcas = length(mrca_nodes);
	if(is.character(mrca_nodes)){
		# tips are given as tip names, not indices
		name2node = 1:Nnodes; 
		names(name2node) = tree$node.label;
		mrca_nodes = name2node[mrca_nodes]; 
		if(check_input && any(is.na(mrca_nodes))) stop(sprintf("ERROR: Unknown node name '%s'",tree$node.label[which(is.na(mrca_nodes))[1]]))
	}else{
		if(!is.numeric(mrca_nodes)) stop("ERROR: mrca_nodes must be a character or integer vector");
		if(check_input){
			min_mrca_node = min(mrca_nodes);
			max_mrca_node = max(mrca_nodes);
			if(min_mrca_node<1 || max_mrca_node>Nnodes) stop(sprintf("ERROR: mrca_nodes must be between 1 and Nnodes (%d); instead, found values from %d to %d",Nnodes,min_mrca_node,max_mrca_node))
		}
	}

	results = get_mrca_defining_tips_CPP(	Ntips 			= Ntips,
											Nnodes			= Nnodes,
											Nedges			= nrow(tree$edge),
											tree_edge 		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
											mrcas			= Ntips+mrca_nodes-1,
											verbose			= FALSE,
											verbose_prefix	= "");

	# reformat results
	mrca2first_tip 	= results$mrca2first_tip + 1;
	mrca2last_tip 	= results$mrca2last_tip + 1;
	mrca_tips 		= results$mrca_tips + 1;	
	return(lapply(1:Nmrcas, function(n) mrca_tips[mrca2first_tip[n]:mrca2last_tip[n]]));
}