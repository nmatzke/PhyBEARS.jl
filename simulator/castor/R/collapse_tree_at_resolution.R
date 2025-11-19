# Collapse tree nodes (and their descending subtrees) into tips, whenever all descending tips have a distance from a node below a certain phylogenetic resolution threshold (but see option criterion)
# Nodes are traversed root-->tips and collapsed as soon as the criterion falls below the threshold (=resolution)
# If criterion=='max_tip_depth': Any node whose distance to all its descending tips is <=resolution, will be collapsed into a single tip
# If criterion=='sum_tip_paths': Any node whose sum-of-descending-edge-lengths is <=resolution, will be collapsed into a single tip
# If criterion=='max_tip_pair_dist': Any node of which all pairs of descending tips have distance <=resolution, will be collapsed into a single tip
# If shorten==TRUE, then collapsed nodes are turned into tips at the same location, thus potentially shortening the tree
# If shorten==FALSE, then collapsed nodes are turned into tips, while their incoming edge is extended by some distance L, where L is the distance to the farthest descending tip
# This function can be used to get the "coarse structure" of a tree
collapse_tree_at_resolution = function(	tree, 
										resolution				= 0, 
										by_edge_count			= FALSE, 
										shorten					= TRUE, 
										rename_collapsed_nodes	= FALSE, 
										criterion				= 'max_tip_depth'){		# (character) criterion to use for collapsing nodes (i.e. how to interpret resolution parameter).
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	Nedges = nrow(tree$edge);
	
	 #basic error checking
	 if(!(criterion %in% c('max_tip_depth','sum_tip_paths','max_tip_pair_dist'))) stop(sprintf("criterion must be one of 'max_tip_depth', 'sum_tip_paths', 'max_tip_pair_dist' (got '%s' instead)",criterion));
	
	# collapse
	results = collapse_tree_at_resolution_CPP(	Ntips			= Ntips,
												Nnodes			= Nnodes,
												Nedges 			= Nedges,
												tree_edge 		= as.vector(t(tree$edge)) - 1,
												edge_length		= (if(by_edge_count || is.null(tree$edge.length)) numeric() else tree$edge.length),
												resolution		= resolution,
												shorten			= shorten,
												criterion		= criterion);
	
	# reformat results into a valid "phylo" object
	# note that some of the old nodes may have turned into new tips
	Ntips_new  		= results$Ntips_new
	Nnodes_new	 	= results$Nnodes_new
	Nclades_new		= Ntips_new+Nnodes_new
	new2old_clade 	= results$new2old_clade + 1; # switch to 1-based indices
	new2old_edge	= results$new2old_edge + 1;
	clade_labels	= c(tree$tip.label, tree$node.label)
	collapsed_nodes	= results$collapsed_nodes + 1;
	collapsed_tree = list(	Nnode 		= Nnodes_new,
							tip.label 	= clade_labels[new2old_clade[1:Ntips_new]],
							node.label 	= (if(is.null(tree$node.label)) NULL else clade_labels[new2old_clade[(Ntips_new+1):Nclades_new]]),
							edge 		= matrix(results$new_tree_edge,ncol=2,byrow=TRUE) + 1,
							edge.length = (if(is.null(tree$edge.length)) NULL else (if(shorten) tree$edge.length[new2old_edge] else results$new_edge_length)),
							root 		= results$new_root+1)
	if(rename_collapsed_nodes){
		old2new_clade	= results$old2new_clade + 1;
		collapsed_tree$tip.label[old2new_clade[Ntips+results$collapsed_nodes+1]] = tree$tip.label[results$farthest_tips+1];
	}
	class(collapsed_tree) = "phylo";
	attr(collapsed_tree,"order") = NULL

	return(list(tree				= collapsed_tree, 
				root_shift			= results$root_shift, # distance between old & new root (will always be non-negative)
				collapsed_nodes		= collapsed_nodes,
				farthest_tips		= results$farthest_tips+1,
				new2old_clade		= new2old_clade, 
				new2old_edge		= new2old_edge))
}
