# Given a rooted phylogenetic tree and a maximum allowed distance from the root (``height''), 
# remove tips and nodes and shorten the remaining terminal edges so that the tree's height does not exceed the specified threshold. 
# This corresponds to drawing the tree in rectangular layout and trimming everything beyond a specific phylogenetic distance from the root. 
# Tips or nodes at the end of trimmed edges are kept, and the affected edges are shortened. 
trim_tree_at_height = function(	tree, 
								height			= Inf, 	# distance from root at which to trim
								by_edge_count	= FALSE){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	Nedges = nrow(tree$edge);
	
	# trim
	results = trim_tree_at_height_CPP(	Ntips,
										Nnodes,
										Nedges,
										tree_edge 				= as.vector(t(tree$edge)) - 1,
										edge_length				= (if(by_edge_count || is.null(tree$edge.length)) numeric() else tree$edge.length),
										max_distance_from_root	= height)
		
	# reformat results into a valid "phylo" object
	# note that some of the old nodes may have turned into new tips
	Ntips_new  			= results$Ntips_new
	Nnodes_new	 		= results$Nnodes_new
	Nclades_new			= Ntips_new+Nnodes_new
	new2old_clade 		= results$new2old_clade + 1 # switch to 1-based indices
	new2old_edge		= results$new2old_edge + 1
	new_edges_trimmed	= results$new_edges_trimmed + 1
	new_tips_ex_nodes 	= results$new_tips_ex_nodes + 1;
	clade_labels		= c(tree$tip.label, tree$node.label)
	trimmed_tree = list(Nnode 		= Nnodes_new,
						tip.label 	= clade_labels[new2old_clade[1:Ntips_new]],
						node.label 	= (if(is.null(tree$node.label)) NULL else clade_labels[new2old_clade[(Ntips_new+1):Nclades_new]]),
						edge 		= matrix(results$new_tree_edge,ncol=2,byrow=TRUE) + 1,
						edge.length = results$new_edge_length,
						root 		= results$new_root+1)
	class(trimmed_tree) = "phylo";
	attr(trimmed_tree,"order") = NULL

	return(list(tree				= trimmed_tree,
				Nedges_trimmed		= results$Nedges_trimmed,
				Nedges_removed		= (Nedges-results$Nedges_new),
				new2old_clade		= new2old_clade, 
				new2old_edge		= new2old_edge,
				new_edges_trimmed	= new_edges_trimmed));
}
