# eliminate zero-length or very short edges by merging nodes into multifurcations
# This is similar to ape::di2multi
# The total number of nodes/tips may decrease
# Note that monofurcations are kept
merge_short_edges = function(	tree, 
								edge_length_epsilon	= 0,
								force_keep_tips		= TRUE,	# if TRUE, all tips are kept even if their incoming edges are short. Note that tip indices may still change. If FALSE, then some tips may be removed (and thus some nodes may become tips).
								new_tip_prefix		= "ex.node.tip."){	# (character) basename for tips in the new tree that used to be nodes. Can be NULL, in which case the original node labels are used
	Ntips 	= length(tree$tip.label);
	Nnodes	= tree$Nnode;
	Nedges	= nrow(tree$edge);
	
	# basic error checking
	if((!force_keep_tips) && is.null(new_tip_prefix) && is.null(tree$node.label)) stop("Missing new_tip_prefix, required because some nodes may turn into tips, and the input tree has no node labels")
	
	results = merge_short_edges_CPP(Ntips,
									Nnodes,
									Nedges,
									tree_edge 			= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
									edge_length 		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
									edge_length_epsilon	= edge_length_epsilon,
									force_keep_tips		= force_keep_tips);
	
	# reformat results into a valid "phylo" object
	# note that some nodes may have turned into tips
	Nnew_tips		= results$Nnew_tips;
	Nnew_nodes		= results$Nnew_nodes;
	new2old_clade 	= results$new2old_clade + 1
	old_labels 		= c(tree$tip.label, tree$node.label)
	new_tree = list(Nnode 		= results$Nnew_nodes,
					tip.label 	= old_labels[new2old_clade[1:Nnew_tips]],
					node.label	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_clade[(Nnew_tips+1):(Nnew_tips+Nnew_nodes)]-Ntips]),
					edge 		= matrix(results$new_tree_edge,ncol=2,byrow=TRUE) + 1,
					edge.length = results$new_edge_length,
					root 		= results$root+1)
	# ensure all tips have labels (including those tips that used to be nodes)
	if((!force_keep_tips) && (!is.null(new_tip_prefix))){
		ex_node_tips = which(new2old_clade[1:Nnew_tips]>Ntips)
		new_tree$tip.label[ex_node_tips] = paste(new_tip_prefix, 1:length(ex_node_tips), sep = "")
	}
	class(new_tree) = "phylo";
	attr(new_tree,"order") = NULL

	return(list(tree 			= new_tree, 
				new2old_clade 	= new2old_clade,
				new2old_edge 	= results$new2old_edge+1,
				Nedges_removed	= Nedges-results$Nnew_edges));
}