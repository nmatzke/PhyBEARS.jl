# Modify the times of specific tips & nodes, by adding or subtracting specific values
# Here, "time" refers to time since the root
# Excessive shifting can result in negative edge lengths, which are corrected in a variety of alternative ways (see option negative_edge_lengths). 
# However, to avoid changing the overall span of the tree (root age and tip times), you should not shift a clade beyond the boundaries of the tree (i.e., resulting in a negative time or a time beyond its descending tips).
# The input tree must be rooted
# The input tree does not need to be ultrametric, but edge lengths are interpreted as time
shift_clade_times = function(	tree,
								clades_to_shift,
								time_shifts,
								shift_descendants 		= FALSE, # if true, then the subclade descending from a shifted clade will move along with it, thus shifting the time of all descendants by the same amount. If false, the descending tips & nodes retain their original time (unless negative edges are created, see option negative_edge_lengths).
								negative_edge_lengths 	= "error"){ # whether and how to fix negative edge lengths resulting from excessive shifting. See negative_edge_lengths_options[] for possible options
	Ntips 	= length(tree$tip.label)
	Nnodes	= tree$Nnode
	Nedges	= nrow(tree$edge)
	
	# basic input checking
	negative_edge_lengths_options = c("error", "allow", "move_all_descendants", "move_all_ancestors", "move_child", "move_parent")
	if(!(negative_edge_lengths %in% negative_edge_lengths_options)) return(list(success=FALSE, error=sprintf("Invalid choice '%s' for negative_edge_lengths. Expected one of '%s'",negative_edge_lengths,paste(negative_edge_lengths_options,collapse="', '"))))
	clades_to_shift = map_tip_or_node_names_to_indices(tree, A=clades_to_shift, type='clade', list_title='clades_to_shift', check_input=TRUE)
	if(length(clades_to_shift)!=length(time_shifts)) return(list(success=FALSE, error=sprintf("Length of clades_to_shift (%d) does not match length of time_shifts (%d)",length(clades_to_shift),length(time_shifts))))
	
	results = shift_clade_times_CPP(Ntips 					= length(tree$tip.label),
									Nnodes					= tree$Nnode,
									Nedges					= nrow(tree$edge),
									tree_edge				= as.vector(t(tree$edge)) - 1,
									edge_length				= (if(is.null(tree$edge.length)) rep(1,Nedges) else tree$edge.length),
									clades_to_shift			= clades_to_shift-1,
									time_shifts				= time_shifts,
									shift_descendants		= shift_descendants,
									negative_edge_lengths	= negative_edge_lengths)
	if(!results$success) return(list(success=FALSE, error=results$error))
	tree$edge.length = 	results$new_edge_length
	
	return(list(success	= TRUE,
				tree	= tree))
}