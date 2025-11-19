# Squared-change parsimony reconstruction of continuous ancestral states
#    Maddison (1991). Squared-change parsimony reconstructions of ancestral states for continuous-valued characters on a phylogenetic tree. Systematic Zoology. 40:304-314.
# The function reroots the tree at each node, to get the global squared-change-max-parsimony state for that node.
#  Requirements:
# 	Tree can be multifurcating, and can also include nodes with a single child
# 	Tree can also include edges with length zero (will be adjusted internally to some small epsilon if weighted==TRUE).
#	Tree must be rooted.
asr_squared_change_parsimony = function(tree, 
										tip_states, 		# numeric vector of size Ntips
										weighted	= TRUE,
										check_input	= TRUE){
	Ntips  = length(tree$tip.label)

	# basic error checking
	if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
	if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be numeric"))
	if(check_input){
		if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
	}

	results = ASR_via_squared_change_parsimony_CPP(	Ntips		= Ntips,
													Nnodes		= tree$Nnode,
													Nedges		= nrow(tree$edge),
													tree_edge	= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													edge_length	= (if((!weighted) || (is.null(tree$edge.length))) numeric() else tree$edge.length),
													tip_states	= tip_states,
													global		= TRUE);
	return(list(success							= TRUE,
				total_sum_of_squared_changes 	= results$TSS,
				ancestral_states				= results$ancestral_states));
}