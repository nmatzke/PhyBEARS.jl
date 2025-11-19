# Given a state for each clade (tip or node), and a list of clade pairs (A[], B[]),
# count the number of state transitions along the shortest path connecting every clade-pair
count_transitions_between_clades = function(tree, 
											A, 		# integer vector of length NP, listing the first clade (tip or node) of each pair
											B, 		# integer vector of length NP, listing the second clade (tip or node) of each pair
											states,	# integer vector of length Nclades (=Ntips+Nnodes), specifying the state at each clade
											check_input = TRUE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode
	
	# basic input checking
	A = map_tip_or_node_names_to_indices(tree, A, type='both', list_title='A', check_input=check_input)
	B = map_tip_or_node_names_to_indices(tree, B, type='both', list_title='B', check_input=check_input)
	if(length(states)!=Ntips+Nnodes) stop(sprintf("Length of states[] must be equal to the number of tips + nodes (%d)",Ntips+Nnodes))
	if(check_input && (!is.null(names(states)))){
		if(any(names(states)[1:Ntips]!=tree$tip.label)) stop(sprintf("The first %d names of states[] do not match the tree's tip labels (make sure they are in the same order)",Ntips))
		if((!is.null(tree$node.label)) && any(names(states)[(Ntips+1):(Ntips+Nnodes)]!=tree$node.label)) stop(sprintf("The last %d names of states[] do not match the tree's node labels (make sure they are in the same order)",Nnodes))
	}
	
	results = count_transitions_between_clades_CPP(	Ntips 			= Ntips,
													Nnodes 			= Nnodes,
													Nedges 			= nrow(tree$edge),
													tree_edge 		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													clade_states	= states,
													cladesA			= A-1,
													cladesB			= B-1);	
	return(results$pair2Ntransitions)
}
