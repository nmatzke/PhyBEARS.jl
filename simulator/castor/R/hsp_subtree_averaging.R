# Hidden state prediction for a continuous trait via subtree averaging, for tips on a tree.
hsp_subtree_averaging = function(	tree, 
									tip_states, 	# numeric vector of size Ntips. Can include NA.
									check_input = TRUE){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	
	# basic error checking
	if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
	if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be a numeric vector"))
	if(check_input){
		if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
	}
	
	# find known_tips, extract known_subtree and synchronize known_tip_states with known_subtree
	known_tips = which(!is.na(tip_states));
	if(length(known_tips)==0) stop("ERROR: All tip states are hidden");
	extraction	 		= get_subtree_with_tips(tree, only_tips=known_tips, omit_tips=FALSE, collapse_monofurcations=TRUE, force_keep_root=TRUE);
	known_subtree		= extraction$subtree;
	known2all_tips		= extraction$new2old_tip
	known2all_nodes		= extraction$new2old_node
	known_tip_states	= tip_states[known2all_tips]

	# perform ancestral state reconstruction on known_subtree
	asr_results = asr_subtree_averaging(known_subtree, 
										tip_states 	= known_tip_states, 
										check_input = FALSE);
	if(!asr_results$success) return(list(success=FALSE))
													
	# forward-project reconstructions to tips with hidden state
	states 									= rep(0, times=Ntips+Nnodes);
	states[known2all_tips] 					= known_tip_states;
	states[known2all_nodes+Ntips] 			= asr_results$ancestral_states;
	states_known 							= rep(FALSE, times=(Ntips+Nnodes))
	states_known[known2all_tips]			= TRUE;
	states_known[known2all_nodes + Ntips] 	= TRUE;
	states = apply_BM_parsimony_to_missing_clades_CPP(	Ntips			= Ntips,
														Nnodes			= Nnodes,
														Nedges			= nrow(tree$edge),
														tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
														states_known	= states_known,
														states			= states);
		
	return(list(success = TRUE, states = states))
}


