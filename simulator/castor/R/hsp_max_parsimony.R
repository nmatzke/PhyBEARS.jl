# Maximum parsimony hidden state prediction of discrete states for tips on a tree.
hsp_max_parsimony = function(	tree, 
								tip_states, 	# integer vector of size Ntips
								Nstates				= NULL, 
								transition_costs	= "all_equal", 
								edge_exponent		= 0.0, 
								weight_by_scenarios = TRUE,
								check_input 		= TRUE){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	
	# basic error checking
	if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
	if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be integers"))
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
	if(is.null(Nstates)) Nstates = max(known_tip_states);

	# some more error checking
	if(check_input){
		min_tip_state = min(known_tip_states)
		max_tip_state = max(known_tip_states)
		if((min_tip_state<1) || (max_tip_state>Nstates)) stop(sprintf("ERROR: Non-NA tip_states must be integers between 1 and %d, but found values between %d and %d",Nstates,min_tip_state,max_tip_state))
	}

	# perform ancestral state reconstruction on known_subtree
	asr_results = asr_max_parsimony(known_subtree, 
									known_tip_states, 
									Nstates = Nstates, 
									transition_costs = transition_costs, 
									edge_exponent = edge_exponent, 
									weight_by_scenarios = weight_by_scenarios);
	if(!asr_results$success) return(list(success=FALSE, likelihoods=NULL))
	Nstates = ncol(asr_results$ancestral_likelihoods);
													
	# forward-project reconstructions to tips with hidden state
	likelihoods = matrix(0, nrow=(Ntips+Nnodes), ncol=Nstates);
	likelihoods[known2all_tips, ] = 0.0;
	likelihoods[cbind(known2all_tips, known_tip_states)] = 1.0;
	likelihoods[known2all_nodes+Ntips, ] 		= asr_results$ancestral_likelihoods;
	likelihoods_known 							= rep(FALSE, times=(Ntips+Nnodes))
	likelihoods_known[known2all_tips] 			= TRUE;
	likelihoods_known[known2all_nodes + Ntips] 	= TRUE;
	likelihoods = apply_attributes_to_descendants_CPP(	Ntips				= Ntips,
														Nnodes				= Nnodes,
														Nedges				= nrow(tree$edge),
														Nattributes			= Nstates,
														tree_edge			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
														attributes_known	= likelihoods_known,
														attributes			= as.vector(t(likelihoods))); 	# flatten in row-major format	
	likelihoods = matrix(likelihoods, ncol=Nstates, byrow=TRUE); # unflatten returned table
	colnames(likelihoods) = colnames(asr_results$ancestral_likelihoods);
		
	return(list(success		= TRUE,
				likelihoods	= likelihoods, 
				total_cost	= asr_results$total_cost))
}


