# extract the subtree spanning the set of tips specified
# this function guarantees that the root will be kept in the subtree
# if the original tree$edge_length was empty or NULL, then subtree$edge_length[e] will be the number of combined edges making up the new edge e
get_subtree_with_tips = function(tree, only_tips=NULL, omit_tips=NULL, collapse_monofurcations=TRUE, force_keep_root=FALSE){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
		
	# figure out which tips to keep
	if(!is.null(only_tips)){ keep_tip = rep(FALSE,times=Ntips); }
	else{ keep_tip = rep(TRUE,times=Ntips); }
	if(!is.null(only_tips)){
		if(is.character(only_tips)){
			only_tips = match(only_tips, tree$tip.label)
			only_tips = only_tips[!is.na(only_tips)]
			keep_tip[only_tips] = TRUE;
		}else{
			keep_tip[only_tips[only_tips<=Ntips]] = TRUE;
		}
	}
	if(!is.null(omit_tips)){
		if(is.character(omit_tips)){
			omit_tips = match(omit_tips, tree$tip.label)
			omit_tips = omit_tips[!is.na(omit_tips)]
			keep_tip[omit_tips] = FALSE;
		}else{
			keep_tip[omit_tips[omit_tips<=Ntips]] = FALSE;
		}
	
	}
	tips_to_keep = which(keep_tip);
	if(length(tips_to_keep)==0) stop("ERROR: No tips kept (all filtered out)")
	if(length(tips_to_keep)==1) stop("ERROR: Only one tip to be kept, so the resulting tree would not be an actual tree")
	
	if(length(tips_to_keep)==Ntips){
		# special case: keeping full tree
		subtree 		= tree
		root_shift 		= 0
		new2old_clade 	= seq_len(Ntips+Nnodes)
		old2new_clade 	= seq_len(Ntips+Nnodes)
		Ntips_kept 		= Ntips
		Nnodes_kept 	= Nnodes
	}else{
		# extract subtree
		results = get_subtree_with_specific_tips_CPP(	Ntips					= Ntips,
														Nnodes					= Nnodes,
														Nedges 					= nrow(tree$edge),
														tree_edge 				= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
														edge_length				= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
														tips_to_keep			= tips_to_keep-1,
														collapse_monofurcations = collapse_monofurcations,
														force_keep_root			= force_keep_root);
		Ntips_kept  	= results$Ntips_kept
		Nnodes_kept 	= results$Nnodes_kept
		new2old_clade 	= results$new2old_clade + 1 # switch to 1-based indices
		root_shift		= results$root_shift
		subtree = list(	Nnode 		= Nnodes_kept,
						tip.label 	= (if(is.null(tree$tip.label)) NULL else tree$tip.label[new2old_clade[1:Ntips_kept]]),
						node.label 	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips]),
						edge 		= matrix(as.integer(results$new_tree_edge),ncol=2,byrow=TRUE) + 1L,
						edge.length = results$new_edge_length,
						root 		= results$new_root+1L,
						root.edge	= (if(results$old_stem_edge<0) (if(!is.null(tree$root.edge)) tree$root.edge else NULL) else (if(!is.null(tree$edge.length)) tree$edge.length[results$old_stem_edge+1L] else NULL)))
		class(subtree) = "phylo"
		attr(subtree,"order") = NULL
		
		# calculate the inverse mapping old-->new clade. Some entries might be 0, i.e. for old clades that were lost.
		old2new_clade = numeric(Ntips+Nnodes)
		old2new_clade[new2old_clade] = seq_len(Ntips_kept+Nnodes_kept)
	}
	
	return(list(subtree			= subtree, 
				root_shift		= root_shift, # distance between old & new root (will always be non-negative)
				new2old_clade	= new2old_clade,
				new2old_tip		= new2old_clade[1:Ntips_kept], 
				new2old_node	= new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips,
				old2new_clade	= old2new_clade,
				old2new_tip		= old2new_clade[1:Ntips],
				old2new_node	= old2new_clade[(Ntips+1):(Ntips+Nnodes)]-Ntips_kept));
}