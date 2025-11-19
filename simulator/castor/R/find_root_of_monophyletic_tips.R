# find the node or tip which, if it became root, would make a set of target tips monophyletic
# specifically, the target tips would descend from a single child of the new root (if as_MRCA==FALSE) or the new root would be the MRCA of the target tips (if as_MRCA==TRUE)
# the returned value will be an integer within 1,..,(Ntips+Nnodes), or -1 on failure
# the input tree can be unrooted
find_root_of_monophyletic_tips = function(tree, monophyletic_tips, as_MRCA=TRUE, is_rooted=FALSE){ 	
	Ntips = length(tree$tip.label)
	if(is.character(monophyletic_tips)){
		# tips are given as names, not indices
		indices	= match(monophyletic_tips, tree$tip.label)
		if(any(is.na(indices))) stop(sprintf("ERROR: Some monophyletic_tips (e.g. '%s') were not found in the tree",monophyletic_tips[which(is.na(indices))[1]]))
		monophyletic_tips = indices
	}
	
	new_root = find_root_for_monophyletic_clade_CPP(Ntips		= Ntips,
													Nnodes		= tree$Nnode,
													Nedges		= nrow(tree$edge),
													tree_edge	= as.vector(t(tree$edge)) - 1,	# flatten in row-major format and adjust clade indices to 0-based
													is_rooted	= is_rooted,
													target_tips	= monophyletic_tips - 1,
													as_MRCA		= as_MRCA)
													
	return(if(new_root<0) NA else as.integer(new_root+1));
}