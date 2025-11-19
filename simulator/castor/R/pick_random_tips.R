# pick Nrepeats random tip subsets of specific size
pick_random_tips = function(tree, size=1, Nsubsets=1, with_replacement=TRUE, drop_dims=TRUE){ 
	Ntips  = length(tree$tip.label);
	if(size>Ntips) stop(sprintf("ERROR: Requested size (%d) exceeds the number of tips in the tree (%d)",size,Ntips))
		
	random_tips = pick_random_tips_CPP(	Ntips				= Ntips,
										Nnodes				= tree$Nnode,
										Nedges 				= nrow(tree$edge),
										tree_edge 			= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
										Nrandoms			= size,
										Nsubsets			= Nsubsets,
										with_replacement	= with_replacement);
											
	# Note that random_tips is a 2D array of size Nsubsets x size in row-major format
	# So expand to 2D matrix if needed
	tips = 1L + as.integer(if(drop_dims && Nsubsets==1) random_tips else matrix(random_tips, ncol=size, byrow=TRUE));
	return(tips);
}