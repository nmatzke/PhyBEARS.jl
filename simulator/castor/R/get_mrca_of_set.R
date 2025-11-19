# Calculate most-recent-common-ancestor (MRCA) for a set of tips/nodes
get_mrca_of_set = function(tree, descendants){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	if(is.character(descendants)){
		# tips/nodes are given as names, not indices
		name2clade 			= 1:(Ntips+Nnodes); 
		names(name2clade) 	= c(tree$tip.label,tree$node.label);
		descendants 		= name2clade[descendants]; 
		descendants 		= descendants[!is.na(descendants)];
	}
	mrca = get_most_recent_common_ancestor_CPP(	Ntips			= Ntips,
												Nnodes			= Nnodes,
												Nedges			= nrow(tree$edge),
												tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												descendants		= descendants-1);
	if(mrca<0){ return(NA); } # something went wrong
	else{ return(mrca+1); }
}