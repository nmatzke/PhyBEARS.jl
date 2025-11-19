# Given a phylogenetic tree, find the root (node with no incoming edge)
# If no root is found, return NA.
# Requirements:
#   The input tree can be multifurcating and/or monofurcating
find_root = function(tree){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	root = 1+get_root_clade_CPP(Ntips,
								Nnodes,
								Nedges 		= nrow(tree$edge),
								tree_edge	= as.vector(t(tree$edge)) - 1)  # flatten in row-major format and adjust clade indices to 0-based
	if(root<=0){ return(NA); }
	else{ return(root); }
}