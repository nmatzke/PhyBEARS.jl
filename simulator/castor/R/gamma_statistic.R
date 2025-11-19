# Given a rooted ultrametric tree, calculate the gamma statistic
# Note that the tree is assumed to be ultrametric, i.e. all tips are assumed to have effectively age 0 (any deviations from this are ignored)
# Requirements:
#	The input tree must be rooted
#   The input tree can be multifurcating and/or monofurcating
gamma_statistic = function(tree){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	gamma_stat = get_gamma_statistic_CPP(	Ntips,
											Nnodes,
											Nedges 		= nrow(tree$edge),
											tree_edge	= as.vector(t(tree$edge)) - 1,  # flatten in row-major format and adjust clade indices to 0-based
											edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length));
											
	return(gamma_stat)
}