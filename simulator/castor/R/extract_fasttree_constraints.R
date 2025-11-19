# Given a rooted tree, extract binary constraints in FastTree alignment format
# Every internal bifurcating node with more than 2 descending tips will constitute an additional constraint
# Requirements:
#   The tree must be rooted; the root should be the unique node with no parent
#   The tree can include multifurcations as well as monofurcations, however only bifurcations are considered as constraints
extract_fasttree_constraints = function(tree){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	
	results = extract_fasttree_constraints_CPP(	Ntips 		= Ntips,
												Nnodes		= Nnodes,
												Nedges		= nrow(tree$edge),
												tree_edge	= as.vector(t(tree$edge))-1);	# flatten in row-major format and make indices 0-based
	Nconstraints = results$Nconstraints
	constraints	 = matrix(results$constraints, ncol=Nconstraints, byrow=TRUE) # unflatten constraints into a 2D matrix
	return(list(Nconstraints	= Nconstraints,
				constraints 	= constraints,
				constraint2node	= results$constraint2node+1));
}