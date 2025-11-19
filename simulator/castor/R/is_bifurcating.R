# Test whether the tree is strictly bifurcating, i.e. every node has exactly 2 children
# Tree need not be rooted
is_bifurcating = function(tree){
	Nchildren = get_child_count_per_node_CPP(	Ntips			= length(tree$tip.label),
												Nnodes			= tree$Nnode,
												Nedges			= nrow(tree$edge),
												tree_edge		= as.vector(t(tree$edge))-1);
	
	return(all(Nchildren==2));
}