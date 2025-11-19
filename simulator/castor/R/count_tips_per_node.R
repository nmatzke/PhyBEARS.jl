# For each node in a tree, count the number of descending (indirectly or directly) tips for each node
# Requirements:
#   The tree must be rooted; the root should be the unique node with no parent
#   The tree can include multifurcations as well as monofurcations
count_tips_per_node = function(tree){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	
	results = get_total_tip_count_per_node_CPP(	Ntips 			= Ntips,
												Nnodes			= Nnodes,
												Nedges			= nrow(tree$edge),
												tree_edge 		= as.vector(t(tree$edge))-1);	# flatten in row-major format and make indices 0-based

	return(results);
}