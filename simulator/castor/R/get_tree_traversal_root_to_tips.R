# returns a list of all nodes (and optionally tips) of a tree, such that each node appears prior to its children
# also returns a list mapping nodes to their outgoing edges (as listed in tree$edge)
# The tree can include multifurcations as well as monofurcations
# The tree must be rooted (i.e. there must exist a node with no incoming edge)
# Returned values:
#	queue: A 1D vector of integers in 1:(Ntips+Nnodes) if include_tips==TRUE, or (Ntips+1):(Ntips+Nnodes) if include_tips==FALSE
#	node2first_edge[p] will be an index pointing node p (p=1:Nnodes) to edges[]
#	node2last_edge[p] will be an index pointing node p (p=1:Nnodes) to edges[]
# 	edges[] will be a list of edge indices (i.e. with values in 1:Nedges), such that edges[node2first_edge[p]],...,edges[node2last_edge[p]] is the set of edges leaving node p
get_tree_traversal_root_to_tips = function(tree, include_tips){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode
	Nedges = nrow(tree$edge)
			
	# use CPP function, don't forget to change indices between 0-based in C++ and 1-based in R
	results = get_tree_traversal_root_to_tips_CPP(	Ntips 			= Ntips,
													Nnodes 			= Nnodes,
													Nedges			= Nedges,
													tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
													include_tips	= include_tips);
	
	# update indices from 0-based to 1-based and wrap results into a list
	return(list(queue=results$queue+1, 
				node2first_edge=results$node2first_edge+1, 
				node2last_edge=results$node2last_edge+1, 
				edges=results$edges+1))	
}