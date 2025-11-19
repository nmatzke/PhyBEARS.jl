# place the root of a tree in the middle of an edge (creating a new node)
root_in_edge = function(tree, root_edge, location=0.5, new_root_name=NULL, collapse_monofurcations=TRUE){
	Ntips 		= length(tree$tip.label)
	Nnodes 		= tree$Nnode
	Nedges 		= nrow(tree$edge)
	new_root 	= Ntips+Nnodes+1 # append new root at the end of the node list
	location	= max(0,min(1,location))
	
	# split edge (downstream half becomes a new edge)
	tree$edge = rbind(tree$edge, c(new_root, tree$edge[root_edge,2]))
	tree$edge[root_edge,2] = new_root
	if(!is.null(tree$edge.length)){
		tree$edge.length = c(tree$edge.length, (1-location)*tree$edge.length[root_edge])
		tree$edge.length[root_edge] = location*tree$edge.length[root_edge]
	}
	
	# add new node
	tree$Nnode = tree$Nnode + 1
	if((!is.null(tree$node.label)) && (!is.null(new_root_name))) tree$node.label = c(tree$node.label,new_root_name)
	
	# re-root at the new root, update indices to meet root indexing convention
	tree = root_at_node(tree, new_root_node=new_root-Ntips, update_indices=TRUE)
	
	# collapse monofurcations (e.g. resulting at old root), if needed
	if(collapse_monofurcations){
		tree = collapse_monofurcations(tree, force_keep_root=FALSE, as_edge_counts=FALSE)$tree
	}
	attr(tree,"order") = NULL
	
	return(tree);
}

