# generate a random tree, based on a list of branching ages
# The oldest branching age will thus be the root age
# Tips are guaranteed to be connected in random order, i.e. this function can also be used to connect a random set of tips into a tree.
# Nodes will be indexed in chronological order (i.e. in order of decreasing age). In particular, node 0 will be the root.
tree_from_branching_ages = function(branching_ages, 			# numeric vector of size Nnodes, listing branching ages (time before present) in ascending order (i.e. root last)
									tip_basename	= "",		# basename for tips (e.g. "tip."). 
									node_basename	= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
									edge_basename	= NULL){	# basename for edge (e.g. "edge."). If NULL, then edges will not have any labels.

	results = get_tree_from_branching_ages_CPP(branching_ages = branching_ages)
	if(!results$success) return(list(success=FALSE, error=results$error)); # something went wrong
	Ntips	= results$Ntips
	Nnodes 	= results$Nnodes
	Nedges 	= results$Nedges
	tree = list(Nnode 		= Nnodes,
				tip.label 	= paste(tip_basename, 1:Ntips, sep=""),
				node.label 	= (if(is.null(node_basename)) NULL else paste(node_basename, 1:Nnodes, sep="")),
				edge.label 	= (if(is.null(edge_basename)) NULL else paste(edge_basename, 1:Nedges, sep="")),
				edge 		= matrix(results$tree_edge,ncol=2,byrow=TRUE) + 1,
				edge.length = results$edge_length,
				root 		= results$root+1)
	class(tree) = "phylo";
	attr(tree,"order") = NULL

	return(list(success = TRUE, tree = tree));
	
}