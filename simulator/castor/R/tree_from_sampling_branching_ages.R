# generate a random bifurcating tree, based on a list of tip-sampling & node-branching ages
# The oldest branching age will be the root age
# Nodes and tips will be indexed in chronological order (i.e. in order of decreasing age). In particular, node 0 will be the root.
tree_from_sampling_branching_ages = function(	sampling_ages,				# numeric vector of size Ntips, listing tip sampling ages (time before present) in ascending order. Note that Ntips must be equal to Nnodes+1.
												branching_ages, 			# numeric vector of size Nnodes, listing branching ages (time before present) in ascending order (i.e. root last)
												tip_basename	= "",		# basename for tips (e.g. "tip."). 
												node_basename	= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
												edge_basename	= NULL){	# basename for edge (e.g. "edge."). If NULL, then edges will not have any labels.
	if(length(sampling_ages)!=length(branching_ages)+1) return(list(success=FALSE, error=sprintf("Number of samplings (%d) is inconsistent with number of branchings (%d); expected Nsamplings = Nbranchines+1",length(sampling_ages),length(branching_ages))))

	results = get_tree_from_sampling_branching_ages_CPP(sampling_ages = sampling_ages, branching_ages = branching_ages)
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