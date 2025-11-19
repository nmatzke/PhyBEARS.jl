# Write tree in Newick format into a string or to a file
# The tree does not need to be rooted. In that case, it is rooted temporarily at the first node.
write_tree = function(	tree, 
						file	= "", 
						append	= FALSE,
						digits	= 10, 
						quoting	= 0,	# whether and how to quote tip & node & edge names. 0:no quoting, 1:single quoting always, 2:double quoting always, -1:quote only when needed and prefer single quotes if possible, -2:quote only when needed and prefer double quotes if possible
						include_edge_labels = FALSE,
						include_edge_numbers = FALSE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	result = tree_to_Newick_string_CPP(	Ntips,
										Nnodes,
										Nedges				= nrow(tree$edge),
										tree_edge			= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
										edge_length			= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
										tip_labels			= (if(is.null(tree$tip.label)) character() else tree$tip.label),
										node_labels			= (if(is.null(tree$node.label)) character() else tree$node.label),
										edge_labels			= (if(is.null(tree$edge.label) || (!include_edge_labels)) character() else tree$edge.label),
										edge_numbers		= (if(is.null(tree$edge.number) || (!include_edge_numbers)) numeric() else tree$edge.number),
										digits				= digits,
										root_edge_length	= (if(is.null(tree$root.edge)) -1 else tree$root.edge),
										quoting				= quoting);
	if(file!=""){
		cat(result, file=file, append=append);
	}else{
		return(result);
	}
}