# Read tree from a Newick-formatted string or file
read_tree = function(	string="",
						file="",
						edge_order				= "cladewise", # how to order edges. Options are "cladewise" (i.e. depth-first-search) or "pruningwise" (i.e. iterating through edge[] leads a post-order traversal) or "none" (unspecified, i.e. depending on how the tree was written in the file)
						include_edge_lengths 	= TRUE, 
						look_for_edge_labels 	= FALSE, 
						look_for_edge_numbers 	= FALSE, 
						include_node_labels 	= TRUE, 
						underscores_as_blanks 	= FALSE, 
						check_label_uniqueness 	= FALSE,
						interpret_quotes		= FALSE, 	# whether to interpret quotes as delimiters of tip/node names, rather than reading quotes just like any other character
						trim_white				= TRUE){	# whether to trim flanking whitespace from tip, node and edge labels
	if(file!=""){
		if(string!="") stop("ERROR: Either string or file must be specified, but not both")
		string = readChar(file, file.info(file)$size)
	}
	if(!(edge_order %in% c("cladewise", "pruningwise"))) stop(sprintf("ERROR: Invalid option '%s' for edge order. Use either 'cladewise' or 'pruningwise'",edge_order))

	results = read_Newick_string_CPP(	input 					= string, 
										underscores_as_blanks 	= underscores_as_blanks,
										interpret_quotes		= interpret_quotes,
										look_for_edge_names		= look_for_edge_labels,
										look_for_edge_numbers	= look_for_edge_numbers);
	if(!results$success) stop(sprintf("ERROR: Could not parse Newick string: %s",results$error))
	if(check_label_uniqueness){
		duplicates = which(duplicated(results$tip_names))
		if(length(duplicates)>0) stop(sprintf("ERROR: Duplicate tip labels (e.g. '%s') found in input tree",results$tip_names[duplicates[1]]))
	}
	
	if(trim_white){
		results$tip_names = trimws(results$tip_names, which="both")
		if(include_node_labels) results$node_names  = trimws(results$node_names, which="both")
		if(look_for_edge_labels) results$edge_names = trimws(results$edge_names, which="both")
	}
	
	tree = list(Nnode 		= results$Nnodes,
				tip.label 	= results$tip_names,
				node.label 	= (if((!include_node_labels) || all(results$node_names=="")) NULL else results$node_names),
				edge 		= matrix(results$tree_edge+1L, ncol=2, byrow=TRUE), # unflatten row-major array
				edge.length = (if((!include_edge_lengths) || all(is.nan(results$edge_lengths))) NULL else results$edge_lengths),
				edge.label 	= (if((!look_for_edge_labels) || all(results$edge_names=="")) NULL else results$edge_names),
				edge.number	= (if((!look_for_edge_numbers) || all(results$edge_numbers<0)) NULL else results$edge_numbers),
				root 		= results$root+1L,
				root.edge	= (if(is.nan(results$root_edge)) NULL else results$root_edge))
	if(!is.null(tree$edge.number)) tree$edge.number[tree$edge.number<0] = NA;
	class(tree) = "phylo";
	attr(tree,"order") = "cladewise"
	
	if(edge_order=="pruningwise"){
		tree = reorder_tree_edges(tree, root_to_tips=FALSE, depth_first_search=TRUE, index_only=FALSE)
	}
	
	return(tree)
}