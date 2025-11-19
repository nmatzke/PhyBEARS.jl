# Test whether a set of tips is monophyletic
# Tree must be rooted
is_monophyletic = function(tree, focal_tips, check_input=TRUE){
	Ntips  		= length(tree$tip.label);
	Nnodes 		= tree$Nnode;
	focal_tips 	= map_tip_or_node_names_to_indices(tree, focal_tips, type="tip", list_title="focal_tips", check_input=check_input);
	monophyletic = is_monophyletic_tip_set_CPP(	Ntips		= Ntips,
												Nnodes		= Nnodes,
												Nedges		= nrow(tree$edge),
												tree_edge	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												focal_tips	= focal_tips-1);
	return(monophyletic);
}