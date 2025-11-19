# For each clade (tip & node) in a tree, find the farthest tip (in terms of cumulative branch length).
# Optionally, the search can be restricted to descending tips.
# Optionally, the search can also be restricted to a subset of target tips.
# Requirements:
#   The input tree must be rooted (root will be determined automatically, as the node that has no incoming edge)
#   The input tree can be multifurcating and/or monofurcating
find_farthest_tips = function(tree, only_descending_tips=FALSE, target_tips=NULL, as_edge_counts=FALSE, check_input=TRUE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	if(!is.null(target_tips)) target_tips = map_tip_or_node_names_to_indices(tree, A=target_tips, type='tip', list_title='target_tips', check_input=TRUE)
	results = get_farthest_tip_per_clade_CPP(	Ntips					= Ntips,
												Nnodes					= Nnodes,
												Nedges					= nrow(tree$edge),
												tree_edge				= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												edge_length				= (if(as_edge_counts || is.null(tree$edge.length)) numeric() else tree$edge.length),
												onlyToTips				= (if(is.null(target_tips)) integer() else (target_tips-1)),
												only_descending_tips	= only_descending_tips,
												verbose					= FALSE,
												verbose_prefix			= "");
	if(only_descending_tips){
		results$farthest_tips 		= results$farthest_descending_tips;
		results$farthest_distances 	= results$farthest_descending_distances;
	}
	results$farthest_tips[results$farthest_tips<0] = NA;
	return(list(farthest_tip_per_tip		= (results$farthest_tips[1:Ntips] + 1), 
				farthest_tip_per_node		= (results$farthest_tips[(Ntips+1):(Ntips+Nnodes)] + 1), 
				farthest_distance_per_tip	= results$farthest_distances[1:Ntips],
				farthest_distance_per_node	= results$farthest_distances[(Ntips+1):(Ntips+Nnodes)]));
}