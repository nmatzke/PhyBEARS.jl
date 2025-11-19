# adjust edge directions such that the new root is at the midpoint (node with smallest max_distance_to_a_tip)
# Note that the number of tips & nodes remains the same
# If update_indices==FALSE, then tip & node indices also remain the same
root_at_midpoint = function(tree, 
							split_edge		= FALSE,
							update_indices	= TRUE,
							as_edge_counts 	= FALSE,	# calculate distances in terms of cumulative edge counts (as as if each edge had length 1)
							is_rooted		= FALSE){ 	# if TRUE, the caller guarantees that the input tree is rooted
	Ntips 	= length(tree$tip.label)
	Nnodes	= tree$Nnode
	Nedges	= nrow(tree$edge)
	
	# root arbitrarily if needed
	if(!is_rooted){
		tree = root_at_node(tree, 1, update_indices=FALSE)
	}
	
	if(!split_edge){
		# figure out midpoint node
		distances = get_farthest_tip_per_clade_CPP(	Ntips					= Ntips,
													Nnodes					= Nnodes,
													Nedges					= nrow(tree$edge),
													tree_edge				= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
													edge_length				= (if(as_edge_counts || is.null(tree$edge.length)) numeric() else tree$edge.length),
													onlyToTips				= integer(),
													only_descending_tips	= FALSE,
													verbose					= FALSE,
													verbose_prefix			= "")
		new_root_node = which.min(distances$farthest_distances[(Ntips+1):(Ntips+Nnodes)])
	
		# place root at midpoint node
		new_edges = root_tree_at_node_CPP(	Ntips			= Ntips,
											Nnodes			= Nnodes,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
											new_root_node	= new_root_node-1);
		new_edges = matrix(new_edges+1, ncol=2, byrow=TRUE); # unflatten returned table and shift clade indices to 1-based
		tree$edge = new_edges;
		tree$root = Ntips + new_root_node;
		tree$root.edge = 0
	
		# update node indices if required
		correct_root_node = 1; # correct index that the root node should have
		if(update_indices && (new_root_node!=correct_root_node)){
			# swap indices with wrong node
			temp_root_index = 0;
			tree$edge[tree$edge==(Ntips+new_root_node)] 	= temp_root_index;
			tree$edge[tree$edge==(Ntips+correct_root_node)] = Ntips+new_root_node;
			tree$edge[tree$edge==temp_root_index] 			= Ntips+correct_root_node;
			tree$root 										= Ntips+correct_root_node;
		
			if(!is.null(tree$node.label)){
				root_label 							= tree$node.label[new_root_node];
				tree$node.label[new_root_node] 		= tree$node.label[correct_root_node];
				tree$node.label[correct_root_node] 	= root_label;
			}
		}
		attr(tree,"order") = NULL

	}else{
		# figure out midpoint edge
		distances = get_farthest_tips_per_edge_CPP(	Ntips		= Ntips,
													Nnodes		= Nnodes,
													Nedges		= nrow(tree$edge),
													tree_edge	= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
													edge_length	= (if(as_edge_counts || is.null(tree$edge.length)) numeric() else tree$edge.length),
													onlyToTips	= integer())
		edge2descending_distance	= distances$farthest_descending_distances 
		edge2upstream_distance		= distances$farthest_upstream_distances
		edge2length					= (if(as_edge_counts || is.null(tree$edge.length)) rep(1,Nedges) else tree$edge.length)
		max_tip_span_per_edge		= edge2descending_distance + edge2upstream_distance + edge2length
		midpoint_edges				= which((edge2descending_distance+edge2length>=edge2upstream_distance) & (edge2upstream_distance+edge2length>=edge2descending_distance))
		new_root_edge 				= midpoint_edges[which.max(max_tip_span_per_edge[midpoint_edges])] # edge within which we should place the new root
		location					= (if(edge2length[new_root_edge]==0) 0.5 else (edge2descending_distance[new_root_edge]+edge2length[new_root_edge]-edge2upstream_distance[new_root_edge])/(2*edge2length[new_root_edge]))
		tree = root_in_edge(tree, 
							root_edge	= new_root_edge, 
							location	= location, 
							collapse_monofurcations = FALSE)
	}
	
	return(tree);
}