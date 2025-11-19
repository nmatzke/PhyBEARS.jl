# get distances between all possible pairs of clades
# returns a matrix of size NC x NC, where NC is the number of clades in the tree (Ntips+Nnodes) or the length of only_clades[] (if specified)
get_all_pairwise_distances = function(	tree, 
										only_clades		= NULL, 	# optional 1D integer vector of clade indices (i.e. in 1,...,Ntips+Nnodes) to which to restrict calculation of pairwise distances. For example, if this is c(1,3,4), then only the distances between clades 1 & 3, 1 & 4 and 3 & 4 are returned. These will correspond to the rows & columns of the returned distance matrix.
										as_edge_counts	= FALSE,
										check_input		= TRUE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	
	# check and reformat input
	if(is.null(only_clades)){
		only_clades = 1:(Ntips+Nnodes)
	}else{
		if((!is.character(only_clades)) && (!is.numeric(only_clades))) stop("ERROR: List of focal clades must be a character or integer vector")
		if(is.character(only_clades)){
			# map names to indices
			name2clade = 1:(Ntips+Nnodes); 
			names(name2clade) = c(tree$tip.label,tree$node.label);
			only_clades_i = name2clade[only_clades]; 
			if(check_input && any(is.na(only_clades_i))) stop(sprintf("ERROR: Unknown tip or node name '%s'",only_clades[which(is.na(only_clades_i))[1]]));
			only_clades = only_clades_i;
		}else if(check_input){
			minC = min(only_clades); maxC = max(only_clades);
			if(minC<1 || maxC>(Ntips+Nnodes)) stop(sprintf("ERROR: only_clades[] must contain values between 1 and Ntips+Nnodes (%d); instead, found values from %d to %d",Ntips+Nnodes,minC,maxC));
		}
	}

	# calculate all pairwise distances													
	distances = get_distance_matrix_between_clades_CPP(	Ntips 			= Ntips,
														Nnodes 			= Nnodes,
														Nedges 			= nrow(tree$edge),
														tree_edge 		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
														edge_length		= (if(as_edge_counts || is.null(tree$edge.length)) numeric() else tree$edge.length),
														focal_clades	= only_clades-1,
														verbose			= FALSE,
														verbose_prefix	= "")
	return(distances)
}
