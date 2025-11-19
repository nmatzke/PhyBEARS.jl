# Calculate most-recent-common-ancestors (MRCA) for pairs of clades (A[] vs B[])
# There's some initial overhead involved with this function, but for large number of clade pairs this becomes more efficient
# Time complexity is O(Ntips+Nnodes).
# Returns an IntegerVector of size Npairs, with each entry being the clade index of the MRCA of the pair
# If one clade is descendant of the other clade, the latter will be returned as MRCA
# Requirements:
#   The input tree must be rooted (root will be determined automatically, as the node that has no incoming edge)
#   The input tree can be multifurcating and/or monofurcating
get_pairwise_mrcas = function(tree, A, B, check_input=TRUE){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	
	# check and reformat input
	if((!is.character(A)) && (!is.numeric(A))) stop("ERROR: List of tips/nodes A must be a character or integer vector")
	if((!is.character(B)) && (!is.numeric(B))) stop("ERROR: List of tips/nodes B must be a character or integer vector")
	if(length(A)!=length(B)) stop(sprintf("ERROR: Lists A & B must have equal lengths; instead, A has length %d and B has length %d",length(A),length(B)))
	if(is.character(A) || is.character(B)){
		name2clade = 1:(Ntips+Nnodes); 
		names(name2clade) = c(tree$tip.label,tree$node.label);
	}
	if(is.character(A)){
		Ai = name2clade[A]; 
		if(check_input && any(is.na(Ai))) stop(sprintf("ERROR: Unknown tip or node name '%s'",A[which(is.na(Ai))[1]]));
		A = Ai;
	}else if(check_input){
		minA = min(A); maxA = max(A);
		if(minA<1 || maxA>(Ntips+Nnodes)) stop(sprintf("ERROR: List A must contain values between 1 and Ntips+Nnodes (%d); instead, found values from %d to %d",Ntips+Nnodes,minA,maxA));
	}
	if(is.character(B)){
		Bi = name2clade[B]; 
		if(check_input && any(is.na(Bi))) stop(sprintf("ERROR: Unknown tip or node name '%s'",B[which(is.na(Bi))[1]]))
		B = Bi;
	}else if(check_input){
		minB = min(B); maxB = max(B);
		if(minB<1 || maxB>(Ntips+Nnodes)) stop(sprintf("ERROR: List B must contain values between 1 and Ntips+Nnodes (%d); instead, found values from %d to %d",Ntips+Nnodes,minB,maxB))
	}
	
	mrcas = get_most_recent_common_ancestors_CPP(Ntips			= Ntips,
												Nnodes			= tree$Nnode,
												Nedges			= nrow(tree$edge),
												tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												cladesA			= A-1,
												cladesB			= B-1,
												verbose			= FALSE,
												verbose_prefix	= "");
# 	mrcas = get_MRCAs_CPP(	Ntips			= Ntips,
# 							Nnodes			= tree$Nnode,
# 							Nedges			= nrow(tree$edge),
# 							tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
# 							cladesA			= A-1,
# 							cladesB			= B-1);
	return(as.integer(mrcas+1));
}