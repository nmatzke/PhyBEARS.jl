# Given a phylogenetic tree in standard representation, calculate an alternative representation of the tree structure by enumerating clades and their properties
# Main output:
#   clades[]: 2D integer matrix of size Nclades x (Nsplits+1), with the following columns:
#     Column 1: Parent clade index
#     Columns 2-(Nsplits+1): Child clade indinces
#   lengths[]: 1D numeric array of size Nclades, listing incoming edge lengths for each clade (will be negative for the root)
#	Each row in clades and lengths[] corresponds to a unique clade (tip or node) in the tree
#   Rows in clades[] and lengths[] can be ordered according to clade indices, or according to post-order traversal.
# Negative values indicate missing values
# This function is loosely analogous to the function phybase::read.tree.nodes
# Requirements:
#    The tree can include mono- and multifurcations
#    The tree must be rooted if postorder==true
get_clade_list = function(tree, postorder = FALSE, missing_value=NA){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	results = tree_to_clade_list_CPP(	Ntips			= Ntips,
										Nnodes			= Nnodes,
										Nedges			= nrow(tree$edge),
										tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
										edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
										postorder		= postorder);
	if(!results$success) return(list(success=FALSE, error=results$error));
	clades				= results$clades
	clades[clades>=0]	= clades[clades>=0] + 1 # make indices 1-based
	clades[clades<0] 	= missing_value;
	clades  			= matrix(clades, ncol=(results$Nsplits+1), byrow=TRUE) # unflatten matrix

	if(is.null(tree$edge.length)){
		lengths = NULL
	}else{
		lengths 			= results$lengths
		lengths[lengths<0] 	= missing_value;
	}
	
	return(list(success			= TRUE,
				Nsplits			= results$Nsplits,
				clades			= clades,
				lengths 		= lengths,
				old2new_clade 	= results$old2new_clade+1))
}