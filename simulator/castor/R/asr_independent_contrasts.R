# Phylogenetic Independent Contrasts (PIC) reconstruction of continuous ancestral states
# PIC ASR for a node only takes into account the subtree descending from that node (hence represents a "local" optimization).
# This corresponds to the "local" squared-change parsimony reconstructions by Maddison (1991), i.e. without rerooting at each node.
#    Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.
#    Maddison (1991). Squared-change parsimony reconstructions of ancestral states for continuous-valued characters on a phylogenetic tree. Systematic Zoology. 40:304-314.
#    Garland et al. (1999). An introduction to phylogenetically based statistical methods, with a new method for confidence intervals on ancestral values. American Zoologist. 39:374-388.
#  Requirements:
# 	Tree can be multifurcating, and can also include nodes with a single child
# 	Tree can also include edges with length zero (will be adjusted internally to some small epsilon if weighted==TRUE).
#	Tree must be rooted.
asr_independent_contrasts = function(	tree, 
										tip_states, 	# numeric vector of size Ntips
										weighted	= TRUE,
										include_CI	= FALSE,
										check_input	= TRUE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;

	# basic error checking
	if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
	if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be numeric"))
	if(check_input){
		if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
	}

	results = ASR_via_independent_contrasts_CPP(	Ntips					= Ntips,
													Nnodes					= tree$Nnode,
													Nedges					= nrow(tree$edge),
													tree_edge				= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													edge_length				= (if((!weighted) || (is.null(tree$edge.length))) numeric() else tree$edge.length),
													tip_states				= tip_states,
													include_standard_errors = include_CI);
	return(list(success				= TRUE,
				ancestral_states	= results$node_states,
				standard_errors		= (if(include_CI) results$node_standard_errors else NULL),
				CI95				= (if(include_CI) results$node_CI95s else NULL)));
}