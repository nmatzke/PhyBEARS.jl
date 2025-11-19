# Given a reference tree ("R") and a target tree ("T"), map T-nodes to concordant R-nodes when possible, and extract divergence times from the R-tree
# This only makes sense if the R-tree is time-calibrated
# A provided mapping specifies which tips in target correspond to which tips in the reference tree
# This function can be used to define 2-ndary dating constraints for a larger target tree, based on a time-calibrated smaller reference tree
#
# Reference: [Eastman et al (2013). Congruification: support for time scaling large phylogenetic trees. Methods in Ecology and Evolution. 4:688-691]
#
# mapping must be one of the following:
#	A 2D integer array of size NM x 2 (with NM<=TNtips), listing Ttips mapped to Rtips (mapping[m,1] --> mapping[m,2])
#	A 2D character array of size NM x 2 (with NM<=TNtips), listing Ttip names mapped to Rtip names (mapping[m,1] --> mapping[m,2])
#	A data frame of size NM x 1, whose row names are target tip labels and whose entries are either integers (R tip indices) or strings (R tip labels). This is the format used by geiger::congruify.phylo
#	A vector of size NM, whose names are target tip labels and whose entries are either integers (R tip indices) or strings (R tip labels).
congruent_divergence_times = function(reference_tree, target_tree, mapping){
	RNtips = length(reference_tree$tip.label)
	congruification = congruify_trees(reference_tree, target_tree, mapping)
	Rages  = get_tree_span(reference_tree)$max_distance - castor::get_all_distances_to_root(reference_tree)
	return(list(reference_nodes = congruification$reference_nodes,
				target_nodes 	= congruification$target_nodes,
				ages			= Rages[congruification$reference_nodes+RNtips]))
}
