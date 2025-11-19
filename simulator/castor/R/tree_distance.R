# Calculate distance (single number) between two trees with identical tips.
# The trees may include multifurcations and monofurcations.
# If normalized==TRUE, the distance will always be within 0 and 1.
# tipsA2B[] should be a 1D array mapping A-tip indices to B-tip indices (one-to-one mapping). If NULL, it is determined by matching tip labels.
tree_distance = function(	treeA, 
							treeB, 
							tipsA2B			= NULL,
							metric			= "RFrooted",	# which distance metric to use
							normalized		= FALSE,
							NLeigenvalues	= 10){	# number of top eigenvalues to consider of the Laplacian Spectrum (e.g., for the metric "WassersteinLaplacianSpectrum"). If <=0, all eigenvalues are considered, which can substantially increase computation time for large trees.
	# basic error checking
	needs_same_tip_count = (metric %in% c("RFrooted", "MeanPathLengthDifference"))
	needs_tip_matching	 = (metric %in% c("RFrooted", "MeanPathLengthDifference"))
	NtipsA  = length(treeA$tip.label)
	NtipsB  = length(treeB$tip.label)
	if(needs_same_tip_count && (NtipsA!=NtipsB)) stop(sprintf("Tip counts don't match, but needed for metric '%s': TreeA has %d tips, treeB has %d tips",metric,NtipsA,NtipsB))
	if(needs_tip_matching){
		if(is.null(tipsA2B)){
			tipsA2B = match(treeA$tip.label, treeB$tip.label)
			if(any(is.na(tipsA2B))) stop(sprintf("Some tip labels in treeA don't match any tip labels in treeB, but needed for metric '%s'",metric))
		}else{
			if(NtipsA!=length(tipsA2B)) stop(sprintf("The length of tipsA2B (%d) differs from the number of tips in treeA (%d)",length(tipsA2B),NtipsA))
		}
	}
	NcladesA  = NtipsA + treeA$Nnode
	NcladesB  = NtipsB + treeB$Nnode
	
	if(metric=="RFrooted"){
		# Robinson-Foulds distance between rooted trees (i.e., outcome depends on rooting)
		# Only considers tree topology, but not branch lengths
		results = get_Robinson_Foulds_distance_CPP(	Ntips		= NtipsA,
													NnodesA		= treeA$Nnode,
													NedgesA		= nrow(treeA$edge),
													tree_edgeA	= as.vector(t(treeA$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													NnodesB		= treeB$Nnode,
													NedgesB		= nrow(treeB$edge),
													tree_edgeB	= as.vector(t(treeB$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													tipsA2B		= tipsA2B-1);
		distance = treeA$Nnode+treeB$Nnode - 2*results$Nmatches
		if(normalized) distance = distance/(treeA$Nnode+treeB$Nnode)
	}else if(metric=="MeanPathLengthDifference"){
		# path-difference metric between unrooted trees (outcome does not depend on rooting)
		# [Steel and Penny, 1993, Systematic Biology. 42:126-141.]
		tip_distancesA = as.vector(get_all_pairwise_distances(treeA, only_clades = seq_len(NtipsA), as_edge_counts = FALSE, check_input = FALSE))
		tip_distancesB = as.vector(get_all_pairwise_distances(treeB, only_clades = tipsA2B, as_edge_counts	= FALSE, check_input = FALSE))
		if(normalized){
			# square root of mean squared difference of normalized patristic distances between all tip pairs
			# normalization of distances is done by dividing by the maximum of the two distances being compared
			both_zero = ((tip_distancesA==0) & (tip_distancesB==0))
			tip_distancesA[both_zero] = 1 # avoid comparison when distance is zero in both trees, due to division by zero
			tip_distancesB[both_zero] = 1 # avoid comparison when distance is zero in both trees, due to division by zero
			distance = sqrt(sum(((tip_distancesA-tip_distancesB)/pmax(tip_distancesA,tip_distancesB))**2)/length(tip_distancesA))
		}else{
			# square root of mean squared difference of patristic distances between all tip pairs
			distance = sqrt(sum((tip_distancesA-tip_distancesB)**2)/length(tip_distancesA))
		}
	}else if(metric=="WassersteinNodeAges"){
		# First Wasserstein ("1-Wasserstein") distance between the empirical distributions of node ages
		# This distance depends on branch lengths and the root placement, but not on tip labeling nor topology (only on the ages of the nodes)
		# Hence, this is strictly speaking only a pseudometric, even in the space of unlabeled trees. It is a metric in the space of tree equivalence classes, where two trees are equivalent iff they have the same node ages.
		root_ageA 	= get_tree_span(treeA)$max_distance
		root_ageB 	= get_tree_span(treeB)$max_distance
		clade_agesA = root_ageA - get_all_distances_to_root(treeA)
		clade_agesB = root_ageB - get_all_distances_to_root(treeB)
		if(normalized){
			# rescale time so that it is relative to the oldest root_age, and thus all clade ages are numbers between 0 and 1
			# this also means that the 1st-Wasserstein distance will be between 0 and 1
			max_root_age 	= max(root_ageA,root_ageB)
			clade_agesA 	= clade_agesA/max_root_age
			clade_agesB		= clade_agesB/max_root_age
		}
		distance = first_Wasserstein_distance(clade_agesA[(NtipsA+1):length(clade_agesA)], clade_agesB[(NtipsB+1):length(clade_agesB)])
	}else if(metric=="WassersteinLaplacianSpectrum"){
		# First Wasserstein ("1-Wasserstein") distance between the eigenspectra of the modified graph Laplacians
		# This distance depends on topology and branch lengths, but not on tip labeling nor on the rooting
		# Hence, this is strictly speaking a metric in the space of unrooted unlabeled trees, but not a metric in the space of labeled trees.
		# This distance is similar to that proposed by Lewitus and Morlon (2016, Systematic Biology. 65:495-507), with the difference that the latter calculates a smoothened version of the spectrum (via a Gaussian convolution kernel) and then calculates the Kullback-Leibler divergence between the smoothened densities.
		# Note that if not all eigenvalues are used (i.e. NLeigenvalues>0 and NLeigenvalues<max(NcladesA,NcladesB)), then this is not even a metric on the space of unlabeled unrooted trees. 
		LaplacianA 	= weighted_graph_Laplacian_of_tree(treeA, sparse=(if((NLeigenvalues>0) && (NLeigenvalues<NcladesA)) TRUE else FALSE))
		LaplacianB 	= weighted_graph_Laplacian_of_tree(treeB, sparse=(if((NLeigenvalues>0) && (NLeigenvalues<NcladesB)) TRUE else FALSE))
		# calculate spectra of the Laplacians. Note that a tree's Laplacian is a symmetric matrix.
		if((NLeigenvalues>0) && (NLeigenvalues<NcladesA)){
			spectrumA = RSpectra::eigs(A=LaplacianA, k=NLeigenvalues, which="LM", opts=list(retvec=FALSE))$values
		}else{
			spectrumA = eigen(x=LaplacianA, symmetric=TRUE, only.values=TRUE)$values
		}
		if((NLeigenvalues>0) && (NLeigenvalues<NcladesB)){
			spectrumB = RSpectra::eigs(A=LaplacianB, k=NLeigenvalues, which="LM", opts=list(retvec=FALSE))$values
		}else{
			spectrumB = eigen(x=LaplacianB, symmetric=TRUE, only.values=TRUE)$values
		}
		if(normalized){
			# rescale eigenvalues so that they are relative to the largest eigenvalue of both trees; hence, all rescaled eigenvalues will be within [0,1]
			# this also means that the 1st-Wasserstein distance will be between 0 and 1
			largest_eigenvalue = max(spectrumA, spectrumB)
			spectrumA = spectrumA / largest_eigenvalue
			spectrumB = spectrumB / largest_eigenvalue
		}
		# calculate first Wasserstein distance between the two spectra
		distance = first_Wasserstein_distance(spectrumA, spectrumB)
	}else{
		stop(sprintf("Unknown metric '%s'",metric))
	}
	return(distance)
}
