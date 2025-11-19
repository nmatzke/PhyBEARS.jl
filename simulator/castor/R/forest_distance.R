# Calculate distance (single number) between two sets of exchangeable unlabeled trees.
# Hence, the distance calculated does not depend in the ordering/labeling of tips in each tree, nor on the order of trees in each tree set.
# The trees may include multifurcations and monofurcations.
# If normalized==TRUE, the distance will always be within 0 and 1. However, it is no longer guaranteed that the normalized distance function satisfies the triangle inequality, as required for actual metrics.
forest_distance = function(	treesA, 		# list of phylo trees
							treesB, 		# list of phylo trees
							metric			= "WassersteinNodeAges",	# which distance function to use
							combine			= "mean_pairwise",			# how the pairwise distances between trees should be combined to form a single distance between the forests.
							normalized		= FALSE,
							NLeigenvalues	= 10){	# number of top eigenvalues to consider of the Laplacian Spectrum (e.g., for the metric "WassersteinLaplacianSpectrum"). If <=0, all eigenvalues are considered, which can substantially increase computation time for large trees.
	NtreesA = length(treesA)
	NtreesB	= length(treesB)
	
	distances = matrix(0, ncol=NtreesA, nrow=NtreesB)
	intermediatesA = vector(mode="list", NtreesA) # temporary storage for any intermediate results for treesA (will be filled as we go)
	intermediatesB = vector(mode="list", NtreesB) # temporary storage for any intermediate results for treesB (will be filled as we go)
	for(a in seq_len(NtreesA)){
		treeA = treesA[[a]]
		NcladesA = length(treeA$tip.label) + treeA$Nnode
		for(b in seq_len(NtreesB)){
			treeB = treesB[[b]]
			NcladesB = length(treeB$tip.label) + treeB$Nnode
			if(metric=="WassersteinNodeAges"){
				# First Wasserstein ("1-Wasserstein") distance between the empirical distributions of node ages
				# This distance depends on branch lengths, but not on tip labeling nor topology (only on the ages of the nodes)
				# Hence, this is strictly speaking only a pseudometric, even in the space of unlabeled trees. It is a metric in the space of tree equivalence classes, where two trees are equivalent iff they have the same node ages.
				if(is.null(intermediatesA[[a]])){
					root_ageA 	= get_tree_span(treeA)$max_distance
					clade_agesA = root_ageA - get_all_distances_to_root(treeA)
					intermediatesA[[a]] = list(root_age=root_ageA, clade_ages=clade_agesA)
				}else{
					root_ageA   = intermediatesA[[a]]$root_age
					clade_agesA = intermediatesA[[a]]$clade_ages
				}
				if(is.null(intermediatesB[[b]])){
					root_ageB 	= get_tree_span(treeB)$max_distance
					clade_agesB = root_ageB - get_all_distances_to_root(treeB)
					intermediatesB[[b]] = list(root_age=root_ageB, clade_ages=clade_agesB)
				}else{
					root_ageB   = intermediatesB[[b]]$root_age
					clade_agesB = intermediatesB[[b]]$clade_ages
				}
				if(normalized){
					# rescale time so that it is relative to the oldest root_age, and thus all clade ages are numbers between 0 and 1
					# this also means that the 1st-Wasserstein distance will be between 0 and 1
					max_root_age 	= max(root_ageA,root_ageB)
					clade_agesA 	= clade_agesA/max_root_age
					clade_agesB		= clade_agesB/max_root_age
				}
				distances[a,b] = first_Wasserstein_distance(clade_agesA[(length(treeA$tip.label)+1):length(clade_agesA)], clade_agesB[(length(treeB$tip.label)+1):length(clade_agesB)])
			}else if(metric=="WassersteinLaplacianSpectrum"){
				# First Wasserstein ("1-Wasserstein") distance between the eigenspectra of the modified graph Laplacians
				# This distance depends on topology and branch lengths, but not on tip labeling nor on the rooting
				# Hence, this is strictly speaking a metric in the space of unrooted unlabeled trees, but not a metric in the space of labeled trees.
				# This metric is similar to that proposed by Lewitus and Morlon (2016, Systematic Biology. 65:495-507), with the difference that the latter calculates a smoothened version of the spectrum (via a Gaussian convolution kernel) and then calculates the Kullback-Leibler divergence between the smoothened densities.
				# Note that if not all eigenvalues are used (i.e. NLeigenvalues>0 and NLeigenvalues<max(NcladesA,NcladesB)), then this is not even a metric on the space of unlabeled unrooted trees. 
				if(is.null(intermediatesA[[a]])){
					LaplacianA = weighted_graph_Laplacian_of_tree(treeA, sparse=((NLeigenvalues>0) && (NLeigenvalues<NcladesA)))
					if((NLeigenvalues>0) && (NLeigenvalues<NcladesA)){
						spectrumA = RSpectra::eigs(A=LaplacianA, k=NLeigenvalues, which="LM", opts=list(retvec=FALSE))$values
					}else{
						spectrumA = eigen(x=LaplacianA, symmetric=TRUE, only.values=TRUE)$values
					}
					intermediatesA[[a]] = list(spectrum=spectrumA)
				}else{
					spectrumA = intermediatesA[[a]]$spectrum
				}
				if(is.null(intermediatesB[[b]])){
					LaplacianB = weighted_graph_Laplacian_of_tree(treeB, sparse=((NLeigenvalues>0) && (NLeigenvalues<NcladesB)))
					if((NLeigenvalues>0) && (NLeigenvalues<NcladesB)){
						spectrumB = RSpectra::eigs(A=LaplacianB, k=NLeigenvalues, which="LM", opts=list(retvec=FALSE))$values
					}else{
						spectrumB = eigen(x=LaplacianB, symmetric=TRUE, only.values=TRUE)$values
					}
					intermediatesB[[b]] = list(spectrum=spectrumB)
				}else{
					spectrumB = intermediatesB[[b]]$spectrum
				}
				if(normalized){
					# rescale eigenvalues so that they are relative to the largest eigenvalue of both trees; hence, all rescaled eigenvalues will be within [0,1]
					# this also means that the 1st-Wasserstein distance will be between 0 and 1
					largest_eigenvalue = max(spectrumA, spectrumB)
					spectrumA = spectrumA / largest_eigenvalue
					spectrumB = spectrumB / largest_eigenvalue
				}
				# calculate first Wasserstein distance between the two spectra
				distances[a,b] = first_Wasserstein_distance(spectrumA, spectrumB)
			}else{
				stop(sprintf("Unknown metric '%s'",metric))
			}
		}
	}
	
	# combine pairwise distances
	if(combine=="mean_pairwise"){
		distance = mean(distances)
	}else if(combine=="max_pairwise"){
		distance = max(distances)
	}else{
		return(list(succes=FALSE, error=sprintf("forest_distance: Unknown option '%s' for combine",combine)))
	}
	return(list(success		= TRUE,
				distance 	= distance,
				distances	= distances))
}
