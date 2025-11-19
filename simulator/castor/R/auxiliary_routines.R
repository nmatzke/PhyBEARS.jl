# convert a list to a vector using unlist, after converting NULLs to NAs
# this function may be needed if your list contains NULLs, in which case the standard unlist removes those entries from the returned vector
unlist_with_nulls = function(L){
	L[sapply(L, is.null)] = NA
	return(unlist(L))
}


get_adjacent_edges_per_edge = function(tree){
	Nedges = nrow(tree$edge)
	adjacents = get_adjacent_edges_per_edge_CPP(	Ntips 		= length(tree$tip.label),
													Nnodes		= tree$Nnode,
													Nedges		= Nedges,
													tree_edge	= as.vector(t(tree$edge)) - 1);
	
	# update indices from 0-based to 1-based
	return(lapply(1:Nedges,FUN = function(edge){ adjacents[[edge]]+1 }))
}



get_outgoing_edges_per_clade = function(tree){
	Nclades = length(tree$tip.label) + tree$Nnode;
	outgoing_edges = get_outgoing_edges_per_clade_CPP(	Ntips 		= length(tree$tip.label),
														Nnodes		= tree$Nnode,
														Nedges		= nrow(tree$edge),
														tree_edge	= as.vector(t(tree$edge)) - 1);
	
	# update indices from 0-based to 1-based
	return(lapply(1:Nclades,FUN = function(clade){ outgoing_edges[[clade]]+1 }))
}


get_incoming_edges_per_clade = function(tree){
	Nclades = length(tree$tip.label) + tree$Nnode;
	incoming_edges = get_incoming_edges_per_clade_CPP(	Ntips 		= length(tree$tip.label),
														Nnodes		= tree$Nnode,
														Nedges		= nrow(tree$edge),
														tree_edge	= as.vector(t(tree$edge)) - 1);
	
	# update indices from 0-based to 1-based
	return(lapply(1:Nclades,FUN = function(clade){ incoming_edges[[clade]]+1 }))
}



get_paths_root_to_tips = function(tree){
	Ntips = length(tree$tip.label)
	paths = get_paths_root_to_tips_CPP(	Ntips 		= Ntips,
										Nnodes		= tree$Nnode,
										Nedges		= nrow(tree$edge),
										tree_edge	= as.vector(t(tree$edge)) - 1);
	
	# update indices from 0-based to 1-based
	return(lapply(1:Ntips,FUN = function(tip){ paths[[tip]]+1 }))
}



# convert a name list to a list of tip or node indices, or retain the original index list if already in integer form
# type can be 'tip', 'node' or 'both'
map_tip_or_node_names_to_indices = function(tree, A, type, list_title, check_input=TRUE){
	Ntips = length(tree$tip.label)
	if(type=='tip'){
		item_title 	= 'tip'
		Nmax_title 	= 'Ntips'
		Nmax 		= Ntips
		if(is.character(A)) name_pool = tree$tip.label;
	}else if(type=='node'){
		item_title 	= 'node';
		Nmax_title 	= 'Nnodes'
		Nmax 		= tree$Nnode;
		if(is.character(A)){
			if(is.null(tree$node.label)) stop(sprintf("ERROR: Tree must include node labels, for mapping %s to node indices",list_title))
			name_pool = tree$node.label
		}
	}else{
		item_title 	= 'tip or node'
		Nmax_title 	= 'Ntips+Nnodes'
		Nmax 		= Ntips+tree$Nnode;
		if(is.character(A)) name_pool = c(tree$tip.label,tree$node.label);
	}
	if((!is.character(A)) && (!is.numeric(A)) && (!is.integer(A))) stop(sprintf("ERROR: %s must be a character or integer vector",list_title))
	if(is.character(A)){
		name2index = 1:Nmax;
		names(name2index) = name_pool;
		Ai = name2index[A]; 
		if(check_input && any(is.na(Ai))) stop(sprintf("ERROR: Unknown %s name '%s'",item_title,A[which(is.na(Ai))[1]]));
		A = Ai;
	}else if(check_input && (length(A)>0)){
		minA = min(A); maxA = max(A);
		if(minA<1 || maxA>Nmax) stop(sprintf("ERROR: %s must contain values between 1 and %s (%d); instead found values from %d to %d",list_title,Nmax_title,Nmax,minA,maxA));
	}
	return(A);
}




# given a Markov transition rate matrix, calculate the transition probability matrix conditional upon a single transition occurring
# input: Q[i,j] is probability rate of transition i-->j
# output: P[i,j] will be the probability of transition i-->j, provided that a single transition i-->* occurred
get_conditional_transition_probabilities = function(Q){
	S = rowSums(Q)-diag(Q)
	S[S<=0] = 1
	P = Q/S
	diag(P) = 0
	return(P)
}



find_edge_splitting_tree = function(tree, target_tips, is_rooted=FALSE){
	Ntips 	= length(tree$tip.label)
	Nnodes 	= tree$Nnode
	if(is.character(target_tips)){
		# target tips are given as names, not indices
		indices	= match(target_tips, tree$tip.label)
		if(any(is.na(indices))) stop(sprintf("ERROR: Some target_tips (e.g. '%s') were not found in the tree",target_tips[which(is.na(indices))[1]]))
		target_tips = indices
	}
	
	results = find_edge_splitting_tree_CPP(	Ntips				= Ntips,
											Nnodes				= tree$Nnode,
											Nedges				= nrow(tree$edge),
											tree_edge			= as.vector(t(tree$edge)) - 1,	# flatten in row-major format and adjust clade indices to 0-based
											is_rooted			= is_rooted,
											target_tips			= target_tips - 1,
											include_misplaced 	= TRUE)
													
	return(list(edge					= (if(results$edge<0) NA else as.integer(results$edge+1)),
				Nmisplaced_targets		= results$Nmisplaced_targets,
				Nmisplaced_nontargets	= results$Nmisplaced_nontargets,
				Ntargets_upstream 		= results$Ntargets_upstream,
				Ntargets_downstream 	= results$Ntargets_downstream,
				misplaced_targets		= results$misplaced_targets,
				misplaced_nontargets	= results$misplaced_nontargets));
				
}




get_subtree_with_clades = function(	tree, 
									clades_to_keep = NULL, 	# integer vector listing tip/node indices to keep
									collapse_monofurcations = TRUE, 
									force_keep_root = FALSE, 
									keep_all_children_of_explicit_clades_to_keep = FALSE, 
									keep_all_tips_of_explicit_clades_to_keep = FALSE){ 
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	
	results = get_subtree_with_specific_clades_CPP(	Ntips 					= Ntips,
													Nnodes 					= Nnodes,
													Nedges					= nrow(tree$edge),
													tree_edge				= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													edge_length				= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													clades_to_keep			= clades_to_keep-1,
													collapse_monofurcations = collapse_monofurcations,
													force_keep_root			= force_keep_root,
													keep_all_children_of_explicit_clades_to_keep 	= keep_all_children_of_explicit_clades_to_keep,
													keep_all_tips_of_explicit_clades_to_keep 		= keep_all_tips_of_explicit_clades_to_keep)
	Ntips_kept  	= results$Ntips_kept
	Nnodes_kept 	= results$Nnodes_kept
	new2old_clade 	= results$new2old_clade + 1L # switch to 1-based indices
	subtree = list(	Nnode 		= Nnodes_kept,
					tip.label 	= (if(is.null(tree$tip.label)) NULL else tree$tip.label[new2old_clade[1:Ntips_kept]]),
					node.label 	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips]),
					edge 		= matrix(as.integer(results$new_tree_edge),ncol=2,byrow=TRUE) + 1L,
					edge.length = results$new_edge_length,
					root 		= results$new_root+1L,
					root.edge	= (if(force_keep_root && (!is.null(tree$root.edge))) tree$root.edge else NULL));
	class(subtree) = "phylo";
	attr(subtree,"order") = NULL
	
	return(list(tree 			= subtree,
				root_shift		= results$root_shift, # distance between old & new root (will always be non-negative)
				new2old_tip		= new2old_clade[1:Ntips_kept], 
				new2old_node	= new2old_clade[(Ntips_kept+1):(Ntips_kept+Nnodes_kept)]-Ntips));
}



# calculate the geometric placement of tips & nodes for plotting a tree as a phylogram (with the root on the left and tips on the right end, edges extend horizontally left to right)
get_phylogram_geometry = function(tree){
	results = get_phylogram_geometry_CPP(	Ntips 					= length(tree$tip.label),
											Nnodes 					= tree$Nnode,
											Nedges					= nrow(tree$edge),
											tree_edge				= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
											edge_length				= (if(is.null(tree$edge.length)) numeric() else tree$edge.length));
	return(list(clade_x	= results$clade_x,		# x-coordinates of tips & nodes
				clade_y	= results$clade_y,		# y-coordinates of tips & nodes
				root_y	= results$root_y,		# the y-coordinate of the root
				min_x	= results$min_x,		# the minimum x-coordinate of any clade (normally this is 0)
				max_x	= results$max_x,		# the maximum x-coordinate of any tip
				min_y	= results$min_y,		# the minimum y-coordinate of any tip
				max_y	= results$max_y));		# the maximum y-coordinate of any tip
}




# assign tips & nodes of a tree to groups, such that each group is monophyletic (a "taxon") represented by exactly one of given representative tips
# each representative tip is taken to represent a different taxon
# tip2taxon[n] or node2taxon[n] will be -1 if the tip/node could not be unambiguously assigned to a taxon (e.g., it contains multiple descending representatives)
assign_clades_to_taxa = function(tree, representative_tips){
	Ntips  = length(tree$tip.label);
	Nnodes = tree$Nnode;
	results = assign_clades_to_taxa_CPP(	Ntips 			= Ntips,
											Nnodes 			= tree$Nnode,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
											representatives	= representative_tips-1);
	return(list(tip2taxon 	= results$clade2taxon[1:Ntips]+1L,
				node2taxon	= results$clade2taxon[(Ntips+1):(Ntips+Nnodes)]+1L));
}




# congruify trees (map nodes in the target tree to "equivalent" nodes in the reference tree)
# [Eastman et al (2013). Congruification: support for time scaling large phylogenetic trees. Methods in Ecology and Evolution. 4:688-691]
# mapping must be one of the following:
#	A 2D integer array of size NM x 2 (with NM<=TNtips), listing Ttips mapped to Rtips (mapping[m,1] --> mapping[m,2])
#	A 2D character array of size NM x 2 (with NM<=TNtips), listing Ttip names mapped to Rtip names (mapping[m,1] --> mapping[m,2])
#	A data frame of size NM x 1, whose row names are target tip labels and whose entries are either integers (R tip indices) or strings (R tip labels). This is the format used by geiger::congruify.phylo
#	A vector of size NM, whose names are target tip labels and whose entries are either integers (R tip indices) or strings (R tip labels).
congruify_trees = function(reference_tree, target_tree, mapping){
	TNtips = length(target_tree$tip.label)
	RNtips = length(reference_tree$tip.label)

	# re-format mapping if needed
	if(is.data.frame(mapping)){
		mapped_Ttips = rownames(mapping)
		mapped_Rtips = (if(is.numeric(mapping[,1])) mapping[,1] else as.character(mapping[,1]))
	}else if(is.vector(mapping)){
		mapped_Ttips = names(mapping)
		mapped_Rtips = (if(is.numeric(mapping)) mapping else as.character(mapping))
	}else{
		# assuming mapping is a 2D array of size NM x 2
		mapped_Ttips = mapping[,1]
		mapped_Rtips = mapping[,2]
	}
	if(is.character(mapped_Ttips)){
		# mapped target tips given as names, not indices
		indices = match(mapped_Ttips, target_tree$tip.label)
		if(any(is.na(indices))) stop(sprintf("ERROR: Some mapped target tips (e.g. '%s') were not found in the target tree",mapped_Ttips[which(is.na(indices))[1]]))
		mapped_Ttips = indices
	}
	if(is.character(mapped_Rtips)){
		# mapped reference tips given as names, not integers
		indices = match(mapped_Rtips, reference_tree$tip.label)
		if(any(is.na(indices))) stop(sprintf("ERROR: Some mapped reference tips (e.g. '%s') were not found in the reference tree",mapped_Rtips[which(is.na(indices))[1]]))
		mapped_Rtips = indices
	}
	mapping = matrix(c(mapped_Ttips,mapped_Rtips),ncol=2,byrow=FALSE)
		
	# congruify
	results = congruify_trees_CPP(	RNtips		= RNtips,
									RNnodes		= reference_tree$Nnode,
									RNedges		= nrow(reference_tree$edge),
									Rtree_edge	= as.vector(t(reference_tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
									TNtips		= TNtips,
									TNnodes		= target_tree$Nnode,
									TNedges		= nrow(target_tree$edge),
									Ttree_edge	= as.vector(t(target_tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
									mapping		= as.vector(t(mapping))-1) # flatten in row-major format and adjust tip indices to 0-based
	return(list(reference_nodes = results$mapped_Rnodes+1L,
				target_nodes 	= results$mapped_Tnodes+1L))
}



# map nodes in tree A to nodes in tree B, assuming both trees have the same topologies (but are potentially indexed differently)
# if tipsA2B is NULL, tips are matched by name
# This function returns an error (success=FALSE) if the trees don't have equivalent topologies, so it can also be used as a simple equivalence test
match_tree_nodes = function(treeA, treeB, tipsA2B=NULL){
	Ntips  = length(treeA$tip.label)
	Nnodes = treeA$Nnode
	if((Ntips!=length(treeB$tip.label)) || (Nnodes!=treeB$Nnode)) return(list(success=FALSE, error=sprintf("Tree sizes don't match: TreeA has %d tips and %d nodes, treeB has %d tips and %d nodes",Ntips,Nnodes,length(treeB$tip.label),treeB$Nnode)))
	if(is.null(tipsA2B)){
		tipsA2B = match(treeA$tip.label, treeB$tip.label)
		if(any(is.na(tipsA2B))) return(list(success=FALSE, error=sprintf("Tip labels in treeA don't match tip labels in treeB")))
	}
	results = match_tree_nodes_CPP(	Ntips		= Ntips,
									Nnodes		= Nnodes,
									Nedges		= nrow(treeA$edge),
									tree_edgeA	= as.vector(t(treeA$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
									tree_edgeB	= as.vector(t(treeB$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
									tipsA2B		= tipsA2B-1);
	if(!results$success) return(list(success=FALSE, error=results$error));
	return(list(success 	= TRUE,
				rootA		= results$rootA,
				rootB		= results$rootB,
				nodesA2B	= results$nodesA2B+1L));
}


get_complement = function(N, indices){
	include 		 = rep(TRUE, times=N)
	include[indices] = FALSE
	return(which(include))
}

get_intersection = function(N, indicesA, indicesB){
	includeA 		 	= rep(FALSE, times=N)
	includeB 		 	= rep(FALSE, times=N)
	includeA[indicesA] 	= TRUE
	includeB[indicesB] 	= TRUE
	return(which(includeA & includeB))
}

# extract the values of independent rates from a transition_matrix, based on a provided index_matrix
# The index matrix should be as generated by get_transition_index_matrix, i.e. a 2D matrix of size Nstates x Nstates, and values in 1,..,Nrates, where Nrates is the number of independent rate variables
extract_independent_rates_from_transition_matrix = function(transition_matrix, index_matrix){
	flattened_index_matrix 	= as.vector(index_matrix)
	flattened_transition_matrix = as.vector(transition_matrix)
	independents 			= seq_along(flattened_index_matrix)[!duplicated(flattened_index_matrix)]
	independents 			= independents[flattened_index_matrix[independents]>0] # omit any zeros from index_matrix
	independent_rates 		= rep(NA,length(independents))
	independent_rates[flattened_index_matrix[independents]] = flattened_transition_matrix[independents]
	return(independent_rates)
}



# get a 1D lookup matrix of size Nstates, mapping birth-rates to indices of a subset of independent birth-rates
# model can be "ER" or "ARD" or a custom index_matrix as if it was generated by this function (in which case it is merely used to determine Nrates)
# This function is analogous to get_transition_index_matrix(..), but for 1D vectors
get_rate_index_vector = function(Nstates, rate_model){
	if (is.character(rate_model)) {
		if(rate_model == "ER"){
			Nrates = 1;
			index_vector = rep(1,Nstates)
		}else if(rate_model == "ARD"){
			Nrates = Nstates;
			index_vector = 1:Nrates;		
		}else{
			stop(sprintf("ERROR: Unknown rate_model '%s'",rate_model))
		}
	}else{
		if(length(rate_model)!=Nstates) stop(sprintf("ERROR: Wrong number of elements in rate model (expected %d, found %d)",Nstates,length(rate_model)));
		index_vector = rate_model
		Nrates = max(rate_model)
	}
	return(list(index_vector=index_vector, Nrates=Nrates))
}


extract_independent_rates_from_vector = function(rates, index_vector){
	independents 			= seq_along(index_vector)[!duplicated(index_vector)]
	independent_rates 		= rep(NA,length(independents))
	independent_rates[index_vector[independents]] = rates[independents]
	return(independent_rates)
}




# guesstimate an Mk transition matrix Q based on transitions along edges, as inferred via max-parsimony ASR
# Convention: Q[i,j] will be an estimate for the probability rate of the transition i-->j
# at least one of tip_states[] or tip_priors[] must be given; tip_states[] is given priority
guesstimate_Mk_transition_rates_via_max_parsimony_ASR = function(	tree, 
																	tip_states			= NULL,	# 1D vector of size Ntips, or NULL
																	tip_priors			= NULL,	# 2D array of size Ntips x Nstates, or NULL 
																	Nstates			 	= NULL, 
																	transition_costs 	= "all_equal",
																	allow_ties			= FALSE){	# logical, specifying whether to include tips with ties in the tip priors, i.e. tips that have the same tip_prior for various states. If TRUE, then ties are broken at random. If FALSE, then tips with ties are omitted from the analysis.
	# basic error checking & input formatting
	if(is.null(tip_states)){
		if(is.null(tip_priors)) return(list(success=FALSE, error="Missing tip_states or tip_priors"));
		if(allow_ties){
			tip_states = max.col(tip_priors, ties.method="random")
		}else{
			tip_states  = max.col(tip_priors, ties.method="first")
			tip_states2 = max.col(tip_priors, ties.method="last")
			tip_states[tip_states!=tip_states2] = NA # in case of ambiguity, assign NA
		}
	}
	if(length(tip_states)!=length(tree$tip.label)) return(list(success=FALSE, error=sprintf("Number of provided tip states (%d) does not match number of tips in the tree (%d)",length(tip_states),length(tree$tip.label))))
	
	# only consider subtree with known tip states
	known_tips = which(!is.na(tip_states));
	if(length(known_tips)<=1) return(list(success=FALSE, error=sprintf("All or almost all tips have unknown or ambiguous state")))
	if(length(known_tips)<length(tip_states)){
		extraction	= get_subtree_with_tips(tree, only_tips=known_tips, omit_tips=FALSE, collapse_monofurcations=TRUE, force_keep_root=TRUE);
		tree		= extraction$subtree;
		tip_states	= tip_states[extraction$new2old_tip]
	}
	Ntips  = length(tree$tip.label)
	Nedges = nrow(tree$edge)
	
	# perform ASR max-parsimony on known subtree
	asr = asr_max_parsimony(	tree				= tree, 
								tip_states			= tip_states, 		
								Nstates				= Nstates, 
								transition_costs	= transition_costs);
	if(!asr$success) return(list(success=FALSE, error="ASR max-parsimony failed"))
	Nstates = ncol(asr$ancestral_likelihoods);
								
	# determine most likely node states
	node_states 	 	= max.col(asr$ancestral_likelihoods, ties.method="first")
	clade_states 		= c(tip_states, node_states)
	state_transitions 	= cbind(clade_states[tree$edge[,1]],clade_states[tree$edge[,2]]) # state_transitions[e,1]-->state_transitions[e,2] is the transition along edge e.
	
	# determine entries Q[i,j] based on transitions i-->j among all edges
	Q 			 = matrix(0,nrow=Nstates,ncol=Nstates)
	Ntransitions = matrix(0,nrow=Nstates,ncol=Nstates)
	edge_lengths = (if(is.null(tree$edge.length)) rep(1.0,Nedges) else tree$edge.length)
	for(i in 1:Nstates){
		for(j in 1:Nstates){
			transitions = which((state_transitions[,1]==i) & (state_transitions[,2]==j));
			Ntransitions[i,j] = length(transitions)
			if(i!=j){ # keep diagonal zero, for correct normalization afterwards
				if(length(transitions)>0){
					mean_transition_time = mean(edge_lengths[transitions])
					Q[i,j] = 1.0/mean_transition_time
				}else{
					Q[i,j] = 0;
				}
			}
		}
	}
	
	# make sure Q has zero sum in each row
	diag(Q) = -rowSums(Q, na.rm = TRUE);
	return(list(success=TRUE, Q=Q, Ntransitions=Ntransitions));
}



# guesstimate the scaled of Mk transition rates based on transitions along independent contrasts
# Attention: Strictly speaking, this function can only estimate the single rate of an "ER" Mk model
# at least one of tip_states[] or tip_priors[] must be given; tip_states[] is given priority
guesstimate_Mk_rate_via_independent_contrasts = function(	tree, 					# rooted phylogenetic tree of class "phylo"
															Nstates,
															tip_states	= NULL,		# 1D vector of size Ntips, or NULL
															tip_priors	= NULL,		# 2D array of size Ntips x Nstates, or NULL 
															allow_ties	= FALSE){	# logical, specifying whether to include tips with ties in the tip priors, i.e. tips that have the same tip_prior for various states. If TRUE, then ties are broken at random. If FALSE, then tips with ties are omitted from the analysis.
	# basic error checking & input formatting
	if(is.null(tip_states)){
		if(is.null(tip_priors)) return(list(success=FALSE, error="Missing tip_states or tip_priors"));
		if(allow_ties){
			tip_states = max.col(tip_priors, ties.method="random")
		}else{
			tip_states  = max.col(tip_priors, ties.method="first")
			tip_states2 = max.col(tip_priors, ties.method="last")
			tip_states[tip_states!=tip_states2] = NA # in case of ambiguity, assign NA
		}
	}
	if(length(tip_states)!=length(tree$tip.label)) return(list(success=FALSE, error=sprintf("Number of provided tip states (%d) does not match number of tips in the tree (%d)",length(tip_states),length(tree$tip.label))))
	
	# only consider subtree with known tip states
	known_tips = which(!is.na(tip_states));
	if(length(known_tips)<=1) return(list(success=FALSE, error=sprintf("All or almost all tips have unknown or ambiguous state")))
	if(length(known_tips)<length(tip_states)){
		extraction	= get_subtree_with_tips(tree, only_tips=known_tips, omit_tips=FALSE, collapse_monofurcations=TRUE, force_keep_root=TRUE);
		tree		= extraction$subtree;
		tip_states	= tip_states[extraction$new2old_tip]
	}
	Ntips  = length(tree$tip.label)
	Nedges = nrow(tree$edge)
	
	# fit a symmetric Mk model to the data
	fit = fit_symmetric_mk(tree, Nstates = Nstates, tip_states	= tip_states, rate_model = "ER", Ntrials = 1)
	if(fit$success){	
		rate = fit$transition_matrix[1,2]
	}else if(!is.null(fit$guess_rate)){
		# fitting failed, so use first-fuess rate
		rate = fit$guess_rate
	}else{
		return(list(success=FALSE, error=sprintf("Unable to even guess the Mk rate: %s",fit$error)))
	}
	
	return(list(success=TRUE, rate=rate));
}



# return the ages of all branching events in an ultrametric timetree, i.e. all nodes accounting for multifurcations
# if the tree is purely bifurcating, this is the same as getting all node ages
# However, if the tree includes multifurcations, these are counted multiple times (since they represent multiple nearby bifurcations)
# Monofurcations are not returned
# Assumes that the tree is ultrametric
get_all_branching_ages = function(tree){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode;
	depths = get_mean_depth_per_node_CPP(	Ntips			= Ntips,
											Nnodes			= Nnodes,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
											edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length));
	Nchildren = get_child_count_per_node_CPP(	Ntips			= Ntips,
												Nnodes			= Nnodes,
												Nedges			= nrow(tree$edge),
												tree_edge		= as.vector(t(tree$edge))-1);
	branching_ages = rep(depths,times=Nchildren-1);
	return(branching_ages);
}


get_child_count_per_node = function(tree){
	Nchildren = get_child_count_per_node_CPP(	Ntips			= length(tree$tip.label),
												Nnodes			= tree$Nnode,
												Nedges			= nrow(tree$edge),
												tree_edge		= as.vector(t(tree$edge))-1);
	return(Nchildren)
}


# given a piecewise polynomial (splines) function f(x), defined as a time series on some x-grid, calculate its antiderivative A(x):=\int_{Xstart}^x f(u) du for an arbitrary number of target x values
# this function is most efficient when the requested target x-values are monotonically increasing or decreasing
get_antiderivative_of_splines_function = function(	Xgrid,			# numeric vector of size NG, listing x-values in ascending order
													Xstart,			# numeric, lower end of the integration, i.e. x-value where antiderivative is set to zero
													Ygrid,			# numeric vector of size NG, listing y-values along Xgrid
													splines_degree,	# integer, either 0,1,2 or 3, specifying the splines degree assumed for Y between grid points
													Xtarget){		# numeric vector of size N, specifying the target x values on which to evaluate the antiderivative. The function is most efficient if Xtarget are in ascending or descending order.
	A = get_antiderivative_CPP(	Xgrid, Xstart, Ygrid, splines_degree, Xtarget);
	return(A);
}


# given a piecewise polynomial (splines) function f(x), defined as a time series on some x-grid, calculate its first derivative A(x):=df(x)/dx at an arbitrary number of target x values
# this function is most efficient when the requested target x-values are monotonically increasing or decreasing
get_derivative_of_splines_function = function(	Xgrid,			# numeric vector of size NG, listing x-values in ascending order
												Ygrid,			# numeric vector of size NG, listing y-values along Xgrid
												splines_degree,	# integer, either 0,1,2 or 3, specifying the splines degree assumed for Y between grid points
												Xtarget){		# numeric vector of size N, specifying the target x values on which to evaluate the derivative. The function is most efficient if Xtarget are in ascending or descending order.
	D = get_derivative_CPP(Xgrid, Ygrid, splines_degree, Xtarget);
	return(D);
}




# given a lineages-through-time curve, defined as a time series on some discrete age grid, extract the branching ages that would have generated that LTT
# ages[] should be a 1D vector of ages (time before present) in ascending order
# LTT[] should be a 1D vector of the same size as ages[], listing the number of lineages at each age
# the LTT is assumed to be linear between adjacent age grid points
# branching points will be associated with those times where the LTT passes through an integer value
get_branching_ages_from_LTT = function(ages, LTT){
	results = get_branching_ages_from_LTT_CPP(ages, LTT);
	if(!results$success){
		return(list(success=FALSE, error=results$error));
	}else{
		return(list(success=TRUE, branching_ages = results$branching_ages));
	}
}


# given some density curve on an X-interval, define a non-uniform X-grid on that interval so that the density of grid points reflects the requested density
# this can be used for example to define an age grid, with the grid density reflecting the number of lineages in a timetree at any given age, e.g. for fitting purposes
# the density curve is specified as a piecewise linear function. The density must be non-negative, and have non-zero total area under the curve.
get_inhomogeneous_grid_1D = function(	Xstart,
										Xend, 
										Ngrid, 		# integer, number of grid points to return, including the edges Xstart & Xend
										densityX, 	# numeric vector of size ND, listing X-values for defining the density, in ascending order
										densityY,	# numeric vector of size ND, listing density values at densityX. 
										extrapolate = FALSE){	# extrapolate density grid as needed, to cover Xstart & Xend. The density will be extrapolated as a constant.
	if(Ngrid<2){
		stop(sprintf("Ngrid must be at least 2"));
	}else if(densityX[1]>=tail(densityX,1)){
		stop(sprintf("Values in densityX must be strictly increasing"))
	}
	if(Xstart<densityX[1]){
		if(extrapolate){
			densityX = c(Xstart,densityX);
			densityY = c(densityY[1], densityY);
		}else{
			stop(sprintf("Xstart (%g) is not covered by the density grid (which starts at %g). Consider setting extrapolate=TRUE.",Xstart,densityX[1]));
		}
	}else if(Xend>tail(densityX,1)){
		if(extrapolate){
			densityX = c(densityX,Xend);
			densityY = c(densityY,tail(densityY,1));
		}else{
			stop(sprintf("Xend (%g) is not covered by the density grid (which ends at %g). Consider setting extrapolate=TRUE.",Xend,tail(densityX,1)));
		}
	}
	return(get_inhomogeneous_grid_1D_CPP(	Xstart		= Xstart, 
											Xend		= Xend, 
											Ngrid		= Ngrid, 
											densityX	= densityX, 
											densityY	= densityY, 
											xepsilon	= 0.0000001*(Xend-Xstart)/Ngrid));
}



# given some distribution of numbers, define a non-uniform X-grid so that the density of grid points reflects the density of the provided random numbers
# This can be used for example to define an age grid, with the grid density reflecting the number of nodes in a timetree at any given age interval, e.g. for fitting purposes
# The returned Xgrid is guaranteed to include Xstart and Xend, to have size Ngrid, and to be non-decreasing.
# Note that if randomX does not cover the full requested interval [Xstart,Xend], then the first and last grid cells might deviate from the intended density
# Also note that if the randomX data is has Dirac-distribution peaks (i.e., many values are equal), then the returned grid might also include multiple grid points at the same location (i.e., the grid is not strictly monotonic)
get_inhomogeneous_grid_1D_from_samples = function(	Xstart,
													Xend, 
													Ngrid, 		# integer, number of grid points to return, including the edges Xstart & Xend
													randomX){ 	# numeric vector, listing random values defining the target density. Must be sorted in ascending order. Typically randomX should roughly cover the full requested interval [Xstart,Xend]
	if(Ngrid<2){
		return(list(success=FALSE, error=sprintf("Ngrid must be at least 2")))
	}else if(randomX[1]>=tail(randomX,1)){
		return(list(success=FALSE, error=sprintf("Values in randomX must be strictly increasing")))
	}
	randomX 		= randomX[(randomX>=Xstart) & (randomX<=Xend)] # ignore samples outside of the requested range
	if(length(randomX)==0) return(list(success=FALSE, error=sprintf("No randomX values found within requested interval")))
	Nrandom			= length(randomX)
	unique_randoms	= 1 + Nrandom - rev(c(1,1+which(diff(rev(randomX))<0))) # indices of unique values in randomX, keeping always the last instance of each value. Assuming that randomX is monotonically increasing.
	randomCDF_X 	= randomX[unique_randoms]
	randomCDF_P		= unique_randoms/Nrandom
	if(Xstart<randomCDF_X[1]){
		randomCDF_X = c(Xstart,randomCDF_X)
		randomCDF_P = c(0,randomCDF_P)
	}
	if(Xend>randomCDF_X[length(randomCDF_X)]){
		randomCDF_X = c(randomCDF_X,Xend)
		randomCDF_P = c(randomCDF_P,1.0001)
	}
	Xgrid = approx(x=randomCDF_P, y=randomCDF_X, xout=seq(from=0,to=1,length.out=Ngrid), yleft=Xstart, yright=Xend)$y
	return(list(success=TRUE,grid=Xgrid))
}



# calculate the pulled speciation rate (PSR) of an HBD congruence class for a given pulled diversification rate (PDR) and product rho*lambda(0)
get_PSR_from_PDR_HBD = function(oldest_age,
								age_grid,				# numeric vector of size NG, listing grid ages in ascending order. Must cover at least age0 and oldest_age.
								PDR				= 0,	# numeric vector of size NG, listing PDRs on the corresponding grid points
								age0			= 0,	# non-negative numeric, specifying the age at which rholambda0 is given, i.e. rholambda0=rho(age0)*lambda(age0)
								rholambda0		= 1,	# positive numeric, product rho(age0)*lambda(age0), where rho is the sampling fraction and lambda is the speciation rate
								splines_degree	= 1,	# either 1, 2 or 3, specifying the degree of the splines defined by the PDR on the age grid.
								relative_dt		= 1e-3,	# numeric, maximum relative time step allowed for integration. Smaller values increase integration accuracy. Typical values are 0.0001-0.001.
								include_nLTT0	= FALSE){	# (logical) whether to also calculate the ratio nLTT0:=LTT(age0)/LTT(present-day)
	# basic error checking
	if(is.null(PDR)) stop("Missing PDR")
	if(is.null(age_grid) || (length(age_grid)<=1)){
		if((!is.null(PDR)) && (length(PDR)!=1)) return(list(success=FALSE, error=sprintf("Invalid number of PDR values; since no age grid was provided, you must provide a single (constant) PDR")))
		# create dummy age grid
		NG 			= 2;
		age_grid	= seq(from=0,to=oldest_age,length.out=NG)
		if(!is.null(PDR)) PDR = rep(PDR,times=NG);
	}else{
		NG = length(age_grid);
		if((age_grid[1]>oldest_age) || (age_grid[NG]<oldest_age)) return(list(success=FALSE, error=sprintf("Age grid must cover the entire requested age interval, including oldest_age (%g)",oldest_age)))
		if((age_grid[1]>age0) || (age_grid[NG]<age0)) return(list(success=FALSE, error=sprintf("Age grid must cover the entire requested age interval, including age0 (%g)",age0)))
		if(include_nLTT0 && (age_grid[1]>0) || (age_grid[NG]<0)) return(list(success=FALSE, error=sprintf("Age grid must cover the present-day age (0) in order to calculate nLTT0")))
		if((!is.null(PDR)) && (length(PDR)!=1) && (length(PDR)!=NG)) return(list(success=FALSE, error=sprintf("Invalid number of PDR values; since an age grid of size %d was provided, you must either provide one or %d PDR",NG,NG)))
		if((!is.null(PDR)) && (length(PDR)==1)) PDR = rep(PDR,times=NG);
	}
	if(rholambda0<=0) return(list(success=FALSE, error=sprintf("rholambda0 must be strictly positive; instead, got %g",rholambda0)))
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success=FALSE, error=sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
	if(age_grid[1]>tail(age_grid,1)) return(list(success=FALSE, error=sprintf("Values in age_grid must be strictly increasing")))

	# calculate PSR
	results = get_PSR_from_PDR_HBD_CPP(	age0			= age0,
										oldest_age 		= oldest_age,
										age_grid		= age_grid,
										PDR				= PDR,
										rholambda0		= rholambda0,
										splines_degree	= splines_degree,
										relative_dt		= relative_dt,
										include_nLTT0	= include_nLTT0)
	if(results$success){
		return(list(success	= TRUE, 
					ages	= results$refined_age_grid,	# numeric vector listing (potentially refined) grid ages, spanning [max(0,age_grid[1]), oldest_age]
					PSR		= results$PSR, 				# numeric vector of the same size as ages[], listing the PSR on the refined grid
					nLTT0	= (if(include_nLTT0) results$nLTT0 else NULL)))
	}else{
		return(list(success=FALSE, error=results$error))
	}
}





# calculate the pulled speciation rate (PSR) of an HBD model for a given speciation rate (lambda), extinction rate (mu) and the sampling fraction (rho0) at some age0>=0
get_PSR_of_HBD_model = function(oldest_age,					# oldest age until which to calculate the PSR
								age_grid		= NULL,		# numeric vector of size NG, listing grid ages in ascending order. Must cover at least age0 and oldest_age. Can also be NULL, in which case the same lambda & mu are assumed everywhere.
								lambda			= 0,		# numeric vector of size NG, listing speciation rates on the corresponding grid points. Can also be a single constant.
								mu				= 0,		# numeric vector of size NG, listing extinction rates on the corresponding grid points. Can also be a single constant.
								age0			= 0,		# numeric, age (time before present) at which the sampling fraction (rho) is specified
								rho0			= 1,		# positive numeric, sampling fraction at age0
								splines_degree	= 1,		# either 1, 2 or 3, specifying the degree of the splines defined by the PDR on the age grid.
								relative_dt		= 1e-3){	# numeric, maximum relative time step allowed for integration. Smaller values increase integration accuracy. Typical values are 0.0001-0.001.
	# basic error checking
	if(rho0<=0) return(list(success=FALSE, error=sprintf("rho must be strictly positive; instead, got %g",rho0)))
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success=FALSE, error=sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
	if(is.null(age_grid) || (length(age_grid)<=1)){
		# this model is time-independent, i.e. lambda & mu are constant over time
		if(length(lambda)!=1) return(list(success = FALSE, error = sprintf("Invalid number of lambda values (%d); since no age grid was provided, you must either provide a single (constant) lambda or none",length(lambda))))
		if(length(mu)!=1) return(list(success = FALSE, error = sprintf("Invalid number of mu values (%d); since no age grid was provided, you must provide a single (constant) mu",length(mu))))
		NG = 1
		constant_rates = TRUE
	}else{
		NG = length(age_grid);
		if((age_grid[1]>oldest_age) || (age_grid[NG]<oldest_age)) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including oldest_age (%g)",oldest_age)))
		if((age_grid[1]>age0) || (age_grid[NG]<age0)) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including age0 (%g)",age0)))
		if((length(lambda)!=1) && (length(lambda)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of lambda values (%d); since an age grid of size %d was provided, you must either provide one or %d lambdas",length(lambda),NG,NG)))
		if((length(mu)!=1) && (length(mu)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of mu values (%d); since an age grid of size %d was provided, you must either provide one or %d mus",length(mu),NG,NG)))
		constant_rates = ((length(lambda)==1) || all(diff(lambda)==0)) && ((length(mu)==1) || all(diff(mu)==0)) # whether this model is in fact a constant-rates model
		if((length(lambda)==1)) lambda = rep(lambda,times=NG)
		if((length(mu)==1)) mu = rep(mu,times=NG)
		if(age_grid[1]>tail(age_grid,1)) return(list(success=FALSE, error=sprintf("Values in age_grid must be strictly increasing")))
	}

	# calculate PSR from lambda, mu & rho0=rho(age0)
	if(constant_rates){
		results = get_PSR_of_CR_HBD_model_CPP(	age0			= age0,
												oldest_age 		= oldest_age,
												lambda			= lambda[1],
												mu				= mu[1],
												rho0			= rho0,
												relative_dt		= relative_dt)
	}else{
		results = get_PSR_of_HBD_model_CPP(	age0			= age0,
											oldest_age 		= oldest_age,
											age_grid		= age_grid,
											lambda			= lambda,
											mu				= mu,
											rho0			= rho0,
											splines_degree	= splines_degree,
											relative_dt		= relative_dt)
		
	}
	if(results$success){
		return(list(success	= TRUE, 
					ages	= results$ages, # numeric vector listing (potentially refined) grid ages in ascending order, spanning [max(0,age_grid[1]), oldest_age]
					PSR		= results$PSR))	# numeric vector of the same size as ages[], listing the PSR on the refined grid
	}else{
		return(list(success=FALSE, error=results$error))
	}
}




# generate a random variable bounded from below but not from above
# used to pick random start params for fitting trials
random_semiboxed_left = function(lower_bound, default, typical_scale, orders_of_magnitude){
	if((default==0) && (lower_bound==0)){
		if(typical_scale==0){
			return(runif(n=1,min=0,max=1))
		}else{
			return(abs(typical_scale) * 10**runif(n=1,min=0,max=orders_of_magnitude))
		}
	}else if((default>lower_bound) && (default>0)){
		return(if(rbinom(n=1,size=1,prob=0.5)==1) (default * 10**runif(n=1,min=0,max=orders_of_magnitude)) else (lower_bound + (default-lower_bound)*runif(n=1,min=0,max=1)))
	}else if(default>lower_bound){
		return(lower_bound + (default-lower_bound) * 10**runif(n=1,min=-orders_of_magnitude/2,max=orders_of_magnitude/2))
	}else{
		return(lower_bound + (0-lower_bound) * 10**runif(n=1,min=-orders_of_magnitude/2,max=orders_of_magnitude/2))
	}
}


# choose random parameter values within boxed constraints
# each lower & upper bound may be Inf
# defaults[], lower_bounds[], upper_bounds[] and scales[] must be 1D numeric vectors of the same length, and must not include NaN or NA 
# lower_bounds[] and upper_bounds[] may include +Inf and -Inf
get_random_params = function(defaults, lower_bounds, upper_bounds, scales, orders_of_magnitude){
	start_values = defaults
	boxed_left	 = which((!is.infinite(lower_bounds)) & is.infinite(upper_bounds))
	boxed_right	 = which((!is.infinite(upper_bounds)) & is.infinite(lower_bounds))
	boxed_dual   = which(!(is.infinite(lower_bounds) | is.infinite(upper_bounds))); # determine fitted params that are boxed, i.e. constrained to within finite lower & upper bounds
	unboxed 	 = which(is.infinite(lower_bounds) & is.infinite(upper_bounds))
	if(length(boxed_dual)>0) 	start_values[boxed_dual] = lower_bounds[boxed_dual] + (upper_bounds[boxed_dual]-lower_bounds[boxed_dual]) * runif(n=length(boxed_dual),min=0,max=1)
	if(length(unboxed)>0) 	 	start_values[unboxed]	 = 10**runif(n=length(unboxed), min=-orders_of_magnitude/2.0, max=orders_of_magnitude/2.0) * start_values[unboxed]
	if(length(boxed_left)>0) 	start_values[boxed_left] = sapply(boxed_left, FUN=function(fp){ random_semiboxed_left(lower_bound=lower_bounds[fp], default=start_values[fp], typical_scale=scales[fp], orders_of_magnitude=orders_of_magnitude) })
	if(length(boxed_right)>0) 	start_values[boxed_right]= sapply(boxed_right, FUN=function(fp){ -random_semiboxed_left(lower_bound=-upper_bounds[fp], default=-start_values[fp], typical_scale=scales[fp], orders_of_magnitude=orders_of_magnitude) })
	start_values = pmax(lower_bounds,pmin(upper_bounds,start_values))
	return(start_values)
}


# check validity and sanitize (standardize format) of model parameters for fitting
# used by various fitting routines
sanitize_parameters_for_fitting = function(	param_values,				# numeric vector of size NP, specifying fixed values for some or all parameters. For fitted (i.e. non-fixed) parameters, use NaN or NA.
											param_guess		= NULL,		# numeric vector of size NP, listing an initial guess for each parameter. For fixed parameters, guess values are ignored.
											param_min		= -Inf,		# numeric vector of size NP, specifying lower bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same lower bound is assumed for all params.
											param_max		= +Inf,		# numeric vector of size NP, specifying upper bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same upper bound is assumed for all params.
											param_scale		= NULL){	# numeric vector of size NP, specifying typical scales for the model parameters. For fixed parameters, scales are ignored. If NULL, scales are automatically estimated from other information (such as provided guess and bounds). May also be a single scalar, in which case the same scale is assumed for all params.
	NP 			= length(param_values);
	param_names	= names(param_values);
	if(is.null(param_guess)){
		if(any(!is.finite(param_values))){
			return(list(success=FALSE, error=sprintf("Missing guessed parameter values")))
		}else{
			param_guess = rep(NA, times=NP);
		}
	}
	if(length(param_guess)!=NP){
		return(list(success=FALSE, error=sprintf("Number of guessed parameters (%d) differs from number of model parameters (%d)",length(param_guess),NP)))
	}else if(!is.null(param_names)){
		names(param_guess) = param_names;
	}
	if((!is.null(param_names)) && (length(param_names)!=NP)){
		return(list(success=FALSE, error=sprintf("Number of parameter names (%d) differs from number of model parameters (%d)",length(param_names),NP)))
	}
	if(is.null(param_min)){
		param_min = rep(-Inf,times=NP);
	}else if(length(param_min)==1){
		param_min = rep(param_min,times=NP);
	}else if(length(param_min)!=NP){
		return(list(success=FALSE, error=sprintf("Length of param_min[] (%d) differs from number of model parameters (%d)",length(param_min),NP)))
	}
	if(is.null(param_max)){
		param_max = rep(+Inf,times=NP);
	}else if(length(param_max)==1){
		param_max = rep(param_max,times=NP);
	}else if(length(param_max)!=NP){
		return(list(success=FALSE, error=sprintf("Length of param_max[] (%d) differs from number of model parameters (%d)",length(param_max),NP)))
	}
	if(is.null(param_scale)){
		param_scale = rep(NA,times=NP);
	}else if(length(param_scale)==1){
		param_scale = rep(param_scale,times=NP);
	}else if(length(param_scale)!=NP){
		return(list(success=FALSE, error=sprintf("Length of param_scale[] (%d) differs from number of model parameters (%d)",length(param_scale),NP)))
	}
	param_values[is.nan(param_values)] = NA # standardize representation of non-fixed params
	param_scale[is.nan(param_scale)] = NA	# standardize representation of unknown param scales		
		
	# determine which parameters are to be fitted
	fitted_params	= which(is.na(param_values))
	fixed_params	= which(!is.na(param_values))
	NFP				= length(fitted_params);
	param_guess[fixed_params] = param_values[fixed_params] # make sure guessed param values are consistent with fixed param values

	if(any(!is.finite(param_guess))) return(list(success=FALSE, error=sprintf("Some guessed parameter values are NA or NaN or Inf; you must specify a valid guess for each non-fixed model parameter")));
	if(any((!is.na(param_scale)) & (param_scale==0))) return(list(success=FALSE, error=sprintf("Some provided parameter scales are zero; expecting non-zero scale for each parameter")));
	
	# determine typical parameter scales, whenever these are not provided
	for(p in fitted_params){
		if(is.na(param_scale[p])){
			if(param_guess[p]!=0){
				param_scale[p] = abs(param_guess[p]);
			}else if((is.finite(param_min[p]) && (param_min[p]!=0)) || (is.finite(param_max[p]) && (param_max[p]!=0))){
				param_scale[p] = mean(abs(c((if(is.finite(param_min[p]) && (param_min[p]!=0)) param_min[p] else NULL), (if(is.finite(param_max[p]) && (param_max[p]!=0)) param_max[p] else NULL))));
			}else{
				param_scale[p] = 1;
			}
		}
	}
	
	return(list(success			= TRUE,
				NP				= NP,
				NFP				= NFP,
				param_names		= param_names,
				param_values	= param_values,
				param_guess		= param_guess,
				param_min		= param_min,
				param_max		= param_max,
				param_scale		= param_scale,
				fitted_params	= fitted_params,
				fixed_params	= fixed_params))
}
	



# function for reformatting/sanitizing generic fit parameters provided by the user
# used by various fitting routines
# prepare_generic_fit_params = function(	param_values,					# numeric vector of size NP, specifying fixed values for a some or all parameters. For fitted (i.e. non-fixed) parameters, use NaN or NA.
# 										param_guess			= NULL,		# numeric vector of size NP, listing an initial guess for each parameter. For fixed parameters, guess values are ignored.
# 										param_min			= -Inf,		# numeric vector of size NP, specifying lower bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same lower bound is assumed for all params.
# 										param_max			= +Inf,		# numeric vector of size NP, specifying upper bounds for the model parameters. For fixed parameters, bounds are ignored. May also be a single scalar, in which case the same upper bound is assumed for all params.
# 										param_scale			= NULL){	# numeric vector of size NP, specifying typical scales for the model parameters. For fixed parameters, scales are ignored. If NULL, scales are automatically estimated from other information (such as provided guess and bounds). May also be a single scalar, in which case the same scale is assumed for all params.
# 	NP = length(param_values)
# 	param_names = names(param_values);
# 	if(is.null(param_guess)){
# 		if(any(is.finite(param_values))){
# 			return(list(success=FALSE, error=sprintf("Missing guessed parameter values")))
# 		}else{
# 			param_guess = rep(NA, times=NP);
# 		}
# 	}
# 	if(length(param_guess)!=NP){
# 		return(list(success=FALSE, error=sprintf("Number of guessed parameters (%d) differs from number of model parameters (%d)",length(param_guess),NP)))
# 	}else if(!is.null(param_names)){
# 		names(param_guess) = param_names;
# 	}
# 	if((!is.null(param_names)) && (length(param_names)!=NP)){
# 		return(list(success=FALSE, error=sprintf("Number of parameter names (%d) differs from number of model parameters (%d)",length(param_names),NP)))
# 	}
# 	if(is.null(param_min)){
# 		param_min = rep(-Inf,times=NP);
# 	}else if(length(param_min)==1){
# 		param_min = rep(param_min,times=NP);
# 	}else if(length(param_min)!=NP){
# 		return(list(success=FALSE, error=sprintf("Length of param_min[] (%d) differs from number of model parameters (%d)",length(param_min),NP)))
# 	}
# 	if(is.null(param_max)){
# 		param_max = rep(+Inf,times=NP);
# 	}else if(length(param_max)==1){
# 		param_max = rep(param_max,times=NP);
# 	}else if(length(param_max)!=NP){
# 		return(list(success=FALSE, error=sprintf("Length of param_max[] (%d) differs from number of model parameters (%d)",length(param_max),NP)))
# 	}
# 	if(is.null(param_scale)){
# 		param_scale = rep(NA,times=NP);
# 	}else if(length(param_scale)==1){
# 		param_scale = rep(param_scale,times=NP);
# 	}else if(length(param_scale)!=NP){
# 		return(list(success=FALSE, error=sprintf("Length of param_scale[] (%d) differs from number of model parameters (%d)",length(param_scale),NP)))
# 	}
# 	return(list(success=TRUE, param_values=param_values, param_guess=param_guess, param_min=param_min, param_max=param_max, param_scale))
# }


# given an undirected graph (nodes,edges), find its maximal connected subgraphs
# any two nodes may be connected by zero, one or multiple edges
# edges[] should be a 2D array of size Nedges x 2, listing source & target nodes of the graph
get_connected_subgraphs = function(Nnodes, edges){
	results = split_undirected_graph_CPP(Nnodes=Nnodes, Nedges=nrow(edges), edges = as.vector(t(edges))-1);
	return(list(Nsubgraphs 		= results$Nsubgraphs,
				subgraph2nodes	= lapply(1:results$Nsubgraphs, FUN=function(n){ results$subgraph2nodes[[n]]+1 }),
				subgraph2edges	= lapply(1:results$Nsubgraphs, FUN=function(n){ results$subgraph2edges[[n]]+1 }),
				node2subgraph	= results$node2subgraph+1,
				edge2subgraph	= results$edge2subgraph+1));
}



# Sarting from a mapping pool-->group, calculate the reverse mapping group-->member_list
# this function returns a list of length Ngroups, where each element is a 1D vector listing member indices
# values in pool2group[] must be between 1 and Ngroups, or NA or <=0 (in which case the item is not affiliated with any group)
get_member_lists_from_group_assignments = function(Ngroups, pool2group){
	if(Ngroups==0) return(list())
	pool2group[is.na(pool2group)] = -1;
	results = get_member_lists_from_group_assignments_CPP(Ngroups=Ngroups, pool2group=pool2group-1)
	return(lapply(1:Ngroups,FUN=function(g){ results$group2members[[g]]+1 }))
}



# evaluate a mathematical expression (univariate function of X) for various X-values
# the input X[] can either be a 1D vector or a 2D matrix
evaluate_univariate_expression = function(expression, Xname="x", X){
	if(is.vector(X)){
		results = evaluate_univariate_expression_CPP(expression=expression, Xname=Xname, X=X);
		return(list(success=results$success, error=results$error, Y=results$Y))
	}else if(is.matrix(X)){
		results = evaluate_univariate_expression_CPP(expression=expression, Xname=Xname, X=as.vector(t(X)));
		if(!results$success) return(list(success=FALSE, error=results$error))
		return(list(success=TRUE, Y=matrix(results$Y,ncol=ncol(X),byrow=TRUE)))
	}else{
		return(list(success=FALSE, error="Unknown data format X: Expecting either a vector or a matrix"))
	}
}


# Split tree into pairs of sister-tips, such that the paths within distinct pairs do not overlap
# If the input tree only contains monofurcations and bifurcations (recommended), it is guaranteed that at most one unpaired tip will be left (i.e., if Ntips was odd)
extract_independent_sister_tips = function(tree){
	results = extract_independent_sister_tips_CPP(	Ntips		= length(tree$tip.label),
													Nnodes		= tree$Nnode,
													Nedges		= nrow(tree$edge),
													tree_edge	= as.vector(t(tree$edge))-1);
	tip_pairs = matrix(as.integer(results$tip_pairs),ncol=2,byrow=TRUE) + 1L;
	return(tip_pairs);
}


# calculate geodesic angle (aka. central angle, in radians) between two geographical locations (assuming the Earth is a sphere)
# based on the Vincenty formula with equal major and minor axis
# geographic coordinates should be given in decimal degrees
geodesic_angle = function(latitude1, longitude1, latitude2, longitude2){
	theta1	= pi*latitude1/180;
	theta2	= pi*latitude2/180;
	phi1 	= pi*longitude1/180.0;
	phi2 	= pi*longitude2/180.0;
	delta	= abs(phi1-phi2);
	angle 	= abs(atan2(sqrt((cos(theta2)*sin(delta))^2 + (cos(theta1)*sin(theta2)-sin(theta1)*cos(theta2)*cos(delta))^2), (sin(theta1)*sin(theta2)+cos(theta1)*cos(theta2)*cos(delta))));
	return(angle);
}

# calculate geodesic angles (in radians) between pairs of coordinates
# this function returns a 1D vector of size equal to the size of the input lists latitudes1[], latitudes2[] etc
geodesic_angles = function(latitudes1, longitudes1, latitudes2, longitudes2){
	return(geodesic_angles_CPP(latitudes1, longitudes1, latitudes2, longitudes2));
}

# calculate all geodesic angles (in radians) between two sets of coordinates
# This function returns a 2D matrix of size N1 x N2
all_pairwise_geodesic_angles = function(latitudes1, longitudes1, latitudes2, longitudes2){
	angles = get_all_pairwise_geodesic_angles_CPP(latitudes1, longitudes1, latitudes2, longitudes2);
	return(matrix(angles, ncol=length(latitudes2), byrow=TRUE));
}


# calculate expectation of sin^2(omega) across multiple independent contrasts with non-equal time steps, where omega is the central transition angle, for the spherical Brownian Motion
get_expected_SBM_sinsquare = function(time_steps, diffusivity, radius){
	return((2/3)*(1-mean(exp(-6*(diffusivity/(radius^2))*time_steps))))
}

expected_SBM_distance = function(time_step, diffusivity, radius){
	tD = time_step * diffusivity/radius^2;
	return(radius * SBM_get_average_transition_angle_CPP(tD = tD, max_error = 1e-7, max_Legendre_terms = 200))
}


# autocorrelation function ACF(t) of Spherical Brownian Motion with constant diffusivity D
# ACF(t) := E <n(t),n(0)> = E cos(omega(t))
# where n(t) is the unit vector pointing to the particle's random location on the unit phere at time t, and <,> is the scalar product, and omega is the transition angle.
get_SBM_ACF = function(times, diffusivity, radius){
	return(exp(-2*(diffusivity/radius^2)*times))
}



# given some rooted timetree, extract information about various sampling and branching events, as defined by the homogenous birth-death-sampling model with continuous sampling and concentrated sampling efforts
# monofurcating nodes are interpreted as sampled nodes, i.e. sampled lineages with additional subsequently sampled descendants
# bifurcating and multifurcating nodes are interpreted as branching events, with multifurcations counted multiple times (i.e., as if they are first split into bifurcations)
extract_HBDS_events_from_tree = function(	tree,
											root_age = NULL, # optional numeric, the age of the root. Can be used to define a time offset, e.g. if the last tip was not actually sampled at the present. If NULL, this will be calculated from the treem and it will be assumed that the last tip was sampled at the present
											CSA_ages = numeric(0), # 1D vector listing ages of concentrated sampling efforts, in ascending order
											age_epsilon = 0){
	Ntips  	= length(tree$tip.label)
	Nnodes 	= tree$Nnode
	NCSA		= length(CSA_ages)
	if((NCSA>=2) && (CSA_ages[1]>tail(CSA_ages,1))) return(list(success=FALSE, error="CSA_ages must be in ascending order"))
	
	# determine branching ages & tip ages & sampled node ages
	if(is.null(root_age)) root_age = get_tree_span(tree)$max_distance
	node2Nchildren = get_child_count_per_node_CPP(	Ntips		= Ntips,
													Nnodes		= Nnodes,
													Nedges		= nrow(tree$edge),
													tree_edge	= as.vector(t(tree$edge))-1);
	clade_heights = get_distances_from_root_CPP(Ntips		= Ntips,
												Nnodes		= Nnodes,
												Nedges		= nrow(tree$edge),
												tree_edge	= as.vector(t(tree$edge))-1,
												edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length))
	sampled_nodes		= which(node2Nchildren==1)
	sampled_node_ages 	= pmax(0, root_age - clade_heights[Ntips + sampled_nodes])
	tip_ages 			= pmax(0, root_age - clade_heights[1:Ntips])
	branchings			= which(node2Nchildren>=2) # determine which nodes are branchings
	branching_ages		= pmax(0, root_age - clade_heights[Ntips + branchings])
	branching_ages		= rep(branching_ages,times=node2Nchildren[branchings]-1) # count multifurcations multiple times
	
	if(length(CSA_ages)>0){
		# determine which tips were sampled during concentrated sampling efforts versus due to Poissonian sampling
		sorted2original = order(tip_ages)
		CSAbinning		= place_sorted_values_into_bins_CPP(items=tip_ages[sorted2original], bin_mins=CSA_ages-age_epsilon, bin_maxs=CSA_ages+age_epsilon)
		tip2CSA 		= 1L + CSAbinning$item2bin
		Ptip_ages		= tip_ages[sorted2original[which(CSAbinning$item2bin<0)]]
		CSA2tips		= lapply(1:NCSA, FUN=function(ce){ sorted2original[1L + CSAbinning$bin2items[[ce]]] });
		CSA_tip_counts 	= sapply(1:NCSA, FUN=function(ce){ length(CSA2tips[[ce]]) })

		# determine which sampled nodes were sampled during concentrated sampling efforts versus due to Poissonian sampling
		sorted2original = order(sampled_node_ages)
		CSAbinning 		= place_sorted_values_into_bins_CPP(items=sampled_node_ages[sorted2original], bin_mins=CSA_ages-age_epsilon, bin_maxs=CSA_ages+age_epsilon)
		Pnode_ages		= sampled_node_ages[sorted2original[which(CSAbinning$item2bin<0)]]
		CSA2nodes		= lapply(1:NCSA, FUN=function(ce){ sampled_nodes[sorted2original[1L + CSAbinning$bin2items[[ce]]]] });
		CSA_node_counts = sapply(1:NCSA, FUN=function(ce){ length(CSA2nodes[[ce]]) })
	}else{
		CSA2tips 		= integer(0)
		CSA_tip_counts 	= integer(0)
		CSA2nodes 		= integer(0)
		CSA_node_counts = integer(0)
		Ptip_ages 		= sort(tip_ages)
		Pnode_ages 		= sort(sampled_node_ages)
	}
	
	return(list(concentrated_tips 			= CSA2tips,
				concentrated_tip_counts 	= CSA_tip_counts,
				concentrated_nodes 			= CSA2nodes,
				concentrated_node_counts 	= CSA_node_counts,
				Ptip_ages					= Ptip_ages, 	# ages of Poissonian tip ages, i.e. of Poissonian terminal sampling events, in ascending order
				Pnode_ages					= Pnode_ages,	# ages of Poissonian node ages, i.e. of Poissonian non-terminal sampling events, in ascending order
				branching_ages				= sort(branching_ages)))
}



# create a finer time grid, based on some old grid and a refinement factor, by splitting each interval of the original grid into refinement_factor sub-intervals
# old_grid[] must be non-decreasing
# refinement_factor must be an integer >= 1
refine_time_grid = function(old_grid, refinement_factor){
	old_diffs = diff(old_grid)
	new_diffs = rep(old_diffs/refinement_factor, rep(refinement_factor, length(old_diffs)))
	new_grid  = old_grid[1] + cumsum(new_diffs)
	return(new_grid)
}


# smoothen a time series using qubic splines
smoothen_time_series = function(times, values){
	refinement_factor 	= max(2,min(100,100*max(abs(diff(values)/sd(values)))))
	refined_times  		= refine_time_grid(old_grid=times, refinement_factor=refinement_factor)
	refined_values 		= evaluate_spline(Xgrid=times, Ygrid=values, splines_degree=3, Xtarget=refined_times, extrapolate="const", derivative=0)
	return(list(times = refined_times, values = refined_values))
}



generate_OU_time_series = function(	times,					# numeric vector of size NT, specifying the times at which to evaluate the OU process. A finer time grid means that the resulting time series will have finer structures, and look more "jagged"
									start_value,			# optional numeric. If NA or NaN or NULL, it will be chosen randomly from the stationary distribution
									stationary_mean,		# numeric
									stationary_std,			# non-negative numeric
									decay_rate,				# strictly positive numeric, in units 1/time
									constrain_min=-Inf,		# optional lower bound for the returned values
									constrain_max=+Inf,		# optional upper bound for the returned values
									smoothen = FALSE){		# logical, specifying whether to smoothen the generated time series using qubic splines onto a refined time grid
	values = get_Ornstein_Uhlenbeck_time_series_CPP(times			= times,
													start_value		= (if((!is.null(start_value)) && is.finite(start_value)) start_value else NaN),
													stationary_mean = stationary_mean,
													stationary_std	= stationary_std,
													decay_rate		= decay_rate)$values
	if(smoothen){
		smoothened = smoothen_time_series(times=times, values=values)
		times  = smoothened$times
		values = smoothened$values
	}
	return(list(times=times, values=pmin(constrain_max,pmax(constrain_min,values))))
}


# generate a time series according to a bounded Brownian Motion, i.e. a diffusing particle inside a finite interval with reflecting boundaries
generate_bounded_BM_time_series = function(	times,					# numeric vector of size NT
											diffusivity,			# strictly positive numeric, in units value^2/time
											lower,					# the lower bound of the reflected Brownian Motion. Either a single numeric (constant bound), or a numeric vector of the same length as times[] (time-dependent bound)
											upper,					# the upper bound of the reflected Brownian Motion. Either a single numeric (constant bound), or a numeric vector of the same length as times[] (time-dependent bound)
											start_value = NULL){	# optional numeric. If NA or NaN or NULL, it will be chosen randomly from the stationary distribution
	values = get_bounded_BM_time_series_CPP(times		= times,
											start_value	= (if((!is.null(start_value)) && is.finite(start_value)) start_value else NaN),
											diffusivity = diffusivity,
											lower		= lower,
											upper		= upper)$values
	return(list(times=times, values=values))
}


# estimate the parameters of a 1-dimensional Ornstein-Uhlenbeck process (proceeding over linear time) via maximum-likelihood, given a single univariate time series
# This is not a phylogenetic routine, i.e., it is not suitable for fitting an OU process evolving & branching along a phylogeny
# The OU process is parameterized through its stationary mean ("mean"), its stationary standard deviation ("std") and its decay_rate (aka. "lambda")
fit_OU1_model = function(	times,	# numeric vector of length NT, listing times in strictly ascending order
							values, # numeric vector of length NT, listing the observed values of the OU process
							min_mean			= -Inf,
							max_mean			= Inf,
							min_std				= 0,
							max_std				= Inf,
							min_decay_rate		= 0,
							max_decay_rate		= Inf,
							fixed_mean			= NA,		# optional numeric, specifying the fixed stationary mean. To fit the stationary mean, set this to NA or NaN.
							fixed_std			= NA,		# optional numeric, specifying the fixed stationary std. To fit the stationary std, set this to NA or NaN.
							fixed_decay_rate	= NA,		# optional numeric, specifying the fixed decay_rate. To fit the decay_rate, set this to NA.
							guess_mean			= NA,		# initial guess for the stationary mean of the OU process. If NA, this will be automatically chosen.
							guess_std			= NA,		# initial guess for the stationary standard deviation of the OU process. If NA, this will be automatically chosen.
							guess_decay_rate 	= NA,		# initial guess for the decay rate ("lambda") of the OU process. If NA, this will be automatically chosen.
							detrend				= TRUE,		# detrend the time series prior to fitting the OU process.
							Ntrials				= 1,		# number of fitting trials to perform, to avoid local non-global likelihood maxima
							max_start_attempts	= 1,		# integer, number of times to attempt finding a valid start point (per trial) before giving up. Randomly choosen start parameters may result in Inf/undefined objective, so this option allows the algorithm to keep looking for valid starting points.
							Nthreads			= 1,		# number of parallel threads to use when performing multiple fitting trials
							fit_control			= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
							verbose				= FALSE,	# boolean, specifying whether to print informative messages
							diagnostics			= FALSE,	# boolean, specifying whether to print detailed info (such as log-likelihood) at every iteration of the fitting. For debugging purposes mainly.
							verbose_prefix		= ""){		# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
	NT 		= length(times)
	Ndata 	= NT-1
	if(NT<2) return(list(success=FALSE, error="Not enough data points"))
	time_steps = diff(times)
	
	if(detrend){
		linfit 	= stats::lm(Y ~ X, data=data.frame(X=times, Y=values), na.action=na.omit)
		values	= values - (linfit$coefficients[1] + linfit$coefficients[2] * times) + mean(values, na.rm=TRUE)
		trend 	= linfit$coefficients[2]
	}
	
	# figure out reasonable parameter guesses where needed
	if(!is.finite(guess_mean)) guess_mean = mean(values)
	if(!is.finite(guess_std)) guess_std = sd(values)
	if(!is.finite(guess_decay_rate)) guess_decay_rate = mean(abs(diff(values)/(diff(times) * values[1:(NT-1)])))
	
	# sanitize model parameters
	sanitized_params = sanitize_parameters_for_fitting(	param_values	= c(fixed_mean, fixed_std, fixed_decay_rate),
														param_guess 	= c(guess_mean, guess_std, guess_decay_rate),
														param_min		= c(min_mean, min_std, min_decay_rate),
														param_max 		= c(max_mean, max_std, max_decay_rate),
														param_scale 	= NULL)
	if(!sanitized_params$success) return(list(success=FALSE, error=sanitized_params$error))
	NP				= sanitized_params$NP
	NFP				= sanitized_params$NFP
	param_names		= sanitized_params$param_names
	param_values	= sanitized_params$param_values
	param_guess		= sanitized_params$param_guess
	param_min		= sanitized_params$param_min
	param_max		= sanitized_params$param_max
	param_scale		= sanitized_params$param_scale
	fitted_params	= sanitized_params$fitted_params
	fixed_params	= sanitized_params$fixed_params

	# set fit-control options, unless provided by the caller
	if(is.null(fit_control)) fit_control = list()
	if(is.null(fit_control$step.min)) fit_control$step.min = 0.001
	if(is.null(fit_control$x.tol)) fit_control$x.tol = 1e-8
	if(is.null(fit_control$iter.max)) fit_control$iter.max = 1000
	if(is.null(fit_control$eval.max)) fit_control$eval.max = 2 * fit_control$iter.max * NFP
							
	# objective function: negated log-likelihood
	# input argument is the subset of fitted parameters, rescaled according to param_scale
	objective_function = function(fparam_values,trial){
		params = param_values; params[fitted_params] = fparam_values * param_scale[fitted_params];
		if(any(is.nan(params)) || any(is.infinite(params))) return(Inf); # catch weird cases where params become NaN
		# calculate joint likelihood of transitions, based on the fact that transitions are independent normally distributed variables
		stationary_mean = params[1]
		stationary_std 	= params[2]
		decay_rate 		= params[3]
		next_stds  		= stationary_std * sqrt(1-exp(-2*time_steps*decay_rate))
		next_means 		= stationary_mean + (values[1:(NT-1)]-stationary_mean) * exp(-time_steps*decay_rate)
		loglikelihood	= -sum(log(next_stds * sqrt(2*pi)) + ((values[2:NT]-next_means)^2)/(2*next_stds^2))
		if(!is.finite(loglikelihood)) loglikelihood = -Inf
		if(diagnostics) cat(sprintf("%s  Trial %s: loglikelihood %.10g\n",verbose_prefix,as.character(trial),loglikelihood))
		return(-loglikelihood)
	}
	
	# fit with various starting points
	fit_single_trial = function(trial){
		scales		 = param_scale[fitted_params]
		lower_bounds = param_min[fitted_params]
		upper_bounds = param_max[fitted_params]
		# randomly choose start values for fitted params (keep trying up to max_start_attempts times)
		Nstart_attempts = 0
		while(Nstart_attempts<max_start_attempts){
			# randomly choose start values for fitted params
			if(trial==1){
				start_values = param_guess[fitted_params]		
			}else{
				start_values = get_random_params(defaults=param_guess[fitted_params], lower_bounds=lower_bounds, upper_bounds=upper_bounds, scales=scales, orders_of_magnitude=4)
			}
			# check if start values yield NaN
			start_objective = objective_function(start_values/scales, trial)
			Nstart_attempts = Nstart_attempts + 1
			if(is.finite(start_objective)) break;
		}
		# run fit with chosen start_values
		if(is.finite(start_objective)){
			fit = stats::nlminb(start		= start_values/scales, 
								objective	= function(pars){ objective_function(pars, trial) },
								lower		= lower_bounds/scales, 
								upper		= upper_bounds/scales, 
								control		= fit_control)
			results = list(objective_value=fit$objective, fparam_values = fit$par*scales, converged=(fit$convergence==0), Niterations=fit$iterations, Nevaluations=fit$evaluations[[1]], Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective)
			if(diagnostics) cat(sprintf("%s  Trial %d: Final loglikelihood %.10g, Niterations %d, Nevaluations %d, converged = %d\n",verbose_prefix,trial,-results$objective_value,results$Niterations, results$Nevaluations, results$converged))
		}else{
			results = list(objective_value=NA, fparam_values = NA, converged=FALSE, Niterations=0, Nevaluations=0, Nstart_attempts=Nstart_attempts, start_values=start_values, start_objective=start_objective)
			if(diagnostics) cat(sprintf("%s  Trial %d: Start objective is non-finite. Skipping trial\n",verbose_prefix,trial))
		}
		return(results)
	}

	# run one or more independent fitting trials
    if((Ntrials>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
		# run trials in parallel using multiple forks
		# Note: Forks (and hence shared memory) are not available on Windows
		if(verbose) cat(sprintf("%sFitting %d model parameters (%d trials, parallelized)..\n",verbose_prefix,NFP,Ntrials))
		fits = parallel::mclapply(	1:Ntrials, 
									FUN = function(trial) fit_single_trial(trial), 
									mc.cores = min(Nthreads, Ntrials), 
									mc.preschedule = FALSE,
									mc.cleanup = TRUE)
	}else{
		# run in serial mode
		if(verbose) cat(sprintf("%sFitting %d model parameters (%s)..\n",verbose_prefix,NFP,(if(Ntrials==1) "1 trial" else sprintf("%d trials",Ntrials))))
		fits = sapply(1:Ntrials,function(x) NULL)
		for(trial in 1:Ntrials){
			fits[[trial]] = fit_single_trial(trial)
		}
	}	

	# extract information from best fit (note that some fits may have LL=NaN or NA)
	objective_values	= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$objective_value))
	valids				= which((!is.na(objective_values)) & (!is.nan(objective_values)) & (!is.null(objective_values)) & (!is.infinite(objective_values)))
	if(length(valids)==0) return(list(success=FALSE, error=sprintf("Fitting failed for all trials")));
	best 				= valids[which.min(sapply(valids, function(i) objective_values[i]))]
	objective_value		= -fits[[best]]$objective_value
	loglikelihood		= objective_value
	fitted_param_values = param_values; fitted_param_values[fitted_params] = fits[[best]]$fparam_values;
	if(is.null(objective_value) || any(is.na(fitted_param_values)) || any(is.nan(fitted_param_values))) return(list(success=FALSE, error=sprintf("Some fitted parameters are NaN")));

	return(list(success					= TRUE,
				objective_value			= objective_value,
				objective_name			= "loglikelihood",
				loglikelihood			= loglikelihood,
				fitted_mean				= fitted_param_values[1],
				fitted_std				= fitted_param_values[2],
				fitted_decay_rate		= fitted_param_values[3],
				linear_trend			= trend,
				NFP						= NFP,
				Ndata					= Ndata,
				AIC						= 2*NFP - 2*loglikelihood,
				BIC						= log(Ndata)*NFP - 2*loglikelihood))
}




# fit a homogenous birth-death model on a grid to a given extant timetree, choosing the "best" grid size according to AIC or BIC
fit_hbd_model_on_best_grid_size = function(	tree, 
											oldest_age			= NULL,		# numeric, the oldest age to consider in the evaluation of the likelihood as well as for defining the age grid. Typically set to the stem age or root age. Can be NULL (equivalent to the root age).
											age0				= 0,		# non-negative numeric, youngest age (time before present) to consider when fitting and with respect to which rho is defined (rho(age0) is the fraction of lineages extant at age0 that are included in the tree)
											grid_sizes			= c(1,10),	# integer vector, listing the grid sizes to consider
											uniform_grid		= TRUE,		# logical, specifying whether the age grid should be uniform (equidistant age intervals). If FALSE, then the grid point density is chosen proportional to the square root of the LTT, hence resulting in higher resolution grid near the present.
											criterion			= "AIC",	# character, how to choose the optimal grid point. Options are "AIC" or "BIC".
											exhaustive			= TRUE,		# logical, whether to try all grid sizes for choosing the "best" one. If FALSE, the grid size is gradually increased until the selectin criterio (e.g., AIC) starts becoming worse, at which point the search is halted.
											min_lambda			= 0,		# numeric, lower bound for the fitted lambdas (applying to all grid points).
											max_lambda			= +Inf,		# numeric, upper bound for the fitted lambdas (applying to all grid points).
											min_mu				= 0,		# numeric, lower bound for the fitted mus (applying to all grid points).
											max_mu				= +Inf,		# numeric, upper bound for the fitted mus (applying to all grid points).
											min_rho0			= 1e-10,	# numeric, lower bound for the fitted rho. Note that rho is always within (0,1]
											max_rho0			= 1,		# numeric, upper bound for the fitted rho.
											guess_lambda		= NULL,		# initial guess for the lambda. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant lambda at all ages), or a function handle (for generating guesses at each grid point; this function may also return NA at some time points).
											guess_mu			= NULL,		# initial guess for the mu. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant mu at all ages), or a function handle (for generating guesses at each grid point; this function may also return NA at some time points).
											guess_rho0			= 1,		# initial guess for rho. Either NULL (an initial guess will be computed automatically) or a single strictly-positive numeric.
											fixed_lambda		= NULL,		# optional fixed lambda value. Either NULL (none of the lambdas are fixed), or a single scalar (all lambdas are fixed to this value) or a function handle specifying lambda for any arbitrary age (lambdas will be fixed at any age for which this function returns a finite number). The function lambda() need not return finite values for all times.
											fixed_mu			= NULL,		# optional fixed mu value. Either NULL (none of the mus are fixed), or a single scalar (all mus are fixed to this value) or a function handle specifying mu for any arbitrary age (mu will be fixed at any age for which this function returns a finite number). The function mu() need not return finite values for all times.
											fixed_rho0			= NULL,		# optional fixed value for rho. If non-NULL and non-NA, then rho is not fitted. 
											const_lambda		= FALSE,	# logical, whether to enforce a constant (time-independent) fitted speciation rate. Only relevant if lambdas are non-fixed.
											const_mu			= FALSE,	# logical, whether to enforce a constant (time-independent) fitted extinction rate. Only relevant if mus are non-fixed.
											splines_degree		= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by lambda and mu on the age grid.
											condition			= "auto",	# one of "crown" or "stem" or "none" or "auto" or "stem2" (or "stem3" etc) or "crown3" (or "crown4" etc), specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age!=root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
											relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
											Ntrials				= 1,
											Nthreads			= 1,
											max_model_runtime	= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
											fit_control			= list(),
											verbose				= FALSE,
											verbose_prefix		= ""){
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input parameters..\n",verbose_prefix))
	root_age = get_tree_span(tree)$max_distance
	if(is.null(oldest_age)) oldest_age = root_age
	if(!is.null(guess_lambda)){
		if(class(guess_lambda) != "function"){
			if(length(guess_lambda)!=1){
				return(list(success=FALSE, error="Expecting either exactly one guess_lambda, or NULL, or a function handle"))
			}else{
				# convert guess_lambda to a function handle
				guess_lambda_value = guess_lambda
				guess_lambda = function(ages){ rep(guess_lambda_value, length(ages)) }
			}
		}
	}else{
		guess_lambda = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(guess_mu)){
		if(class(guess_mu) != "function"){
			if(length(guess_mu)!=1){
				return(list(success=FALSE, error="Expecting either exactly one guess_mu, or NULL, or a function handle"))
			}else{
				# convert guess_mu to a function handle
				guess_mu_value = guess_mu
				guess_mu = function(ages){ rep(guess_mu_value, length(ages)) }
			}
		}
	}else{
		guess_mu = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(fixed_lambda)){
		if(class(fixed_lambda) != "function"){
			if(length(fixed_lambda)!=1){
				return(list(success=FALSE, error="Expecting either exactly one fixed_lambda, or NULL, or a function handle"))
			}else{
				# convert fixed_lambda to a function handle
				fixed_lambda_value = fixed_lambda
				fixed_lambda = function(ages){ rep(fixed_lambda_value, length(ages)) }
			}
		}
	}else{
		fixed_lambda = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(fixed_mu)){
		if(class(fixed_mu) != "function"){
			if(length(fixed_mu)!=1){
				return(list(success=FALSE, error="Expecting either exactly one fixed_mu, or NULL, or a function handle"))
			}else{
				# convert fixed_mu to a function handle
				fixed_mu_value = fixed_mu
				fixed_mu = function(ages){ rep(fixed_mu_value, length(ages)) }
			}
		}
	}else{
		fixed_mu = function(ages){ rep(NA, length(ages)) }
	}
	if(length(min_lambda)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one min_lambda; instead, received %d",length(min_lambda))))
	if(length(max_lambda)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one max_lambda; instead, received %d",length(max_lambda))))
	if(!(criterion %in% c("AIC", "BIC"))) return(list(success=FALSE, error=sprintf("Invalid model selection criterion '%s'. Expected 'AIC' or 'BIC'",criterion)))
	Nmodels  = length(grid_sizes)
	
	# calculate tree LTT if needed
	if(!uniform_grid){
		LTT = count_lineages_through_time(tree=tree, Ntimes = max(100,10*max(grid_sizes)), regular_grid = TRUE, ultrametric=TRUE)
		LTT$ages = root_age - LTT$times
	}
	
	# determine order in which to examine models
	if(exhaustive){
		model_order = seq_len(Nmodels)
	}else{
		# examine models in the order of increasing grid sizes
		model_order = order(grid_sizes)
	}
	
	# fit HBD model on various grid sizes, keeping track of the "best" Ngrid
	if(verbose) cat(sprintf("%sFitting models with %s%d different grid sizes..\n",verbose_prefix,(if(exhaustive) "" else "up to "),Nmodels))
	AICs 		= rep(NA, times=Nmodels)
	BICs 		= rep(NA, times=Nmodels)
	best_fit	= NULL
	for(m in model_order){
		Ngrid = grid_sizes[m]
		if(uniform_grid || (Ngrid==1)){
			age_grid = seq(from=age0, to=oldest_age, length.out=Ngrid)
		}else{
			age_grid = get_inhomogeneous_grid_1D(Xstart = age0, Xend = oldest_age, Ngrid = Ngrid, densityX = rev(LTT$ages), densityY=sqrt(rev(LTT$lineages)), extrapolate=TRUE)
		}
		if(verbose) cat(sprintf("%s  Fitting model with grid size %d..\n",verbose_prefix,Ngrid))
		fit = fit_hbd_model_on_grid(tree				= tree, 
									oldest_age			= oldest_age,
									age0				= age0,
									age_grid			= age_grid,
									min_lambda			= min_lambda,
									max_lambda			= max_lambda,
									min_mu				= min_mu,
									max_mu				= max_mu,
									min_rho0			= min_rho0,
									max_rho0			= max_rho0,
									guess_lambda		= guess_lambda(age_grid),
									guess_mu			= guess_mu(age_grid),
									guess_rho0			= guess_rho0,
									fixed_lambda		= fixed_lambda(age_grid),
									fixed_mu			= fixed_mu(age_grid),
									fixed_rho0			= fixed_rho0,
									const_lambda		= const_lambda,
									const_mu			= const_mu,
									splines_degree		= splines_degree,
									condition			= condition,
									relative_dt			= relative_dt,
									Ntrials				= Ntrials,
									Nthreads			= Nthreads,
									max_model_runtime	= max_model_runtime,
									fit_control			= fit_control)
		if(!fit$success) return(list(success=FALSE, error=sprintf("Fitting model with grid size %d failed: %s",Ngrid,fit$error)))
		criterion_value = fit[[criterion]]
		if(is.null(best_fit)){
			best_fit = fit
			worsened = FALSE
		}else if(criterion_value<best_fit[[criterion]]){
			best_fit = fit
			worsened = FALSE
		}else{
			worsened = TRUE
		}
		AICs[m] = fit$AIC
		BICs[m] = fit$BIC
		if(verbose) cat(sprintf("%s  --> %s=%.10g. Best grid size so far: %d\n",verbose_prefix,criterion,criterion_value,length(best_fit$age_grid)))
		if((!exhaustive) && worsened) break; # model selection criterion became worse compared to the previous grid size, so stop search and keep best model found so far
	}
	
	return(list(success	 	= TRUE,
				best_fit 	= best_fit,
				grid_sizes	= grid_sizes,
				AICs		= AICs,
				BICs		= BICs))
}



# fit a homogenous birth-death model on a grid to a given extant timetree, choosing the "best" grid size according to AIC or BIC
fit_hbds_model_on_best_grid_size = function(tree, 
											root_age 			= NULL,		# optional numeric, the age of the root. Can be used to define a time offset, e.g. if the last tip was not actually sampled at the present. If NULL, this will be calculated from the tree and it will be assumed that the last tip was sampled at the present (age 0).
											oldest_age			= NULL,		# numeric, specifying the oldest age to consider. Can also be NULL, in which case this is set to the root age
											grid_sizes			= c(1,10),	# integer vector, listing the grid sizes to consider
											uniform_grid		= TRUE,		# logical, specifying whether the age grid should be uniform (equidistant age intervals). If FALSE, then the grid point density is chosen proportional to the square root of the LTT
											criterion			= "AIC",	# character, how to choose the optimal grid point. Options are "AIC" or "BIC".
											exhaustive			= TRUE,		# logical, whether to try all grid sizes for choosing the "best" one. If FALSE, the grid size is gradually increased until the selectin criterio (e.g., AIC) starts becoming worse, at which point the search is halted.
											min_lambda			= 0,		# numeric, lower bound for the fitted lambdas (applying to all grid points).
											max_lambda			= +Inf,		# numeric, upper bound for the fitted lambdas (applying to all grid points).
											min_mu				= 0,		# numeric, lower bound for the fitted mus (applying to all grid points).
											max_mu				= +Inf,		# numeric, upper bound for the fitted mus (applying to all grid points).
											min_psi				= 0,		# numeric, lower bound for the fitted psis (applying to all grid points).
											max_psi				= +Inf,		# numeric, upper bound for the fitted psis (applying to all grid points).
											min_kappa			= 0,		# numeric, lower bound for the fitted kappas (applying to all grid points).
											max_kappa			= +Inf,		# numeric, upper bound for the fitted kappas (applying to all grid points).
											guess_lambda		= NULL,		# initial guess for the lambda. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant lambda at all ages), or a function handle (for generating guesses at each grid point; this function may also return NA at some time points).
											guess_mu			= NULL,		# initial guess for the mu. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant mu at all ages), or a function handle (for generating guesses at each grid point; this function may also return NA at some time points).
											guess_psi			= NULL,		# initial guess for the psi. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant psi at all ages), or a function handle (for generating guesses at each grid point; this function may also return NA at some time points).
											guess_kappa			= NULL,		# initial guess for the kappa. Either NULL (an initial guess will be computed automatically), or a single numeric (guessing a constant kappa at all ages), or a function handle (for generating guesses at each grid point; this function may also return NA at some time points).
											fixed_lambda		= NULL,		# optional fixed lambda value. Either NULL (none of the lambdas are fixed), or a single scalar (all lambdas are fixed to this value) or a function handle specifying lambda for any arbitrary age (lambdas will be fixed at any age for which this function returns a finite number). The function lambda() need not return finite values for all times.
											fixed_mu			= NULL,		# optional fixed mu value. Either NULL (none of the mus are fixed), or a single scalar (all mus are fixed to this value) or a function handle specifying mu for any arbitrary age (mu will be fixed at any age for which this function returns a finite number). The function mu() need not return finite values for all times.
											fixed_psi			= NULL,		# optional fixed psi value. Either NULL (none of the psis are fixed), or a single scalar (all psis are fixed to this value) or a function handle specifying psi for any arbitrary age (psi will be fixed at any age for which this function returns a finite number). The function psi() need not return finite values for all times.
											fixed_kappa			= NULL,		# optional fixed kappa value. Either NULL (none of the kappas are fixed), or a single scalar (all kappas are fixed to this value) or a function handle specifying kappa for any arbitrary age (kappa will be fixed at any age for which this function returns a finite number). The function kappa() need not return finite values for all times.
											const_lambda		= FALSE,	# logical, whether to enforce a constant (time-independent) fitted speciation rate. Only relevant if lambdas are non-fixed.
											const_mu			= FALSE,	# logical, whether to enforce a constant (time-independent) fitted extinction rate. Only relevant if mus are non-fixed.
											const_psi			= FALSE,	# logical, whether to enforce a constant (time-independent) fitted sampling rate. Only relevant if psis are non-fixed.
											const_kappa			= FALSE,	# logical, whether to enforce a constant (time-independent) fitted retention probability. Only relevant if kappas are non-fixed.
											splines_degree		= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by lambda and mu on the age grid.
											condition			= "auto",	# one of "crown" or "stem" or "none" or "auto", specifying whether to condition the likelihood on the survival of the stem group or the crown group. It is recommended to use "stem" when oldest_age!=root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
											ODE_relative_dt		= 0.001,	# positive unitless number, relative integration time step for the ODE solvers. Relative to the typical time scales of the dynamics, as estimated from the theoretically maximum possible rate of change. Typical values are 0.001 - 0.1.
											ODE_relative_dy		= 1e-3,		# positive unitless mumber, unitless number, relative step for interpolating simulated values over time. So a ODE_relative_dy of 0.001 means that E is recorded and interpolated between points between which E differs by roughy 0.001. Typical values are 0.01-0.0001. A smaller E_value_step increases interpolation accuracy, but also increases memory requirements and adds runtime (scales with the tree's age span, not Ntips).
											Ntrials				= 1,
											max_start_attempts	= 1,		# integer, number of times to attempt finding a valid start point (per trial) before giving up. Randomly choosen start parameters may result in Inf/undefined objective, so this option allows the algorithm to keep looking for valid starting points.
											Nthreads			= 1,
											max_model_runtime	= NULL,		# maximum time (in seconds) to allocate for each likelihood evaluation. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
											fit_control			= list(),	# a named list containing options for the nlminb fitting routine (e.g. iter.max and rel.tol)
											verbose				= FALSE,
											verbose_prefix		= ""){
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input parameters..\n",verbose_prefix))
	if(is.null(root_age)) root_age = get_tree_span(tree)$max_distance
	if(is.null(oldest_age)) oldest_age = root_age
	if(!is.null(guess_lambda)){
		if(class(guess_lambda) != "function"){
			if(length(guess_lambda)!=1){
				return(list(success=FALSE, error="Expecting either exactly one guess_lambda, or NULL, or a function handle"))
			}else{
				# convert guess_lambda to a function handle
				guess_lambda_value = guess_lambda
				guess_lambda = function(ages){ rep(guess_lambda_value, length(ages)) }
			}
		}
	}else{
		guess_lambda = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(guess_mu)){
		if(class(guess_mu) != "function"){
			if(length(guess_mu)!=1){
				return(list(success=FALSE, error="Expecting either exactly one guess_mu, or NULL, or a function handle"))
			}else{
				# convert guess_mu to a function handle
				guess_mu_value = guess_mu
				guess_mu = function(ages){ rep(guess_mu_value, length(ages)) }
			}
		}
	}else{
		guess_mu = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(guess_psi)){
		if(class(guess_psi) != "function"){
			if(length(guess_psi)!=1){
				return(list(success=FALSE, error="Expecting either exactly one guess_psi, or NULL, or a function handle"))
			}else{
				# convert guess_psi to a function handle
				guess_psi_value = guess_psi
				guess_psi = function(ages){ rep(guess_psi_value, length(ages)) }
			}
		}
	}else{
		guess_psi = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(guess_kappa)){
		if(class(guess_kappa) != "function"){
			if(length(guess_kappa)!=1){
				return(list(success=FALSE, error="Expecting either exactly one guess_kappa, or NULL, or a function handle"))
			}else{
				# convert guess_kappa to a function handle
				guess_kappa_value = guess_kappa
				guess_kappa = function(ages){ rep(guess_kappa_value, length(ages)) }
			}
		}
	}else{
		guess_kappa = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(fixed_lambda)){
		if(class(fixed_lambda) != "function"){
			if(length(fixed_lambda)!=1){
				return(list(success=FALSE, error="Expecting either exactly one fixed_lambda, or NULL, or a function handle"))
			}else{
				# convert fixed_lambda to a function handle
				fixed_lambda_value = fixed_lambda
				fixed_lambda = function(ages){ rep(fixed_lambda_value, length(ages)) }
			}
		}
	}else{
		fixed_lambda = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(fixed_mu)){
		if(class(fixed_mu) != "function"){
			if(length(fixed_mu)!=1){
				return(list(success=FALSE, error="Expecting either exactly one fixed_mu, or NULL, or a function handle"))
			}else{
				# convert fixed_mu to a function handle
				fixed_mu_value = fixed_mu
				fixed_mu = function(ages){ rep(fixed_mu_value, length(ages)) }
			}
		}
	}else{
		fixed_mu = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(fixed_psi)){
		if(class(fixed_psi) != "function"){
			if(length(fixed_psi)!=1){
				return(list(success=FALSE, error="Expecting either exactly one fixed_psi, or NULL, or a function handle"))
			}else{
				# convert fixed_psi to a function handle
				fixed_psi_value = fixed_psi
				fixed_psi = function(ages){ rep(fixed_psi_value, length(ages)) }
			}
		}
	}else{
		fixed_psi = function(ages){ rep(NA, length(ages)) }
	}
	if(!is.null(fixed_kappa)){
		if(class(fixed_kappa) != "function"){
			if(length(fixed_kappa)!=1){
				return(list(success=FALSE, error="Expecting either exactly one fixed_kappa, or NULL, or a function handle"))
			}else{
				# convert fixed_kappa to a function handle
				fixed_kappa_value = fixed_kappa
				fixed_kappa = function(ages){ rep(fixed_kappa_value, length(ages)) }
			}
		}
	}else{
		fixed_kappa = function(ages){ rep(NA, length(ages)) }
	}
	if(length(min_lambda)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one min_lambda; instead, received %d",length(min_lambda))))
	if(length(max_lambda)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one max_lambda; instead, received %d",length(max_lambda))))
	if(length(min_mu)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one min_mu; instead, received %d",length(min_mu))))
	if(length(max_mu)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one max_mu; instead, received %d",length(max_mu))))
	if(length(min_psi)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one min_psi; instead, received %d",length(min_psi))))
	if(length(max_psi)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one max_psi; instead, received %d",length(max_psi))))
	if(length(min_kappa)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one min_kappa; instead, received %d",length(min_kappa))))
	if(length(max_kappa)!=1) return(list(success=FALSE, error=sprintf("Expecting exactly one max_kappa; instead, received %d",length(max_kappa))))
	if(!(criterion %in% c("AIC", "BIC"))) return(list(success=FALSE, error=sprintf("Invalid model selection criterion '%s'. Expected 'AIC' or 'BIC'",criterion)))
	Nmodels  = length(grid_sizes)
	
	# calculate tree LTT if needed
	if(!uniform_grid){
		LTT = count_lineages_through_time(tree=tree, Ntimes = max(100,10*max(grid_sizes)), regular_grid = TRUE)
		LTT$ages = root_age - LTT$times
	}
	
	# determine order in which to examine models
	if(exhaustive){
		model_order = seq_len(Nmodels)
	}else{
		# examine models in the order of increasing grid sizes
		model_order = order(grid_sizes)
	}
	
	# fit HBD model on various grid sizes, keeping track of the "best" Ngrid
	if(verbose) cat(sprintf("%sFitting models with %s%d different grid sizes..\n",verbose_prefix,(if(exhaustive) "" else "up to "),Nmodels))
	AICs 		= rep(NA, times=Nmodels)
	BICs 		= rep(NA, times=Nmodels)
	best_fit	= NULL
	for(m in model_order){
		Ngrid = grid_sizes[m]
		if(splines_degree==0){
			if(uniform_grid || (Ngrid==1)){
				#age_grid = c(seq(from=0,to=root_age*(Ngrid-1)/Ngrid,length.out=Ngrid),root_age*2)
				age_grid = seq(from=0,to=root_age*(Ngrid-1)/Ngrid,length.out=Ngrid)
			}else{
				#age_grid = c(get_inhomogeneous_grid_1D(Xstart=0, Xend=root_age*(Ngrid-1)/Ngrid, Ngrid=Ngrid, densityX = rev(LTT$ages), densityY=rev(sqrt(LTT$lineages)), extrapolate=TRUE),root_age*2)
				age_grid = get_inhomogeneous_grid_1D(Xstart=0, Xend=root_age*(Ngrid-1)/Ngrid, Ngrid=Ngrid, densityX = rev(LTT$ages), densityY=rev(sqrt(LTT$lineages)), extrapolate=TRUE)
			}
		}else{
			if(uniform_grid || (Ngrid==1)){
				age_grid = seq(from=0,to=root_age,length.out=Ngrid)
			}else{
				age_grid = get_inhomogeneous_grid_1D(Xstart=0, Xend=root_age, Ngrid=Ngrid, densityX = rev(LTT$ages), densityY=rev(sqrt(LTT$lineages)), extrapolate=TRUE)
			}
		}
		
		if(verbose) cat(sprintf("%s  Fitting model with grid size %d..\n",verbose_prefix,Ngrid))
		fit = fit_hbds_model_on_grid(	tree				= tree, 
										root_age			= root_age,
										oldest_age			= oldest_age,
										age_grid			= age_grid,
										min_lambda			= min_lambda,
										max_lambda			= max_lambda,
										min_mu				= min_mu,
										max_mu				= max_mu,
										min_psi				= min_psi,
										max_psi				= max_psi,
										min_kappa			= min_kappa,
										max_kappa			= max_kappa,
										guess_lambda		= guess_lambda(age_grid),
										guess_mu			= guess_mu(age_grid),
										guess_psi			= guess_psi(age_grid),
										guess_kappa			= guess_kappa(age_grid),
										fixed_lambda		= fixed_lambda(age_grid),
										fixed_mu			= fixed_mu(age_grid),
										fixed_psi			= fixed_psi(age_grid),
										fixed_kappa			= fixed_kappa(age_grid),
										const_lambda		= const_lambda,
										const_mu			= const_mu,
										const_psi			= const_psi,
										const_kappa			= const_kappa,
										splines_degree		= splines_degree,
										condition			= condition,
										ODE_relative_dt		= ODE_relative_dt,
										ODE_relative_dy		= ODE_relative_dy,
										Ntrials				= Ntrials,
										max_start_attempts	= max_start_attempts,
										Nthreads			= Nthreads,
										max_model_runtime	= max_model_runtime,
										fit_control			= fit_control,
										verbose				= FALSE,
										verbose_prefix		= paste0(verbose_prefix,"    "))
		if(!fit$success) return(list(success=FALSE, error=sprintf("Fitting model with grid size %d failed: %s",Ngrid,fit$error)))
		criterion_value = fit[[criterion]]
		if(is.null(best_fit)){
			best_fit = fit
			worsened = FALSE
		}else if(criterion_value<best_fit[[criterion]]){
			best_fit = fit
			worsened = FALSE
		}else{
			worsened = TRUE
		}
		AICs[m] = fit$AIC
		BICs[m] = fit$BIC
		if(verbose) cat(sprintf("%s  --> %s=%.10g. Best grid size so far: %d\n",verbose_prefix,criterion,criterion_value,length(best_fit$age_grid)))
		if((!exhaustive) && worsened) break; # model selection criterion became worse compared to the previous grid size, so stop search and keep best model found so far
	}
	
	return(list(success	 	= TRUE,
				best_fit 	= best_fit,
				grid_sizes	= grid_sizes,
				AICs		= AICs,
				BICs		= BICs))
}




# given a list of objects, unpack every object that is a pure list (non recursively)
# "pure list" means that class(object) = c("list")
# For example, list(a=3, b=list(c=4,d=5), e=list(f=list(6,7), g=8, h=phylo_object, i=list())) will become list(a=3, c=4, d=5, f=list(6,7), g=8, h=phylo_object)
# This function may be needed when returning a large list of lists from an Rcpp function (which may necessitate list nesting due to limitations of Rcpp) and you want to unpack the first-level lists.
flatten_list_first_level = function(L){
	N  = length(L)
	NF = sum(sapply(seq_len(N), FUN=function(i){ (if(all(class(L[[i]])=="list")) length(L[[i]]) else 1) })) # predicted length of the flattened list
	LF = vector(mode="list", NF)
	Lnames = names(L)
	LFnames = character(NF)
	f  = 1 # next index in flattened list at which to place item
	for(i in seq_len(N)){
		if(all(class(L[[i]])=="list")){
			if(length(L[[i]])>0){
				LF[f:(f+length(L[[i]])-1)] = L[[i]]
				if((!is.null(Lnames)) && (!is.null(names(L[[i]])))) LFnames[f:(f+length(L[[i]])-1)] = names(L[[i]])
			}
			f = f + length(L[[i]])
		}else{
			if(!is.null(Lnames)) LFnames[f] = Lnames[i]
			LF[[f]] = L[[i]]
			f = f + 1
		}
	}
	names(LF) = LFnames
	return(LF)
}



# given a list of birth events and death events, calculate the corresponding LTT
# Note that birth & death events need not be consistent with a tree topology, i.e. the LTT simply decreases (over increasing time) with every death event and increases with every birth event
count_lineages_through_time_BD = function(	birth_times, 	# numeric vector of arbitrary length, listing times of birth events. Need not necessarily be sorted.
											death_times, 	# numeric vector of arbitrary length, listing times of death events. Need not necessarily be sorted.
											time_grid, 		# times at which lineage counts are requested. Need not necessarily be sorted.
											initial_count){ # integer, the number of lineages assumed at -infinity
	if(length(time_grid)==0) return(numeric())
	if(any(diff(time_grid)<0)){
		# time grid is not sorted in ascending order, so figure out proper order
		time_grid_order = order(time_grid)
	}else{
		time_grid_order = c(1:length(time_grid))
	}
	lineages = numeric(length(time_grid))
	lineages[time_grid_order] = initial_count + get_LTT_BD_CPP(birth_times = sort(birth_times), death_times = sort(death_times), time_grid=time_grid[time_grid_order]);
}



# extract independent contrasts from a bifurcating tree, for fitting SBM diffusion models
get_SBM_independent_contrasts = function(	tree,
											tip_latitudes, 						# numeric vector of size Ntips, listing geographical latitudes of the tips (in decimal degrees)
											tip_longitudes, 					# numeric vector of size Ntips, listing geographical longitudes of the tips (in decimal degrees)
											radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
											clade_states			= NULL,		# optional, either an integer vector of length Ntips+Nnodes (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), specifying the discrete "state" of each tip and node in each tree. This can be used to limit independent contrasts to tip pairs whose total number of state-transitions (along their shortest path) is zero.
											planar_approximation	= FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
											only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
											only_distant_tip_pairs	= FALSE,	# logical, whether to only consider tip pairs located at distinct geographic locations
											min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
											max_MRCA_age			= Inf,		# numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set max_MRCA_age=Inf to disable this filter.
											max_phylodistance		= Inf,		# numeric, maximum allowed geodistance for an independent contrast to be included
											no_state_transitions	= FALSE,	# if TRUE, only tip pairs without state transitions along their shortest paths are considered. In particular, only tips in the same state are considered. Requires that clade_states[] is provided.
											only_state				= NULL){	# optional integer, specifying the state in which tip pairs (and their connecting ancestors) must be in order to be considered. Requires that clade_states[] is provided.
	Ntips = length(tree$tip.label)
	if(("list" %in% class(tip_latitudes)) && (length(tip_latitudes)==Ntips)){
		tip_latitudes = unlist(tip_latitudes)
	}
	if(("list" %in% class(tip_latitudes)) && (length(tip_longitudes)==Ntips)){
		tip_longitudes = unlist(tip_longitudes)
	}
	if((!is.null(clade_states)) && ("list" %in% class(clade_states)) && (length(clade_states)==Ntips+tree$Nnode)){
		clade_states = unlist(clade_states)
	}
	if((!is.null(only_state)) && is.null(clade_states)) return(list(success=FALSE, error="Missing clade_states[], needed when only_state is specified"))
	if(no_state_transitions && is.null(clade_states)) return(list(success=FALSE, error="Missing clade_states[], needed when no_state_transitions=TRUE"))
		
	# make sure tree does not have multifurcations
	# note that this does not change the tip & node indices, so the returned tip_pairs[] will still refer to the correct tips in the original tree
	if(tree_has_multifurcations_CPP(Ntips = Ntips, Nnodes = tree$Nnode, Nedges = nrow(tree$edge), tree_edge = as.vector(t(tree$edge)) - 1)){
		tree = multifurcations_to_bifurcations(tree)$tree
	}

	# extract independent pairs of sister tips
	tip_pairs = extract_independent_sister_tips(tree)
	if(only_basal_tip_pairs){
		# calculate number of nodes between tip pairs
		edge_counts = get_pairwise_distances(tree, A=tip_pairs[,1], B=tip_pairs[,2], as_edge_counts=TRUE, check_input=FALSE)
		# only keep tip pairs with at most 2 edges connecting them
		keep_pairs 	= which(edge_counts<=2)
		tip_pairs 	= tip_pairs[keep_pairs,,drop=FALSE]
	}
	
	# calculate MRCAs
	MRCAs = get_pairwise_mrcas(tree, tip_pairs[,1], tip_pairs[,2], check_input=FALSE)
	
	# calculate clade times
	clade_times = castor::get_all_distances_to_root(tree)
	
	# filter tip pairs based on MRCA time & age, if needed
	if((min_MRCA_time>0) || (max_MRCA_age<Inf)){
		tree_span 	= max(clade_times)
		keep_pairs	= which((clade_times[MRCAs]>=min_MRCA_time) & (tree_span-clade_times[MRCAs]<=max_MRCA_age))
		tip_pairs	= tip_pairs[keep_pairs,,drop=FALSE]
		MRCAs		= MRCAs[keep_pairs]
	}
	if(nrow(tip_pairs)==0) return(list(success=FALSE, error="No valid tip pairs left for extracting independent contrasts"))
		
	# calculate phylogenetic divergences and geodesic distances between sister tips
	phylodistances 	= get_pairwise_distances(tree, A=tip_pairs[,1], B=tip_pairs[,2], check_input=FALSE)
	geodistances 	= radius * sapply(1:nrow(tip_pairs), FUN=function(p){ geodesic_angle(tip_latitudes[tip_pairs[p,1]],tip_longitudes[tip_pairs[p,1]],tip_latitudes[tip_pairs[p,2]],tip_longitudes[tip_pairs[p,2]]) })

	# omit tip pairs with zero phylogenetic distance, because in that case the likelihood density is pathological
	# also omit tip pairs located at the same geographic location, if requested
	keep_pair = (phylodistances>0)
	if(only_distant_tip_pairs) keep_pair = keep_pair & (geodistances>0)
	if(max_phylodistance<Inf) keep_pair = keep_pair & (phylodistances<=max_phylodistance)
	if(no_state_transitions || (!is.null(only_state))){
		Ntransitions = count_transitions_between_clades(tree=tree, A=tip_pairs[,1], B=tip_pairs[,2], states=clade_states, check_input=TRUE)
		if(no_state_transitions) keep_pair = keep_pair & (Ntransitions==0)
		if(!is.null(only_state)) keep_pair = keep_pair & (Ntransitions==0) & (clade_states[tip_pairs[,1]]==only_state) & (clade_states[tip_pairs[,2]]==only_state)
	}
	tip_pairs		= tip_pairs[keep_pair,,drop=FALSE]
	MRCAs			= MRCAs[keep_pair]
	phylodistances 	= phylodistances[keep_pair]
	geodistances 	= geodistances[keep_pair]
	NC 				= length(phylodistances)
	if(NC==0) return(list(success=FALSE, error="No valid tip pairs left for extracting independent contrasts"))

	# determine MRCA and tip times
	MRCA_times = clade_times[MRCAs]
	tip_times1 = clade_times[tip_pairs[,1]]
	tip_times2 = clade_times[tip_pairs[,2]]

	return(list(success			= TRUE,
				NC				= NC,
				tip_pairs		= tip_pairs,		# numeric matrix of size NC x 2, listing the tip indices used to define the indepent contrasts
				phylodistances	= phylodistances,
				geodistances	= geodistances,
				MRCA_times		= MRCA_times,
				child_times1	= tip_times1,
				child_times2	= tip_times2))
}



# fit a homogenous birth-death model on a grid to a given extant timetree, choosing the "best" grid size according to AIC or BIC
fit_sbm_on_best_grid_size = function(	tree, 
										tip_latitudes, 						# numeric vector of size Ntips, listing geographical latitudes of the tips (in decimal degrees)
										tip_longitudes, 					# numeric vector of size Ntips, listing geographical longitudes of the tips (in decimal degrees)
										radius,								# numeric, radius to assume for the sphere (e.g. Earth). Use this e.g. if you want to hange the units in which diffusivity is estimated. Earth's mean radius is about 6371e3 m.
										clade_states			= NULL,		# optional, either an integer vector of length Ntips+Nnodes (if trees[] is a single tree) or a list of 1D vectors (if trees[] is a list of trees), specifying the discrete "state" of each tip and node in each tree. This can be used to limit independent contrasts to tip pairs whose total number of state-transitions (along their shortest path) is zero.
										planar_approximation	= FALSE,	# logical, specifying whether the estimation formula should be based on a planar approximation of Earth's surface, i.e. geodesic angles are converted to distances and then those are treated as if they were Euclideanon a 2D plane. This approximation substantially increases the speed of computations.
										only_basal_tip_pairs	= FALSE,	# logical, specifying whether only immediate sister tips should be considered, i.e. tip pairs with at most 2 edges between the two tips
										only_distant_tip_pairs	= FALSE,	# logical, whether to only consider tip pairs located at distinct geographic locations
										min_MRCA_time			= 0,		# numeric, specifying the minimum allowed height (distance from root) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set min_MRCA_time=0 to disable this filter.
										max_MRCA_age			= Inf,		# numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set max_MRCA_age=Inf to disable this filter.
										max_phylodistance		= Inf,		# numeric, maximum allowed geodistance for an independent contrast to be included in the SBM fitting
										no_state_transitions	= FALSE,	# if TRUE, only tip pairs without state transitions along their shortest paths are considered. In particular, only tips in the same state are considered. Requires that clade_states[] is provided.
										only_state				= NULL,		# optional integer, specifying the state in which tip pairs (and their connecting ancestors) must be in order to be considered. Requires that clade_states[] is provided.
										grid_sizes				= c(1,10),	# integer vector, listing the grid sizes to consider
										uniform_grid			= TRUE,		# logical, specifying whether the age grid should be uniform (equidistant age intervals). If FALSE, then the grid point density is chosen proportional to the square root of the independent contrasts density, hence resulting in higher resolution grid in areas where there is more data available.
										guess_diffusivity		= NULL,		# optional numeric, first guess for the diffusivity (at all grid points). If NULL, this will be automatically chosen.
										min_diffusivity			= NULL,		# numeric, lower bound for the fitted diffusivity (applying to all grid points). If NULL, this is automatically chosen.
										max_diffusivity			= Inf,		# numeric, upper bound for the fitted diffusivity (applying to all grid points).
										criterion				= "AIC",	# character, how to choose the optimal grid point. Options are "AIC" or "BIC".
										exhaustive				= TRUE,		# logical, whether to try all grid sizes for choosing the "best" one. If FALSE, the grid size is gradually increased until the selectin criterio (e.g., AIC) starts becoming worse, at which point the search is halted.
										Ntrials					= 1,
										Nthreads				= 1,
										Nbootstraps				= 0,		# (integer) optional number of parametric-bootstrap samples for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
										Ntrials_per_bootstrap	= NULL,		# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
										NQQ						= 0,		# (integer) optional number of simulations to perform for creating Q-Q plots of the theoretically expected distribution of geodistances vs the empirical distribution of geodistances (across independent contrasts). The resolution of the returned QQ plot will be equal to the number of independent contrasts used for fitting.
										fit_control				= list(),
										SBM_PD_functor			= NULL,		# internally used SBM probability density functor
										verbose					= FALSE,
										verbose_prefix			= ""){
	# basic error checking
	if(verbose) cat(sprintf("%sChecking input parameters..\n",verbose_prefix))
	if((!is.null(guess_diffusivity)) && (length(guess_diffusivity)!=1)) return(list(success=FALSE, error="Expecting either exactly one guess_diffusivity, or NULL"))									
	if((!is.null(min_diffusivity)) && (length(min_diffusivity)!=1)) return(list(success=FALSE, error="Expecting either exactly one min_diffusivity, or NULL"))									
	if((!is.null(max_diffusivity)) && (length(max_diffusivity)!=1)) return(list(success=FALSE, error="Expecting either exactly one min_diffusivity, or NULL"))									
	if(!(criterion %in% c("AIC", "BIC"))) return(list(success=FALSE, error=sprintf("Invalid model selection criterion '%s'. Expected 'AIC' or 'BIC'",criterion)))
	root_age = get_tree_span(tree)$max_distance
	Nmodels  = length(grid_sizes)
		
	# determine order in which to examine models
	if(exhaustive){
		model_order = c(1:Nmodels)
	}else{
		# examine models in the order of increasing grid sizes
		model_order = order(grid_sizes)
	}
	
	if(!uniform_grid){
		# get independent contrasts density through time, needed for defining non-uniform grid
		if(verbose) cat(sprintf("%sDetermining independent contrasts density..\n",verbose_prefix))
		ICs = get_SBM_independent_contrasts(tree					= tree,
											tip_latitudes			= tip_latitudes,
											tip_longitudes			= tip_longitudes,
											radius					= radius,
											clade_states			= clade_states,
											planar_approximation	= planar_approximation,
											only_basal_tip_pairs	= only_basal_tip_pairs,
											only_distant_tip_pairs	= only_distant_tip_pairs,
											min_MRCA_time			= min_MRCA_time,
											max_MRCA_age			= max_MRCA_age,
											max_phylodistance		= max_phylodistance,
											no_state_transitions	= no_state_transitions,
											only_state				= only_state)
		ICdensity = list(times = seq(from=root_age/1000, to=root_age, length.out=200))
		ICdensity$lineages = count_lineages_through_time_BD(birth_times		= c(ICs$MRCA_times,ICs$MRCA_times),
															death_times		= c(ICs$child_times1,ICs$child_times2),
															time_grid		= ICdensity$times, 
															initial_count	= 0)
		if(!ICs$success) return(list(success=FALSE, error=sprintf("Failed to acquire independent contrasts: %s",ICs$error)))
	}
	
	# fit SBM model on various grid sizes, keeping track of the "best" Ngrid
	if(verbose) cat(sprintf("%sFitting SBM models with %s%d different grid sizes..\n",verbose_prefix,(if(exhaustive) "" else "up to "),Nmodels))
	AICs 		= rep(NA, times=Nmodels)
	BICs 		= rep(NA, times=Nmodels)
	best_fit	= NULL
	for(m in model_order){
		Ngrid = grid_sizes[m]
		if(uniform_grid || (Ngrid==1)){
			time_grid = seq(from=0, to=root_age, length.out=Ngrid)
		}else{
			time_grid = get_inhomogeneous_grid_1D(Xstart = 0, Xend = root_age, Ngrid = Ngrid, densityX = ICdensity$times, densityY=sqrt(ICdensity$lineages), extrapolate=TRUE)
		}
		if(verbose) cat(sprintf("%s  Fitting SBM model with grid size %d..\n",verbose_prefix,Ngrid))
		fit = fit_sbm_on_grid(	tree					= tree, 
								tip_latitudes			= tip_latitudes,
								tip_longitudes			= tip_longitudes,
								radius					= radius,
								clade_states			= clade_states,
								planar_approximation	= planar_approximation,
								only_basal_tip_pairs	= only_basal_tip_pairs,
								only_distant_tip_pairs	= only_distant_tip_pairs,
								min_MRCA_time			= min_MRCA_time,
								max_MRCA_age			= max_MRCA_age,
								max_phylodistance		= max_phylodistance,
								no_state_transitions	= no_state_transitions,
								only_state				= only_state,
								time_grid				= time_grid,
								guess_diffusivity		= guess_diffusivity,
								min_diffusivity			= min_diffusivity,
								max_diffusivity			= max_diffusivity,
								Ntrials					= Ntrials,
								Nthreads				= Nthreads,
								Nbootstraps				= 0,
								Ntrials_per_bootstrap	= Ntrials_per_bootstrap,
								NQQ						= NQQ,
								fit_control				= fit_control,
								SBM_PD_functor			= (if(!is.null(best_fit)) best_fit$SBM_PD_functor else SBM_PD_functor),
								verbose					= verbose,
								verbose_prefix			= paste0(verbose_prefix,"    "))
		if(!fit$success) return(list(success=FALSE, error=sprintf("Fitting SBM model with grid size %d failed: %s",Ngrid,fit$error)))
		criterion_value = fit[[criterion]]
		if(is.null(best_fit)){
			best_fit = fit
			worsened = FALSE
		}else if(criterion_value<best_fit[[criterion]]){
			best_fit = fit
			worsened = FALSE
		}else{
			worsened = TRUE
		}
		AICs[m] = fit$AIC
		BICs[m] = fit$BIC
		if(verbose) cat(sprintf("%s  --> %s=%.10g. Best grid size so far: %d\n",verbose_prefix,criterion,criterion_value,length(best_fit$time_grid)))
		if((!exhaustive) && worsened) break; # model selection criterion became worse compared to the previous grid size, so stop search and keep best model found so far
	}
	
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sRepeating fit of model with grid size %d, for bootstrapping..\n",verbose_prefix,length(best_fit$time_grid)))
		best_fit = fit_sbm_on_grid(	tree					= tree, 
									tip_latitudes			= tip_latitudes,
									tip_longitudes			= tip_longitudes,
									radius					= radius,
									clade_states			= clade_states,
									planar_approximation	= planar_approximation,
									only_basal_tip_pairs	= only_basal_tip_pairs,
									only_distant_tip_pairs	= only_distant_tip_pairs,
									min_MRCA_time			= min_MRCA_time,
									max_MRCA_age			= max_MRCA_age,
									max_phylodistance		= max_phylodistance,
									no_state_transitions	= no_state_transitions,
									only_state				= only_state,
									time_grid				= best_fit$time_grid,
									guess_diffusivity		= guess_diffusivity,
									min_diffusivity			= min_diffusivity,
									max_diffusivity			= max_diffusivity,
									Ntrials					= Ntrials,
									Nthreads				= Nthreads,
									Nbootstraps				= Nbootstraps,
									Ntrials_per_bootstrap	= Ntrials_per_bootstrap,
									NQQ						= NQQ,
									fit_control				= fit_control,
									SBM_PD_functor			= best_fit$SBM_PD_functor,
									verbose					= verbose,
									verbose_prefix			= paste0(verbose_prefix,"  "))

	}
	
	return(list(success	 	= TRUE,
				best_fit 	= best_fit,
				grid_sizes	= grid_sizes,
				AICs		= AICs,
				BICs		= BICs))
}


read_fasta = function(	file,
						include_headers		= TRUE,
						include_sequences	= TRUE,
						truncate_headers_at	= NULL){ # optional needle string, at which to truncate headers (i.e. remove everything at and after the first instance of the needle)
	results = read_fasta_from_file_CPP(	fasta_path			= file,
										include_headers		= include_headers,
										include_sequences	= include_sequences)
	if(!results$success) return(list("success"=FALSE, error=results$error))
	if(include_headers && (!is.null(truncate_headers_at))){
		results$headers = sapply(seq_len(length(results$headers)), FUN=function(h){ strsplit(results$headers[h],split=truncate_headers_at,fixed=TRUE)[[1]][1] })
	}
	return(list(headers		= (if(include_headers) results$headers else NULL),
				sequences	= (if(include_sequences) results$sequences else NULL),
				Nlines		= results$Nlines,
				Nsequences	= results$Nsequences))
}





# make a time series monotonically increasing or decreasing, by setting problematic values to NaN
monotonize_series_by_pruning = function(values,				# numeric vector of size N, listing the scalar time series values. May include NaN.
										increasing, 		# logical, specifying whether the resulting time series should be monotonically increasing (rather than decreasing)
										prefer_later_data){	# logical, specifying whether later data (rather than earlier data) should be kept when resolving monotonicity conflicts
	results = monotonize_series_by_pruning_CPP(values = values, increasing = increasing, prefer_later_data = prefer_later_data)
	return(results)
}


# make a time series monotonically increasing or decreasing, by replacing problematic values with linearly interpolations between valid values.
monotonize_time_series_via_interpolation = function(times,				# numeric vector of size N, listing the "times" of the time series, in ascending order. May not include NaN.
													values,				# numeric vector of size N, listing the scalar time series values. May include NaN.
													increasing, 		# logical, specifying whether the resulting time series should be monotonically increasing (rather than decreasing)
													prefer_later_data){	# logical, specifying whether later data (rather than earlier data) should be kept when resolving monotonicity conflicts
	results = monotonize_series_via_interpolation_CPP(times = times, values = values, increasing = increasing, prefer_later_data = prefer_later_data)
	return(results)
}



bootstraps_to_confidence_intervals = function(bootstrap_samples){		# 2D numeric matrix of size NB x NP, where NB is the number of bootstraps and NP is the number of distinct parameters fitted at each bootstraps
	NP = ncol(bootstrap_samples)
	NB = nrow(bootstrap_samples)
	if(NP==0){
		return(list(mean		= numeric(0),
					median		= numeric(0),
					CI50lower	= numeric(0),
					CI50upper	= numeric(0),
					CI95lower	= numeric(0),
					CI95upper	= numeric(0)))
	}
	quantiles 	= sapply(seq_len(NP), FUN=function(k){ quantile(bootstrap_samples[,k], probs=c(0.5, 0.25, 0.75, 0.025, 0.975), na.rm=TRUE, type=8) })
	if(NP==1) quantiles = matrix(quantiles, ncol=NP) # make sure quantiles[,] is 2D matrix, even if only a single parameter was boostrapped
	return(list(mean 		= colMeans(bootstrap_samples, na.rm=TRUE),
				median 		= quantiles[1,],
				CI50lower 	= quantiles[2,],
				CI50upper 	= quantiles[3,],
				CI95lower 	= quantiles[4,],
				CI95upper 	= quantiles[5,]))
}




# Locally estimate the exponential growth rate of a time series in a sliding window
fit_local_exponential_rate = function(	X,					# 1D numeric vector of length N, listing X values in ascending order
										Y,					# 1D numeric vector of length N, listing non-negative Y-values corresponding to X[]
										ref						= NULL,			# optional numeric vector of length N, listing reference values for judging the adequacy of data in a sliding window.
										window_size				= 2,			# strictly positive integer, specifying the size of the sliding window (number of data points per fitting)
										trim_window_at_bounds 	= TRUE,			# logical, specifying whether to trim the sliding window when hitting the data's X-bounds. If false, then the sliding window always has the specified size, but may not always be centered around the point of evaluation, and toward the left & right bound the fitted params will be constant. If TRUE, the window becomes smaller towards the edges, hence estimates towards the edges become more noisy.
										normalize_by_Xstep		= FALSE,		# logical, specifying whether to divide each Y value by the corresponding X-step (this is only relevant if X is a non-regular grid). For example, if X are times and Y are new infections since the last time point, then the exponential growth rate of the disease should be calculated using the normalized Y.
										model 					= "lognormal",	# character, specifying the stochastic model to assume for the Y values. Available options are "lognormal" and "Poisson" (only suitable for count data).
										min_Npoints 			= 2,			# integer, specifying the minimum number of valid data points for fitting, per sliding window. If a sliding window covers fewer usable points than this threshold, the corresponding estimate is set to NA.
										min_Npositives 			= 2,			# integer, specifying the minimum number of valid positive (Y>0) data points for fitting, per sliding window. If a sliding window covers fewer usable points than this threshold, the corresponding estimate is set to NA. Not relevant for the lognormal model, since zeros are always ignored for that model.
										min_Nref				= 0,			# integer, specifying the minimum sum of reference values in a sliding window, needed for fitting. If the sum of reference values in a sliding window is below this threshold, the corresponding estimate is set to NA.
										Nbootstraps 			= 0,			# integer, specifying the optional number of boostraps for estimating confidence intervals
										Nthreads				= 1){			# integer, number of parallel threads to use for bootstrapping
	NX 				= length(X)
	window_size 	= max(1,window_size)
	min_Npoints 	= min(min_Npoints,window_size)
	min_Npositives 	= min(min_Npositives,window_size)
	
	Ymodified = Y # may be further modified below, prior to estimation of exponential rates
	if(model == "lognormal"){
		if(normalize_by_Xstep && (NX>=2)){
			# divide Y values by the X-step lengths. For Y[1] we don't know the X-step, so we can't normalize it, hence we set Y[1]=NA
			Ymodified = c(NA,Y[2:NX]/diff(X))
		}
		fit = fit_exp_LeastLogSquares_moving_window_CPP(X = X, Y = Ymodified, window_size = window_size, trim_window_at_bounds = trim_window_at_bounds)
		fit$rate[fit$Npoints<min_Npoints] = NA
		if(!is.null(ref)){
			Nrefs_per_window = sapply(seq_len(length(fit$window_starts)), FUN=function(w) sum(ref[c(fit$window_starts[w]:fit$window_ends[w])],na.rm=TRUE))
			fit$rate[Nrefs_per_window<min_Nref] = NA
		}
		valid_points_for_bootstrap = which(is.finite(fit$predicted_logY) & is.finite(fit$log_variance) & is.finite(X) & is.finite(Y) & (Y>0))
	}else if(model=="Poisson"){
		if(normalize_by_Xstep && (NX>=2)){
			# the appropriate scalings must be passed to the max-likelihood fitting routine, i.e. we can't just rescale the Y-values
			scalings = c(1, diff(X))
			Ymodified = c(NA, Y[2:NX])
		}else{
			scalings = rep(1, NX)
		}
		fit = fit_exp_Poisson_moving_window_CPP(X = X, Y = Ymodified, scalings = scalings, window_size = window_size, trim_window_at_bounds = trim_window_at_bounds)	
		fit$rate[(fit$Npoints<min_Npoints) | (fit$Npositives<min_Npositives)] = NA
		if(!is.null(ref)){
			Nrefs_per_window = sapply(seq_len(length(fit$window_starts)), FUN=function(w) sum(ref[c(fit$window_starts[w]:fit$window_ends[w])],na.rm=TRUE))
			fit$rate[Nrefs_per_window<min_Nref] = NA
		}
		valid_points_for_bootstrap = which(is.finite(X) & is.finite(Y) & (Y>=0))
	}else{
		stop(sprintf("Invalid model '%s'",model))
	}
	
	if(Nbootstraps>0){
		boostrap_rates = matrix(NA, nrow=Nbootstraps, ncol=NX)
		if(length(valid_points_for_bootstrap)>0){
			aux_bootstrap = function(b){
				# generate random time series bootstrap_Y, according to the fitted model params 
				bootstrap_Y = Y
				if(model=="lognormal"){
					# assuming a log-normal model, i.e. where log(Y) ~ normal(mean=fit$predicted_logY, variance=fit$log_variance)
					bootstrap_Y[valid_points_for_bootstrap] = exp(rnorm(n=length(valid_points_for_bootstrap), mean=fit$predicted_logY[valid_points_for_bootstrap], sd=sqrt(pmax(0,fit$log_variance[valid_points_for_bootstrap]))))
				}else if(model=="Poisson"){
					bootstrap_Y[valid_points_for_bootstrap] = rpois(n=length(valid_points_for_bootstrap), lambda=Y[valid_points_for_bootstrap])				
				}
				bfit = fit_local_exponential_rate(	X						= X,
													Y						= bootstrap_Y,
													ref						= ref,
													window_size				= window_size,
													trim_window_at_bounds	= trim_window_at_bounds,
													model					= model,
													min_Npoints				= min_Npoints,
													min_Npositives			= min_Npositives,
													min_Nref				= min_Nref,
													Nbootstraps				= 0,
													Nthreads				= 1)
				return(bfit$rate)
			}
			if((Nthreads>1) && (Nbootstraps>1) && (.Platform$OS.type!="windows")){
				boostrap_rates = matrix(unlist(parallel::mclapply(seq_len(Nbootstraps), FUN = function(b) aux_bootstrap(b), mc.cores=min(Nthreads, Nbootstraps), mc.preschedule=TRUE, mc.cleanup=TRUE)), nrow=Nbootstraps, byrow=TRUE)
			}else{
				boostrap_rates = t(sapply(seq_len(Nbootstraps), FUN = function(b) aux_bootstrap(b)))
			}
		}
		rate_CI = bootstraps_to_confidence_intervals(bootstrap_samples=boostrap_rates)
		for(i in seq_len(length(rate_CI))) rate_CI[[i]][!is.finite(fit$rate)] = NA
	}
	
	return(list(rate 	= fit$rate, 
				rate_CI = (if(Nbootstraps>0) rate_CI else NULL)))
}


# create a latitude x longitude grid, such that every grid cell has the same surface area
# returns a vector of latitudes (of size Nlat+1) and longitudes (of size Nlon), such that latitudes[1]=-90, latitudes[end]=+90, longitues[1]=-180 and longitudes[end]=+180
split_sphere_into_equisized_tiles = function(Nlat, Nlon){
	latitudes  = asin(2*seq(from=0,to=1,length.out=(Nlat+1)) - 1) * 180/pi # split sphere into Nlat latitudinal segments ("onion rings")
	longitudes = seq(from=-pi,to=pi,length.out=(Nlon+1)) * 180/pi
	return(list(latitudes=latitudes, longitudes=longitudes))
}


assign_points_to_tiles_on_sphere = function(point_latitudes, 
											point_longitudes,
											tile_latitudes,		# numeric vector of size Nlat+1, listing latitudes in ascending order from -90 up to 90.
											tile_longitudes){	# numeric vector of size Nlon+1, listing longitudes in ascending order from -180 to 180.
	Nlat = length(tile_latitudes)-1
	Nlon = length(tile_longitudes)-1
	NP	 = length(point_latitudes)
	
	# first bin by latitude
	# point2lat[p] will be the latitudinal tile index (from 1 to Nlat) for point p
	lat_order = order(point_latitudes)
	lat_binning = place_sorted_values_into_bins_CPP(items		= point_latitudes[lat_order],
													bin_mins 	= tile_latitudes[1:Nlat],
													bin_maxs	= tile_latitudes[2:(Nlat+1)])
	point2lat = integer(NP)
	point2lat[lat_order] = lat_binning$item2bin+1
	
	# next bin by longitude
	# point2lon[p] will be the longitudinal tile index (from 1 to Nlon) for point p
	lon_order = order(point_longitudes)
	lon_binning = place_sorted_values_into_bins_CPP(items		= point_longitudes[lon_order],
													bin_mins 	= tile_longitudes[1:Nlon],
													bin_maxs	= tile_longitudes[2:(Nlon+1)])
	point2lon = integer(NP)
	point2lon[lon_order] = lon_binning$item2bin+1
	
	return(list(point2lat=point2lat, point2lon=point2lon))
}



# generate random birth-death trees, simulate a constant-diffusivity SBM on each tree, and subsample each tree with geographic biases
# The geographically variable sampling probability can either be provided (see tile_counts[]) or determined based on a set of provided points on the sphere, by splitting the sphere into a number of equally sized tiles
# Note that the returned trees might have somewhat different tip counts than requested, depending on whether some simulated tips land in zero-probability tiles
simulate_geobiased_sbm = function(	Nsims,
									Ntips,
									radius,
									diffusivity,
									lambda					= 1,
									mu						= 0,
									rarefaction 			= 1,
									crown_age				= NULL,
									stem_age				= NULL,
									root_latitude			= NULL,
									root_longitude			= NULL,
									Nthreads				= 1,
									omit_failed_sims		= FALSE,
									Ntiles					= NULL,	# number of tiles into which to split sphere, for calculating sampling density. Only needed if tile_counts==NULL.
									reference_latitudes		= NULL,	# optional numeric vector of size NR, specifying the latitudes of the reference points on the sphere, i.e. based on which the sampling probability density will be determined. Only needed if tile_counts==NULL.
									reference_longitudes	= NULL,	# optional numeric vector of size NR, specifying the latitudes of the reference points on the sphere, i.e. based on which the sampling probability density will be determined. Only needed if tile_counts==NULL.
									tile_counts				= NULL,	# optional 2D numeric matrix listing reference counts (~density) for each spherical tile. The normalization of tile_counts[] does not matter, i.e. only relative values matter.
									tile_latitudes			= NULL,	# optional numeric vector of size nrow(tile_counts)+1, listing tile latitudes, ranging from -90 to +90. Must be provided iff tile_counts is provided.
									tile_longitudes			= NULL){# optional numeric vector of size ncol(tile_counts)+1, listing tile longitudes, ranging from -180 to +180. Must be provided iff tile_counts is provided.

	# determine spherical tiles
	if(!is.null(tile_counts)){
		if(is.null(tile_latitudes)) return(list(success=FALSE, error="Missing tile_latitudes"))
		if(is.null(tile_longitudes)) return(list(success=FALSE, error="Missing tile_longitudes"))
		Ntiles = nrow(tile_counts) * ncol(tile_counts)
		Nlat   = nrow(tile_counts)
		Nlon   = ncol(tile_counts)
		if(length(tile_latitudes)!=Nlat+1) return(list(success=FALSE, error=sprintf("Number of tile_latitudes (%d) differs from expectation (%d)",length(tile_latitudes),Nlat+1)))
		if(length(tile_longitudes)!=Nlon+1) return(list(success=FALSE, error=sprintf("Number of tile_longitudes (%d) differs from expectation (%d)",length(tile_longitudes),Nlon+1)))
	}else{
		if(is.null(reference_latitudes)) return(list(success=FALSE, error="Missing reference_latitudes"))
		if(is.null(reference_longitudes)) return(list(success=FALSE, error="Missing reference_longitudes"))
		if(is.null(Ntiles)) Ntiles = max(8,length(reference_latitudes)/10)
		Nlat 	= max(1,as.integer(round(sqrt(Ntiles/2))))
		Nlon 	= max(1,as.integer(round(Ntiles/Nlat)))
		tiles 	= split_sphere_into_equisized_tiles(Nlat=Nlat, Nlon=Nlon)
		tile_latitudes  = tiles$latitudes
		tile_longitudes = tiles$longitudes

		# determine number of reference points in each tile (i.e. geographic sampling density)
		tile_counts = matrix(0, nrow=Nlat, ncol=Nlon)
		for(r in seq_len(Nlat)){
			for(k in seq_len(Nlon)){
				tile_counts[r,k] = sum((reference_latitudes<tile_latitudes[r+1]) & (reference_latitudes>=tile_latitudes[r]) & (reference_longitudes<tile_longitudes[k+1]) & (reference_longitudes>=tile_longitudes[k]))
			}
		}
	}
		
	# core function for running a single simulation
	single_simulation = function(r){
		# generate random tree
		tree_sim = generate_tree_hbd_reverse(	Ntips		= Ntips/rarefaction,
												lambda		= lambda,
												mu			= mu,
												rho			= 1,
												crown_age	= crown_age,
												stem_age	= stem_age,
												relative_dt	= 0.01)
		if(!tree_sim$success) return(list(success=FALSE, error=sprintf("Simulation #%d failed: %s",r,tree_sim$error)))
		tree = tree_sim$trees[[1]]
		# simulate SBM on tree
		sbm_sim = simulate_sbm(	tree			= tree, 
								radius			= radius,
								diffusivity		= diffusivity,
								root_latitude	= root_latitude,
								root_longitude	= root_longitude)
		if(!sbm_sim$success) return(list(success=FALSE, error=sprintf("SBM simulation #%d failed: %s",r,sbm_sim$error)))
		# rarefy tree, picking tips according to the previously determined spherical probability density
		tip2tile 		= assign_points_to_tiles_on_sphere(point_latitudes = sbm_sim$tip_latitudes, point_longitudes = sbm_sim$tip_longitudes, tile_latitudes = tile_latitudes, tile_longitudes = tile_longitudes)
		tip2tile_count  = tile_counts[cbind(tip2tile$point2lat,tip2tile$point2lon)] # tip2tile_count[i] will be the tile_count value for tip i, i.e. the number of reference points found in the same tile as this tip
		if(sum(tip2tile_count>0)<2) return(list(success=FALSE, error=sprintf("SBM simulation #%d failed: Nearly all simulated tips are in empty tiles",r)))
		tip2prob 		= tip2tile_count/sum(tip2tile_count)
		keep_tips		= sample(x=seq_len(length(tip2prob)), size=min(Ntips,sum(tip2prob>0)), replace=FALSE, prob=tip2prob)
		if(length(keep_tips)<2) return(list(success=FALSE, error=sprintf("SBM simulation #%d failed: Insufficient number of tips sampled",r)))
		rarefying		 = get_subtree_with_tips(tree, only_tips=keep_tips)
		return(list(success=TRUE, tree=rarefying$subtree, latitudes=sbm_sim$tip_latitudes[rarefying$new2old_tip], longitudes=sbm_sim$tip_longitudes[rarefying$new2old_tip]))
	}

	
	# simulate multiple trees & SBMs
    if((Nsims>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
		# simulate trees in parallel using multiple forks
		# Note: Forks (and hence shared memory) are not available on Windows
		sims = parallel::mclapply(	seq_len(Nsims), 
									FUN = function(r){ single_simulation(r) }, 
									mc.cores = min(Nthreads, Nsims), 
									mc.preschedule = TRUE, 
									mc.cleanup = TRUE)
	}else{
		# run in serial mode
		sims = vector(mode="list", Nsims)
		for(r in 1:Nsims){
			sims[[r]] = single_simulation(r)
		}
	}

	# remove failed sims if requested
	if(omit_failed_sims){
		sims = sims[which(sapply(seq_len(Nsims), FUN=function(r) sims[[r]]$success))]
	}

	return(list(success			= TRUE,
				sims			= sims,
				Nlat			= Nlat,
				Nlon			= Nlon,
				tile_latitudes 	= tile_latitudes,
				tile_longitudes	= tile_longitudes,
				tile_counts		= tile_counts))
}


# remove values in the left- and right-tail of a distribution of numbers
# For example, if outlier_prob=0.5, then only the percentile 0.025 - 0.975 will be kept.
remove_outliers = function(X, outlier_prob){
	Q = quantile(X, probs=c(outlier_prob/2,1-outlier_prob/2), na.rm=TRUE, type=8)
	return(X[(X>Q[1]) & (X<Q[2])])
}



# simulate a single trajectory on a sphere according to Spherical Brownian Motion
simulate_SBM_trajectory = function(	times,						# numeric vector, listing times in ascending order, at which to evaluate the SBM trajectory
									radius, 					# numeric, radius of the sphere
									diffusivity,				# diffusivity, in distance units^2 per time
									start_latitude	= NULL,		# initial latitude decimal degrees, from -90 to 90. If NULL, it will be randomly generated
									start_longitude	= NULL){	# initial longitude decimal degrees, from -180 to 180. If NULL, it will be randomly generated
	# determine initial coordinates if needed
	if(is.null(start_latitude) && is.null(start_longitude)){
		# randomly choose root latitude & longitude
		start_longitude = runif(1, min=-180, max=180)
		start_latitude  = (180/pi) * asin(2*runif(1, min=0, max=1)-1) # randomly drawn latitude = arcsin(2*U-1), where U is uniformly distributed in [0,1]
	}else if(is.null(start_longitude)){
		# root latitude is given, so randomly choose root longitude
		start_longitude = runif(1, min=-180, max=180)
	}else if(is.null(start_latitude)){
		# root longitude is given, so randomly choose root latitude
		start_latitude = (180/pi) * asin(2*runif(1, min=0, max=1)-1) # randomly drawn latitude = arcsin(2*U-1), where U is uniformly distributed in [0,1]
	}

	# simulate trajectory on sphere
	sim = simulate_SBM_trajectory_CPP(	times		= times,
										radius		= radius,
										diffusivity	= diffusivity,
										start_theta	= pi*start_latitude/180,
										start_phi	= pi*start_longitude/180)
	return(list(success		= TRUE,
				latitudes 	= sim$thetas * 180/pi,
				longitudes	= sim$phis * 180/pi,
				distances	= sim$omegas * radius)) # geodesic distances between successive time points. So distances[t] is the geodesic distance between times[t-1] and times[t].
}


# estimate empirical cumulative distribution function from random samples
# this function can accommodate discrete as well as continuous distributions
CDF_from_samples = function(samplesX){
	samplesX 	= sort(samplesX[is.finite(samplesX)])
	N 			= length(samplesX)
	
	unique_samples	= 1 + N - rev(c(1,1+which(diff(rev(samplesX))<0))) # indices of unique values in samplesX, keeping always the last instance of each value. Assuming that samplesX is monotonically increasing.
	CDF_grid 		= samplesX[unique_samples]
	CDF_values		= unique_samples/N

	return(list(CDF_grid=CDF_grid, CDF_values=CDF_values))
}


# run a  Kolmogorov-Smirnov test to check whether an empirical set of samples was generated by some null model (or some fitted model), using parametric boostraps of that model
# For example, empirical[] may be branch lengths of some empirical tree, while bootstraps[] are branch lengths of trees simulated under a fitted birth-death model
bootstrap_Kolmogorov_Smirnov_test = function(	bootstraps,	# list of length NB, listing samples from NB bootstraps (e.g., simulations of some stochastic null model or fitted model). Hence, bootstraps[[b]] is a numeric vector, listing values generated during the b-th bootstrap.
												empirical){	# numeric vector of arbitrary size, listing samples whose deviation from the null/fitted model is to be evaluated
	NB = length(bootstraps)
	
	# calculate empirical CDF
	empirical_CDF 		 = CDF_from_samples(empirical)
	empirical_CDF_grid 	 = empirical_CDF$CDF_grid
	empirical_CDF_values = empirical_CDF$CDF_values

	# calculate bootstrap CDFs and mean bootstrap CDF (will be used as "reference"), evaluated at the same points as the empirical CDF
	bootstrap_CDF_values = matrix(NA, nrow=NB, ncol=length(empirical_CDF_grid))
	for(b in seq_len(NB)){
		NBG = length(bootstraps[[b]])
		bootstrap_CDF = CDF_from_samples(bootstraps[[b]])
		bootstrap_CDF_values[b,] = evaluate_spline(Xgrid=bootstrap_CDF$CDF_grid, Ygrid=bootstrap_CDF$CDF_values, splines_degree=0, Xtarget=empirical_CDF_grid, extrapolate="const", derivative=0)
	}
	reference_CDF_values = colMeans(bootstrap_CDF_values)	

	# calculate KS distance (i.e. the Kolmogorov-Smirnov test statistic) of empirical data from the reference_CDF
	empirical_KS = max(abs(empirical_CDF_values-reference_CDF_values))
	
	# calculate KS distance of each boostrap to the reference CDF
	boostrap_KSs = rep(NA, NB)
	for(b in seq_len(NB)){
		boostrap_KSs[b] = max(abs(bootstrap_CDF_values[b,]-reference_CDF_values))
	}

	# compare empirical_KS to bootstrap_KSs
	mean_bootstrap_KS 	= mean(boostrap_KSs)
	median_bootstrap_KS = median(boostrap_KSs)
	Pvalue 				= mean(boostrap_KSs>=empirical_KS)
	return(list(empirical_KS		= empirical_KS,
				mean_bootstrap_KS	= mean_bootstrap_KS, 
				median_bootstrap_KS	= median_bootstrap_KS, 
				Pvalue				= Pvalue,
				CDF_grid			= empirical_CDF_grid,
				reference_CDF_values= reference_CDF_values,
				empirical_CDF_values= empirical_CDF_values))
}



# calculate equal-tailed credible intervals of a collection of curves defined on a time grid of size NG (or more abstractly, of NG-dimensional numeric vectors)
# each curve is considered a "replicate", and the CIs are computed at each grid point based on all replicates
calculate_equal_tailed_CIs_of_curves = function(curves){ 	# 2D numeric matrix of size NB x NG, listing the values of NB replicate curves on NG grid points
	NB 			= nrow(curves)
	NG			= ncol(curves)
	CI50lower 	= numeric(NG)
	CI50upper 	= numeric(NG)
	CI95lower 	= numeric(NG)
	CI95upper 	= numeric(NG)
	medians 	= numeric(NG)
	means	 	= colMeans(curves, na.rm=TRUE)
	for(g in seq_len(NG)){
		quantiles = quantile(curves[,g], probs=c(0.25, 0.75, 0.025, 0.975, 0.5), na.rm=TRUE)
		CI50lower[g] = quantiles[1]
		CI50upper[g] = quantiles[2]
		CI95lower[g] = quantiles[3]
		CI95upper[g] = quantiles[4]
		medians[g] 	 = quantiles[5]	
	}
	return(list(means		= means,
				medians 	= medians,
				CI50lower 	= CI50lower,
				CI50upper 	= CI50upper,
				CI95lower	= CI95lower,
				CI95upper	= CI95upper))
}



# draw random numbers from the 4-parameter beta distribution, i.e. a beta distribution scaled and shifted to an arbitrary interval [minx, maxx]
# The parameterization used here is in terms of the mean (mu), standard deviation (sigma) and the interval bounds minx & maxx
# Note that mu & sigma must satisfy the condition sigma^2 < (mu-minx)*(maxx-mu).
rbeta4 = function(n, mu, sigma, minx=0, maxx=1){
	# determine the corresponding parameterization to use for the standard (2-parameter) beta distribution in [0,1]
	mu0 	= (mu-minx)/(maxx-minx)
	sigma0 	= sigma/(maxx-minx) # note that sigma0 and mean0 must satisfy the condition sigma0^2 < mu0*(1-mu0), which is equivalent to sigma^2 < (mu-minx)*(maxx-mu)
	nu0		= mu0*(1-mu0)/(sigma0^2) - 1
	alpha0	= mu0 * nu0
	beta0	= (1-mu0)*nu0
	# draw from the standard beta distribution, then rescale & shift
	x0 = rbeta(n=n, shape1=alpha0, shape2=beta0)
	return(minx + (maxx-minx)*x0)
}


# calculate the probability density of the 4-parameter beta distribution, i.e. a beta distribution scaled and shifted to an arbitrary interval [minx, maxx]
# The parameterization used here is in terms of the mean (mu), standard deviation (sigma) and the interval bounds minx & maxx
# Note that mu & sigma must satisfy the condition sigma^2 < (mu-minx)*(maxx-mu).
dbeta4 = function(x, mu, sigma, minx=0, maxx=1){
	# determine the corresponding parameterization to use for the standard (2-parameter) beta distribution in [0,1]
	mu0 	= (mu-minx)/(maxx-minx)
	sigma0 	= sigma/(maxx-minx) # note that sigma0 and mean0 must satisfy the condition sigma0^2 < mu0*(1-mu0), which is equivalent to sigma^2 < (mu-minx)*(maxx-mu)
	nu0		= mu0*(1-mu0)/(sigma0^2) - 1
	alpha0	= mu0 * nu0
	beta0	= (1-mu0)*nu0
	# get density of standard beta distribution
	x0	= (x-minx)/(maxx-minx)
	rho0 = dbeta(x=x0, shape1=alpha0, shape2=beta0)
	return(rho0/(maxx-minx))
}


# calculate the variance of entries in each row, thus returning a numeric vector of length nrow(A)
rowVars = function(A, d=1, na.rm=FALSE){
	means = rowMeans(A, na.rm=na.rm)
	vars  = rowSums((A-means)^2, na.rm=na.rm)/(ncol(A)-d)
	return(vars)
}


# calculate the covvariance of entries in each row between two equally-sized matrixes A & B, thus returning a numeric vector of length nrow(A) = nrow(B)
rowCovs = function(A, B, d=1, na.rm=FALSE){
	meansA = rowMeans(A, na.rm=na.rm)
	meansB = rowMeans(B, na.rm=na.rm)
	covs   = rowSums((A-meansA)*(B-meansB), na.rm=na.rm)/(ncol(A)-d)
	return(covs)
}



# a basic ABC-MCMC Metropolis-Hastings sampler for Bayesian parameter inference when the likelihood is intractable but the model itself can be efficiently simulated to generate random samples for any given choice of parameters
# ABC-MCMC: Approximate Bayesian Computation Markov Chain Monte Carlo
# Literature:
#   First introduced by: Marjoram et a. (2003). Markov chain Monte Carlo without likelihoods. PNAS. 100:15324-15328.
#   Reviewed by: Beaumont (2019). Approximate Bayesian Computation. Annual Review of Statistics and Its Application. 6:379-403. Section 4.1.
ABC_MCMC = function(observation,					# a data structure representing the original observations to which the model should be fitted, or summary statistics thereof. Must be of the same type as generated by the model.
					model,							# functional, generating data (e.g., a gene tree in the case of a Multispecies Coalescent model) according to specific input parameters, or returning low-dimensional summary statistics of generated data. Thus, model(theta) returns data or a low-dimensional summary statistic thereof. May also occasionally return NULL (e.g. model failure for some parameter values).
					metric,							# functional, returning the distance between two data points (e.g., between two gene trees), or between low-dimensional summary statistics of two data points. Hence, metric(observation, model(theta)) should return a numeric scalar.
					prior_sampler,					# functional, generating random draws from the prior parameter distribution. prior_sampler() should return a random choice of parameters. If numeric_theta==TRUE, then prior_sampler() must not generate thetas outside of the box [min_theta,max_theta].
					prior_density,					# functional, returning the prior probability density at a specific location in parameter space (need not be normalized). prior_density(theta) should be a numeric scalar.
					proposal_sampler		= "uniform",# functional, generating new proposals conditional on a specific current location in parameter space. proposal_sampler(theta_old) should return a new random choice of parameters. If numeric_theta==TRUE, then proposal_sampler() must not generate thetas outside of the box [min_theta,max_theta]. May also be "beta" or "uniform" (if numeric_theta==TRUE).
					proposal_density		= NULL,		# functional, returning the density of new proposals conditional on a specific current location in parameter space. proposal_density(theta_new, theta_old) should be a numeric scalar, corresponding to the probability density P(theta_new | theta_old). May also be NULL (if numeric_theta==TRUE), in which case a default proposal is used.
					bandwidth				= NULL,		# numeric, acceptance threshold for generated samples ("epsilon"). Either bandwidth or relative_bandwidth must be provided.
					relative_bandwidth		= NULL,		# numeric, specifying the relative acceptance threshold for generated samples (epsilon relative to the typical distances of samples drawn from the prior). Typical suitable values are 0.01 - 0.1, but this also strongly depends on your metric.
					numeric_theta			= TRUE,		# logical, specifying whether it can be assumed that the parameters theta are always a numeric of a specific size. If FALSE, theta is considered an abstract data structure and thus computation is somewhat less efficient (also, statistical analyses of theta will be omitted).
					min_theta				= -Inf,		# either a single numeric, or a numeric vector of length NP, specifying the lower bounds for the parameters theta. Only relevant if numeric_theta==TRUE. May include -Inf.
					max_theta				= +Inf,		# either a single numeric, or a numeric vector of length NP, specifying the upper bounds for the parameters theta. Only relevant if numeric_theta==TRUE. May include +Inf.
					proposal_rescaling		= 1,		# singla numeric or a numeric vector of size NP, specifying a factor for the scales (e.g., standard deviation) of the proposal. Only relevant if numeric_theta==TRUE and proposal_sampler==NULL.
					chain_length			= 10000,	# integer, the length of each MCMC chain (including burnin and prior to thinning)
					burnin					= 1000,		# integer, how many initial samples should be discarded as burnin from each MCMC chain
					thinning_step			= 1,		# integer, thinning step. A thinning step of 10 means that every 10-th sample is kept from the chain (after burnin)
					max_ACF_lag				= 1000,		# integer, maximum time lag to consider when calculating ACFs. Only relevant if numeric_theta=TRUE.
					Nchains					= 1,		# integer, number of independent MCMC chains to run. If Nthreads>1, each chain is run on a separate thread.
					Nthreads				= 1,		# integer, number of parallel threads to use for computing. Only relevant if Nchains>1.
					max_start_attempts		= 100,		# integer, maximum number of attempts to start a chain (i.e., initiate a jump)
					start_bandwidth_factor	= NULL,		# optional numeric, specifying the how much larger the initial bandwidth of each chain should be relative to the final bandwidth, during the early phase of "burnin". If NULL, this is automatically chosen.
					verbose_prefix			= "",		# character, line prefix to use before any console output
					verbose					= FALSE,	# logical, whether to print progress reports to standard output
					diagnostics				= FALSE){	# logical, whether to print technical details at each iteration. For debugging purposes mainly.
	# basic error checking and preparations
	if(chain_length<=burnin) return(list(success=FALSE, error="The chain length cannot be shorter than the burnin"))
	if(chain_length-burnin<=thinning_step) return(list(success=FALSE, error="The chain length minus the burnin cannot be shorter than the thinning_step"))
	thinning_step = max(1,thinning_step)
	adaptive_proposal = FALSE
	if(is.null(bandwidth) && is.null(relative_bandwidth)) return(list(success=FALSE, error="Either bandwidth or relative_bandwidth must be specified"))
	if((!is.null(bandwidth)) && (!is.null(relative_bandwidth))) return(list(success=FALSE, error="Only one of bandwidth or relative_bandwidth must be specified"))

	# check if prior sampler and density are functional
	if(verbose) cat(sprintf("%sChecking functionality of prior..\n",verbose_prefix))
	example_thetas = lapply(seq_len(10), FUN=function(k) prior_sampler())
	example_prior_densities = unlist(lapply(example_thetas, FUN=function(theta) prior_density(theta)))
	if(any(!is.finite(example_prior_densities))) return(list(success=FALSE, error=sprintf("%d out of %d randomly drawn prior thetas have a non-finite probability density",sum(!is.finite(example_prior_densities)),length(example_prior_densities))))

	if(numeric_theta){
		example_theta = prior_sampler()
		if(!is.numeric(example_theta)) return(list(success=FALSE, error="Example theta (drawn from the prior) is not numeric. Consider setting numeric_theta=FALSE"))
		NP = length(example_theta) # dimensionality of parameter space
		if(length(min_theta)==1) min_theta = rep(min_theta, NP)
		if(length(min_theta)!=NP) return(list(success=FALSE, error=sprintf("min_theta must either be of length 1 or %d (NP); instead, got length %d",NP,length(min_theta))))
		if(length(max_theta)==1) max_theta = rep(max_theta, NP)
		if(length(max_theta)!=NP) return(list(success=FALSE, error=sprintf("max_theta must either be of length 1 or %d (NP); instead, got length %d",NP,length(max_theta))))
		if(is.character(proposal_sampler) && (!is.null(proposal_density))) return(list(success=FALSE, error=sprintf("Since proposal sampler is set to '%s', the proposal density must be left at NULL.",proposal_sampler)))
		if((!is.character(proposal_sampler)) && is.null(proposal_density)) return(list(success=FALSE, error=sprintf("Missing proposal density")))
		if(is.character(proposal_sampler)){
			# define a reasonable proposal sampler and its corresponding probability density, accounting for potential lower & upper bounds
			if(length(proposal_rescaling)==1) proposal_rescaling = rep(proposal_rescaling,NP)
			if(length(proposal_rescaling)!=NP) return(list(success=FALSE, error=sprintf("proposal_rescaling must either be of length 1 or %d (NP); instead, got length %d",NP,length(proposal_rescaling))))
			typical_thetas 	  = sapply(seq_len(1000), FUN=function(k) prior_sampler())
			proposal_scales	  = proposal_rescaling*0.1*rowMeans(abs(typical_thetas), na.rm=TRUE)
			adaptive_proposal = TRUE
			if(proposal_sampler=="beta"){
				# Beta distribution
				if(verbose) cat(sprintf("%sSetting proposal sampler to beta distribution\n",verbose_prefix))
				proposal_sampler = function(old, scales){ rbeta4(n=NP, mu=old, sigma=pmin(scales,0.5*(old-min_theta),0.5*(max_theta-old)), minx=min_theta, maxx=max_theta) }
				proposal_density = function(theta, old, scales){ dbeta4(x=theta, mu=old, sigma=pmin(scales,0.5*(old-min_theta),0.5*(max_theta-old)), minx=min_theta, maxx=max_theta) } # note that the density must be consistent with proposal_sampler()
			}else if(proposal_sampler=="uniform"){
				# uniform distribution within a small box
				if(verbose) cat(sprintf("%sSetting proposal sampler to uniform distribution in a box\n",verbose_prefix))
				proposal_sampler = function(old, scales){ return(runif(n=NP, min=pmax(min_theta,old-scales), max=pmin(max_theta,old+scales))) }
				proposal_density = function(theta, old, scales){
					box_mins = pmax(min_theta,old-scales)
					box_maxs = pmin(max_theta,old+scales)
					return(1.0/prod(box_maxs-box_mins))
				}
			}else{
				return(list(success=FALSE, error=sprintf("Invalid proposal_sampler '%s'",proposal_sampler)))
			}
		}
	}else{
		if(is.null(proposal_sampler)) return(list(success=FALSE, error="Proposal sampler must be non-NULL when numeric_theta=FALSE"))
		if(is.null(proposal_density)) return(list(success=FALSE, error="Proposal density must be non-NULL when numeric_theta=FALSE"))
	}
	
	# determine bandwidths if needed
	if(is.null(bandwidth) || is.null(start_bandwidth_factor)){
		# determine typical distances based on prior
		if(verbose) cat(sprintf("%sDetermining typical distances from observation based on prior models..\n",verbose_prefix))
		Ncalibrations = 1000
		get_typical_distance = function(k){ theta = prior_sampler(); sam = model(theta); return(if(is.null(sam)) NA else metric(observation, sam)); }
		if((Ncalibrations>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
			typical_distances = parallel::mcmapply(seq_len(Ncalibrations), FUN=get_typical_distance, mc.cores=min(Nthreads, Ncalibrations), mc.preschedule=TRUE, mc.cleanup=TRUE, SIMPLIFY=TRUE, USE.NAMES=FALSE)
		}else{
			typical_distances = sapply(seq_len(Ncalibrations), FUN=get_typical_distance)
		}
		mean_typical_distance = mean(typical_distances, na.rm=TRUE)
		if(verbose) cat(sprintf("%s  Note: Mean prior distance = %g\n",verbose_prefix,mean_typical_distance))
	}
	if(is.null(bandwidth)){
		bandwidth = relative_bandwidth * mean_typical_distance
		if(verbose) cat(sprintf("%sSetting long-term bandwidth to %g (corresponding to quantile %g among prior distances)\n",verbose_prefix,bandwidth,mean(typical_distances<=bandwidth,na.rm=TRUE)))
	}
	if(is.null(start_bandwidth_factor) || (!is.finite(start_bandwidth_factor))){
		if(burnin<=1){
			start_bandwidth = bandwidth
		}else{
			# figure out a reasonable start bandwidth, based on the typical distances of randomly generated data from the observation
			start_bandwidth = quantile(typical_distances, probs=0.5, na.rm=TRUE) # choose initial bandwidth so that half of the typical samples would be kept
			if(verbose) cat(sprintf("%sSetting start (early burnin) bandwidth to %g (median of prior distances)\n",verbose_prefix,start_bandwidth))
		}
	}else{
		start_bandwidth = start_bandwidth_factor * bandwidth
	}	
	
	aux_single_ABC_MCMC_chain = function(n){
		if(verbose) cat(sprintf("%s  Running chain %d..\n",verbose_prefix,n))
		if(numeric_theta){
			thetas = matrix(NA, nrow=NP, ncol=floor((chain_length-burnin)/thinning_step))
		}else{
			thetas = vector(mode="list", floor((chain_length-burnin)/thinning_step)) # pre-allocate space for storing samples of this chain (after burnin and thinning)
		}
		if(numeric_theta && adaptive_proposal){
			# reserve space for storing burnin thetas whose sample ended up within the desired bandwidth. Only store up to 1000 thetas.
			burnin_good_thetas  = matrix(NA, nrow=NP, ncol=min(burnin,1000))
			Nburnin_good_thetas = 0
		}
		distances				= rep(NA,ncol(thetas))
		Naccepts				= 0
		Naccepts_after_burnin 	= 0
		Nproposals_after_burnin	= 0
		Nstart_attempts 		= 1
		next_slot 				= 1 # next slot in thetas[] for storing a sample
		theta_old 				= prior_sampler() # draw random starting point from prior
		k						= 1
		while(k<=chain_length){
			current_bandwidth = (if(k>burnin) bandwidth else (bandwidth*(k/burnin) + start_bandwidth*(1-k/burnin))) # if within the burnin phase, gradually reduce the bandwidth to its final value
			if(numeric_theta && adaptive_proposal && (k==burnin+1) && (Nburnin_good_thetas>=10)){
				# finalize proposal scales based on good thetas seen so far falling within the bandwidth
				burnin_theta_stds = sqrt(rowVars(burnin_good_thetas[,seq_len(Nburnin_good_thetas)], na.rm=TRUE))
				if(diagnostics) cat(sprintf("%s  Chain %d, end of burnin: Standard deviations of thetas within bandwidth (N=%d): %s\n",verbose_prefix,n,Nburnin_good_thetas,paste(sprintf("%.3g",burnin_theta_stds),collapse=", ")))
				proposal_scales = proposal_rescaling * burnin_theta_stds # For classical Metropolis-Hastings MCMC, [Yang and Rodriguez, 2013, PNAS. 110:19307-19312] show that for a Gaussian proposal and a Gaussian target distribution, the optimal proposal standard deviation should be 2.5 x the target standard deviation. A similar result is reported by [Thawornwattana et al. 2018. Bayesian Analysis. 13:1033-1059]
			}
			theta_new = (if(adaptive_proposal) proposal_sampler(theta_old, proposal_scales) else proposal_sampler(theta_old))
			if(k>burnin) Nproposals_after_burnin = Nproposals_after_burnin + 1
			#  decide whether to accept this proposal
			acceptance_probability = min(1, (prior_density(theta_new)/prior_density(theta_old)) * (if(adaptive_proposal) (proposal_density(theta_old, theta_new, proposal_scales)/proposal_density(theta_new, theta_old, proposal_scales)) else (proposal_density(theta_old, theta_new)/proposal_density(theta_new, theta_old))))
			accept 	 = TRUE
			distance = NA
			if(runif(n=1, min=0, max=1)>acceptance_probability){
				accept = FALSE
			}else{
				if(diagnostics) cat(sprintf("%s  Chain %d, iteration %d: simulating model for theta=%s..\n",verbose_prefix,n,k,paste(sprintf("%.3g",theta_new),collapse=", ")))
				sample_new = model(theta_new) # simulate data from the model based on the proposed parameters
				if(is.null(sample_new)){
					# the model failed to generate data for this specific choice of parameters, so don't accept it
					accept = FALSE
				}else{
					if(diagnostics) cat(sprintf("%s  Chain %d, iteration %d: Calculating distance to observation..\n",verbose_prefix,n,k))
					distance = metric(observation, sample_new)
					if(is.null(distance) || (!is.finite(distance)) || (distance>current_bandwidth)){
						accept = FALSE
					}
					# keep record of this theta if it is a "good" one (i.e., within the desired bandwidth) and if we are in the burnin phase
					if(numeric_theta && (k<=burnin) && (Nburnin_good_thetas<ncol(burnin_good_thetas)) && (!is.null(distance)) && is.finite(distance) && (distance<=bandwidth)){
						Nburnin_good_thetas = Nburnin_good_thetas + 1
						burnin_good_thetas[,Nburnin_good_thetas] = theta_new
					}
				}
			}
			if(diagnostics) cat(sprintf("%s  Chain %d, iteration %d: distance=%g, accept=%d (proposed theta=%s), Naccepts=%d\n",verbose_prefix,n,k,distance,accept,paste(sprintf("%.3g",theta_new),collapse=", "),Naccepts))
			if(accept){
				Naccepts = Naccepts + 1
				if(k>burnin) Naccepts_after_burnin = Naccepts_after_burnin + 1
			}else{
				theta_new = theta_old
				if((Naccepts==0) && (Nstart_attempts<max_start_attempts)){
					# try to start the chain again
					Nstart_attempts = Nstart_attempts + 1
					if(diagnostics) cat(sprintf("%s  Chain %d: Attempting restart (attempt %d)\n",verbose_prefix,n,Nstart_attempts))
					theta_old = prior_sampler()
					next
				}
			}
			if((k>burnin) && (((k-burnin-1) %% thinning_step) == 0)){
				# store this sample
				if(numeric_theta){
					thetas[,next_slot] = theta_new
				}else{
					thetas[[next_slot]] = theta_new
				}
				distances[next_slot] = distance
				next_slot = next_slot + 1
			}
			theta_old = theta_new
			k = k + 1
			if((k %%100) ==0) gc()
		}
		if(verbose) cat(sprintf("%s  Finished chain %d\n",verbose_prefix,n))
		return(list(success=TRUE, thetas=thetas, Naccepts=Naccepts, Naccepts_after_burnin=Naccepts_after_burnin, acceptance_rate_after_burnin=Naccepts_after_burnin/Nproposals_after_burnin, distances=distances, Nstart_attempts=Nstart_attempts))
	}
	if((Nchains>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
		if(verbose) cat(sprintf("%sRunning %d MCMC chains (parallelized)..\n",verbose_prefix,Nchains))
		chains = parallel::mclapply(seq_len(Nchains), 
									FUN 			= function(n){ aux_single_ABC_MCMC_chain(n) }, 
									mc.cores 		= min(Nthreads, Nchains), 
									mc.preschedule 	= FALSE, 
									mc.cleanup 		= TRUE)
	}else{
		if(verbose) cat(sprintf("%sRunning %d MCMC chains (sequential)..\n",verbose_prefix,Nchains))
		chains = vector(mode="list", Nchains)
		for(n in seq_len(Nchains)){
			chains[[n]] = aux_single_ABC_MCMC_chain(n)
		}
	}
	
	# discard failed or stuck chains
	valid_chains = which(sapply(seq_len(Nchains), FUN=function(n) chains[[n]]$success))
	if(length(valid_chains)<Nchains){
		if(length(valid_chains)==0) return(list(success=FALSE, error=sprintf("All MCMC chains failed: %s",chains[[1]]$error)))
		chains  = chains[valid_chains]
		Nchains = length(valid_chains)
	}
	stuck_chains = which(sapply(seq_len(Nchains), FUN=function(n) (chains[[n]]$Naccepts_after_burnin==0)))
	if(length(stuck_chains)==Nchains) return(list(success=FALSE, error=sprintf("All MCMC chains got stuck at a single value after burnin")))
	if(length(stuck_chains)>0){
		if(verbose) cat(sprintf("%sWARNING: %d out of %d MCMC chains stayed at a single value after burnin. These chains will be discarded.\n",verbose_prefix,length(stuck_chains),Nchains))
		chains  = chains[get_complement(Nchains,stuck_chains)]
		Nchains = length(chains)
	}
		
	if(numeric_theta){
		# calculate autocorrelation function (ACF) and other statistics, separately for each chain and each component of theta
		# ACFs[[n]] will be a 2D numeric matrix listing the ACF for the n-th chain, with ACFs[[n]][p,l] being the autocorrelation of the p-th component at time-lag l in chain n
		if(verbose) cat(sprintf("%sCalculating autocorrelations..\n",verbose_prefix))
		ACFs 	= vector(mode="list", Nchains)
		ESS 	= rep(NA, Nchains) # effective sample size per chain
		for(n in seq_len(Nchains)){
			thetas		= chains[[n]]$thetas
			max_lag		= min(max_ACF_lag, ncol(thetas)-2)
			ACFs[[n]] 	= matrix(NA, nrow=NP, ncol=max_lag)
			theta_mean 	= rowMeans(thetas)
			theta_var	= rowMeans((thetas-theta_mean)**2)
			for(l in seq_len(max_lag)){
				ACFs[[n]][,l] = rowMeans((thetas[,1:(ncol(thetas)-l),drop=FALSE]-theta_mean)*(thetas[,(1+l):ncol(thetas),drop=FALSE]-theta_mean),na.rm=TRUE)/theta_var
			}
			# determine "effective" number of samples, based on how much lag is needed between two samples to effectively become uncorrelated
			# we use a definition similar to that used by BEAST: https://jrnold.github.io/bayesian_notes/mcmc-diagnostics.html
			decoherence_lags = sapply(seq_len(NP), FUN=function(p) which(ACFs[[n]][p,]<=0)[1])
			if(all(is.finite(decoherence_lags))){
				ESS[n] = min(sapply(seq_len(NP), FUN=function(p) ncol(thetas)/(1 + 2*sum(ACFs[[n]][p,seq_len(max(1,decoherence_lags[p]-1))],na.rm=TRUE))))
			}else{
				# for at least one component of theta the ACF never drops to zero no matter how large the time-lag, so effectively we only have one sample (being conservative here)
				ESS[n] = 1
			}
		}
		
		# calculate the (original) Gelman-Rubin diagnostic, otherwise known as potential scale reduction factor (PSRF), denoted \hat{R}, separately for each parameter component
		# A PSRF < 1.1 is often used as a rough criterion for terminating the chain, although for critical applications the threshold should be closer to 1 (Vats and Knudson, 2020, Revisiting the Gelman-Rubin Diagnostic)
		# Because the commonly used PSRF formula (Brooks and Gelman, 1998) assumes equal-length chains, we only use the last N samples from every chain, where N is the minimum chain length available
		# The formula used here is also the same as used by the R function coda::gelman.diag.
		# Reference:
		#	Gelman and Rubin (1992). Inference from iterative simulation using multiple sequences. Statistical science. 7:457-472.
		#	Brooks and Gelman (1998). General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics. 7:434-455.
		if(verbose) cat(sprintf("%sCalculating Gelman-Rubin diagnostic R (for each model parameter)..\n",verbose_prefix))
		chain_means = matrix(NA,nrow=NP, ncol=Nchains)
		chain_vars  = matrix(NA,nrow=NP, ncol=Nchains)
		min_chain_length = min(sapply(seq_len(Nchains), FUN=function(n) ncol(chains[[n]]$thetas)))
		for(n in seq_len(Nchains)){
			N = ncol(chains[[n]]$thetas)
			chain_means[,n] = rowMeans(chains[[n]]$thetas[,(N-min_chain_length):N], na.rm=TRUE)
			chain_vars[,n]  = rowSums((chains[[n]]$thetas[,(N-min_chain_length):N]-chain_means[,n])^2, na.rm=TRUE)/(min_chain_length-1)
		}
		pooled_means = rowMeans(chain_means) # average of the chain means, separately for each parameter component
		Bn 			 = rowVars(chain_means,d=1)
		W  			 = rowMeans(chain_vars)
		sigma2 		 = ((min_chain_length-1)/min_chain_length) * W + Bn # Equation (3) in Gelman and Rubin (1992)
		Vhat		 = sigma2 + Bn/Nchains
		varVhat		 = ((min_chain_length-1)/min_chain_length)^2 * (1/Nchains) * rowVars(chain_vars,d=1) + ((Nchains+1)/Nchains)^2*(2/(Nchains-1))*(Bn^2) + (2*(Nchains+1)*(min_chain_length-1)/(Nchains*Nchains*min_chain_length)) * (rowCovs(chain_vars,chain_means^2,d=1) - 2*pooled_means*rowCovs(chain_vars,chain_means,d=1)) # Equation 4 in Gelman and Rubin (1992).
		degrees		 = 2*(Vhat^2)/varVhat
		PSRF 		 = sqrt(((degrees+3)/(degrees+1)) * Vhat/W) # \hat{R}_c on page 438 in Brooks and Gelman (1998)
		if(verbose) cat(sprintf("%s  Note: max R = %g..\n",verbose_prefix,max(PSRF)))
		
		# calculate equal-tailed credible intervals of parameters (using the pooled chains)
		if(verbose) cat(sprintf("%sCalculating posteror credible intervals..\n",verbose_prefix))
		thetas 			= sapply(seq_len(Nchains), FUN=function(n) chains[[n]]$thetas)
		theta_CI50lower = numeric(NP)
		theta_CI50upper = numeric(NP)
		theta_CI95lower = numeric(NP)
		theta_CI95upper = numeric(NP)
		theta_median	= numeric(NP)
		theta_mean 		= rowMeans(thetas, na.rm=TRUE)
		for(p in seq_len(NP)){
			quantiles 			= quantile(thetas[p,], probs=c(0.25, 0.75, 0.025, 0.975, 0.5), na.rm=TRUE)
			theta_CI50lower[p] 	= quantiles[1]
			theta_CI50upper[p] 	= quantiles[2]
			theta_CI95lower[p] 	= quantiles[3]
			theta_CI95upper[p] 	= quantiles[4]
			theta_median[p] 	= quantiles[5]	
		}
	}
	
	return(list(success					= TRUE,
				Nchains					= Nchains,
				chain_thetas			= lapply(seq_len(Nchains), FUN=function(n) chains[[n]]$thetas),
				chain_distances			= lapply(seq_len(Nchains), FUN=function(n) chains[[n]]$distances),
				acceptance_rate			= sapply(seq_len(Nchains), FUN=function(n) chains[[n]]$acceptance_rate_after_burnin),
				ESS 					= (if(numeric_theta) ESS else NULL), # numeric vector of size Nchains
				PSRF					= (if(numeric_theta) PSRF else NULL), # numeric vector of size NP, Gelman-Rubin convergence diagnostic, aka. potential scale reduction factor, per parameter.
				Nstart_attempts			= sapply(seq_len(Nchains), FUN=function(n) chains[[n]]$Nstart_attempts),
				ACF						= (if(numeric_theta) ACFs else NULL),
				theta_mean				= (if(numeric_theta) theta_mean else NULL),
				theta_median			= (if(numeric_theta) theta_median else NULL),
				theta_CI50lower			= (if(numeric_theta) theta_CI50lower else NULL),
				theta_CI50upper			= (if(numeric_theta) theta_CI50upper else NULL),
				theta_CI95lower			= (if(numeric_theta) theta_CI95lower else NULL),
				theta_CI95upper			= (if(numeric_theta) theta_CI95upper else NULL)))
	
}


# calculate the first Wasserstein distance between two sets of numbers, X & Y
# This is the 1-Wasserstein distance between two probability distributions, each consisting of a sum of Dirac distributions
# As shown by Ramdas et al. (2017), the p-Wasserstein distance between two probability distributions can be written as:
#	D(X,Y) = \int_0^1 |g(u) - f(u)| du
# where g and f are the quantile functions of the two distributions. For the special case where p=1 (first Wasserstein distance), this is equivalent to:
#	D(X,Y) = \int_{-\infty}^\infty |G(t) - F(t)| dt
# where G and F are the cumulative distribution functions (CDFs). For finite sets of numbers X and Y, generated by the two distributions, one replaces G and F by their empirical CDFs.
# This approach is used in python's function scipy.stats.wasserstein_distance
# References:
#   Ramdas et al. (2017). On Wasserstein two-sample testing and related families of nonparametric tests. Entropy. 19(2):47. Proposition 1.
first_Wasserstein_distance = function(X, Y){
	return(first_Wasserstein_distance_CPP(sort(X), sort(Y)))
}


# calculate the weighted Graph Laplacian of a phylogenetic tree (Lewitus and Morlon, 2016, Systematic Biology. 65:495-507)
# If normalized=TRUE, then the Laplacian is ensured to have eigenvalues in [0,2], see for example: http://www2.cs.cas.cz/semincm/lectures/2010-04-13-Hall.pdf
weighted_graph_Laplacian_of_tree = function(tree, normalized=FALSE, sparse=FALSE){
	Nclades = length(tree$tip.label) + tree$Nnode
	edge_length_sums = get_sum_of_edge_lengths_per_clade_CPP(	Ntips 		= length(tree$tip.label),
																Nnodes		= tree$Nnode,
																Nedges		= nrow(tree$edge),
																tree_edge	= as.vector(t(tree$edge)) - 1,
																edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length))
	bidirectional_edges   = rbind(tree$edge,tree$edge[,c(2,1)]) # duplicate list of edges, to include both directions
	bidirectional_lengths = c(tree$edge.length, tree$edge.length)
	if(normalized){
		diagonal = rep(1, Nclades)
		bidirectional_weights = -bidirectional_lengths/sqrt(edge_length_sums[bidirectional_edges[,1]]*edge_length_sums[bidirectional_edges[,2]])
	}else{
		diagonal = edge_length_sums
		bidirectional_weights = -bidirectional_lengths
	}
	if(sparse){
		# construct matrix in sparse symmetric format
		upper_triangulars = which(bidirectional_edges[,1]>=bidirectional_edges[,2])
		L = Matrix::sparseMatrix(i=c(seq_len(Nclades),bidirectional_edges[upper_triangulars,1]), j=c(seq_len(Nclades),bidirectional_edges[upper_triangulars,2]), x=c(diagonal,bidirectional_weights[upper_triangulars]), dims=c(Nclades,Nclades), symmetric=TRUE)
	}else{
		L 			 = matrix(0, ncol=Nclades, nrow=Nclades)
		diag(L) 	 = diagonal
		L[bidirectional_edges] = bidirectional_weights
	}
	return(L)
}


# Given a phylogenetic tree, extract a subset of pairwise tip distances (with replacement)
sample_pairwise_tip_distances = function(	tree, 
											Npairs,	# number of pairwise distances to consider
											as_edge_counts	= FALSE){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode
	distances = get_distances_between_clades_CPP(	Ntips 			= Ntips,
													Nnodes 			= Nnodes,
													Nedges 			= nrow(tree$edge),
													tree_edge 		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
													edge_length		= (if(as_edge_counts || is.null(tree$edge.length)) numeric() else tree$edge.length),
													cladesA			= sample.int(n=Ntips,size=Npairs, replace=TRUE)-1,
													cladesB			= sample.int(n=Ntips,size=Npairs, replace=TRUE)-1,
													verbose			= FALSE,
													verbose_prefix	= "")
	return(distances)
}




eliminate_bifurcating_root = function(tree){
	Ntips  = length(tree$tip.label)
	Nnodes = tree$Nnode
	if(Ntips<=2) return(list(success=FALSE, error="Input tree must have at least 3 tips"))
	results = eliminate_bifurcating_root_CPP(	Ntips 			= Ntips,
												Nnodes 			= Nnodes,
												Nedges 			= nrow(tree$edge),
												tree_edge 		= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
												edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length));
	if(!results$changed){
		return(list(success=TRUE, changed=FALSE, tree=tree))
	}else{
		new2old_clade = results$new2old_clade+1
		new_tree = list(Nnode 		= results$Nnodes,
						tip.label 	= (if(is.null(tree$tip.label)) NULL else tree$tip.label[new2old_clade[1:Ntips]]),
						node.label 	= (if(is.null(tree$node.label)) NULL else tree$node.label[new2old_clade[(Ntips+1):results$Nclades]-Ntips]),
						edge 		= matrix(as.integer(results$new_tree_edge),ncol=2,byrow=TRUE) + 1L,
						edge.length = results$new_edge_length)
		class(new_tree) = "phylo"
		attr(new_tree,"order") = NULL
		return(list(success=TRUE, changed=TRUE, tree=new_tree, new2old_clade=new2old_clade))
	}
}


