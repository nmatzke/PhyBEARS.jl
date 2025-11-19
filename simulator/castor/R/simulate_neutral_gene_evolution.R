# WARNING: THIS FUNCTION IS UNDER DEVELOPMENT AND NOT DEBUGGED YET. 
# IT LIKELY WORKS OK.
#
# Perform random simulation of neutral gene evolution, moving from root to tips.
# The root's state is explicitly specified at each simulation.
# Optionally, multiple independent simulations can be performed using the same model (e.g. as part of some Monte Carlo integration)
# Each site of the gene is subject to a fixed mutation rate (mutations per site per edge_length_unit).
# Mutations are assumed to be independent for each site, and happen at an exponential rate, at equal probability to each state (all-to-all).
# Optionally, distances between alleles can be used to calculate new edge lengths (edge length = number of site differences between parent & child allele)
simulate_neutral_gene_evolution = function(	tree, 
											Nsites,								# number of sites (e.g. nucleotides) at which the gene can vary neutrally
											Nstates,							# number of states each site can take (e.g. 4 for nucleotides). States are indexed 0,..,Nstates-1
											mutation_rate			= NULL,		# mutation probability rate (mutations per site per edge_length_unit). 
											include_tips			= TRUE, 
											include_nodes			= TRUE, 
											include_gene_distances	= FALSE,
											root_states 			= NULL, 	# 2D numeric matrix of size NR x Nsite (where NR can be arbitrary), specifying states for the root. If NR is smaller than Nsimulations, then values are recycled. If NULL, zero is used as root state for all sites.
											Nsimulations			= 1,
											drop_dims				= TRUE){
	Ntips  	= length(tree$tip.label);
	Nnodes  = tree$Nnode;
	Nedges	= nrow(tree$edge);
	if(is.null(root_states)) root_states = numeric()
	
	# run simulation
	results = simulate_neutral_gene_evolution_CPP(	Ntips					= Ntips,
													Nnodes					= Nnodes,
													Nedges					= Nedges,
													Nsites					= Nsites,
													Nstates					= Nstates,
													tree_edge 				= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
													edge_length				= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
													root_states				= as.vector(t(root_states)),
													mutation_rate			= mutation_rate,
													include_tips			= include_tips,
													include_nodes			= include_nodes,
													include_gene_distances	= include_gene_distances,
													Nsimulations			= Nsimulations);

	# unflatten returned arrays
	tip_states  = NULL
	node_states = NULL
	if(include_tips){
		if(drop_dims && (Nsimulations==1) && (Nsites==1)){ 
			tip_states = results$tip_states;
		}else if(drop_dims && (Nsites==1)){
			tip_states = matrix(results$tip_states, ncol=Ntips, byrow=TRUE)
		}else if(drop_dims && (Nsimulations==1)){
			tip_states = matrix(results$tip_states, ncol=Nsites, byrow=TRUE)
		}else{ 
			tip_states = aperm(array(results$tip_states,dim=c(Nsites,Ntips,Nsimulations)),c(3,2,1))
		}
	}
	if(include_nodes){
		if(drop_dims && (Nsimulations==1) && (Nsites==1)){ 
			node_states = results$node_states;
		}else if(drop_dims && (Nsites==1)){
			node_states = matrix(results$node_states, ncol=Nnodes, byrow=TRUE)
		}else if(drop_dims && (Nsimulations==1)){
			node_states = matrix(results$node_states, ncol=Nsites, byrow=TRUE)
		}else{ 
			node_states = aperm(array(results$node_states,dim=c(Nsites,Nnodes,Nsimulations)),c(3,2,1))
		}
	}
	if(include_gene_distances){
		if(drop_dims && (Nsimulations==1)){
			gene_distances = results$gene_distances;
		}else{
			gene_distances = matrix(results$gene_distances,ncol=Nedges,byrow=TRUE);
		}
	}
	return(list(tip_states		= tip_states, 		# 3D matrix of size Nsimulations x Ntips x Nsites
				node_states		= node_states,		# 3D matrix of size Nsimulations x Nnodes x Nsites
				gene_distances	= gene_distances));	# 2D matrix of size Nsimulations x Nedges
}
