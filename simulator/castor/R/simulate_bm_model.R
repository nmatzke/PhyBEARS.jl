# Perform random simulation of a Brownian Motion model of continuous multivariate trait evolution, moving from root to tips.
# The root's state is explicitly specified at each simulation.
# Optionally, multiple independent simulations can be performed using the same model (e.g. as part of some Monte Carlo integration)
# The diffusivity matrix D is a non-negative definite symmetric matrix of size Ntraits x Ntraits, such that exp(-X^T*D^{-1}*X/(4*L))/sqrt(det(2*pi*D)) is the probability density for the multidimensional trait vector X after phylogenetic distance L, if initially located at the origin.
# Alternatively, the noise-amplitude matrix sigma is a Ntraits x Ndegrees matrix, such that dX = sigma * dB is the SDE of the model
simulate_bm_model = function(	tree, 
								diffusivity		= NULL,		# either a single number, or a 2D array of size Ntraits x Ntraits, diffusivity matrix of the model (in units trait_units^2/edge_length)
								sigma			= NULL,		# either a single number, or a 2D array of size Ntraits x Ndegrees, noise-amplitude matrix of the model (in units trait_units/sqrt(edge_length)). Can be provided alternatively to diffusivity.
								include_tips	= TRUE, 
								include_nodes	= TRUE, 
								root_states 	= NULL, 	# 2D numeric matrix of size NR x Ntrait (where NR can be arbitrary), specifying states for the root. If NR is smaller than Nsimulations, then values are recycled. If NULL, zero is used as root state for all traits.
								Nsimulations	= 1,
								drop_dims		= TRUE){
	Ntips  	 	= length(tree$tip.label);
	Nnodes  	= tree$Nnode;
	if(is.null(root_states)) root_states = numeric()

	# check model structure
	if(is.null(sigma) && is.null(diffusivity)) stop("ERROR: Incomplete model specification: Both sigma and diffusivity are NULL")
	if((!is.null(sigma)) && (!is.null(diffusivity))) stop("ERROR: Redundant model specification: Both sigma and diffusivity are non-NULL, but only one should be specified")
	if(!is.null(sigma)){
		if(is.vector(sigma)){
			Ntraits  = 1; 
			Ndegrees = 1;
		}else{
			Ntraits  = nrow(sigma);
			Ndegrees = ncol(sigma);
		}
		diffusivity = sigma %*% t(sigma);
	}else{
		if(is.vector(diffusivity)){
			Ntraits = 1;
		}else{
			Ntraits  = nrow(diffusivity)
			if(nrow(diffusivity)!=ncol(diffusivity)) stop(sprintf("ERROR: Diffusivity matrix must be quadratic (instead, it has dimensions %d x %d)",nrow(diffusivity),ncol(diffusivity)))
		}
		Ndegrees = Ntraits;
		sigma 	 = sqrt(2) * t(chol(diffusivity)); # Cholesky decomposition, i.e. lower-triangular matrix such that diffusivity = 0.5 * cholesky * cholesky^T
	}
	
	# run simulation
	if((Ntraits==1) && (Ndegrees==1)){
		results = simulate_scalar_Brownian_motion_model_CPP(Ntips				= Ntips,
															Nnodes				= Nnodes,
															Nedges				= nrow(tree$edge),
															tree_edge 			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
															edge_length		 	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
															root_states			= root_states,
															diffusivity			= diffusivity,
															include_tips		= include_tips,
															include_nodes		= include_nodes,
															Nsimulations		= Nsimulations);	
	}else{
		results = simulate_multivariate_Brownian_motion_model_CPP(	Ntips			= Ntips,
																	Nnodes			= Nnodes,
																	Nedges			= nrow(tree$edge),
																	Ntraits			= Ntraits,
																	Ndegrees		= Ndegrees,
																	tree_edge 		= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
																	edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
																	root_states		= as.vector(t(root_states)),
																	sigma			= as.vector(t(sigma)), 		# flatten in row-major format
																	include_tips	= include_tips,
																	include_nodes	= include_nodes,
																	Nsimulations	= Nsimulations);
	}
	# unflatten returned arrays
	tip_states  = NULL
	node_states = NULL
	if(include_tips){
		if(drop_dims && (Nsimulations==1) && (Ntraits==1)){ 
			tip_states = results$tip_states;
		}else if(drop_dims && (Ntraits==1)){
			tip_states = matrix(results$tip_states, ncol=Ntips, byrow=TRUE)
		}else if(drop_dims && (Nsimulations==1)){
			tip_states = matrix(results$tip_states, ncol=Ntraits, byrow=TRUE)
		}else{ 
			tip_states = aperm(array(results$tip_states,dim=c(Ntraits,Ntips,Nsimulations)),c(3,2,1))
		}
	}
	if(include_nodes){
		if(drop_dims && (Nsimulations==1) && (Ntraits==1)){ 
			node_states = results$node_states;
		}else if(drop_dims && (Ntraits==1)){
			node_states = matrix(results$node_states, ncol=Nnodes, byrow=TRUE)
		}else if(drop_dims && (Nsimulations==1)){
			node_states = matrix(results$node_states, ncol=Ntraits, byrow=TRUE)
		}else{ 
			node_states = aperm(array(results$node_states,dim=c(Ntraits,Nnodes,Nsimulations)),c(3,2,1))
		}
	}	
	return(list(tip_states	= tip_states, 	# 3D matrix of size Nsimulations x Ntips x Ntraits
				node_states	= node_states));# 3D matrix of size Nsimulations x Nnodes x Ntraits
}
