# simulate Spherical Brownian Motion along a rooted phylogenetic tree, from root to tips
simulate_sbm = function(	tree, 
							radius,						# sphere radius. For example Earth's mean radius is 6371 km.
							diffusivity,				# diffusion coefficient in distance_units^2/time_unit. Either a single numeric, or a numeric vector of length NT (for time-dependent diffusivities), corresponding to the times in time_grid[].
							time_grid 		= NULL,		# numeric vector of length NT, listing times (distance from root) in ascending order. Can also be NULL, in which case diffusivity[] must be a single numeric. time_grid should cover at least the interval [0,root_age], otherwise diffusivity is extrapolated as a constant if needed.
							splines_degree	= 1,		# integer, either 1 or 2 or 3, specifying the splines degree assumed for the diffusivity on the time_grid
							root_latitude	= NULL, 	# latitude of the tree root, in decimal degrees. If NULL, it will be chosen randomly.
							root_longitude	= NULL){	# longitude of the tree root, in decimal degrees. If NULL, it will be chosen randomly.
	Ntips	= length(tree$tip.label)					
	Nnodes	= tree$Nnode
	
	# basic error checking
	if(!("numeric" %in% class(diffusivity))) return(list(success=FALSE, error="Expected numeric vector for diffusivity"))
	diffusivity = pmax(diffusivity,0);
	if(is.null(time_grid) || (length(time_grid)<=1)){
		if(length(diffusivity)!=1){
			return(list(success=FALSE, error=sprintf("Expecting exactly one diffusivity, since time_grid was NULL; instead, got %d diffusivities",length(diffusivity))))
		}else{
			NT = 1
		}
	}else if(!is.null(time_grid)){
		NT = length(time_grid)
		if(length(diffusivity)==1){
			diffusivity = rep(diffusivity, NT)
		}else if(length(diffusivity)!=NT){
			return(list(success=FALSE, error=sprintf("Expecting either one diffusivity or exactly %d diffusivities, since time_grid had length %d; instead, got %d diffusivities",NT,NT,length(diffusivity))))
		}
		if(time_grid[1]>time_grid[NT]) return(list(success=FALSE, error=sprintf("time_grid must be in ascending order")))
	}
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
	
	# determine root coordinates if needed
	if(is.null(root_latitude) && is.null(root_longitude)){
		# randomly choose root latitude & longitude
		root_longitude = runif(1, min=-180, max=180)
		root_latitude  = (180/pi) * asin(2*runif(1, min=0, max=1)-1) # randomly drawn latitude = arcsin(2*U-1), where U is uniformly distributed in [0,1]
	}else if(is.null(root_longitude)){
		# root latitude is given, so randomly choose root longitude
		root_longitude = runif(1, min=-180, max=180)
	}else if(is.null(root_latitude)){
		# root longitude is given, so randomly choose root latitude
		root_latitude = (180/pi) * asin(2*runif(1, min=0, max=1)-1) # randomly drawn latitude = arcsin(2*U-1), where U is uniformly distributed in [0,1]
	}
	
	if(NT==1){
		# simulate time-independent SBM
		results = simulate_SBM_on_tree_CPP(	Ntips		= Ntips,
											Nnodes		= Nnodes,
											Nedges		= nrow(tree$edge),
											tree_edge	= as.vector(t(tree$edge))-1,
											edge_length	= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
											radius		= radius,
											diffusivity	= diffusivity,
											root_theta	= pi*root_latitude/180,
											root_phi	= pi*root_longitude/180)
	}else{
		# simulate time-variable SBM
		results = simulate_TSBM_on_tree_CPP(Ntips			= Ntips,
											Nnodes			= Nnodes,
											Nedges			= nrow(tree$edge),
											tree_edge		= as.vector(t(tree$edge))-1,
											edge_length		= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
											radius			= radius,
											time_grid		= time_grid,
											diffusivities	= diffusivity,
											splines_degree	= splines_degree,
											root_theta		= pi*root_latitude/180,
											root_phi		= pi*root_longitude/180)
		
	}
	return(list(success			= TRUE,
				tip_latitudes 	= 180 * results$clade_theta[1:Ntips]/pi,
				tip_longitudes	= 180 * results$clade_phi[1:Ntips]/pi,
				node_latitudes	= 180 * results$clade_theta[(Ntips+1):(Ntips+Nnodes)]/pi,
				node_longitudes	= 180 * results$clade_phi[(Ntips+1):(Ntips+Nnodes)]/pi));
}
