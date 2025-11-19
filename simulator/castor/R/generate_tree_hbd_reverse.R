# generate a random ultrametric timetree of extant species according to the homogenous birth-death model, in backward time, based on either:
#	A provided pulled speciation rate (PSR)
#	A provided pulled diversification rate (PDR) and product rholambda0:=rho*lambda(0)
#	A provided speciation rate (lambda), extinction rate (mu) and present-day sampling fraction (rho)
# The tree will be conditioned on the number of sampled extant species at present day (i.e. the number of tips)
# Optionally, the tree can be conditioned on either:
#	The crown age
#	The stem age
# If the tree is neither conditioned on a crown age nor stem age, the PDR or PSR will be extrapolated to older ages as constants if needed, if they do not cover the entire required age interval for tree coalescence.
generate_tree_hbd_reverse = function(	Ntips,							# (integer) Number of tips in the tree
										stem_age			= NULL,		# either NULL (don't condition on the stem age) or a strictly positive numeric, specifying the stem age of the tree. If <=0, this is the same as NULL.
										crown_age			= NULL,		# either NULL (don't condition on the crown age) or a strictly positive numeric, specifying the crown age of the tree. If <=0, this is the same as NULL.
										age_grid			= NULL,		# either NULL, or empty, or a numeric vector of size NG, listing ages in ascending order, on which the PSR is specified. If NULL or empty, then the PSR must be a single scalar.
										lambda				= NULL,		# either a single numeric (constant lambda over time), or a numeric vector of size NG (listing lambda at each age in age_grid[]). Can also be NULL.
										mu					= NULL,		# either a single numeric (constant mu over time), or a numeric vector of size NG (listing mu at each age in age_grid[]). Can also be NULL.
										rho					= NULL,		# numeric, species sampling fraction at present day. Can also be NULL.
										PSR					= NULL,		# either a single numeric (constant PSR over time), or a numeric vector of size NG (listing PSR at each age in age_grid[]). Can also be NULL.
										PDR					= NULL,		# either a single numeric (constant PDR over time), or a numeric vector of size NG (listing PDR at each age in age_grid[]). Can also be NULL.
										rholambda0			= NULL,		# numeric, specifying the product rho*lambda(0). Can also be NULL, if PSR is provided.
										force_max_age		= Inf,		# numeric, optional maximum allowed age for the tree's root. If the tree ends up ranging beyond that age, all remaining lineages are forced to coalesce at that age. This is not statistically consistent; it's just a way to prevent excessively large trees if the PSR is close to zero. To disable this feature set it to Inf.
										splines_degree		= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by PSR on the age grid.
										relative_dt			= 1e-3,		# maximum relative time step allowed for integration. Smaller values increase the accuracy of the computed likelihoods, but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
										Ntrees				= 1,		# positive integer, the number of trees to generate
										tip_basename		= "",		# (character) basename for tips (e.g. "tip."). 
										node_basename		= NULL,		# (character) basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
										edge_basename		= NULL){	# (character) basename for edge (e.g. "edge."). If NULL, then edges will not have any labels.
	# basic error checking
	if(Ntrees<=0) return(list(success=TRUE, trees=list())) # nothing to be done
	if(is.null(PSR) && (is.null(PDR) || is.null(rholambda0)) && (is.null(lambda) || is.null(mu) || is.null(rho))) return(list(success = FALSE, error = sprintf("Insufficient information: Expecting either PSR, or PDR & rholambda0, or lambda & mu & rho.")))
	if(!is.null(PSR)){
		if(!is.null(PDR)) return(list(success = FALSE, error = sprintf("PDR and PSR cannot both be specified at the same time; must set at least one of them to NULL")))
		if(!is.null(rholambda0)) return(list(success = FALSE, error = sprintf("rholambda0 and PSR cannot both be specified at the same time; must set at least one of them to NULL")))
		if(!is.null(lambda)) return(list(success = FALSE, error = sprintf("lambda and PSR cannot both be specified at the same time; must set at least one of them to NULL")))
		if(!is.null(mu)) return(list(success = FALSE, error = sprintf("mu and PSR cannot both be specified at the same time; must set at least one of them to NULL")))
		if(!is.null(rho)) return(list(success = FALSE, error = sprintf("rho and PSR cannot both be specified at the same time; must set at least one of them to NULL")))
	}else if(!is.null(PDR)){
		if(is.null(rholambda0)) return(list(success = FALSE, error = sprintf("Need rholambda0 in addition to PDR")))
		if(rholambda0<=0) return(list(success = FALSE, error = sprintf("Invalid rholambda0 (%g); must be strictly positive",rholambda0)))
	}else if(!is.null(rholambda0)){
		if(is.null(PDR)) return(list(success = FALSE, error = sprintf("Need PDR in addition to rholambda0")))
	}else if(!is.null(lambda)){
		if(is.null(mu)) return(list(success = FALSE, error = sprintf("Need mu in addition to lambda")))
		if(is.null(rho)) return(list(success = FALSE, error = sprintf("Need rho in addition to lambda")))
	}else if(!is.null(mu)){
		if(is.null(lambda)) return(list(success = FALSE, error = sprintf("Need lambda in addition to mu")))
		if(is.null(rho)) return(list(success = FALSE, error = sprintf("Need rho in addition to mu")))
	}else if(!is.null(rho)){
		if(is.null(lambda)) return(list(success = FALSE, error = sprintf("Need lambda in addition to rho")))
		if(is.null(mu)) return(list(success = FALSE, error = sprintf("Need mu in addition to rho")))
	}
	if(is.null(age_grid) || (length(age_grid)<=1)){
		if((!is.null(PSR)) && (length(PSR)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of PSR values (%d); since no age grid was provided, you must either provide a single (constant) PSR or none",length(PSR))))
		if((!is.null(PDR)) && (length(PDR)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of PDR values (%d); since no age grid was provided, you must either provide a single (constant) PDR or none",length(PDR))))
		if((!is.null(lambda)) && (length(lambda)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of lambda values (%d); since no age grid was provided, you must either provide a single (constant) lambda or none",length(lambda))))
		if((!is.null(mu)) && (length(mu)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of mu values (%d); since no age grid was provided, you must either provide a single (constant) mu or none",length(mu))))
		# create dummy age grid
		NG 		 = 2;
		if((!is.null(crown_age)) || (!is.null(stem_age))){
			oldest_age = min(crown_age,stem_age)
		}else{
			# estimate how far back we might need to go, in order to cover the likeky root ages
			r = lambda-mu
			if(lambda>mu){
				oldest_age = 10*log(Ntips)/(lambda-mu)
			}else{
				oldest_age = log(Ntips)/(1000*lambda) # really a pathological case that should not be requested by the user, so this is just a wild guess
			}
		}
		age_grid = seq(from=0,to=oldest_age,length.out=NG)
		if(!is.null(PSR)) PSR = rep(PSR,times=NG);
		if(!is.null(PDR)) PSR = rep(PSR,times=NG);
	}else{
		NG = length(age_grid);
		if((!is.null(PSR)) && (length(PSR)==1)) PSR = rep(PSR,times=NG);
		if((!is.null(PDR)) && (length(PDR)==1)) PDR = rep(PDR,times=NG);
		if((!is.null(lambda)) && (length(lambda)==1)) lambda = rep(lambda,times=NG);
		if((!is.null(mu)) && (length(mu)==1)) mu = rep(mu,times=NG);
		if(any(diff(age_grid)<=0)) return(list(success = FALSE, error = sprintf("Values in age_grid must be strictly increasing")))
	}
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
	if(age_grid[1]>0) return(list(success = FALSE, error = sprintf("Age grid must cover the present-day (age 0)")))
	if((!is.null(crown_age)) && (crown_age>0) && (crown_age>tail(age_grid,1))) return(list(success = FALSE, error = sprintf("age_grid does not cover the crown age (%g)", crown_age)))
	if((!is.null(stem_age)) && (stem_age>0) && (stem_age>tail(age_grid,1))) return(list(success = FALSE, error = sprintf("age_grid does not cover the stem age (%g)", stem_age)))
	if((!is.null(crown_age)) && (crown_age>0) && (!is.null(stem_age)) && (stem_age>0)){
		if(crown_age>stem_age) return(list(success = FALSE, error = sprintf("Requested stem age (%g) is smaller than requested crown age (%g)", stem_age, crown_age)))
		stem_age = NULL; # will be conditioning on the crown age anyway, so set stem age to NULL
	}
		
	if(!is.null(lambda)){
		# calculate the PSR from lambda & mu & rho
		sim = get_PSR_of_HBD_model(	oldest_age		= min(crown_age,stem_age,tail(age_grid,1)),
									age_grid		= age_grid,
									lambda			= lambda,
									mu				= mu,
									age0			= 0,
									rho0			= rho,
									splines_degree	= splines_degree,
									relative_dt		= relative_dt)
		if(!sim$success) return(list(success = FALSE, error = sprintf("Could not calculate PSR from lambda & mu: %s",sim$error)))
		age_grid = sim$ages
		PSR 	 = sim$PSR;
		NG		 = length(age_grid)
	}else if(!is.null(PDR)){
		# calculate the PSR from the PDR & rholambda0
		sim = get_PSR_from_PDR_HBD(	oldest_age 		= min(crown_age,stem_age,tail(age_grid,1)),
									age_grid		= age_grid,
									PDR				= PDR,
									age0			= 0,
									rholambda0		= rholambda0,
									splines_degree	= splines_degree,
									relative_dt		= relative_dt,
									include_nLTT0	= FALSE);
		if(!sim$success) return(list(success = FALSE, error = sprintf("Could not calculate PSR from PDR & rholambda0: %s",sim$error)))
		age_grid = sim$ages
		PSR 	 = sim$PSR;
		NG		 = length(age_grid)
	}
		
	# generate tree based on the PSR
	results = generate_tree_from_PSR_CPP(	age_grid		= age_grid, 
											PSR				= PSR, 
											splines_degree	= splines_degree, 
											Ntips			= Ntips, 
											stem_age		= (if(is.null(stem_age)) -1 else stem_age),
											crown_age		= (if(is.null(crown_age)) -1 else crown_age),
											relative_dt		= relative_dt,
											force_max_age	= force_max_age,
											Ntrees			= Ntrees);
	if(!results$success) return(list(success = FALSE, error = results$error))

	# parse results into phylo trees
	# note that all trees will have the same number of tips, nodes and edges, but different topologies, edge lengths and potentially root indices
	Ntips	= results$Ntips
	Nnodes 	= results$Nnodes
	Nedges 	= results$Nedges
	trees = vector(Ntrees,mode="list")
	for(n in 1:Ntrees){
		trees[[n]] = list(	Nnode 		= Nnodes,
							tip.label 	= paste(tip_basename, 1:Ntips, sep=""),
							node.label 	= (if(is.null(node_basename)) NULL else paste(node_basename, 1:Nnodes, sep="")),
							edge.label 	= (if(is.null(edge_basename)) NULL else paste(edge_basename, 1:Nedges, sep="")),
							edge 		= matrix(results$tree_edge[[n]],ncol=2,byrow=TRUE) + 1,
							edge.length = results$edge_length[[n]],
							root 		= results$root[[n]]+1)
		class(trees[[n]]) = "phylo";
		attr(trees[[n]],"order") = NULL
	}
	return(list(success	= TRUE,
				trees	= trees));
}

