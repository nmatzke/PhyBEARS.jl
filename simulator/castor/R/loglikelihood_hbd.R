# Calculate the loglikelihood of a homogenous speciation-extinction (HSE) cladogenic model, for a given ultrametric timetree and given one of the two combinations of input variables:
#	1. speciation/extinction rates over time (lambda & mu), and rho of the tree (rho)
#	2. pulled diversification rate over time (PDR), and the product rho*lambda0, where rho is the rho and lambda0 is the present-day speciation rate
#	3. pulled speciation rate over time (PSR)
# The speciation rate lambda, extinction rate mu, PDR and PSR must be specified on a discrete age-grid, and are assumed to vary linearly (or polynomially, as splines) between grid points (see "degree" argument).
#
# References:
#	Morlon et al. (2011). Reconciling molecular phylogenies with the fossil record. PNAS 108:16327-16332
#
# Tree requirements:
#   Tree can include multi- and mono-furcations.
#   Tree must be rooted. Root will be determined automatically as the node with no parent.
#   Tree must be ultrametric (e.g. a timetree of extant species). In particular, all tips are assumed to have age 0.
loglikelihood_hbd = function(	tree, 
								oldest_age			= NULL,			# either a numeric specifying the oldest age to consider or NULL (equivalent to the root age). This is similar to the "tot_time" option in the R function RPANDA::likelihood_bd. If oldest_age>root_age, this is assumed to be the stem age. If oldest_age<root_age, the tree is "split" into multiple subtrees at that age, and each subtree is considered an independent realization of the HBD model stemming at that age.
								age0				= 0,			# non-negative numeric, youngest age (time before present) to consider when cancluating the loglikelihood and with respect to which rholambda0 is defined (i.e. rholambda0 = rho(age0)*lambda(age0))
								rho0				= NULL,			# numeric within (0,1], specifying the fraction of extant diversity at age0 that is represented in the tree.
								rholambda0			= NULL,			# either NULL or numeric, specifying the product between the sampling fraction at age0 and the speciation rate at age0
								age_grid			= NULL,			# either NULL, or empty, or a numeric vector of size NG, listing ages in ascending order, on which birth/mu are specified. If NULL or empty, then lambda and mu must be a single scalar.
								lambda				= NULL,			# either NULL, or a single scalar (constant speciation rate over time), or a numeric vector of size NG (listing speciation rates at each age in grid_ages[]).
								mu					= NULL,			# either NULL, or a single scalar (constant extinction rate over time), or a numeric vector of size NG (listing extinction rates at each age in grid_ages[]).
								PDR					= NULL,			# either NULL, or a single scalar (constant PDR over time), or a numeric vector of size NG (listing PDR at each age in grid_ages[]).
								PSR					= NULL,			# either NULL, or a single scalar (constant PSR over time), or a numeric vector of size NG (listing PSR at each age in grid_ages[]).
								splines_degree		= 1,			# integer, either 1 or 2 or 3, specifying the degree for the splines defined by lambda, mu and PDR on the age grid.
								condition			= "auto",		# one of "crown", "stem", "auto" or "none" (or FALSE), specifying whether to condition the likelihood on the survival of the stem group, the crown group or none (not recommended, and only available when lambda/mu are provided). It is recommended to use "stem" when oldest_age>root_age, and "crown" when oldest_age==root_age. This argument is similar to the "cond" argument in the R function RPANDA::likelihood_bd. Note that "crown" really only makes sense when oldest_age==root_age.
								max_model_runtime	= -1,			# maximum time (in seconds) to allocate for the likelihood evaluation. If negative, this option is ignored.
								relative_dt			= 1e-3){		# maximum relative time step allowed for integration. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
	# basic input error checking
	if(tree$Nnode<2) return(list(success = FALSE, error="Input tree is too small"));
	if(age0<0) return(list(success = FALSE, error="age0 must be non-negative"));

	# trim tree at age0 if needed, while shifting time for the subsequent analyses (i.e. new ages will start counting at age0)
	if(age0>0){
		root_age = get_tree_span(tree)$max_distance
		if(root_age<age0) return(list(success=FALSE, error=sprintf("age0 is older than the root age (%g)",root_age)));
		tree = trim_tree_at_height(tree,height=root_age-age0)$tree
		if(tree$Nnode<2) return(list(success = FALSE, error=sprintf("Tree is too small after trimming at age0 (%g)",age0)));
		if(!is.null(oldest_age)) oldest_age	= oldest_age - age0	
		if(!is.null(age_grid)) age_grid 	= age_grid - age0
	}

	PDR_based = (!is.null(PDR))
	PSR_based = (!is.null(PSR))
	
	# get branching ages (=node ages) in ascending order
	# branching ages must be in ascending order when provided to the C++ routines below
	sorted_node_ages	= sort(get_all_branching_ages(tree));
	root_age 		 	= tail(sorted_node_ages,1);

	# check validity of input variables
	if(PDR_based){
		if(!is.null(lambda)) return(list(success = FALSE, error = sprintf("lambda must be NULL when PDR is provided")))
		if(!is.null(mu)) return(list(success = FALSE, error = sprintf("mu must be NULL when PDR is provided")))
		if(!is.null(rho0)) return(list(success = FALSE, error = sprintf("rho0 must be NULL when PDR is provided")))
		if(is.null(rholambda0)) return(list(success = FALSE, error = sprintf("rholambda0 must be non-NULL when PDR is provided")))
		if(!is.null(PSR)) return(list(success = FALSE, error = sprintf("PSR must be NULL when PDR is provided")))
	}else if(PSR_based){
		if(!is.null(lambda)) return(list(success = FALSE, error = sprintf("lambda must be NULL when PSR is provided")))
		if(!is.null(mu)) return(list(success = FALSE, error = sprintf("mu must be NULL when PSR is provided")))
		if(!is.null(rho0)) return(list(success = FALSE, error = sprintf("rho0 must be NULL when PSR is provided")))
		if(!is.null(rholambda0)) return(list(success = FALSE, error = sprintf("rholambda0 must be NULL when PSR is provided")))
		if(!is.null(PDR)) return(list(success = FALSE, error = sprintf("PDR must be NULL when PSR is provided")))
	}else{
		if(is.null(lambda) || is.null(mu)) return(list(success = FALSE, error = sprintf("Either lambda/mu, or PDR, must be NULL, but not both")))
		if(!is.null(rholambda0)) return(list(success = FALSE, error = sprintf("rholambda0 must be NULL when lambda and mu are provided")))
		if(is.null(rho0)) return(list(success = FALSE, error = sprintf("rho must be non-NULL when lambda and mu are provided")))
	}
	if(condition==FALSE) condition = "none"
	if(PDR_based || PSR_based){
		if(!(condition %in% c("stem","crown","auto"))) return(list(success = FALSE, error = sprintf("Invalid condition option '%s'; expected either 'crown', 'stem' or 'auto'",condition)))
	}else{
		if(!(condition %in% c("stem","crown","auto","none"))) return(list(success = FALSE, error = sprintf("Invalid condition option '%s'; expected either 'crown', 'stem', 'auto' or 'none' (not recommended)",condition)))
	}
	if(is.null(oldest_age)) oldest_age = root_age;
	if(condition=="auto") condition = (if(abs(oldest_age-root_age)<=1e-10*root_age) "crown" else "stem")
	if(is.null(age_grid) || (length(age_grid)==0) || (length(age_grid)==1)){
		constant_rates = TRUE
		if(PDR_based){
			if(length(PDR)!=1) return(list(success = FALSE, error = sprintf("Invalid number of PDR; since no non-trivial age grid was provided, you must provide a single (constant) PDR")))
		}else if(PSR_based){
			if(length(PSR)!=1) return(list(success = FALSE, error = sprintf("Invalid number of PSR; since no non-trivial age grid was provided, you must provide a single (constant) PDR")))
		}else{
			if(length(lambda)!=1) return(list(success = FALSE, error = sprintf("Invalid number of lambda; since no non-trivial age grid was provided, you must provide a single (constant) lambda")))
			if(length(mu)!=1) return(list(success = FALSE, error = sprintf("Invalid number of mu; since no non-trivial age grid was provided, you must provide a single (constant) mu")))
		}
		# create dummy age grid
		NG		 = 2;
		age_grid = seq(from=0,to=oldest_age,length.out=NG)
		if(PDR_based){
			PDR = rep(PDR,times=NG);
		}else if(PSR_based){
			PSR = rep(PSR,times=NG);
		}else{
			lambda = rep(lambda,times=NG);
			mu = rep(mu,times=NG);
		}
	}else{
		constant_rates = FALSE
		NG = length(age_grid);
		if((age_grid[1]>0) || (age_grid[NG]<oldest_age)) return(list(success = FALSE, error = sprintf("Age grid must cover all ages from age0 (%g) until oldest_age (%g)",age0,oldest_age+age0)))
		if(PDR_based){
			if((length(PDR)!=1) && (length(PDR)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of PDR (%d); since an age grid of size %d was provided, you must either provide one or %d PDR",length(PDR),NG,NG)));
			if(length(PDR)==1) PDR = rep(PDR,times=NG);
		}else if(PSR_based){
			if((length(PSR)!=1) && (length(PSR)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of PSR (%d); since an age grid of size %d was provided, you must either provide one or %d PSR",length(PSR),NG,NG)));
			if(length(PSR)==1) PSR = rep(PSR,times=NG);
		}else{
			if((length(lambda)!=1) && (length(lambda)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of lambda (%d); since an age grid of size %d was provided, you must either provide one or %d lambda",length(lambda),NG,NG)));
			if((length(mu)!=1) && (length(mu)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of mu (%d); since an age grid of size %d was provided, you must either provide one or %d mu",length(mu),NG,NG)));
			if(length(lambda)==1) lambda = rep(lambda,times=NG);
			if(length(mu)==1) mu = rep(mu,times=NG);
			constant_rates = (length(unique(lambda))==1) && (length(unique(mu))==1)
		}
	}
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree: Extected one of 0,1,2,3.")));
		
	# calculate log-likelihood
	if(PDR_based){
		results = HBD_PDR_loglikelihood_CPP(branching_ages		= sorted_node_ages,
											oldest_age			= oldest_age,
											rholambda0 			= rholambda0,
											age_grid 			= age_grid,
											PDRs 				= PDR,
											splines_degree		= splines_degree,
											condition			= condition,
											relative_dt			= relative_dt,
											runtime_out_seconds	= max_model_runtime,
											diff_PDR			= numeric(),
											diff_PDR_degree		= 0);
	}else if(PSR_based){
		results = HBD_PSR_loglikelihood_CPP(branching_ages		= sorted_node_ages,
											oldest_age			= oldest_age,
											age_grid 			= age_grid,
											PSRs 				= PSR,
											splines_degree		= splines_degree,
											condition			= condition,
											relative_dt			= relative_dt,
											runtime_out_seconds	= max_model_runtime);
	}else{
		if(constant_rates){
			# This is a constant-rates model (i.e., lambda & mu are constant over time), so use more specialized (efficient) loglikelihood routine
			results = CR_HBD_model_loglikelihood_CPP(	branching_ages		= sorted_node_ages,
														oldest_age			= oldest_age,
														rarefaction			= rho0,
														lambda	 			= lambda[1],
														mu 					= mu[1],
														condition			= condition);
		}else{
			results = HBD_model_loglikelihood_CPP(	branching_ages		= sorted_node_ages,
													oldest_age			= oldest_age,
													rarefaction			= rho0,
													age_grid 			= age_grid,
													lambdas 			= lambda,
													mus					= mu,
													splines_degree		= splines_degree,
													condition			= condition,
													relative_dt			= relative_dt,
													runtime_out_seconds	= max_model_runtime);
		}
	}
	if(!results$success) return(list(success = FALSE, error = sprintf("Could not calculate loglikelihood: %s",results$error)))
	loglikelihood = results$loglikelihood;
		
	return(list(success			= TRUE,
				loglikelihood	= loglikelihood));
}


