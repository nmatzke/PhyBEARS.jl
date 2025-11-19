# generate a random phylogenetic tree according to a homogenous birth-death-sampling process
# the speciation, extinction and continuous (Poissonian) sampling rate can each be time-dependent, and there may be additional discrete sampling times included
# The simulation proceeds in forward time (starting from the root) until one of the stopping criteria are met, OR once all lineages are extinct.
generate_tree_hbds = function(	max_sampled_tips		= NULL, 
								max_sampled_nodes		= NULL, 
								max_extant_tips			= NULL,
								max_extinct_tips		= NULL,
								max_tips				= NULL, 	# integer, max number of tips of any type (extant + extinct + sampled). The simulation is halted once this number is reached.
								max_time				= NULL,
								include_extant			= FALSE,	# logical, whether to include extant non-sampled tips in the final tree
								include_extinct			= FALSE,	# logical, whether to include extinct non-sampled tips in the final tree
								as_generations			= FALSE,	# if FALSE, then edge lengths correspond to time. If TRUE, then edge lengths correspond to generations (hence if include_extant==true and include_extinct==true, all edges will have unit length).
								time_grid				= NULL,		# numeric vector listing grid times in ascending order. The time grid should generally cover the maximum possible simulation time, otherwise everything is polynomially extrapolated (according to splines_degree).
								lambda					= NULL,		# numeric vector of the same length as time_grid[], listing per-capita birth rates (speciation rates) at each time_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
								mu						= NULL,		# numeric vector of the same length as time_grid[], listing per-capita death rates (extinction rates) at each time_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
								psi						= NULL,		# numeric vector of the same length as time_grid[], listing per-capita sampling rates (Poissonian detection rates) at each time_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
								kappa					= NULL,		# numeric vector of the same length as time_grid[], listing the retention probability (upon sampling) at each time_grid point, i.e. the probability that a sampled lineage remains in the pool. If 0, then every sampled lineage becomes a tip.
								splines_degree			= 1,		# polynomial degree of time-dependent model parameters (lambda, mu, psi, kappa) between time-grid points
								CSA_times				= NULL,		# optional numeric vector listing concentrated sampling times, in ascending order
								CSA_probs				= NULL,		# optional numeric vector listing concentrated sampling probabilities, corresponding to CSA_times[]
								CSA_kappas				= NULL,		# optional numeric vector listing retention probabilities during concentrated sampling attempts, corresponding to CSA_times[]
								no_full_extinction		= FALSE,	# if true, then extinction of the entire tree is prevented. This is done by temporarily disabling extinctions when the number of extant tips is 1.
								max_runtime				= NULL,		# maximum time (in seconds) to allow for the computation; if the computation roughly exceeds this threshold, it is aborted. Use this as protection against badly parameterized models. If NULL or <=0, this option is ignored.
								tip_basename			= "",		# basename for tips (e.g. "tip."). 
								node_basename			= NULL,		# basename for nodes (e.g. "node."). If NULL, then nodes will not have any labels.
								edge_basename			= NULL,		# basename for edge (e.g. "edge."). If NULL, then edges will not have any labels.
								include_birth_times		= FALSE,
								include_death_times		= FALSE){
	# basic input checking
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
	if(is.null(max_tips) && is.null(max_sampled_tips) && is.null(max_extant_tips) && is.null(max_extinct_tips) && is.null(max_sampled_nodes) && is.null(max_time)) return(list(success=FALSE, error="ERROR: At least one of {max_tips, max_sampled_tips, max_sampled_nodes, max_extant_tips, max_extinct_tips, max_time} must be non-NULL"));
	if(is.null(time_grid) || (length(time_grid)<=1)){
		if((!is.null(lambda)) && (length(lambda)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of lambda values (%d); since no time grid was provided, you must either provide a single (constant) birth rate or NULL",length(lambda))))
		if((!is.null(mu)) && (length(mu)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of mu values (%d); since no time grid was provided, you must either provide a single (constant) death rate or none",length(mu))))
		if((!is.null(psi)) && (length(psi)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of psi values (%d); since no time grid was provided, you must either provide a single (constant) sampling rate or none",length(psi))))
		if((!is.null(kappa)) && (length(kappa)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of kappa values (%d); since no time grid was provided, you must either provide a single (constant) retention probability or none",length(kappa))))
		# create dummy time grid
		NG = 2;
		time_grid = seq(from=0,to=1,length.out=NG)
		if(!is.null(lambda)){
			lambda = rep(lambda,times=NG)
		}else{
			lambda = rep(0,times=NG)
		}
		if(!is.null(mu)){
			mu = rep(mu,times=NG)
		}else{
			mu = rep(0,times=NG)
		}
		if(!is.null(psi)){
			psi = rep(psi,times=NG)
		}else{
			psi = rep(0,times=NG)
		}
		if(!is.null(kappa)){
			kappa = rep(kappa,times=NG)
		}else{
			kappa = rep(0,times=NG)
		}
	}else{
		NG = length(time_grid);
		if(is.null(lambda)){
			lambda = rep(0,times=NG)
		}else if(length(lambda)==1){
			lambda = rep(lambda,times=NG)
		}else if(length(lambda)!=NG){
			return(list(success=FALSE, error=sprintf("Expected either a single birth-rate lambda or exactly %d birth-rates (=time_grid length), but instead got %d",NG,length(lambda))))
		}
		if(is.null(mu)){
			mu = rep(0,times=NG)
		}else if(length(mu)==1){
			mu = rep(mu,times=NG)
		}else if(length(mu)!=NG){
			return(list(success=FALSE, error=sprintf("Expected either a single death-rate mu or exactly %d death-rates (=time_grid length), but instead got %d",NG,length(mu))))
		}
		if(is.null(psi)){
			psi = rep(0,times=NG)
		}else if(length(psi)==1){
			psi = rep(psi,times=NG)
		}else if(length(psi)!=NG){
			return(list(success=FALSE, error=sprintf("Expected either a single sampling-rate psi or exactly %d sampling-rates (=time_grid length), but instead got %d",NG,length(psi))))
		}
		if(is.null(kappa)){
			kappa = rep(0,times=NG)
		}else if(length(kappa)==1){
			kappa = rep(kappa,times=NG)
		}else if(length(kappa)!=NG){
			return(list(success=FALSE, error=sprintf("Expected either a single retention probability kappa or exactly %d probabilities (=time_grid length), but instead got %d",NG,length(kappa))))
		}
		if(any(diff(time_grid)<=0)) return(list(success = FALSE, error = sprintf("Values in time_grid must be strictly increasing")))
	}
	NCSA = (if(is.null(CSA_times)) 0 else length(CSA_times))
	if((NCSA==0) && (!is.null(CSA_probs)) && (length(CSA_probs)>0)) return(list(success=FALSE, error="CSA_times is missing while CSA_probs was provided; either provide both or none"))
	if((NCSA==0) && (!is.null(CSA_kappas)) && (length(CSA_probs)>0)) return(list(success=FALSE, error="CSA_times is missing while CSA_kappas was provided; either provide both or none"))
	if((NCSA>0) && is.null(CSA_probs)) return(list(success=FALSE, error="CSA_probs is missing while CSA_times was provided; either provide both or none"))
	if((NCSA>0) && is.null(CSA_kappas)) return(list(success=FALSE, error="CSA_kappas is missing while CSA_times was provided; either provide both or none"))
	if((NCSA>0) && (!is.null(CSA_probs)) && (!is.null(CSA_kappas))){
		if(length(CSA_probs)==1) CSA_probs 	 = rep(CSA_probs, times=NCSA)
		if(length(CSA_kappas)==1) CSA_kappas = rep(CSA_kappas, times=NCSA)
		if(length(CSA_times)!=length(CSA_probs)) return(list(success=FALSE, error="Number of CSA_times (%d) differs from number of CSA_probs (%d)",length(CSA_times),length(CSA_probs)))
		if(length(CSA_times)!=length(CSA_kappas)) return(list(success=FALSE, error="Number of CSA_times (%d) differs from number of CSA_kappas (%d)",length(CSA_times),length(CSA_kappas)))
		if(any(diff(CSA_times)<=0)) return(list(success=FALSE, error="CSA_times must be in strictly increasing order"))
		if(any(CSA_probs<0) || any(CSA_probs>1)) return(list(success=FALSE, error="CSA_probs must be true probabilities, and thus between 0 and 1"))
		if(any(CSA_kappas<0) || any(CSA_kappas>1)) return(list(success=FALSE, error="CSA_kappas must be true probabilities, and thus between 0 and 1"))
	}
	if(is.null(max_runtime)) max_runtime = 0
	
	# generate tree
	results = generate_random_tree_HBDS_CPP(max_sampled_tips		= (if(is.null(max_sampled_tips)) -1 else max_sampled_tips),
											max_sampled_nodes		= (if(is.null(max_sampled_nodes)) -1 else max_sampled_nodes),
											max_extant_tips			= (if(is.null(max_extant_tips)) -1 else max_extant_tips),
											max_extinct_tips		= (if(is.null(max_extinct_tips)) -1 else max_extinct_tips),
											max_tips				= (if(is.null(max_tips)) -1 else max_tips),
											max_time				= (if(is.null(max_time)) -1 else max_time),
											time_grid				= time_grid,
											birth_rates				= lambda,
											death_rates				= mu,
											sampling_rates			= psi,
											retention_probs			= kappa,
											splines_degree			= splines_degree,
											CSA_times				= (if(is.null(CSA_times)) numeric() else CSA_times),
											CSA_probs				= (if(is.null(CSA_probs)) numeric() else CSA_probs),
											CSA_kappas				= (if(is.null(CSA_kappas)) numeric() else CSA_kappas),
											as_generations			= as_generations,
											no_full_extinction		= no_full_extinction,
											runtime_out_seconds		= max_runtime,
											include_extant			= include_extant,
											include_extinct			= include_extinct,
											include_birth_times		= include_birth_times,
											include_death_times		= include_death_times)
	if(!results$success) return(list(success=FALSE, error=results$error)); # something went wrong
	Ntips	= results$Ntips
	Nnodes 	= results$Nnodes
	Nedges 	= results$Nedges
	tree = list(Nnode 		= Nnodes,
				tip.label 	= paste(tip_basename, 1:Ntips, sep=""),
				node.label 	= (if(is.null(node_basename)) NULL else paste(node_basename, 1:Nnodes, sep="")),
				edge.label 	= (if(is.null(edge_basename)) NULL else paste(edge_basename, 1:Nedges, sep="")),
				edge 		= matrix(results$tree_edge,ncol=2,byrow=TRUE) + 1,
				edge.length = results$edge_length,
				root 		= results$root+1)
	class(tree) = "phylo";
	attr(tree,"order") = NULL

	return(list(success				= TRUE,
				tree				= tree,
				root_time			= results$root_time,
				final_time			= results$final_time,
				root_age			= results$final_time - results$root_time,
				Nbirths		 		= results$Nbirths,
				Ndeaths				= results$Ndeaths,
				Nsamplings			= results$Nsamplings,
				Nretentions			= results$Nretentions,
				sampled_clades		= results$sampled_clades+1,
				retained_clades		= results$retained_clades+1,
				extant_tips			= (if(include_extant) results$extant_tips+1 else integer()),
				extinct_tips		= (if(include_extinct) results$extinct_tips+1 else integer())));
	
}