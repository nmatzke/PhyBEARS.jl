# assess the adequacy of a (presumably fitted) homogenous-birth-death model to explaining a given (presumably real) timetree
# Each provided model should be a named list with the following elements:
#	ages	: numeric vector of size NG, listing ages in ascending order, covering at least 0 until the root
#	PSR		: numeric vector of size NG, listing the pulled speciation rate at each age-grid point
model_adequacy_hbd = function(	tree, 						# ultrametric timetree of class "phylo"
								models,						# named list specifying a single HBD model, or a list of such lists specifying multiple HBD models. Models are sampled randomly during bootstrapping.
								splines_degree	= 1,		# integer
								extrapolate		= FALSE,	# boolean, specifying whether to extrapolate model variables beyond their original time grid all the way to root_age and present-day if needed (as constants).
								Nbootstraps 	= 1000,		# integer, number of bootstraps to perform e.g. for checking statistical significances
								Nthreads		= 1){		# integer, number of parallel threads to use
	Nbootstraps 	= max(2,Nbootstraps)
	Nthreads 		= pmax(1,Nthreads)
	tree			= castor::multifurcations_to_bifurcations(tree)$tree # make sure tree is bifurcating, so that node & edge counts are the same as for the bootstraps
	Ntips			= length(tree$tip.label)
	Nnodes			= tree$Nnode
	root_age 		= get_tree_span(tree)$max_distance
	
	# check models
	if("list" %in% class(models[[1]])){
		Nmodels = length(models)
	}else{
		models  = list(models)
		Nmodels = 1
	}
	for(m in seq_len(Nmodels)){
		if(!is.null(models[[m]]$ages)){
			if(any(diff(models[[m]]$ages)<0)) return(list(success=FALSE, error=sprintf("Age grid  of model %d is not in increasing order",m)))
			if(models[[m]]$ages[1]>0){
				# the model's age grid does not cover the present-day
				if(extrapolate){
					models[[m]]$ages = c(0,models[[m]]$ages)
					models[[m]]$PSR  = (if(is.null(models[[m]]$PSR)) NULL else c(models[[m]]$PSR[1],models[[m]]$PSR))
				}else{
					return(list(success=FALSE, error=sprintf("Age grid must cover the present-day (age 0), but only reaches %g in model %d. Consider setting extrapolate=TRUE.",models[[m]]$ages[1],m)))
				}
			}
			if(tail(models[[m]]$ages,1)<root_age){
				# the model's age grid does not cover the root
				if(extrapolate){
 					models[[m]]$ages = c(models[[m]]$ages,root_age*1.001)
 					models[[m]]$PSR  = (if(is.null(models[[m]]$PSR)) NULL else c(models[[m]]$PSR,1.001*tail(models[[m]]$PSR,1)))
				}else{
					return(list(success=FALSE, error=sprintf("Age grid must cover the root age (%g), but only reaches %g in model %d. Consider setting extrapolate=TRUE.",root_age,tail(models[[m]]$ages,1),m)))
				}
			}
		}
	}	
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))

	# calculate various statistics for the provided tree
	tree_gamma 		= gamma_statistic(tree)
	tree_Colless	= tree_imbalance(tree, type="Colless")
	tree_Sackin		= tree_imbalance(tree, type="Sackin")
	tree_LTT		= count_lineages_through_time(tree,Ntimes=10*log2(Ntips));
	tree_LTT$ages 	= pmax(0,root_age-tree_LTT$times)
	tree_node_ages	= root_age - get_all_distances_to_root(tree)[(1+Ntips):(2*Ntips-1)]
	tree_node_sizes	= as.numeric(count_tips_per_node(tree))
	tree_Blum		= sum(log(tree_node_sizes))
		
	# simulate bootstrap trees
	# for every unique model, perform all replicate simulations at once for efficiency
	bootstrap_trees = vector(mode="list", Nbootstraps)
	bootstrap2model = sample.int(n=Nmodels, size=Nbootstraps, replace=TRUE)
	for(m in seq_len(Nmodels)){
		bootstraps_this_model = which(bootstrap2model==m)
		if(length(bootstraps_this_model)==0) next
		sims = generate_tree_hbd_reverse(	Ntips			= Ntips,
											crown_age		= root_age,
											age_grid		= models[[m]]$ages,
											PSR				= models[[m]]$PSR,
											splines_degree	= splines_degree,
											Ntrees			= length(bootstraps_this_model))
		if(!sims$success) return(list(success=FALSE, error=sprintf("Failed to simulate bootstrap trees%s: %s",(if(Nmodels==1) "" else sprintf(" (model %d)",m)),sims$error)))
		bootstrap_trees[bootstraps_this_model] = sims$trees
	}
	
		
	# analyze bootstrap trees in parallel										
	analyze_bootstrap = function(b){
		bootstrap_tree = bootstrap_trees[[b]]	
		node_sizes = count_tips_per_node(bootstrap_tree)
		return(list(edge_lengths	= bootstrap_tree$edge.length,
					node_ages		= root_age - get_all_distances_to_root(bootstrap_tree)[(1+Ntips):(2*Ntips-1)],
					node_sizes		= as.numeric(node_sizes),
					Blum			= sum(log(node_sizes)), # Blum and Francois (2006). Which random processes describe the Tree of Life? A large-scale study of phylogenetic tree imbalance. Systematic Biology. 55:685-691.
					gamma 			= gamma_statistic(bootstrap_tree),
					Colless			= tree_imbalance(bootstrap_tree, type="Colless"),
					Sackin			= tree_imbalance(bootstrap_tree, type="Sackin"),
					LTT				= count_lineages_through_time(bootstrap_tree,times=tree_LTT$times)$lineages))
	}
	if((Nbootstraps>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
		# run trials in parallel using multiple forks
		# Note: Forks (and hence shared memory) are not available on Windows
		bootstraps = parallel::mclapply(seq_len(Nbootstraps), 
										FUN = analyze_bootstrap, 
										mc.cores = min(Nthreads, Nbootstraps), 
										mc.preschedule = TRUE, 
										mc.cleanup = TRUE)
	}else{
		# run in serial mode
		bootstraps = vector(mode="list", Nbootstraps)
		for(b in seq_len(Nbootstraps)){
			bootstraps[[b]] = analyze_bootstrap(b)
		}
	}
	if(is.null(bootstraps)) return(list(success=FALSE, error="An unknown error caused all bootstraps to fail"))
	
	# calculate P-value of gamma statistic
	bootstrap_gammas		= sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$gamma)
	bootstrap_mean_gamma 	= mean(bootstrap_gammas, na.rm=TRUE)
	Pgamma 					= mean(abs(bootstrap_gammas-bootstrap_mean_gamma)>=abs(tree_gamma-bootstrap_mean_gamma),na.rm=TRUE)
	
	# calculate P-value of Colless imbalance statistic
	bootstrap_Colless		= sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$Colless)
	bootstrap_mean_Colless	= mean(bootstrap_Colless, na.rm=TRUE)
	PColless				= mean(abs(bootstrap_Colless-bootstrap_mean_Colless)>=abs(tree_Colless-bootstrap_mean_Colless),na.rm=TRUE)

	# calculate P-value of Sackin imbalance statistic
	bootstrap_Sackin		= sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$Sackin)
	bootstrap_mean_Sackin	= mean(bootstrap_Sackin, na.rm=TRUE)
	PSackin					= mean(abs(bootstrap_Sackin-bootstrap_mean_Sackin)>=abs(tree_Sackin-bootstrap_mean_Sackin),na.rm=TRUE)

	# calculate P-value of Blum imbalance statistic
	bootstrap_Blum		= sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$Blum)
	bootstrap_mean_Blum	= mean(bootstrap_Blum, na.rm=TRUE)
	PBlum				= mean(abs(bootstrap_Blum-bootstrap_mean_Blum)>=abs(tree_Blum-bootstrap_mean_Blum),na.rm=TRUE)
	
	# calculate confidence intervals of bootstrap LTTs
	bootstrap_LTTs		= t(sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$LTT))
	bootstrap_LTT_CIs 	= calculate_equal_tailed_CIs_of_curves(bootstrap_LTTs)
	fraction_LTT_in_CI95= mean((tree_LTT$lineages<=bootstrap_LTT_CIs$CI95upper) & (tree_LTT$lineages>=bootstrap_LTT_CIs$CI95lower)) # fraction of the tree's LTT contained in the bootstrap CI95 interval
	
	# compare edge lengths of tree to bootstraps using Kolmogorov-Smirnov test
	bootstrap_edge_lengths 	= lapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$edge_lengths)
	edge_Kolmogorov_Smirnov = bootstrap_Kolmogorov_Smirnov_test(bootstraps=bootstrap_edge_lengths, empirical=tree$edge.length)

	# compare node ages of tree to bootstraps using Kolmogorov-Smirnov test
	# this is essentially a comparison of LTTs, since the LTT is the (non-normalized) cumulative distribution function of node ages
	bootstrap_node_ages 	= lapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$node_ages)
	node_Kolmogorov_Smirnov = bootstrap_Kolmogorov_Smirnov_test(bootstraps=bootstrap_node_ages, empirical=tree_node_ages)

	# compare node sizes (Ntips per node) of tree to bootstraps using Kolmogorov-Smirnov test
	bootstrap_node_sizes 	= lapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$node_sizes)
	size_Kolmogorov_Smirnov = bootstrap_Kolmogorov_Smirnov_test(bootstraps=bootstrap_node_sizes, empirical=tree_node_sizes)
	
	# wrap all results into a "report" list
	return(list(success						= TRUE,
				Nbootstraps					= Nbootstraps,
				tree_gamma					= tree_gamma,
				bootstrap_mean_gamma 		= bootstrap_mean_gamma,
				Pgamma						= Pgamma,
				tree_Colless				= tree_Colless,
				bootstrap_mean_Colless 		= bootstrap_mean_Colless,
				PColless					= PColless,
				tree_Sackin					= tree_Sackin,
				bootstrap_mean_Sackin 		= bootstrap_mean_Sackin,
				PSackin						= PSackin,
				tree_Blum					= tree_Blum,
				bootstrap_mean_Blum 		= bootstrap_mean_Blum,
				PBlum						= PBlum,
				tree_edgeKS					= edge_Kolmogorov_Smirnov$empirical_KS,
				bootstrap_mean_edgeKS 		= edge_Kolmogorov_Smirnov$mean_bootstrap_KS,
				PedgeKS						= edge_Kolmogorov_Smirnov$Pvalue,
				tree_nodeKS					= node_Kolmogorov_Smirnov$empirical_KS,
				bootstrap_mean_nodeKS 		= node_Kolmogorov_Smirnov$mean_bootstrap_KS,
				PnodeKS						= node_Kolmogorov_Smirnov$Pvalue,
				tree_sizeKS					= size_Kolmogorov_Smirnov$empirical_KS,
				bootstrap_mean_sizeKS 		= size_Kolmogorov_Smirnov$mean_bootstrap_KS,
				PsizeKS						= size_Kolmogorov_Smirnov$Pvalue,
				statistical_tests			= data.frame(	statistic			= c("gamma", "Colless", "Sackin", "Blum", "edge-lengths Kolmogorov-Smirnov", "node-ages Kolmogorov-Smirnov", "node-sizes Kolmogorov-Smirnov"),
															tree_value	 		= c(tree_gamma, tree_Colless, tree_Sackin, tree_Blum, edge_Kolmogorov_Smirnov$empirical_KS, node_Kolmogorov_Smirnov$empirical_KS, size_Kolmogorov_Smirnov$empirical_KS),
															bootstrap_mean		= c(bootstrap_mean_gamma, bootstrap_mean_Colless, bootstrap_mean_Sackin, bootstrap_mean_Blum, edge_Kolmogorov_Smirnov$mean_bootstrap_KS, node_Kolmogorov_Smirnov$mean_bootstrap_KS, size_Kolmogorov_Smirnov$mean_bootstrap_KS),
															Pvalue	 			= c(Pgamma, PColless, PSackin, PBlum, edge_Kolmogorov_Smirnov$Pvalue, node_Kolmogorov_Smirnov$Pvalue, size_Kolmogorov_Smirnov$Pvalue)),
				LTT_ages				= tree_LTT$ages,
				tree_LTT				= tree_LTT$lineages,
				bootstrap_LTT_CI		= bootstrap_LTT_CIs,
				fraction_LTT_in_CI95	= fraction_LTT_in_CI95))
}
