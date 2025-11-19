# assess the adequacy of a (presumably fitted) homogenous-birth-death-sampling model to explaining a given (presumably real) timetree (e.g., a viral phylogeny sampled through time)
# "Age" is time-before-present, measured from tips to root, with the input tree's youngest tip being at age 0.
# Each provided model should be a named list with the following elements:
#	ages 		: numeric vector listing grid ages in ascending order. The age grid must cover the present day all the way back to stem_time.
#	lambda		: numeric vector of the same length as age_grid[], listing per-capita birth rates (speciation rates) at each age_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
#	mu			: numeric vector of the same length as age_grid[], listing per-capita death rates (extinction rates) at each age_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
#	psi			: numeric vector of the same length as age_grid[], listing per-capita sampling rates (Poissonian detection rates) at each age_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
#	kappa		: numeric vector of the same length as age_grid[], listing the retention probability (upon sampling) at each age_grid point, i.e. the probability that a sampled lineage remains in the pool. If 0, then every sampled lineage becomes a tip.
#	CSA_ages	: optional numeric vector listing concentrated sampling ages, in ascending order
#	CSA_probs	: optional numeric vector listing concentrated sampling probabilities, corresponding to CSA_ages[]
#	CSA_kappas	: optional numeric vector listing retention probabilities during concentrated sampling attempts, corresponding to CSA_ages[]
#	stem_age	: numeric, time-before-present at which the HBDS process starts. If NULL, this is set to the root age of the input tree.
#	end_age		: numeric, time-before-present at which the HBDS process halts. Typically this will be at the input tree's youngest tip (i.e. age 0) or shortly after (i.e. slightly negative age).
model_adequacy_hbds = function(	tree,
								models,							# named list specifying a single HBDS model, or a list of such lists specifying multiple HBDS models. Models are sampled randomly during bootstrapping.
								splines_degree		= 1,		# polynomial degree of time-dependent model parameters (lambda, mu, psi, kappa) between age_grid points
								extrapolate			= FALSE,	# boolean, specifying whether to extrapolate model variables beyond the provided age grid all the way to stem_age and end_age if needed (as constants).
								Nbootstraps 		= 1000,		# integer, number of bootstraps to perform e.g. for checking statistical significances
								max_sim_attempts 	= 1000,		# integer, maximum number of simulation attempts per bootstrap, before giving up. Multiple attempts may be needed if the HBDS model has a high probability of leading to extinction early on
								Nthreads			= 1,		# integer, number of parallel threads to use
								max_extant_tips		= NULL,		# optional maximum number of extant tips per simulation. A simulation is aborted (and that bootstrap iteration skipped) if the number of extant tips exceeds this threshold. Use this to avoid occasional explosions of runtimes (e.g., due to very large generated trees).
								max_model_runtime	= NULL){	# optional maximum time (in seconds) to allocate for each HBDS model simulation (i.e. each bootstrap simulation). Use this to avoid occasional explosions of runtimes (e.g., due to very large generated trees).
	Nbootstraps 	= max(2,Nbootstraps)
	Nthreads 		= max(1,Nthreads)
	max_sim_attempts= max(1,max_sim_attempts)
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
		if(is.null(models[[m]]$stem_age)) models[[m]]$stem_age = root_age
		if(is.null(models[[m]]$end_age)) models[[m]]$end_age = 0
		if(!is.null(models[[m]]$ages)){
			if(any(diff(models[[m]]$ages)<0)) return(list(success=FALSE, error=sprintf("Age grid  of model %d is not in increasing order",m)))
			if(models[[m]]$ages[1]>models[[m]]$end_age){
				# the model's age grid does not cover the end_age
				if(extrapolate){
					models[[m]]$ages 	= c(models[[m]]$end_age,models[[m]]$ages)
					models[[m]]$lambda 	= (if(is.null(models[[m]]$lambda)) NULL else c(models[[m]]$lambda[1],models[[m]]$lambda))
					models[[m]]$mu 		= (if(is.null(models[[m]]$mu)) NULL else c(models[[m]]$mu[1],models[[m]]$mu))
					models[[m]]$psi 	= (if(is.null(models[[m]]$psi)) NULL else c(models[[m]]$psi[1],models[[m]]$psi))
					models[[m]]$kappa 	= (if(is.null(models[[m]]$kappa)) NULL else c(models[[m]]$kappa[1],models[[m]]$kappa))
				}else{
					return(list(success=FALSE, error=sprintf("Age grid must cover the end age (%g), but only reaches %g in model %d. Consider setting extrapolate=TRUE.",models[[m]]$end_age,models[[m]]$ages[1],m)))
				}
			}
			if(tail(models[[m]]$ages,1)<models[[m]]$stem_age){
				# the model's age grid does not cover the stem_age
				if(extrapolate){
					models[[m]]$ages	= c(models[[m]]$ages,models[[m]]$stem_age)
					models[[m]]$lambda 	= (if(is.null(models[[m]]$lambda)) NULL else c(models[[m]]$lambda,tail(models[[m]]$lambda,1)))
					models[[m]]$mu 		= (if(is.null(models[[m]]$mu)) NULL else c(models[[m]]$mu,tail(models[[m]]$mu,1)))
					models[[m]]$psi 	= (if(is.null(models[[m]]$psi)) NULL else c(models[[m]]$psi,tail(models[[m]]$psi,1)))
					models[[m]]$kappa 	= (if(is.null(models[[m]]$kappa)) NULL else c(models[[m]]$kappa,tail(models[[m]]$kappa,1)))
				}else{
					return(list(success=FALSE, error=sprintf("Age grid must cover the stem age (%g), but only reaches %g in model %d. Consider setting extrapolate=TRUE.",models[[m]]$stem_age,tail(models[[m]]$ages,1),m)))
				}
			}
		}
	}
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
	
	# calculate various statistics for the provided tree
	tree_Colless	= tree_imbalance(tree, type="Colless")
	tree_Sackin		= tree_imbalance(tree, type="Sackin")
	tree_LTT		= count_lineages_through_time(tree,Ntimes=10*log2(Ntips))
	tree_LTT$ages 	= pmax(0,root_age-tree_LTT$times)
	tree_clade_ages	= root_age - get_all_distances_to_root(tree)
	tree_tip_ages	= tree_clade_ages[1:Ntips]
	tree_node_ages	= tree_clade_ages[(1+Ntips):(2*Ntips-1)]

	# convert ages to simulation times (time 0 is at the stem, i.e. when the HBDS process starts)
	for(m in seq_len(Nmodels)){
		models[[m]]$times	 	= rev(models[[m]]$stem_age - models[[m]]$ages)
		models[[m]]$lambda 		= (if(is.null(models[[m]]$lambda)) NULL else rev(models[[m]]$lambda))
		models[[m]]$mu 			= (if(is.null(models[[m]]$mu)) NULL else rev(models[[m]]$mu))
		models[[m]]$psi 		= (if(is.null(models[[m]]$psi)) NULL else rev(models[[m]]$psi))
		models[[m]]$kappa 		= (if(is.null(models[[m]]$kappa)) NULL else rev(models[[m]]$kappa))
		models[[m]]$CSA_times	= (if(is.null(models[[m]]$CSA_ages)) NULL else rev(models[[m]]$stem_age - models[[m]]$CSA_ages))
		models[[m]]$CSA_probs 	= (if(is.null(models[[m]]$CSA_probs)) NULL else rev(models[[m]]$CSA_probs))
		models[[m]]$CSA_kappas 	= (if(is.null(models[[m]]$CSA_kappas)) NULL else rev(models[[m]]$CSA_kappas))	
	}

	# simulate & analyze bootstrap trees in parallel
	run_bootstrap = function(b){
		m = sample.int(n=Nmodels, size=1) # pick a random model from the pool
		model = models[[m]]
		sim_attempts = 0
		while(TRUE){
			sim_attempts = sim_attempts + 1
			bootstrap = tryCatch({ generate_tree_hbds(	max_time				= model$stem_age - model$end_age,
														max_extant_tips			= (if(is.null(max_extant_tips)) NULL else 1+max_extant_tips),
														include_extant			= FALSE,
														include_extinct			= FALSE,
														time_grid				= model$times,
														lambda					= model$lambda,
														mu						= model$mu,
														psi						= model$psi,
														kappa					= model$kappa,
														splines_degree			= splines_degree,
														CSA_times				= model$CSA_times,
														CSA_probs				= model$CSA_probs,
														CSA_kappas				= model$CSA_kappas,
														no_full_extinction		= FALSE,
														max_runtime				= max_model_runtime,
														include_birth_times		= FALSE,
														include_death_times		= FALSE)
								}, error = function(e){ list(success=FALSE, error="An unknown error occurred") })
			if((!bootstrap$success) || ((!is.null(max_extant_tips)) && (length(bootstrap$tree$tip.label)>max_extant_tips))){
				# this simulation failed
				if(sim_attempts>=max_sim_attempts){
					return(list(success=FALSE, error=sprintf("Failed to simulate bootstrap tree%s after %d attempts: %s",(if(Nmodels>1) sprintf(" (model %d)",m) else ""),sim_attempts,bootstrap$error)))
				}else{
					# try again
					next
				}
			}else{
				# this simulation succeeded
				break
			}
		}
		bootstrap_tree 		= bootstrap$tree
		NBtips				= length(bootstrap_tree$tip.label)
		bootstrap_root_age 	= get_tree_span(bootstrap_tree)$max_distance
		clade_ages			= bootstrap_root_age - get_all_distances_to_root(bootstrap_tree)
		return(list(success			= TRUE,
					Ntips			= NBtips,
					edge_lengths	= bootstrap_tree$edge.length,
					tip_ages		= clade_ages[1:NBtips],
					node_ages		= clade_ages[(1+NBtips):(2*NBtips-1)],
					Colless			= tree_imbalance(bootstrap_tree, type="Colless"),
					Sackin			= tree_imbalance(bootstrap_tree, type="Sackin"),
					LTT				= count_lineages_through_time(bootstrap_tree,times=(tree_LTT$times+model$stem_age-root_age-bootstrap$root_time))$lineages))
	}
	if((Nbootstraps>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
		# run trials in parallel using multiple forks
		# Note: Forks (and hence shared memory) are not available on Windows
		bootstraps = parallel::mclapply(seq_len(Nbootstraps), 
										FUN = run_bootstrap, 
										mc.cores = min(Nthreads, Nbootstraps), 
										mc.preschedule = TRUE, 
										mc.cleanup = TRUE)
	}else{
		# run in serial mode
		bootstraps = vector(mode="list", Nbootstraps)
		for(b in seq_len(Nbootstraps)){
			bootstraps[[b]] = run_bootstrap(b)
		}
	}
	if(is.null(bootstraps)) return(list(success=FALSE, error="An unknown error caused all bootstraps to fail"))
	
	# omit failed bootstraps
	valid_bootstraps = which(sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$success))
	if(length(valid_bootstraps)<2){
		invalid_bootstraps = which(sapply(seq_len(Nbootstraps), FUN=function(b) !bootstraps[[b]]$success))
		return(list(success=FALSE, error=sprintf("Nearly all bootstraps failed: %s",bootstraps[[invalid_bootstraps[1]]]$error)))
	}
	bootstraps  = bootstraps[valid_bootstraps]
	Nbootstraps = length(bootstraps)

	# calculate P-value of Ntips
	bootstrap_Ntips			= sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$Ntips)
	bootstrap_mean_Ntips	= mean(bootstrap_Ntips, na.rm=TRUE)
	PNtips					= mean(abs(bootstrap_Ntips-bootstrap_mean_Ntips)>=abs(Ntips-bootstrap_mean_Ntips),na.rm=TRUE)

	# calculate P-value of Colless imbalance statistic
	bootstrap_Colless		= sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$Colless)
	bootstrap_mean_Colless	= mean(bootstrap_Colless, na.rm=TRUE)
	PColless				= mean(abs(bootstrap_Colless-bootstrap_mean_Colless)>=abs(tree_Colless-bootstrap_mean_Colless),na.rm=TRUE)

	# calculate P-value of Sackin imbalance statistic
	bootstrap_Sackin		= sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$Sackin)
	bootstrap_mean_Sackin	= mean(bootstrap_Sackin, na.rm=TRUE)
	PSackin					= mean(abs(bootstrap_Sackin-bootstrap_mean_Sackin)>=abs(tree_Sackin-bootstrap_mean_Sackin),na.rm=TRUE)
	
	# calculate confidence intervals of bootstrap LTTs
	bootstrap_LTTs		= t(sapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$LTT))
	bootstrap_LTT_CIs 	= calculate_equal_tailed_CIs_of_curves(bootstrap_LTTs)
	fraction_LTT_in_CI95= mean((tree_LTT$lineages<=bootstrap_LTT_CIs$CI95upper) & (tree_LTT$lineages>=bootstrap_LTT_CIs$CI95lower)) # fraction of the tree's LTT contained in the bootstrap CI95 interval
	
	# compare edge lengths of tree to bootstraps using Kolmogorov-Smirnov test
	bootstrap_edge_lengths 	= lapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$edge_lengths)
	edge_Kolmogorov_Smirnov = bootstrap_Kolmogorov_Smirnov_test(bootstraps=bootstrap_edge_lengths, empirical=tree$edge.length)

	# compare tip ages of tree to bootstraps using Kolmogorov-Smirnov test
	bootstrap_tip_ages 		= lapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$tip_ages)
	tip_Kolmogorov_Smirnov 	= bootstrap_Kolmogorov_Smirnov_test(bootstraps=bootstrap_tip_ages, empirical=tree_tip_ages)
	
	# compare node ages of tree to bootstraps using Kolmogorov-Smirnov test
	bootstrap_node_ages 	= lapply(seq_len(Nbootstraps), FUN=function(b) bootstraps[[b]]$node_ages)
	node_Kolmogorov_Smirnov = bootstrap_Kolmogorov_Smirnov_test(bootstraps=bootstrap_node_ages, empirical=tree_node_ages)

	# wrap all results into a "report" list
	return(list(success						= TRUE,
				Nbootstraps					= Nbootstraps,
				tree_Ntips					= Ntips,
				bootstrap_mean_Ntips 		= bootstrap_mean_Ntips,
				PNtips						= PNtips,
				tree_Colless				= tree_Colless,
				bootstrap_mean_Colless	 	= bootstrap_mean_Colless,
				PColless					= PColless,
				tree_Sackin					= tree_Sackin,
				bootstrap_mean_Sackin 		= bootstrap_mean_Sackin,
				PSackin						= PSackin,
				tree_edgeKS					= edge_Kolmogorov_Smirnov$empirical_KS,
				bootstrap_mean_edgeKS 		= edge_Kolmogorov_Smirnov$mean_bootstrap_KS,
				PedgeKS						= edge_Kolmogorov_Smirnov$Pvalue,
				tree_tipKS					= tip_Kolmogorov_Smirnov$empirical_KS,
				bootstrap_mean_tipKS 		= tip_Kolmogorov_Smirnov$mean_bootstrap_KS,
				PtipKS						= tip_Kolmogorov_Smirnov$Pvalue,
				tree_nodeKS					= node_Kolmogorov_Smirnov$empirical_KS,
				bootstrap_mean_nodeKS 		= node_Kolmogorov_Smirnov$mean_bootstrap_KS,
				PnodeKS						= node_Kolmogorov_Smirnov$Pvalue,
				statistical_tests			= data.frame(	statistic			= c("Ntips", "Colless", "Sackin", "edge-lengths Kolmogorov-Smirnov", "tip-ages Kolmogorov-Smirnov", "node-ages Kolmogorov-Smirnov"),
															tree_value	 		= c(Ntips, tree_Colless, tree_Sackin, edge_Kolmogorov_Smirnov$empirical_KS, tip_Kolmogorov_Smirnov$empirical_KS, node_Kolmogorov_Smirnov$empirical_KS),
															bootstrap_mean		= c(bootstrap_mean_Ntips, bootstrap_mean_Colless, bootstrap_mean_Sackin, edge_Kolmogorov_Smirnov$mean_bootstrap_KS, tip_Kolmogorov_Smirnov$mean_bootstrap_KS, node_Kolmogorov_Smirnov$mean_bootstrap_KS),
															Pvalue	 			= c(PNtips, PColless, PSackin, edge_Kolmogorov_Smirnov$Pvalue, tip_Kolmogorov_Smirnov$Pvalue, node_Kolmogorov_Smirnov$Pvalue)),
				LTT_ages				= tree_LTT$ages,
				tree_LTT				= tree_LTT$lineages,
				bootstrap_LTT_CI		= bootstrap_LTT_CIs,
				fraction_LTT_in_CI95	= fraction_LTT_in_CI95))
}
