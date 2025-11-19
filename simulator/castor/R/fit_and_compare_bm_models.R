# Fit a Brownian Motion model of continuous trait evolution to two different sets of trees, and compare the fitted models
# This function can calculate the two-sided statistical significance of the diffusivity ratio (D1/D2), or equivalently, the statistical significance of the log_difference |log(D1)-log(D2)| (separately for each entry in D, in the case of multivariate BM)
# The main inputs are:
# 	Two lists of trees
#	The corresponding numerical tip states
fit_and_compare_bm_models = function(	trees1, 				# either a single tree in phylo format, or a list of trees
										tip_states1,			# Numerical trait values for the tips of each tree in trees1. If trees1[] is a single tree, then tip_states1[] should be either a 1D vector of size Ntips or a 2D matrix of size Ntips x Ntraits. If trees1[] is a list of trees, then tip_states1 should either be a list of 1D vectors or a 2D matrix of size Ntips x Ntraits.
										trees2, 				# either a single tree in phylo format, or a list of trees
										tip_states2,			# Numerical trait values for the tips of each tree in trees2. If trees2[] is a single tree, then tip_states2[] should be either a 1D vector of size Ntips or a 2D matrix of size Ntips x Ntraits. If trees2[] is a list of trees, then tip_states2 should either be a list of 1D vectors or a 2D matrix of size Ntips x Ntraits.
										Nbootstraps		= 0,	# integer, number of bootstraps to perform for estimating confidence intervals. Set to <=0 for no bootstrapping.
										Nsignificance	= 0,	# integer, number of random simulations to perform for testing the statistical significance of the log_difference of diffusivities (|log(D1)-log(D2)|). If <=0, the statistical significance is not calculated
										check_input		= TRUE,
										verbose			= FALSE,
										verbose_prefix 	= ""){
	# fit BM to each set of trees
	if(verbose) cat(sprintf("%sFitting BM to first tree set..\n",verbose_prefix))
	fit1 = fit_bm_model(trees=trees1, tip_states=tip_states1, Nbootstraps=Nbootstraps, check_input=check_input)
	if(!fit1$success) return(list(success=FALSE, error=sprintf("Failed to fit BM to tree set 1: %s",fit1$error)))

	if(verbose) cat(sprintf("%sFitting BM to second tree set..\n",verbose_prefix))
	fit2 = fit_bm_model(trees=trees2, tip_states=tip_states2, Nbootstraps=Nbootstraps, check_input=check_input)
	if(!fit2$success) return(list(success=FALSE, error=sprintf("Failed to fit BM to tree set 2: %s",fit2$error)))
	
	log_difference 	= abs(log(fit1$diffusivity) - log(fit2$diffusivity))
	if(Nsignificance>0){
		# calculate statistical significance of log_difference
		# Note: In the case of multivariate BM (i.e. Ntraits>1), the log_difference is a matrix and hence we are calculating Ntraits x Ntraits statistical significances
		if(verbose) cat(sprintf("%sCalculating statistical significance of ratio D1/D2..\n",verbose_prefix))
		# bring all data into a standard format
		if("phylo" %in% class(trees1)) trees1 = list(trees1)
		if("phylo" %in% class(trees2)) trees2 = list(trees2)
		Ntrees1 = length(trees1)
		Ntrees2 = length(trees2)
		if(!("list" %in% class(tip_states1))) tip_states1 = list(tip_states1)
		if(!("list" %in% class(tip_states2))) tip_states2 = list(tip_states2)
		if(verbose) cat(sprintf("%s  Fitting common BM model to both tree sets..\n",verbose_prefix))
		fit_common = fit_bm_model(trees=c(trees1,trees2), tip_states=c(tip_states1,tip_states2), Nbootstraps=0, check_input=FALSE)
		if(verbose) cat(sprintf("%s  Assessing significance over %d BM simulations..\n",verbose_prefix,Nsignificance))
		random_tip_states1 	= vector(mode="list", Ntrees1)
		random_tip_states2 	= vector(mode="list", Ntrees2)
		Ngreater = 0
		Nsuccess = 0
		for(r in 1:Nsignificance){
			# simulate the same BM model on each of the input trees
			for(tr in 1:Ntrees1){
				random_tip_states1[[tr]] = simulate_bm_model(trees1[[tr]], diffusivity=fit_common$diffusivity, include_tips=TRUE, include_nodes=FALSE, drop_dims=TRUE)$tip_states
			}
			for(tr in 1:Ntrees2){
				random_tip_states2[[tr]] = simulate_bm_model(trees2[[tr]], diffusivity=fit_common$diffusivity, include_tips=TRUE, include_nodes=FALSE, drop_dims=TRUE)$tip_states
			}
			# fit BM separately to each tree set
			random_fit1 = fit_bm_model(trees=trees1, tip_states=random_tip_states1, Nbootstraps=0, check_input=FALSE)
			if(!random_fit1$success){
				# fitting failed
				if(verbose) cat(sprintf("%s  WARNING: BM fitting failed for random simulation #%d, for tree set 1\n",verbose_prefix,r))
				next;
			}
			random_fit2 = fit_bm_model(trees=trees2, tip_states=random_tip_states2, Nbootstraps=0, check_input=FALSE)
			if(!random_fit2$success){
				# fitting failed
				if(verbose) cat(sprintf("%s  WARNING: BM fitting failed for random simulation #%d, for tree set 2\n",verbose_prefix,r))
				next;
			}
			Nsuccess = Nsuccess + 1
			# compare the obtained log_difference to the original one
			random_log_difference = abs(log(random_fit1$diffusivity) - log(random_fit2$diffusivity))
			Ngreater = Ngreater + (random_log_difference>=log_difference)
		}
		significance = Ngreater / Nsuccess
	}
		
	return(list(success			= TRUE,
				fit1			= fit1,
				fit2			= fit2,
				log_difference 	= log_difference,
				significance	= (if(Nsignificance>0) significance else NULL),
				fit_common		= (if(Nsignificance>0) fit_common else NULL)))

}