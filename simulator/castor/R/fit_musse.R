# Fit MuSSE (Multiple State Speciation Extinction) model for a single evolving discrete trait on a tree, based on the observed states at the tips
#    MuSSE is an extension of the BiSSE model by Maddison (2007), allowing for more than 2 states.
#    The model can account for tips with unknown state as well as for sampling biases (species are included in the tree non-uniformly) and reveal biases (tips with known state are not uniformly chosen).
#    This implementation also extends the classical MuSSE model to allow for Poissonian sampling through time, assuming that sampled tips are removed from the pool of extant tips.
#
# The function also supports an extension of MuSSE, called Hidden State Speciation and Extinction (HiSSE) model, which distinguishes between latent and observed ("proxy") states. [Beaulieu and Meara, 2016]
#   Birth and death rates are determined by the latent states, and a transition matrix defines the Markov transition rates between latent states, but only proxy states are known/observed for the tips. 
#   Each proxy state represents a distinct set of latent states. For example, proxy state 1 may represent latent states A,B,C, and proxy state 2 may represent latent states D,E.
#   The terminology used below is that "states" refers to the trait actually influencing birth/death rates (in a HiSSE or MuSSE model), while "proxy-states" refers to the trait actually observed for the tips.
#	For the special case of MuSSE (i.e. non-HiSSE) models, "proxy-states" and "states" are exactly the same and no distinction is made between them.
#   To define a HiSSE model, use the NPstates and proxy_map options below. 
#   If NPstates==Nstates, each proxy-state is assumed to correspond to a state, and thus a plain MuSSE model is used instead.
#
# The function can also estimate confidence intervals using parametric bootstrapping.
#
# References:
#   Maddison et al (2007). Estimating a binary character's effect on speciation and extinction. Systematic Biology. 56:701-710
#   FitzJohn et al. (2009). Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology. 58:595-611
#   FitzJohn (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution. 3:1084-1092
#	Beaulieu and Meara (2016). Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic Biology. 65:583-601
#	Louca and Pennell (2020). A general and efficient algorithm for the likelihood of diversification and discrete-trait evolutionary models. Systematic Biology. 69:545-556
#	Kuhnert, Stadler et al. (2016). Phylodynamics with migration: A computational framework to quantify population structure from genomic data. Molecular Biology and Evolution. 33:2102-2116
#
# Requirements:
#   Tree can include multi- and mono-furcations.
#   Tree must be rooted. Root will be determined automatically as the node with no parent.
#   Tree must be ultrametric (e.g. a timetree of extant species). In particular, all tips are assumed to have age 0.
fit_musse = function(	tree, 
						Nstates,								# number of discrete states, influencing birth and death rates.
						NPstates				= NULL,			# optional number of discrete "proxy" states. Each proxy-state may correspond to one or more distinct states. Proxy-states do not directly influence birth/death rates, but each proxy-state represents a set of states. Proxy-states are only needed when defining a "hidden state speciation and extinction" (HiSSE) model [Beaulieu and Meara, 2016]. If this is NULL or 0, then proxy-states are equal to states, and NPstates=Nstates; this is the original MuSSE model (i.e. no latent states).
						proxy_map				= NULL,			# optional 1D integer vector of size Nstates, mapping states to proxy-states. Hence, proxy_map[s] is an integer in 1:NPstates, specifying which proxy-state the state s belongs to. Only relevant if NPstates!=NULL and NPstates!=Nstates.
						state_names				= NULL,			# optional 1D character vector of size Nstates, specifying a name/description for each state. This does not influence any of the calculations, and is merely used to add human-readable names (rather than integers) to the returned vectors/matrices.
						tip_pstates				= NULL,			# 1D integer array of size Ntips, listing integers in 1:NPstates, specifying the proxy state at each tip. Can also be NULL. May include NA (unknown tip states).
						tip_priors 				= NULL, 		# 2D numerical array of either size Ntips x Nstates or size Ntips x NPstates, listing the prior likelihood at each state (or proxy-state) at each tip. Specifically, tip_priors[i,s] is the likelihood of observing the data at tip i (sampling the tip and observing the observed state), if tip i had state (or proxy-state) s. Can be provided alternatively to tip_pstates.
						sampling_fractions		= 1,			# optional numerical vector of size NPstates, indicating the present-day sampling fractions (fraction of extant species included as tips in the tree) conditional upon each species' proxy state. This can be used to incorporate detection biases for species, depending on their proxy state. Can also be NULL or a single number (in which case sampling fractions are assumed to be independent of proxy state).
						reveal_fractions		= 1,			# optional numerical vector of size NPstates, indicating the resolved fractions (fraction of tree-tips with revealed, i.e. non hidden, state) conditional upon each tip's proxy state. Hence, reveal_fractions[s] is the probability that a species with proxy state s will have revealed state, conditional upon being included in the tree. This can be used to incorporate reveal biases for tips, depending on their proxy state. Can also be NULL or a single number (in which case reveal fractions are assumed to be independent of proxy state).
						sampling_rates			= 0,			# optional numerical vector of size NPstates, indicating the Poissonian lineage sampling rates, conditional upon each species' proxy state. This can be used to account for tips sampled continuously through time rather than at present-day. Can also be a single number (in which case Poissonian sampling rates are assumed to be independent of proxy state). Can also be NULL, in which case Poissonian sampling is assumed to not have occurred.
						transition_rate_model	= "ARD",		# either "ER" or "SYM" or "ARD" or "SUEDE" or "SRD" or a 2D integer matrix of size Nstates x Nstates mapping entries of the transition matrix to a set of independent transition-rate parameters. The format and interpretation is the same as for index_matrix generated by the function get_transition_index_matrix(..).
						birth_rate_model		= "ARD",		# either "ER" (equal rates) or "ARD" (all rates different) or an integer vector of length Nstates, mapping birth-rates to a set of indepedent birth-rate parameters (indices 1,2,..). For example, the model vector c(1,2,1) means that the birth-rates lambda1 & lambda3 must be the same, but lambda2 is independent.
						death_rate_model		= "ARD",		# either "ER" (equal rates) or "ARD" (all rates different) or an integer vector of length Nstates, mapping death-rates to a set of indepedent death-rate parameters (indices 1,2,..). For example, the model vector c(1,2,1) means that the death-rates mu1 & mu3 must be the same, but mu2 is independent.
						transition_matrix		= NULL,			# either NULL, or a transition matrix of size Nstates x Nstates, such that transition_matrix[r,c] is the transition rate r-->c. If NULL, the entire transition matrix will be fitted via maximum-likelihood. NA-values are allowed in the matrix, corresponding to rates that are to be fitted.
						birth_rates				= NULL,			# either NULL, or a single scalar, or a numeric vector of size Nstates, listing speciation rates at each state. Can include NA, corresponding to rates that are to be fitted. Non-NA values are assumed to be fixed.
						death_rates				= NULL,			# either NULL, or a single scalar, or a numeric vector of size Nstates, listing extinction rates at each state. Can include NA, corresponding to rates that are to be fitted. Non-NA values are assumed to be fixed.
						first_guess				= NULL,			# either NULL, or a named list containing optional entries 'transition_matrix', 'birth_rates' and/or 'death_rates'. Each entry specifies start values (first guess) for the particular model parameter. Some start values may include NA, corresponding to unknowns which are to be guessed internally. Note that any provided lower & upper bounds overrule the first_guess, i.e. first_guess is forced to be within any specified bounds.
						lower					= NULL,			# either NULL, or a named list containing optional entries 'transition_matrix', 'birth_rates' and/or 'death_rates'. Each entry specifies lower bounds for model parameters (only relevant if a parameter is not fixed). Parameters not listed in lower[] are assumed not explicitly bounded (except for >=0 when reasonable). To keep a parameter unbounded, you can also use NA, -Inf or NaN.
						upper					= NULL,			# either NULL, or a named list containing optional entries 'transition_matrix', 'birth_rates' and/or 'death_rates'. Each entry specifies upper bounds for model parameters (only relevant if a parameter is not fixed). Parameters not listed in upper[] are assumed not explicitly bounded. To keep a parameter unbounded, you can also use NA, +Inf or NaN.
						root_prior 				= "auto",		# can be 'flat', 'empirical', 'likelihoods', 'max_likelihood', 'auto' or a numeric vector of size Nstates. Specifies the weights for averaging the root's state likelihoods in order to obtain an overall model likelihood. 'empirical' means the root's prior is set to the proportions of (estimated) extant species in each state (correcting for sampling_fractions and reveal_fractions). 'likelihoods' means that the computed state-likelihoods of the root (i.e. D_s) are used, after normalizing to obtain an probability distribution; this is the approach used in hisse::hisse if root.p=NULL, and the approach in diversitree when root=ROOT.OBS.
						root_conditioning		= "auto",		# can be 'none', 'madfitz', 'herr_als', 'stem', 'crown' or 'auto', specifying how to condition the state-likelihoods at the root prior to averaging. "none" (no conditioning) uses the raw state-likelihoods, as originally described by Maddison (2007). "madfitz" and "herr_als" are the options implemented in the hisse package, conditioning the root's state-likelihoods based on the birth-rates and the computed extinction probability (after or before averaging, respectively). For typical use cases, we suspect that this option does not matter much.
						oldest_age				= NULL,			# oldest age to consider. Typically this will be <=root_age. If NULL, this will be set to root_age.
						Ntrials 				= 1,			# (integer) number of trials (starting points) for fitting the model. If 0, then no fitting is performed, and only the a first-guess (start params) is evaluated and returned.
						max_start_attempts		= 10,			# integer, number of times to attempt finding a valid start point (per trial) before giving up. Randomly choosen start parameters may result in Inf/undefined objective, so this option allows the algorithm to keep looking for valid starting points.
						optim_algorithm 		= "nlminb",		# either "optim", "nlminb" or "subplex". What algorithm to use for fitting.
						optim_max_iterations	= 10000,		# maximum number of iterations of the optimization algorithm (per trial)
						optim_max_evaluations	= NULL,			# maximum number of evaluations of the objective function (per trial). If left unspecified, it is chosen automatically based on the number of iterations.
						optim_rel_tol			= 1e-6,			# relative tolerance when optimizing the objective function
						check_input 			= TRUE,			# (bool) perform some basic sanity checks on the input data. Set this to FALSE if you're certain your input data is valid.
						include_ancestral_likelihoods = FALSE,	# (boolean) whether to also calculate ancestral state likelihoods at each node, based on their descending subtree. These are the D values from the MuSSE model, and may be used as "local" ancestral state estimates.
						Nthreads 				= 1,			# (integer) number of threads for running multiple fitting trials in parallel
						Nbootstraps				= 0,			# (integer) optional number of parametric-bootstrap samples (random simulations of the ML-fitted model) for estimating confidence intervals of fitted parameters. If 0, no parametric bootstrapping is performed. Typical values are 10-100.
						Ntrials_per_bootstrap	= NULL,			# (integer) optional number of fitting trials for each bootstrap sampling. If NULL, this is set equal to Ntrials. A smaller Ntrials_per_bootstrap will reduce computation, at the expense of increasing the estimated confidence intervals (i.e. yielding more conservative estimates of confidence).
						max_condition_number	= 1e4,			# (positive unitless number) the maximum acceptable condition number for Dmap (as estimated from the linearized dynamics), when choosing the integration interval size. A larger max_condition number leads to fewer age-splits, thus faster computation but also lower accuracy. Hence, this number controls the trade-off between speed and accuracy. Typical values are 1e4 (slower, more accurate) up to 1e8 (faster, less accurate).
						relative_ODE_step		= 0.1,			# (positive unitless number) relative integration time step for the ODE solvers. Relative to the typical time scales of the dynamics, as estimated from the theoretically maximum possible rate of change of D or E. Typical values are 0.001 - 0.1.
						E_value_step			= 1e-4,			# (positive unitless mumber) unitless number, relative step for interpolating E over time. So a E_value_step of 0.001 means that E is recorded and interpolated between points between which E differs by roughy 0.001. Typical values are 0.01-0.0001. A smaller E_value_step increases interpolation accuracy, but also increases memory requirements and adds runtime (scales with the tree's age span, not Ntips).
						D_temporal_resolution	= 100,			# (positive unitless number) relative resolution for interpolating Dmap over time. This is relative to the "typical" time scales at which Dmap varies. So a resolution of 10 means for every typical time scale there will be 10 interpolation points. Typical values are 1-100. A greater resolution increases interpolation accuracy, but also increases memory requirements and adds runtime (scales with the tree's age span, not Ntips).
						max_model_runtime		= NULL,			# maximum time (in seconds) to allocate for each evaluation of a model. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If NULL or <=0, this option is ignored.
						verbose					= TRUE,			# boolean, specifying whether to print informative messages
						diagnostics				= FALSE,		# boolean, specifying whether to print detailed info (such as log-likelihood) at every iteration of the fitting. For debugging purposes mainly.
						verbose_prefix			= ""){			# string, specifying the line prefix when printing messages. Only relevant if verbose==TRUE.
    Ntips 					= length(tree$tip.label)
    Nnodes 					= tree$Nnode
    Nedges 					= nrow(tree$edge)
	loglikelihood 			= NULL # value will be calculated as we go
	ancestral_likelihoods 	= NULL # values will be calculated as we go, if needed
	return_value_on_failure = list(success=FALSE, Nstates = NULL, loglikelihood=NULL, parameters=NULL)
	has_sampling_fractions  = ((!is.null(sampling_fractions)) && (length(sampling_fractions)>0))
	has_reveal_fractions  	= ((!is.null(reveal_fractions)) && (length(reveal_fractions)>1) && (!all(reveal_fractions<=0)) && (length(unique(reveal_fractions))>1))
	has_sampling_rates  	= ((!is.null(sampling_rates)) && (length(sampling_rates)>0))
	birth_rate_indices		= get_rate_index_vector(Nstates, birth_rate_model)$index_vector
	death_rate_indices		= get_rate_index_vector(Nstates, death_rate_model)$index_vector
	is_hisse_model			= !(is.null(NPstates) || (NPstates==0) || (NPstates==Nstates))
	root_age 				= get_tree_span(tree)$max_distance
	if(!is_hisse_model) NPstates = Nstates
	if(is.null(oldest_age)) oldest_age = root_age
			
	# check validity of input variables, and determine known tips
	if(verbose) cat(sprintf("%sChecking input variables..\n",verbose_prefix))
	if(!(optim_algorithm %in% c("subplex", "optim", "nlminb"))) stop(sprintf("ERROR: Invalid optimization algorithm '%s'",optim_algorithm))
	if((!is.null(tip_pstates)) && (!is.null(tip_priors))) stop("ERROR: tip_pstates and tip_priors are both non-NULL, but exactly one of them should be NULL")
	else if(is.null(tip_pstates) && is.null(tip_priors))  stop("ERROR: tip_pstates and tip_priors are both NULL, but exactly one of them should be non-NULL")
	if(! (root_conditioning %in% c("none","madfitz","herr_als","stem","crown","auto"))) stop("ERROR: root_conditioning must be one of 'none', 'madfitz', 'herr_als', 'stem', 'crown' or 'auto'")
	if(root_conditioning=="auto") root_conditioning = (if(abs(oldest_age-root_age)<=1e-10*root_age) "crown" else "stem");
	if(is.null(Nstates) || is.na(Nstates) || (Nstates==0)) stop("ERROR: Nstates must be a strictly positive integer")
	if((!is.null(state_names)) && (length(state_names)!=Nstates)) stop(sprintf("ERROR: Number of provided state_names (%d) differs from Nstates (%d)",length(state_names),Nstates))
	if(is_hisse_model && is.null(proxy_map)) stop("ERROR: Missing proxy_map, needed for HiSSE model")
	if(is_hisse_model && (length(proxy_map)!=Nstates)) stop("ERROR: proxy_map has length %d, but should have length %d (Nstates)",length(proxy_map),Nstates)
	if(is_hisse_model && (length(unique(proxy_map))!=NPstates)) stop("ERROR: Not all %d proxy states are represented in proxy_map",NPstates)
	if((!is_hisse_model) && (!is.null(proxy_map)) & ((length(proxy_map)!=Nstates) || (any(proxy_map!=(1:Nstates))))) stop("ERROR: Non-trivial proxy_map contradicts non-HiSSE model")
	if(!is_hisse_model) proxy_map = c(1:Nstates)
	if(is.character(transition_rate_model)){
		if(!(transition_rate_model %in% c("ER", "SYM", "ARD", "SUEDE", "SRD"))) stop(sprintf("ERROR: Unknown transition_rate_model '%s'",transition_rate_model))
		transition_indices = get_transition_index_matrix(Nstates,transition_rate_model)$index_matrix
	}else{
		# transition_rate_model seems to be (or should be) a 2D matrix of size Nstates x Nstates
		if((nrow(transition_rate_model)!=Nstates) || (ncol(transition_rate_model)!=Nstates)) stop(sprintf("ERROR: transition_rate_mode, which was provided as 2D matrix, must have size %d x %d, but instead got %d x %d",Nstates,Nstates,nrow(transition_rate_model),ncol(transition_rate_model)))
		transition_indices = transition_rate_model
	}
	if(!is.null(transition_matrix)){
		if((nrow(transition_matrix)!=Nstates) || (ncol(transition_matrix)!=Nstates)) stop(sprintf("ERROR: transition_matrix must have size %d x %d, but instead got %d x %d",Nstates,Nstates,nrow(transition_matrix),ncol(transition_matrix)))
		if(check_input){
			# make sure entries in transition_matrix are consistent with the transition_rate_model, i.e. check that whenever transition_indices[i,j]==transition_indices[k,l] for some i,j,k,l one also has transition_matrix[i,j]=transition_matrix[k,l]
			transition_indices_flat = as.vector(transition_indices)
			transition_matrix_flat 	= as.vector(transition_matrix)
			unique_indices 			= unique(transition_indices_flat)
			for(index in unique_indices){
				if(index==0) next;
				unique_rates = unique(transition_matrix_flat[transition_indices_flat==index])
				if(length(unique_rates)>1) stop(sprintf("ERROR: Entries in transition_matrix are inconsistent with transition_rate_model: Found inconsistent rates (e.g. %g and %g) for the same rate index %d",unique_rates[1],unique_rates[2],index))
			}
		}
	}
	max_start_attempts 	= (if(is.null(max_start_attempts)) 1 else max(1,max_start_attempts))
	if((Nbootstraps>0) && (!is.null(Ntrials_per_bootstrap)) && (Ntrials_per_bootstrap<=0)) stop(sprintf("ERROR: Ntrials_per_bootstrap must be strictly positive, if bootstrapping is requested"))
	if((!is.null(birth_rates)) && (length(birth_rates)!=1) && (length(birth_rates)!=Nstates)) stop(sprintf("ERROR: Invalid number of birth-rates; got %d, but expected 1 or %d (Nstates)",length(birth_rates),Nstates))
	if((!is.null(death_rates)) && (length(death_rates)!=1) && (length(death_rates)!=Nstates)) stop(sprintf("ERROR: Invalid number of death-rates; got %d, but expected 1 or %d (Nstates)",length(death_rates),Nstates))
	if(is.null(max_model_runtime)) max_model_runtime = 0
	if(!is.null(tip_pstates)){
		if(!is.numeric(tip_pstates)) stop(sprintf("ERROR: tip_pstates must be integers"))
		if(length(tip_pstates)==0) stop("ERROR: tip_pstates is non-NULL but empty")
		if(length(tip_pstates)!=Ntips) stop(sprintf("ERROR: Length of tip_pstates (%d) is not the same as the number of tips in the tree (%d)",length(tip_pstates),Ntips));
		known_tips  = which(!is.na(tip_pstates));
		Nknown_tips = length(known_tips)
		if(Nknown_tips==0) stop("ERROR: All tip states are unspecified");
		max_tip_state = max(tip_pstates[known_tips])
		if(check_input){
			min_tip_state = min(tip_pstates[known_tips])
			if((min_tip_state<1) || (max_tip_state>NPstates)) stop(sprintf("ERROR: Known tip states must be integers between 1 and %d, but found values between %d and %d",NPstates,min_tip_state,max_tip_state))
			if((!is.null(names(tip_pstates))) && any(names(tip_pstates)!=tree$tip.label)) stop("ERROR: Names in tip_pstates and tip labels in tree don't match (must be in the same order).")
		}
		unknown_tips = get_complement(Ntips, known_tips)
		# count number of tips in tree (known to be) in each state
		Nknown_tips_per_pstate = sapply(1:NPstates, function(state) sum(tip_pstates[known_tips]==state));
		if(verbose && (Nknown_tips<Ntips)) cat(sprintf("%sNote: %d out of %d tips have unspecified state\n",verbose_prefix,Ntips-Nknown_tips,Ntips))

	}else{
		if(nrow(tip_priors)==0) stop("ERROR: tip_priors is non-NULL but has zero rows")
		if((Nstates != ncol(tip_priors)) && (NPstates != ncol(tip_priors))) stop(sprintf("ERROR: The number of columns in tip_priors (%d) must be equal to Nstates (%d) or NPstates (%d)",ncol(tip_priors),Nstates,NPstates))
		if(Ntips != nrow(tip_priors)) stop(sprintf("ERROR: Number of tips (%d) differs from the number of rows in tip_priors (%d)",Ntips,nrow(tip_priors)))
		known_priors = which(rowSums(is.na(tip_priors))==0);
		if(length(known_priors)==0) stop("ERROR: All tip priors are unspecified");
		if(check_input){
			max_tip_prior = max(tip_priors[known_priors,])
			if(max_tip_prior>1.0) stop(sprintf("ERROR: Some tip_priors are larger than 1.0 (max was %g)",max_tip_prior))
			if((!is.null(rownames(tip_priors))) && (!is.null(tree$tip.label)) && (rownames(tip_priors)!=tree$tip.label)) stop("ERROR: Row names in tip_priors and tip labels in tree don't match")
		}
		Nknown_priors  = length(known_priors)
		unknown_priors = get_complement(Ntips, known_priors)
		if(verbose && (Nknown_priors<Ntips)) cat(sprintf("%sNote: %d out of %d tips have unspecified priors\n",verbose_prefix,Ntips-Nknown_priors,Ntips))
		# tip states are not explicitly given, but instead we got prior state likelihoods, so calculate number of known tips per pstate using a probabilistic averaging approach
		# only tips with specified priors are counted; sum(Nknown_tips_per_pstate) will be equal to Nknown_priors
		known_tip_prior_row_sums = rowSums(tip_priors[known_priors,])
		if(check_input && any(known_tip_prior_row_sums==0)) stop(sprintf("ERROR: Some specified tip priors are all zero (e.g., for tip %d)",known_priors[which(known_tip_prior_row_sums==0)[1]]))
		if(ncol(tip_priors)==NPstates){
			Nknown_tips_per_pstate = colSums(tip_priors[known_priors,]/known_tip_prior_row_sums)
		}else{
			Nknown_tips_per_state  = colSums(tip_priors[known_priors,]/known_tip_prior_row_sums)
			Nknown_tips_per_pstate = sapply(1:NPstates, FUN=function(p) sum(Nknown_tips_per_state[proxy_map==p]))
		}
	}
	Nstates_per_pstate = sapply(1:NPstates, FUN=function(p) sum(proxy_map==p)) # number of states per proxy-state
	states_per_pstate  = lapply(1:NPstates, FUN=function(p) which(proxy_map==p)) # states_per_pstate[p] is a 1D vector, listing the states represented by the proxy-state p
		
	# pre-process sampling_fractions & reveal_fractions & sampling_rates
	if(!has_sampling_fractions){
		sampling_fractions = rep(1, times=NPstates) # sampling fractions missing, so assume that the phylogeny is complete
	}else if(length(sampling_fractions)==1){
		sampling_fractions = rep(sampling_fractions[1], times=NPstates); # assume the same sampling fraction for all proxy-states
	}else if(length(sampling_fractions)!=NPstates){
		stop(sprintf("ERROR: Expected 1 or %d sampling fractions, but got %d",NPstates,length(sampling_fractions)))
	}
	if(has_reveal_fractions){
		if(length(reveal_fractions)!=NPstates) stop(sprintf("ERROR: Expected %d reveal fractions, but got %d",NPstates,length(reveal_fractions)))
		if(any((Nknown_tips_per_pstate>0) & (reveal_fractions==0))) stop("ERROR: Reveal fractions for some states are zero, despite the presence of tips with that known state")
		# rescale reveal fractions to be consistent with the number of known tips in each state within the considered tree
		reveal_fractions = reveal_fractions*sum(Nknown_tips_per_pstate[Nknown_tips_per_pstate>0]/reveal_fractions[Nknown_tips_per_pstate>0])/Ntips
	}else{
		# reveal fractions missing, so assume tips are revealed at the same probabilities (i.e. regardless of state)
		if(!is.null(tip_pstates)){
			reveal_fractions = rep(Nknown_tips/Ntips, times=NPstates);
		}else{
			reveal_fractions = rep(Nknown_priors/Ntips, times=NPstates);
		}
	}
	if(!has_sampling_rates){
		sampling_rates = rep(0, times=NPstates) # sampling rates missing, so assume that Poissonian sampling did not occur
	}else if(length(sampling_rates)==1){
		sampling_rates = rep(sampling_rates[1], times=NPstates); # assume the same sampling rates for all proxy-states
	}else if(length(sampling_rates)!=NPstates){
		stop(sprintf("ERROR: Expected 1 or %d sampling rates, but got %d",NPstates,length(sampling_rates)))
	}
	if(all(sampling_fractions==0) && all(sampling_rates==0)) return(list(success=FALSE, error="At least one of sampling_fractions[] or sampling_rates[] must contain non-zeros"))
	
	# determine Poissonian-sampled tips and extant tips (i.e. sampled at present-day)
	# also determine number of extant tips per pstate (i.e., omitting Poissonian sampled tips)
	if(all(sampling_rates==0)){
		# no Poissonian sampling
		Psampled_tips 	= integer(0)
		extant_tips		= seq_len(Ntips)
		NextantKnown_tips_per_pstate = Nknown_tips_per_pstate
	}else if(all(sampling_fractions==0)){
		# no present-day sampling
		Psampled_tips 	= seq_len(Ntips)
		extant_tips		= integer(0)
		NextantKnown_tips_per_pstate = rep(0, times=NPstates)
	}else{
		# tips are probably a mix of Poissonian-sampled and present-day-sampled
		tree_events 	= extract_HBDS_events_from_tree(tree, root_age = root_age, CSA_ages = 0, age_epsilon = root_age/1000)
		Psampled_tips	= tree_events$concentrated_tips[[1]]
		extant_tips		= get_complement(Ntips, Psampled_tips)
		if(verbose) cat(sprintf("%sNote: Found %d Poissonian-sampled tips and %d present-day tips\n",verbose_prefix,length(Psampled_tips),length(extant_tips)))
		if(!is.null(tip_pstates)){
			known_extant_tips = get_intersection(Ntips,known_tips,extant_tips)
			NextantKnown_tips_per_pstate = sapply(1:NPstates, function(p) sum(tip_pstates[known_extant_tips]==p));
		}else{
			known_extant_priors = get_intersection(Ntips,known_priors,extant_tips)
			if(ncol(tip_priors)==NPstates){
				NextantKnown_tips_per_pstate = colSums(tip_priors[known_extant_priors,]/rowSums(tip_priors[known_extant_priors,]))
			}else{
				NextantKnown_tips_per_state  = colSums(tip_priors[known_extant_priors,]/rowSums(tip_priors[known_extant_priors,]))
				NextantKnown_tips_per_pstate = sapply(1:NPstates, FUN=function(p) sum(NextantKnown_tips_per_state[proxy_map==p]))
			}

		}
	}
	
	
	# Estimate proportions of species per pstate
	# Also calculate overall present-day rarefaction of the tree, i.e. fraction of extant species included in the tree
	if(sum(NextantKnown_tips_per_pstate)>0){
		# Some tips are extant, i.e. sampled at present-day. Use these to estimate species proportions per pstate
		#   For any given proxy state, the number of tips in the tree with that proxy state is: NextantKnown_tips_per_pstate[state]/reveal_fractions[state]
		#   Hence, for any given proxy state, the number of extant species with that proxy state is: NextantKnown_tips_per_pstate[state]/(reveal_fractions[state]*sampling_fractions[state])
		NextantSpecies_per_pstate = NextantKnown_tips_per_pstate; 
		NextantSpecies_per_pstate[NextantKnown_tips_per_pstate>0] = NextantKnown_tips_per_pstate[NextantKnown_tips_per_pstate>0]/(reveal_fractions[NextantKnown_tips_per_pstate>0]*sampling_fractions[NextantKnown_tips_per_pstate>0]);
		overall_rarefaction = length(extant_tips)/sum(NextantSpecies_per_pstate)	
		species_proportions_per_pstate = NextantSpecies_per_pstate/sum(NextantSpecies_per_pstate)
	}else{
		# no extant known tips available, i.e., all known tips are Poissonian-sampled; so use Poisson-sampled known tips instead
		overall_rarefaction = 1 # no real information available for estimating this, since there's no extant known tips in the tree
		Nspecies_per_pstate = Nknown_tips_per_pstate; 
		Nspecies_per_pstate[Nknown_tips_per_pstate>0] = Nknown_tips_per_pstate[Nknown_tips_per_pstate>0]/(reveal_fractions[Nknown_tips_per_pstate>0]*sampling_rates[Nknown_tips_per_pstate>0]);
		species_proportions_per_pstate = Nspecies_per_pstate/sum(Nspecies_per_pstate)
	}

	
	# if some model parameters were passed as NaN or NULL, replace with NAs for consistency
	# also, if some birth/death/sampling rates were passed as single scalars, multiply to Nstates
	if((!is.null(birth_rates)) && (length(birth_rates)==1)) birth_rates = rep(birth_rates,Nstates)
	if((!is.null(death_rates)) && (length(death_rates)==1)) death_rates = rep(death_rates,Nstates)
	if(!is.null(transition_matrix)) transition_matrix[is.nan(transition_matrix)] = NA
	else transition_matrix = matrix(rep(NA,Nstates*Nstates),nrow=Nstates);
	if(!is.null(birth_rates))	birth_rates[is.nan(birth_rates)] = NA
	else birth_rates = rep(NA, Nstates)
	if(!is.null(death_rates))	death_rates[is.nan(death_rates)] = NA
	else death_rates = rep(NA, Nstates)
	
    # figure out probability distribution for root if needed
    if(root_prior[1]=="auto"){
    	if(oldest_age>=(1+1e-10)*root_age){
    		root_prior = "likelihoods"
    	}else{
    		root_prior = "max_likelihood"
    	}
    }
    if(root_prior[1]=="flat"){
    	root_prior_type = "custom"
    	root_prior_probabilities = rep(1.0/Nstates, times=Nstates);
    }else if(root_prior[1]=="empirical"){
    	# if is_hisse_model, then the empirical probability of each proxy-state is subdivided equally among all states represented by that proxy state
    	root_prior_type = "custom"
		root_prior_probabilities = (species_proportions_per_pstate/Nstates_per_pstate)[proxy_map];
	}else if(root_prior[1]=="likelihoods"){
		# root prior will be set to the computed state-likelihoods at the root, and thus depend on the particular model parameters. 
    	root_prior_type = "likelihoods"
		root_prior_probabilities = numeric(0)
	}else if(root_prior[1]=="max_likelihood"){
		# root prior will be set to a Dirac distribution, localized at the state with highest likelihood
    	root_prior_type = "max_likelihood"
		root_prior_probabilities = numeric(0)
	}else if((length(root_prior)==1) && (is.character(root_prior[[1]]))){
		stop(sprintf("ERROR: Unknown root_prior '%s'",root_prior[[1]]))
	}else{
		# basic error checking
		if(length(root_prior)!=Nstates) stop(sprintf("ERROR: root_prior has length %d, expected %d (Nstates)",length(root_prior),Nstates))
		if(check_input){
			if(any(root_prior<0)) stop(sprintf("ERROR: root_prior contains negative values (down to %g)",min(root_prior)))
			if(abs(1.0-sum(root_prior))>1e-6) stop(sprintf("ERROR: Entries in root prior do not sum up to 1 (sum=%.10g)",sum(root_prior)))
		}
    	root_prior_type = "custom"
		root_prior_probabilities = root_prior
	}
	
	# define initial conditions for E & D (tip_priors), depending on sampling fractions & reveal fractions
	if(is.null(tip_priors)){
		if(verbose) cat(sprintf("%sPreparing tip priors based on tip pstates..\n",verbose_prefix))
		tip_priors = matrix(0, nrow=Ntips, ncol=Nstates);
		if(Nknown_tips>0){
			if(length(Psampled_tips)>0){
				# for each known Psampled tip i, with proxy-state p=tip_pstates[i], set: tip_priors[i,states_per_pstate[p]] = sampling_rates[p] * reveal_fractions[p]
				known_Psampled_tips				= get_intersection(Ntips,known_tips,Psampled_tips)
				known_Psampled_tip_pstates 	 	= tip_pstates[known_Psampled_tips]
				Nstates_per_known_Psampled_tip 	= Nstates_per_pstate[known_Psampled_tip_pstates]
				tip_priors[cbind(rep(known_Psampled_tips,times=Nstates_per_known_Psampled_tip), unlist(states_per_pstate[known_Psampled_tip_pstates]))] = rep(sampling_rates[known_Psampled_tip_pstates] * reveal_fractions[known_Psampled_tip_pstates], times=Nstates_per_known_Psampled_tip)
			}
			if(length(extant_tips)>0){
				# for each known extant tip i, with proxy-state p=tip_pstates[i], set: tip_priors[i,states_per_pstate[p]] = sampling_fractions[p] * reveal_fractions[p]
				known_extant_tips				= get_intersection(Ntips,known_tips,extant_tips)
				known_extant_tip_pstates 	 	= tip_pstates[known_extant_tips]
				Nstates_per_known_extant_tip 	= Nstates_per_pstate[known_extant_tip_pstates]
				tip_priors[cbind(rep(known_extant_tips,times=Nstates_per_known_extant_tip), unlist(states_per_pstate[known_extant_tip_pstates]))] = rep(sampling_fractions[known_extant_tip_pstates] * reveal_fractions[known_extant_tip_pstates], times=Nstates_per_known_extant_tip)
			}
		}
		if(length(unknown_tips)>0){
			if(length(Psampled_tips)>0){
				unknown_Psampled_tips 				= get_intersection(Ntips,unknown_tips,Psampled_tips)
				likelihood_per_pstate 				= sampling_rates*(1-reveal_fractions)
				tip_priors[unknown_Psampled_tips,] 	= matrix(rep(likelihood_per_pstate[proxy_map],times=length(unknown_Psampled_tips)),nrow=length(unknown_Psampled_tips),ncol=Nstates,byrow=TRUE)
			}
			if(length(extant_tips)>0){
				unknown_extant_tips				 = get_intersection(Ntips,unknown_tips,extant_tips)
				likelihood_per_pstate 			 = sampling_fractions*(1-reveal_fractions)
				tip_priors[unknown_extant_tips,] = matrix(rep(likelihood_per_pstate[proxy_map],times=length(unknown_extant_tips)),nrow=length(unknown_extant_tips),ncol=Nstates,byrow=TRUE)
			}
		}
	}else{
		if(length(unknown_priors)>0){
			if(verbose) cat(sprintf("%sSetting unspecified priors to default..\n",verbose_prefix))
			# some tips have unspecified priors, so set to default
			if(length(Psampled_tips)>0){
				unknown_Psampled_priors					= get_intersection(Ntips,unknown_priors,Psampled_tips)
				likelihood_per_pstate					= sampling_rates*(1-reveal_fractions)
				tip_priors[unknown_Psampled_priors,] 	= matrix(rep(likelihood_per_pstate[proxy_map],times=length(unknown_Psampled_priors)),nrow=length(unknown_Psampled_priors),ncol=NPstates,byrow=TRUE)	
			}
			if(length(extant_tips)>0){
				unknown_extant_priors 				= get_intersection(Ntips,unknown_priors,extant_tips)
				likelihood_per_pstate 				= sampling_fractions*(1-reveal_fractions)
				tip_priors[unknown_extant_priors,] 	= matrix(rep(likelihood_per_pstate[proxy_map],times=length(unknown_extant_priors)),nrow=length(unknown_extant_priors),ncol=NPstates,byrow=TRUE)	
			}
		}
		if(ncol(tip_priors)==NPstates){
			# tip_priors are given as Ntips x NPstates, so reformat to make sure they are of dimensions Ntips x Nstates
			# basically, for every state s, set new_tip_priors[,s] = tip_priors[,proxy_map[s]]
			tip_priors = tip_priors[,proxy_map];
		}
	}
	initial_extinction_probabilities = (1-sampling_fractions[proxy_map]);

	# figure out which parameters are fitted vs fixed
	# get set of independent parameters (as a condensed vector)
	# For example, if the birth_rate_model is "ER", only one birth rate is considered an independent model parameter (regardless of whether it is fixed or fitted)
	provided_param_values 	= compress_musse_params(transition_matrix, birth_rates, death_rates, Nstates, transition_indices, birth_rate_indices, death_rate_indices);
	fitted_params			= which(is.na(provided_param_values))
	fixed_params			= which(!is.na(provided_param_values))
	NIP						= length(provided_param_values) # total number of independent model parameters (fixed+fitted)
	NFP						= length(fitted_params) # number of fitted (non-fixed) parameters
	NP						= Nstates*Nstates + 2*Nstates # total number of parameters (including redundancies): transition-rates, birth-rates, death-rates
	
	# determine lower & upper bounds for parameters
	if(is.null(lower)) lower = list();
	if(is.null(upper)) upper = list();
	param_mins = compress_musse_params(lower$transition_matrix, lower$birth_rates, lower$death_rates, Nstates, transition_indices, birth_rate_indices, death_rate_indices)
	param_mins[is.na(param_mins) | is.nan(param_mins) | (param_mins<0)] = 0;
	param_maxs = compress_musse_params(upper$transition_matrix, upper$birth_rates, upper$death_rates, Nstates, transition_indices, birth_rate_indices, death_rate_indices)
	param_maxs[is.na(param_maxs) | is.nan(param_maxs)] = Inf
	
	# determine optimization options if not already provided
	if(is.null(optim_max_iterations) || is.na(optim_max_iterations)) optim_max_iterations = 10000
	if(is.null(optim_max_evaluations) || is.na(optim_max_evaluations) || is.nan(optim_max_evaluations)) optim_max_evaluations = optim_max_iterations*NFP*10

	# determine start parameters
	if(is.null(first_guess)) first_guess = list();
	if((!is.null(first_guess$birth_rates)) && (length(first_guess$birth_rates)==1)) first_guess$birth_rates = rep(first_guess$birth_rates,Nstates)
	if((!is.null(first_guess$death_rates)) && (length(first_guess$death_rates)==1)) first_guess$death_rates = rep(first_guess$death_rates,Nstates)
	first_guess_compr = compress_musse_params(first_guess$transition_matrix, first_guess$birth_rates, first_guess$death_rates, Nstates, transition_indices, birth_rate_indices, death_rate_indices)
	first_guess_compr[!is.na(provided_param_values)] = provided_param_values[!is.na(provided_param_values)] # incorporate fixed values into start values
	first_guess_compr = pmin(pmax(first_guess_compr, param_mins), param_maxs) # make sure provisionary start params are within bounds
	first_guess = uncompress_musse_params(first_guess_compr, Nstates, transition_indices, birth_rate_indices, death_rate_indices, NULL)
	if(all(is.na(first_guess$birth_rates)) || all(is.na(first_guess$death_rates))){
		# Either all birth_rates and/or all death_rates are non-fixed and have unknown start values, so guesstimate by fitting a birth-death-sampling model
		if(verbose) cat(sprintf("%sGuesstimating start speciation/extinction rates..\n",verbose_prefix))
		if(!all(is.na(first_guess$birth_rates))) hbds_guess_lambda = mean(first_guess$birth_rates, na.rm=TRUE)
		else hbds_guess_lambda = NULL
		if(!all(is.na(first_guess$death_rates))) hbds_guess_mu = mean(first_guess$death_rates, na.rm=TRUE)
		else hbds_guess_mu = NULL
		if(!all(is.na(birth_rates))) hbds_fixed_lambda = mean(birth_rates, na.rm=TRUE)
		else hbds_fixed_lambda = NULL
		if(!all(is.na(death_rates))) hbds_fixed_mu = mean(death_rates, na.rm=TRUE)
		else hbds_fixed_mu = NULL
		
		if(all(sampling_rates==0)){
			# faster implementation when Poissonian sampling is zero
			fit = fit_tree_model(	tree, 
									parameters			= list(birth_rate_intercept=0, birth_rate_factor=hbds_fixed_lambda, birth_rate_exponent=1, death_rate_intercept=0, death_rate_factor=hbds_fixed_mu, death_rate_exponent=1, rarefaction=overall_rarefaction, resolution=0),
									first_guess 		= list(birth_rate_factor=hbds_guess_lambda, death_rate_factor=hbds_guess_mu),
									min_age				= 0,
									max_age	 			= 0,
									age_centile			= NULL,
									Ntrials				= max(1,Ntrials/5),
									Nthreads			= Nthreads,
									coalescent			= TRUE,
									discovery_fraction	= NULL,
									fit_control			= list(iter.max=optim_max_iterations, rel.tol=optim_rel_tol),
									min_R2				= 0.5,
									min_wR2				= 0.5,
									grid_size			= 100,
									max_model_runtime	= max_model_runtime,
									objective			= 'MRD');
			if(fit$success){
				# fitting succeeded, so use fitted speciation/extinction parameters
				first_guess$birth_rates[is.na(first_guess$birth_rates)] = fit$parameters$birth_rate_factor
				first_guess$death_rates[is.na(first_guess$death_rates)] = fit$parameters$death_rate_factor
			}else{
				# fitting failed, so use the birth-death-model start params as a last resort (these are probably reasonable by order of magnitude)
				first_guess$birth_rates[is.na(first_guess$birth_rates)] = fit$start_parameters$birth_rate_factor
				first_guess$death_rates[is.na(first_guess$death_rates)] = fit$start_parameters$death_rate_factor		
			}
		}else{
			fit = fit_hbds_model_on_grid(	tree, 
											oldest_age				= oldest_age,
											age_grid				= NULL,	# fit a constant-rates model, i.e. with a single grid point
											CSA_ages				= 0,
											guess_lambda			= hbds_guess_lambda,
											guess_mu				= hbds_guess_mu,
											fixed_lambda			= hbds_fixed_lambda,
											fixed_mu				= hbds_fixed_mu,
											fixed_psi				= mean(sampling_rates),
											fixed_kappa				= 0,
											fixed_CSA_probs			= overall_rarefaction,
											fixed_CSA_kappas		= 0,
											condition				= "auto",
											ODE_relative_dt			= 0.001,
											ODE_relative_dy			= 1e-3,
											Ntrials					= max(1,Ntrials/5),
											max_start_attempts		= 10,
											Nthreads				= Nthreads,
											max_model_runtime		= max_model_runtime,
											verbose					= FALSE,
											diagnostics				= diagnostics,
											verbose_prefix			= paste0(verbose_prefix,"  "))
			if(fit$success){
				# fitting succeeded, so use fitted speciation/extinction parameters
				first_guess$birth_rates[is.na(first_guess$birth_rates)] = fit$param_fitted$lambda
				first_guess$death_rates[is.na(first_guess$death_rates)] = fit$param_fitted$mu
			}else{
				# fitting failed, so use the birth-death-sampling-model start params as a last resort (these are probably reasonable by order of magnitude)
				first_guess$birth_rates[is.na(first_guess$birth_rates)] = fit$param_guess$lambda
				first_guess$death_rates[is.na(first_guess$death_rates)] = fit$param_guess$mu		
			}
		}
	}else{
		if(any(is.na(first_guess$birth_rates))){
			# some but not all guess birth rates are specified
			first_guess$birth_rates[is.na(first_guess$birth_rates)] = mean(first_guess$birth_rates, na.rm=TRUE)
		}
		if(any(is.na(first_guess$death_rates))){
			# some but not all guess death rates are specified
			first_guess$death_rates[is.na(first_guess$death_rates)] = mean(first_guess$death_rates, na.rm=TRUE)
		}
	}
	if(any(is.na(first_guess$transition_matrix))){
		# some transition rates are non-fixed and have unknown start values, so guesstimate by fitting an Mk model
		if(verbose) cat(sprintf("%sGuesstimating transition rate scale based on independent contrasts..\n",verbose_prefix))
		guessQ = guesstimate_Mk_rate_via_independent_contrasts(tree, Nstates=Nstates, tip_states=NULL, tip_priors=tip_priors, allow_ties=TRUE)
		if(guessQ$success){
			guessQ$Q = matrix(guessQ$rate, nrow=Nstates, ncol=Nstates)
		}else{
			if(verbose) cat(sprintf("%s  WARNING: Independent-contrasts guesstimate of Q failed: %s\n",verbose_prefix,guessQ$error))
		}
		if(!guessQ$success){
			if(verbose) cat(sprintf("%sIndependent-contrasts guesstimate of Q failed, so guesstimating via max-parsimony..\n",verbose_prefix))
			guessQ = guesstimate_Mk_transition_rates_via_max_parsimony_ASR(tree, tip_states=NULL, tip_priors=tip_priors, Nstates=Nstates, transition_costs = "all_equal", allow_ties=TRUE)
			if((!guessQ$success) && verbose) cat(sprintf("%s  WARNING: Max-parsimony guesstimate of Q failed: %s\n",verbose_prefix,guessQ$error))
		}
		if(!guessQ$success){
			if(verbose) cat(sprintf("%sMax-parsimony guesstimate of Q failed, so fitting Mk model..\n",verbose_prefix))
			fit = fit_mk(	tree,
							Nstates					= Nstates,
							tip_states				= NULL,
							tip_priors 				= tip_priors,
							rate_model 				= transition_rate_model,
							root_prior 				= (if((root_prior=="flat") || (root_prior=="empirical")) root_prior else "flat"),
							oldest_ages				= oldest_age,
							Ntrials 				= max(1,Ntrials/5),
							max_model_runtime		= max_model_runtime,
							check_input 			= TRUE,
							Nthreads 				= Nthreads,
							Nbootstraps				= 0,
							optim_max_iterations	= optim_max_iterations,
							optim_rel_tol			= optim_rel_tol,
							verbose					= FALSE,
							verbose_prefix			= paste0(verbose_prefix,"  "))
			if(fit$success){
				# fitting succeeded
				guessQ = list(success=TRUE, Q=fit$transition_matrix)
			}else{
				# fitting failed, so use a rough guesstimate based on the number of tips
				guessQ = list(success=TRUE, rate = Nstates/((if(is.null(tree$edge.length)) 1 else mean(tree$edge.length))*log(Nknown_tips)/log(2.0)))
				guessQ$Q = matrix(guessQ$rate, nrow=Nstates, ncol=Nstates)
			}
		}
		# adopt guesstimated transition rates. Note that the guesstimated Q may not be a fully valid transition rate matrix, for example the rowSums may not be 0
		first_guess$transition_matrix[is.na(first_guess$transition_matrix)] = guessQ$Q[is.na(first_guess$transition_matrix)]
		# make sure first-guess transition matrix is a valid transition matrix
		diag(first_guess$transition_matrix) = 0
		diag(first_guess$transition_matrix) = -rowSums(first_guess$transition_matrix)
	}
	first_guess_compr = compress_musse_params(first_guess$transition_matrix, first_guess$birth_rates, first_guess$death_rates, Nstates, transition_indices, birth_rate_indices, death_rate_indices)
	first_guess_compr = pmin(pmax(first_guess_compr, param_mins), param_maxs) # make sure finalized start params are within bounds
	first_guess 	  = uncompress_musse_params(first_guess_compr, Nstates, transition_indices, birth_rate_indices, death_rate_indices, state_names)
	return_value_on_failure$start_parameters = first_guess
	
	# at this point, all start parameters (first_guess/first_guess_compr) are assumed to be well defined (non-NA), either based on provided/fixed start values and/or based on some guesstimates
	# now proceed with the fitting trials

	if(verbose) cat(sprintf("%sPreparing fitting of full model..\n",verbose_prefix))

	# pre-calculate clade ages, used in loglikelihood function
	clade_ages = root_age - castor::get_all_distances_to_root(tree)
	node_ages  = clade_ages[(Ntips+1):(Ntips+Nnodes)]
			
	# determine typical parameter scales
	param_scales = lapply(first_guess,function(x) abs(x))
	if(all(param_scales$birth_rates==0)){
		# all birth-rate scales are zero, so determine scales based on a simple Yule model
		param_scales$birth_rates[] = log(Ntips)/root_age;
	}else{
		# some (but not all) birth-rate scales are zero, so determine scale based on mean birth-rate scale
		param_scales$birth_rates[param_scales$birth_rates==0] = mean(param_scales$birth_rates)
	}
	if(all(param_scales$death_rates==0)){
		# all death-rate scales are zero, so determine scale based on birth-rates
		param_scales$death_rates[] = 0.1*mean(param_scales$birth_rates)
	}else{
		# some (but not all) death-rate scales are zero, so determine scale based on mean death-rate scale
		param_scales$death_rates[param_scales$death_rates==0] = mean(param_scales$death_rates)
	}
	if(all(param_scales$transition_matrix==0)){
		# all transition rate scales are zero, so set transition rates equal to mean birth-rate
		param_scales$transition_matrix = matrix(mean(param_scales$birth_rates), nrow=Nstates, ncol=Nstates)
	}else{
		# some (but not all) transition rate scales are zero, so determine scale based on mean absolute transition-rate scale
		param_scales$transition_matrix[param_scales$transition_matrix==0] = 0.1*mean(param_scales$transition_matrix)
	}
	param_scales = compress_musse_params(param_scales$transition_matrix, param_scales$birth_rates, param_scales$death_rates, Nstates, transition_indices, birth_rate_indices, death_rate_indices)
							
	# define objective function to be minimized (negated log-likelihood)
	# the input to the objective function must be scaled and free (non-fixed) independent parametes
	objective_function = function(fparam_values, trial){
		# # reverse parameter transformation. fparam_values are shifted with regards to their lower bounds, i.e. the true parameter values are: SQ(fparam_values)+param_mins
		if(any(is.nan(fparam_values)) || any(is.infinite(fparam_values))) return(if(optim_algorithm == "optim") 1e100 else Inf)
		fparam_values 	= unscale_params(fparam_values, param_scales[fitted_params], param_mins[fitted_params])
		if(any(fparam_values<param_mins[fitted_params]) || any(fparam_values>param_maxs[fitted_params])) return(Inf)
		param_values 	= provided_param_values; param_values[fitted_params] = fparam_values; # merge fixed & fitted parameter values
		param_values 	= uncompress_musse_params(param_values, Nstates, transition_indices, birth_rate_indices, death_rate_indices, NULL)
		results = tryCatch({ get_MuSSE_loglikelihood_CPP(	Ntips 							= Ntips,
												Nnodes							= Nnodes,
												Nedges							= Nedges,
												Nstates							= Nstates,
												oldest_age						= oldest_age,
												tree_edge 						= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
												clade_ages 						= clade_ages,
												transition_rates 				= as.vector(t(param_values$transition_matrix)),	# flatten in row-major format
												speciation_rates 				= param_values$birth_rates,
												extinction_rates 				= param_values$death_rates,
												sampling_rates 					= sampling_rates[proxy_map],
												initial_D_per_tip 				= as.vector(t(tip_priors)), 		# flatten in row-major format
												initial_E_per_state				= initial_extinction_probabilities,
												root_prior_type					= root_prior_type,
												root_prior 						= root_prior_probabilities,
												root_conditioning				= root_conditioning,
												include_ancestral_likelihoods 	= FALSE,
												include_warnings				= FALSE,
												max_condition_number			= max_condition_number,
												relative_ODE_step				= relative_ODE_step,
												E_value_step					= E_value_step,
												D_temporal_resolution			= D_temporal_resolution,
												runtime_out_seconds				= max_model_runtime)
							}, error = function(e){ list(loglikelihood=NaN, success=FALSE, error="Unknown runtime exception") })
		loglikelihood = if((!results$success) || is.na(results$loglikelihood) || is.nan(results$loglikelihood)) (if(optim_algorithm == "optim") -1e100 else -Inf) else results$loglikelihood
		if(diagnostics){
			if(results$success){ cat(sprintf("%s  Trial %d: loglikelihood %.10g, model runtime %.5g sec\n",verbose_prefix,trial,loglikelihood,results$runtime)) }
			else{ cat(sprintf("%s  Trial %d: Model evaluation failed: %s\n",verbose_prefix,trial,results$error)) }
		}
		return(-loglikelihood);
	}

	# fit starting with various start-params, keep track of best fit
	scaled_param_mins = scale_params(param_mins,param_scales,param_mins);
	scaled_param_maxs = scale_params(param_maxs,param_scales,param_mins);
	fit_single_trial = function(trial){
		# randomly choose start values for fitted params (keep trying up to max_start_attempts times)
		Nstart_attempts = 0
		while(Nstart_attempts<max_start_attempts){
			# randomly choose start values for fitted params
			if(trial==1){
				initial_param_values = first_guess_compr
			}else{
				initial_param_values = get_random_params(defaults=first_guess_compr, lower_bounds=param_mins, upper_bounds=param_maxs, scales=param_scales, orders_of_magnitude=max(1,min(6,log2(10.0)*sqrt(Ntrials))))
			}
			scaled_initial_fparam_values = scale_params(initial_param_values,param_scales,param_mins)[fitted_params]
			# check if start values yield NaN
			start_LL = -objective_function(fparam_values=scaled_initial_fparam_values, trial=0)
			Nstart_attempts = Nstart_attempts + 1
			if(is.finite(start_LL)) break;
		}
		# run fit with chose start values
		if(is.finite(start_LL)){
			if(optim_algorithm == "optim"){
				fit = stats::optim(	par		= scaled_initial_fparam_values, 
									fn		= function(pars){ objective_function(pars, trial) }, 
									method 	= "L-BFGS-B", 
									lower 	= scaled_param_mins[fitted_params],
									upper 	= scaled_param_maxs[fitted_params],
									control = list(maxit=optim_max_iterations, reltol=optim_rel_tol))
				LL 				= -fit$value
				Nevaluations 	= fit$counts
				Niterations		= NA
				converged 		= (fit$convergence==0)
			}else if(optim_algorithm == "nlminb"){
				fit = stats::nlminb(start		= scaled_initial_fparam_values, 
									objective	= function(pars){ objective_function(pars, trial) }, 
									lower 		= scaled_param_mins[fitted_params], 
									upper 		= scaled_param_maxs[fitted_params],
									control 	= list(iter.max = optim_max_iterations, eval.max = optim_max_evaluations, rel.tol = optim_rel_tol))
				LL 				= -fit$objective;
				Nevaluations 	= fit$evaluations[1]
				Niterations		= fit$iterations
				converged		= (fit$convergence==0)
			}else if(optim_algorithm == "subplex"){
				fit = nloptr::nloptr(	x0 		= scaled_initial_fparam_values, 
										eval_f 	= function(pars){ objective_function(pars, trial) },
										lb 		= scaled_param_mins[fitted_params], 
										ub 		= scaled_param_maxs[fitted_params], 
										opts 	= list(algorithm = "NLOPT_LN_SBPLX", maxeval = optim_max_evaluations, ftol_rel = optim_rel_tol))
				LL 				= -fit$objective
				Nevaluations 	= NA
				Niterations 	= (fit$iterations)
				converged 		= (fit$status %in% c(0,3,4))
				fit$par 		= fit$solution
				#fit = subplex::subplex(par = initial_param_values[fitted_params], 
				#						fn = function(pars){ objective_function(pars, trial) },
				#						control = list(maxit=optim_max_iterations, reltol = optim_rel_tol, parscale=0.1), 
				#						hessian = FALSE)
				#LL 				= -fit$value
				#Nevaluations 	= fit$count
				#converged 		= (fit$convergence>=0)
			}
			if(is.null(LL)){ LL = NaN; converged = FALSE; }
			if(diagnostics) cat(sprintf("%s  Trial %d: Final loglikelihood %.10g, converged = %d\n",verbose_prefix,trial,LL,converged))
			return(list(LL=LL, start_LL=start_LL, Nevaluations=Nevaluations, Niterations=Niterations, converged=converged, fit=fit, Nstart_attempts=Nstart_attempts))
		}else{
			if(diagnostics) cat(sprintf("%s  Trial %d: Start loglikelihood is non-finite. Skipping trial\n",verbose_prefix,trial))
			return(list(LL=NA, start_LL = NA, Niterations=0, Nevaluations=0, converged=FALSE, Nstart_attempts=Nstart_attempts))
		}
	}
	
	# run one or more independent fitting trials
	if(NFP==0){
		# all parameters are fixed, so no fitting needed
		if(verbose) cat(sprintf("%sAll parameters are fixed; simply calculating loglikelihood of model..\n",verbose_prefix))
		best_param_values 	= uncompress_musse_params(provided_param_values, Nstates, transition_indices, birth_rate_indices, death_rate_indices, state_names)
		Nevaluations 		= NA
		Niterations			= NA
		converged 			= TRUE
		point_name			= "fixed params"
		LLs					= rep(NA, Ntrials)

	}else if(Ntrials==0){
		# fitting not requested, so just evaluate once at start params
		if(verbose) cat(sprintf("%sRequested zero trials; the loglikelihood will only be evaluated once at the starting point..\n",verbose_prefix))
		best_param_values 	= first_guess
		Nevaluations 		= 1
		Niterations			= NA
		converged 			= FALSE
		point_name 			= "start params"
		LLs					= numeric()

	}else{
		if((Ntrials>1) && (Nthreads>1) && (.Platform$OS.type!="windows")){
			if(verbose) cat(sprintf("%sFitting %d model parameters (%d trials, parallelized)..\n",verbose_prefix,NFP,Ntrials))
			# run trials in parallel using multiple forks
			# Note: Forks (and hence shared memory) are not available on Windows
			fits = parallel::mclapply(	1:Ntrials, 
										FUN = function(trial) fit_single_trial(trial), 
										mc.cores = min(Nthreads, Ntrials), 
										mc.preschedule = FALSE, 
										mc.cleanup = TRUE);
		}else{
			# run in serial mode
			if(verbose) cat(sprintf("%sFitting %d model parameters (%s)..\n",verbose_prefix,NFP,(if(Ntrials==1) "1 trial" else sprintf("%d trials",Ntrials))))
			fits = vector("list", Ntrials)
			for(trial in 1:Ntrials){
				fits[[trial]] = fit_single_trial(trial)
			}
		}
	
		# extract information from best fit (note that some fits may have LL=NaN or NA)
		if(verbose && (Ntrials>1)) cat(sprintf("%sChoosing best fit among all %d trials..\n",verbose_prefix,Ntrials))
		LLs		= sapply(1:Ntrials, function(trial) (if(is.null(fits[[trial]]$LL)) NA else fits[[trial]]$LL))
		valids 	= which((!is.na(LLs)) & (!is.nan(LLs)) & (!is.infinite(LLs)) & (!is.null(LLs)) & sapply(1:Ntrials, function(trial) (!any(is.na(fits[[trial]]$fit$par))) && (!any(is.nan(fits[[trial]]$fit$par))) && (!any(is.null(fits[[trial]]$fit$par)))))
		if(length(valids)==0){
			return_value_on_failure$error = "Fitting failed for all trials"
			return(return_value_on_failure);
		}
		best				= valids[which.max(LLs[valids])]
		loglikelihood		= fits[[best]]$LL;
		best_fparams		= unscale_params(fits[[best]]$fit$par, param_scales[fitted_params], param_mins[fitted_params]) # reverse rescaling
		best_param_values 	= provided_param_values; best_param_values[fitted_params] = best_fparams; # merge fixed & fitted parameter values
		if(is.null(loglikelihood) || any(is.na(best_param_values)) || any(is.nan(best_param_values))){
			return_value_on_failure$error = "Fitting yielded NaN loglikelihood and/or rates"
			return(return_value_on_failure);
		}
		best_param_values 	= uncompress_musse_params(best_param_values, Nstates, transition_indices, birth_rate_indices, death_rate_indices, state_names)
		Nevaluations 		= fits[[best]]$Nevaluations
		Niterations 		= fits[[best]]$Niterations
		converged	 		= fits[[best]]$converged
		
		# keep record of these new results, in case we fail in the future
		return_value_on_failure$parameters 		= best_param_values
		return_value_on_failure$Nevaluations	= Nevaluations
		point_name = "final point (best fit)"
	}

	# perform one final evaluation of the model
	if(verbose) cat(sprintf("%sEvaluating model at %s..\n",verbose_prefix,point_name))
	final = get_MuSSE_loglikelihood_CPP(Ntips 							= Ntips,
										Nnodes							= Nnodes,
										Nedges							= Nedges,
										Nstates							= Nstates,
										oldest_age						= oldest_age,
										tree_edge 						= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based,
										clade_ages 						= clade_ages,
										transition_rates 				= as.vector(t(best_param_values$transition_matrix)),	# flatten in row-major format
										speciation_rates 				= best_param_values$birth_rates,
										extinction_rates 				= best_param_values$death_rates,
										sampling_rates 					= sampling_rates,
										initial_D_per_tip 				= as.vector(t(tip_priors)), 		# flatten in row-major format
										initial_E_per_state				= initial_extinction_probabilities,
										root_prior_type					= root_prior_type,
										root_prior 						= root_prior_probabilities,
										root_conditioning				= root_conditioning,
										include_ancestral_likelihoods 	= include_ancestral_likelihoods,
										include_warnings				= TRUE,
										max_condition_number			= max_condition_number,
										relative_ODE_step				= relative_ODE_step,
										E_value_step					= E_value_step,
										D_temporal_resolution			= D_temporal_resolution,
										runtime_out_seconds				= max_model_runtime);
	if(!final$success){
		return_value_on_failure$error = sprintf("Evaluation at %s failed: %s",point_name,final$error)
		return(return_value_on_failure)
	}
	loglikelihood = final$loglikelihood
	if(include_ancestral_likelihoods){
		ancestral_likelihoods = matrix(final$ancestral_likelihoods, ncol=Nstates, byrow=TRUE, dimnames=list(tree$tip.label,NULL)) # unflatten
		if(oldest_age<root_age) ancestral_likelihoods[node_ages>oldest_age,] = NA; # state likelihoods are not calculated for nodes older than oldest_age
	}
	final$ML_subroot_states	= final$ML_subroot_states + 1 # make index 1-based
	final$ML_substem_states	= final$ML_substem_states + 1 # make index 1-based
	final$subroots 			= final$subroots + 1
	
	
	# estimate confidence intervals if needed, using parametric bootstrapping
	if(Nbootstraps>0){
		if(verbose) cat(sprintf("%sEstimating confidence intervals using parametric bootstrapping..\n",verbose_prefix))
		if(is.null(Ntrials_per_bootstrap)) Ntrials_per_bootstrap = max(1,Ntrials)
		best_param_values_flat 	= flatten_musse_params(best_param_values, Nstates);
		bootstrap_param_values 	= matrix(NA,nrow=Nbootstraps,ncol=NP)
		NBsucceeded 			= 0
		for(b in 1:Nbootstraps){
			# simulate model with fitted parameters
			if(verbose) cat(sprintf("%s  Bootstrap #%d..\n",verbose_prefix,b))
			simulation = simulate_musse(Nstates, 
										NPstates 				= NPstates,
										proxy_map 				= proxy_map,
										parameters 				= best_param_values,
										max_extant_tips 		= (if(length(Psampled_tips)<=length(extant_tips)) Ntips/overall_rarefaction else NULL),
										max_Psampled_tips		= (if(length(Psampled_tips)>length(extant_tips)) length(Psampled_tips) else NULL),
										sampling_fractions		= sampling_fractions,
										reveal_fractions		= reveal_fractions,
										sampling_rates			= sampling_rates,
										coalescent 				= TRUE,
										no_full_extinction		= TRUE,	
										include_labels			= FALSE)
			# fit MuSSE model to simulated tree
			# use "true" parameters as first_guess
			fit = fit_musse(tree					= simulation$tree, 
							Nstates					= Nstates,
							NPstates				= NPstates,
							proxy_map				= proxy_map,
							tip_pstates				= (if(is_hisse_model) simulation$tip_proxy_states else simulation$tip_states),
							sampling_fractions		= sampling_fractions,
							reveal_fractions		= reveal_fractions,
							sampling_rates			= sampling_rates,
							transition_rate_model	= transition_rate_model,
							birth_rate_model		= birth_rate_model,
							death_rate_model		= death_rate_model,
							transition_matrix		= transition_matrix,
							birth_rates				= birth_rates,
							death_rates				= death_rates,
							first_guess				= best_param_values,
							lower					= lower,
							upper					= upper,
							root_prior 				= root_prior,
							root_conditioning		= root_conditioning,
							oldest_age				= oldest_age,
							Ntrials 				= Ntrials_per_bootstrap,
							max_start_attempts		= max_start_attempts,
							optim_algorithm 		= optim_algorithm,
							optim_max_iterations	= optim_max_iterations,
							optim_max_evaluations	= optim_max_evaluations,
							optim_rel_tol			= optim_rel_tol,
							check_input 			= FALSE,
							include_ancestral_likelihoods = FALSE,
							Nthreads 				= Nthreads,
							Nbootstraps				= 0,
							max_condition_number	= max_condition_number,
							relative_ODE_step		= relative_ODE_step,
							E_value_step			= E_value_step,
							D_temporal_resolution	= D_temporal_resolution,
							max_model_runtime		= max_model_runtime,
							verbose					= verbose,
							diagnostics				= diagnostics,
							verbose_prefix			= paste0(verbose_prefix,"    "))
			if(!fit$success){
				if(verbose) cat(sprintf("%s  WARNING: Fitting failed for this bootstrap\n",verbose_prefix))
			}else{
				bootstrap_param_values[b,] = flatten_musse_params(fit$parameters, Nstates)
				#means			= means + bootstrapped_param_values_flat
				#standard_errors = standard_errors + bootstrapped_param_values_flat^2
				NBsucceeded = NBsucceeded + 1
			}
		}
		# calculate standard errors and confidence intervals from distribution of bootstrapped parameters
		standard_errors_flat = sqrt(pmax(0, colMeans(bootstrap_param_values^2, na.rm=TRUE) - colMeans(bootstrap_param_values, na.rm=TRUE)^2))
		standard_errors = unflatten_musse_params(standard_errors_flat, Nstates, state_names)
		quantiles = sapply(1:NP, FUN=function(p) quantile(bootstrap_param_values[,p], probs=c(0.25, 0.75, 0.025, 0.975), na.rm=TRUE, type=8))
		CI50lower = unflatten_musse_params(quantiles[1,], Nstates, state_names)
		CI50upper = unflatten_musse_params(quantiles[2,], Nstates, state_names)
		CI95lower = unflatten_musse_params(quantiles[3,], Nstates, state_names)
		CI95upper = unflatten_musse_params(quantiles[4,], Nstates, state_names)
		CI = cbind(best_param_values_flat,standard_errors_flat,t(as.matrix(quantiles)))
		rownames(CI) = get_flat_musse_param_names(Nstates)
		colnames(CI) = c("ML_estimate", "standard_error", "CI50lower", "CI50upper", "CI95lower", "CI95upper")
	}
		
	return(list(success						= TRUE, 
				Nstates						= Nstates,
				NPstates					= NPstates,
				root_prior					= root_prior,
				parameters					= best_param_values, 
				start_parameters			= first_guess,
				loglikelihood				= loglikelihood,
				AIC							= 2*NFP - 2*loglikelihood,
				Niterations					= Niterations, # may be NA, depending on the optimization algorithm
				Nevaluations				= Nevaluations, # may be NA, depending on the optimization algorithm
				converged					= converged,
				warnings					= (if(length(final$warnings)>0) final$warnings else NULL),
				subroots					= final$subroots, # clade indices of subtree roots
				ML_subroot_states			= final$ML_subroot_states, # maximum-likelihood estimate of the state at each subtree's root, based on the computed state-likelihoods
				ML_substem_states			= final$ML_substem_states, # maximum-likelihood estimate of the state at each subtree's stem, based on the computed state-likelihoods
				trial_start_loglikelihoods	= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$start_LL)),
				trial_loglikelihoods		= LLs,
				trial_Nstart_attempts		= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Nstart_attempts)),
				trial_Niterations			= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Niterations)),
				trial_Nevaluations			= unlist_with_nulls(sapply(1:Ntrials, function(trial) fits[[trial]]$Nevaluations)),
				standard_errors				= (if(Nbootstraps>0) standard_errors else NULL),
				CI50lower					= (if(Nbootstraps>0) CI50lower else NULL),
				CI50upper					= (if(Nbootstraps>0) CI50upper else NULL),
				CI95lower					= (if(Nbootstraps>0) CI95lower else NULL),
				CI95upper					= (if(Nbootstraps>0) CI95upper else NULL),
				CI							= (if(Nbootstraps>0) CI else NULL), # 2D matrix of size NP x 4, listing 50% 95% confidence intervals
				ancestral_likelihoods 		= ancestral_likelihoods));
}




#########################################################
# AUXILIARY FUNCTIONS


# generate a compressed flat vector of independent model parameter values, with parameters listed in a specific order
# for example, if the birth-rate model is "ER" (as encoded by birth_rate_indices), then only one birth rate will be included in the flattened param vector
compress_musse_params = function(transition_matrix, birth_rates, death_rates, Nstates, transition_indices, birth_rate_indices, death_rate_indices){
	if(is.null(transition_matrix)){
		free_transition_rates = rep(NA, max(transition_indices))
	}else{
		if(length(transition_matrix)==1){
			# if only one transition rate is provided, assume it applies to all state-pairs
			transition_matrix = matrix(abs(transition_matrix[1]), nrow=Nstates, ncol=Nstates); 
			diag(transition_matrix) = -rowSums(transition_matrix); # make sure this is a valid transition rate matrix (row-sums = 0)
		}
		free_transition_rates = extract_independent_rates_from_transition_matrix(transition_matrix,transition_indices)
	}
	if(is.null(birth_rates)){
		free_birth_rates = rep(NA, max(birth_rate_indices))
	}else{
		if(length(birth_rates)==1) birth_rates = rep(birth_rates, Nstates); # if only one birth rate is provided, assume it applies to all states
		free_birth_rates = extract_independent_rates_from_vector(birth_rates, birth_rate_indices)
	}
	if(is.null(death_rates)){
		free_death_rates = rep(NA, max(death_rate_indices))
	}else{
		if(length(death_rates)==1) death_rates = rep(death_rates, Nstates); # if only one death rate is provided, assume it applies to all states
		free_death_rates = extract_independent_rates_from_vector(death_rates, death_rate_indices)
	}
	params 	= c(free_transition_rates, free_birth_rates, free_death_rates)
	params[is.nan(params)] = NA # replace NaN with NA for uniformity
	return(params)
}


# given a compressed vector of independent model parameters (e.g. generated by compress_musse_params()), return an unflattened named list of model parameters
# the returned list will have entries such as transition_matrix[], birth_rates[] and death_rates[], with the proper dimensionalities and potentially redundant values
# For example, if the birth_rate_model is "ER" (as encoded by birth_rate_indices), then the returned birth_rates will be a vector with Nstates identical entries
uncompress_musse_params = function(params, Nstates, transition_indices, birth_rate_indices, death_rate_indices, state_names){
	NFT = max(as.vector(transition_indices)) # number of free (fitted) transition rates
	NFB = max(birth_rate_indices) # number of free (fitted) birth rates
	NFD = max(death_rate_indices) # number of free (fitted) death rates
	
	transition_matrix = get_transition_matrix_from_rate_vector(params[1:NFT], transition_indices, Nstates)
	rownames(transition_matrix) = state_names
	colnames(transition_matrix) = state_names
	
	free_birth_rates = params[(NFT+1):(NFT+NFB)]
	birth_rates = setNames(free_birth_rates[birth_rate_indices],state_names)
	
	free_death_rates = params[(NFT+NFB+1):(NFT+NFB+NFD)]
	death_rates = setNames(free_death_rates[death_rate_indices],state_names)
	
	return(list(transition_matrix=transition_matrix, birth_rates=birth_rates, death_rates=death_rates))
}


# return flattened list of model parameter values (transition matrix, birth rates, death rates, sampling rates)
# matrices are flattened in row-major format
flatten_musse_params = function(params, Nstates){
	flattened = c(as.vector(t(params$transition_matrix)), params$birth_rates, params$death_rates)
	return(flattened)
}

# return unflattened list of model parameters (transition matrix, birth rates, death rates, sampling rates)
# this is the reverse of flatten_musse_params(..)
unflatten_musse_params = function(params_flat, Nstates, state_names){
	unflattened = list(	transition_matrix 	= matrix(params_flat[1:(Nstates*Nstates)], ncol=Nstates, nrow=Nstates, byrow=TRUE, dimnames=list(state_names,state_names)),
						birth_rates			= setNames(params_flat[((Nstates*Nstates)+(0*Nstates)+1):((Nstates*Nstates)+(1*Nstates))], state_names),
						death_rates			= setNames(params_flat[((Nstates*Nstates)+(1*Nstates)+1):((Nstates*Nstates)+(2*Nstates))], state_names))
	return(unflattened)
}

get_name_series = function(base_name, N){
	return(sapply(1:N, function(n) sprintf("%s%d",base_name,n)));
}

# return parameter names synchronized with the flattened lists generated by flatten_musse_params(..)
get_flat_musse_param_names = function(Nstates){
	Q = matrix(NA,nrow=Nstates,ncol=Nstates) # auxiliary matrix for getting flattened matrix indices (content does not matter)
	return(c(	paste("Q",paste(as.character(as.vector(t(row(Q)))),as.character(as.vector(t(col(Q)))),sep="."),sep=""),
				get_name_series("birth_rate.",Nstates), 
				get_name_series("death_rate.",Nstates)))
}


scale_params = function(params, scales, mins){
	scaled_params = (params-mins)/scales;
	return(scaled_params)
}


unscale_params = function(scaled_params, scales, mins){
	params = scales * scaled_params + mins;
	return(params)
}





