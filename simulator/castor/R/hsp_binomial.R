# Hidden state prediction for a binary trait, by fitting the probability of a binomial distribution to empirical state frequencies.
# For each internal node, the probability P1 (probability of a randomly chosen descending tip being in state 1) is fitted separately for each node based on the observed states in the descending tips via maximum likelihood.
# The function can account for potential state-measurement errors, hidden states and reveal biases (ASR phase).
# Only nodes with a number of non-hidden tips above a certain threshold are included in the estimation phase (ASR-phase). 
# All other nodes and hidden tips are then assigned the probabilities of a sufficiently ancestral node with estimated probabilities (HSP-phase)
# This function is similar to hsp_empirical_probabilities, but it accounts for potential state-measurement errors and reveal biases.
hsp_binomial = function(tree, 
						tip_states, 	# integer vector of size Ntips, listing measured tip states (potentially erroneously). May include NAs for hidden states, or 1 or 2 for measured states.
						reveal_probs	= NULL,	# 2D matrix of size Ntips x 2, listing the reveal probabilities at each tip conditional on the tip's true state. May also be a vector of length 2 (same reveal_probs for all tips) or NULL (unbiased reveal probs)
						state1_probs	= NULL,	# 2D matrix of size Ntips x 2, listing the probability of measuring state 1 (potentially erroneously) at each tip conditional upon its true state and conditional upon revealing (measuring) some state. For example, for an incompletely sequenced genome with completion level C_i and state1=absence and state2=presence of a gene, one has state1_probs[i,1] = 1, state1_probs[i,2]=1 - C_i. May also be a vector of length 2 (same probabilities for all tips) or NULL. If NULL, state measurements are assumed error-free, and hence this is the same as c(1,0).
						min_revealed 	= 1,	# minimum number of descending tips with revealed state for inferring P1 for a clade during the ASR-phase. For clades with too few descending tips with revealed state, the probability P1 will not be estimated in the ASR-phase. It is advised to set this greater than zero. The larger this number, the more nodes will not have their P1 calculated.
						max_STE			= Inf,	# maximum acceptable standard error for the estimated probability P1 for a clade. If the STE for a clade exceeds this threshold, the P1 for that clade is set to the P1 of the nearest ancestor with STE below that threshold. Setting this to Inf disables this functionality.
						check_input		= TRUE){
	Ntips   = length(tree$tip.label);
	Nnodes  = tree$Nnode;
	Nclades = Ntips + Nnodes
	
	# basic error checking
	if(length(tip_states)!=Ntips) return(list(success=FALSE, error=sprintf("Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips)))
	if(check_input){
		if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) return(list(success=FALSE, error="Names in tip_states and tip labels in tree don't match (must be in the same order)."))
		if(any(!(is.na(tip_states) | (tip_states==1) | (tip_states==2)))) return(list(success=FALSE, error="All tip states must be either 1, 2 or NA"))
	}

	if(is.null(reveal_probs)) reveal_probs = c(1,1)*sum(!(is.na(tip_states)))/Ntips; # no reveal biases, i.e. reveal_probs[1]==reveal_probs[1]. In that case the scaling of reveal_probs does not matter, so we can w.l.o.g. set them both to the empirical reveal fraction (fraction of non-hidden tips)
	if(is.vector(reveal_probs)){
		if(length(reveal_probs)==2){
			reveal_probs = matrix(rep(reveal_probs,times=Ntips),ncol=2,byrow=TRUE)
		}else{
			return(list(success=FALSE, error=sprintf("reveal_probs was provided as a vector, in which case it must contain exactly 2 elements; instead, found %d elements",length(reveal_probs))))
		}
	}else if(is.matrix(reveal_probs)){
		if(ncol(reveal_probs)!=2) return(list(success=FALSE, error=sprintf("reveal_probs was provided as a matrix, in which case it must contain exactly 2 columns; instead, found %d columns",ncol(reveal_probs))))
		if(nrow(reveal_probs)==1){
			# assume the same reveal_probs for all tips
			reveal_probs = matrix(rep(reveal_probs,times=Ntips),ncol=2,byrow=TRUE)
		}else if(nrow(reveal_probs)!=Ntips){
			return(list(success=FALSE, error=sprintf("reveal_probs was provided as a matrix, which case it must have %d rows (=Ntips); instead, found %d rows",Ntips,nrow(reveal_probs))))
		}
	}else{
		return(list(success=FALSE, error=sprintf("Unknown data type for reveal_probs: Expected either a vector of length 2 or a matrix with 2 columns and %d rows",Ntips)))
	}
	
	if(is.null(state1_probs)) state1_probs = c(1,0); # error-free state measurements
	if(is.vector(state1_probs)){
		if(length(state1_probs)==2){
			state1_probs = matrix(rep(state1_probs,times=Ntips),ncol=2,byrow=TRUE)
		}else{
			return(list(success=FALSE, error=sprintf("state1_probs was provided as a vector, in which case it must contain exactly 2 elements; instead, found %d elements",length(state1_probs))))
		}
	}else if(is.matrix(state1_probs)){
		if(ncol(state1_probs)!=2) return(list(success=FALSE, error=sprintf("state1_probs was provided as a matrix, in which case it must contain exactly 2 columns; instead, found %d columns",ncol(state1_probs))))
		if(nrow(state1_probs)==1){
			# assume the same state1_probs for all tips
			state1_probs = matrix(rep(state1_probs,times=Ntips),ncol=2,byrow=TRUE)
		}else if(nrow(state1_probs)!=Ntips){
			return(list(success=FALSE, error=sprintf("state1_probs was provided as a matrix, which case it must have %d rows (=Ntips); instead, found %d rows",Ntips,nrow(state1_probs))))
		}
	}else{
		return(list(success=FALSE, error=sprintf("Unknown data type for state1_probs: Expected either a vector of length 2 or a matrix with 2 columns and %d rows",Ntips)))
	}
	
	# estimate empirical probabilities for ancestral nodes (ASR-phase)
	# Note that the C++ routine assumes binary states within {0,1}, whereas in the R routines we assume states within {1,2}.
	input_tip_states = tip_states - 1 # make zero-based for C++ routine
	input_tip_states[is.na(input_tip_states)] = -1;
	asr_results = ASR_binomial_CPP(	Ntips			= Ntips,
									Nnodes			= Nnodes,
									Nedges			= nrow(tree$edge),
									tree_edge		= as.vector(t(tree$edge))-1, # flatten in row-major format
									tip_states		= input_tip_states,
									reveal_probs	= as.vector(t(reveal_probs)), # flatten in row-major format
									state0_probs	= as.vector(t(state1_probs)), # flatten in row-major format
									min_revealed	= min_revealed);
	if(!asr_results$success) return(list(success=FALSE, error=sprintf("ASR failed: %s",asr_results$error)))
	P1s 			= asr_results$P0s
	STEs			= asr_results$STEs # standard error estimates
	P1_known 		= rep(TRUE, times=Nclades)
	P1_known[P1s<0]	= FALSE
	P1_known[STEs>max_STE] = FALSE
	# for tips with definitely known true state, don't use the ML estimates and just fix P1 based on their state
	definites = which(sapply(1:Ntips, FUN=function(tip) return((!is.na(tip_states[tip])) && ((tip_states[tip]==1 && state1_probs[tip,2]==0) || (tip_states[tip]==2 && state1_probs[tip,1]==1)))))
	if(length(definites)>0){
		P1s[definites] 		= (tip_states[definites]==1);
		STEs[definites] 	= 0
		P1_known[definites] = TRUE
	}
	inheritted = rep(TRUE, times=Nclades)
	inheritted[P1_known] = FALSE

	# propagate P1s and STEs to tips with hidden state (HSP-phase)
	P1s = apply_attributes_to_descendants_CPP(	Ntips				= Ntips,
												Nnodes				= Nnodes,
												Nedges				= nrow(tree$edge),
												Nattributes			= 1,
												tree_edge			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												attributes_known	= P1_known,
												attributes			= P1s);
	P1s[P1s<0] = NaN
	STEs = apply_attributes_to_descendants_CPP(	Ntips				= Ntips,
												Nnodes				= Nnodes,
												Nedges				= nrow(tree$edge),
												Nattributes			= 1,
												tree_edge			= as.vector(t(tree$edge))-1,	# flatten in row-major format and make indices 0-based
												attributes_known	= P1_known,
												attributes			= STEs);
	STEs[STEs<0] = NaN
		
	return(list(success			= TRUE, 
				P1				= P1s,
				STE				= STEs, 
				reveal_counts	= asr_results$reveal_counts,
				inheritted		= inheritted))
}


