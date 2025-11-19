# create a random transition rate matrix for a fixed-rates continuous-time Markov model of discrete character evolution
# rate_model can be "ER", "SYM", "ARD" or "SUEDE"
# 0<=min_rate<=max_rate
get_random_mk_transition_matrix = function(Nstates, rate_model, min_rate=0, max_rate=1){
	if((min_rate<0) || (max_rate<min_rate)) stop("ERROR: Invalid min_rate and/or max_rate; must satisfy 0<=min_rate<=max_rate")
	temp = get_transition_index_matrix(Nstates=Nstates, rate_model=rate_model);
	Q 	 = get_transition_matrix_from_rate_vector(stats::runif(n=temp$Nrates,min=min_rate,max=max_rate), index_matrix=temp$index_matrix, Nstates=Nstates)
	return(Q);
}