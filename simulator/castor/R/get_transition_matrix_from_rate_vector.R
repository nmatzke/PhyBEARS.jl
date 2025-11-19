get_transition_matrix_from_rate_vector = function(rates, index_matrix, Nstates){
	transition_matrix = matrix(c(0, rates)[index_matrix + 1], nrow=Nstates, ncol=Nstates, byrow = FALSE);
	diag(transition_matrix) = -rowSums(transition_matrix, na.rm = TRUE);
	return(transition_matrix);
}