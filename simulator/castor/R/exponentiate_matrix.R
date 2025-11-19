# calculate exp(T*A) for some quadratic matrix A and several scalar scalings T
# this function becomes very efficient when the number of scalings is large
# returns the exponentials as a 3D array, with the [r,c,s]-th entry being the (r,c)-th entry in exp(scalings[s]*A)
# max_absolute_error is in terms of the Hilbert-Schmidt L2 norm.
# if enforce_probability_matrix==TRUE, then each column of each exponential is guaranteed to be a valid probability vector
exponentiate_matrix = function(	A, 
								scalings			= 1, 
								max_absolute_error	= 1e-3, 
								min_polynomials		= 1, 
								max_polynomials		= 1000){
	exponentials_flat = exponentiate_matrix_for_multiple_scalings_CPP(NR=nrow(A), A=as.vector(t(A)), scalings=scalings, epsilon=max_absolute_error, NPmin=min_polynomials, NPmax=max_polynomials, enforce_probability_matrix=FALSE);
	return(aperm(array(exponentials_flat,dim=c(ncol(A),nrow(A),length(scalings))),c(2,1,3)))
}

