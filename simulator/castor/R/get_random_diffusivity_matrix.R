# create a random diffusivity matrix for for a multivariate Brownian motion model of continuous trait co-evolution
# The matrix is drawn from the Wishart distribution of symmetric, nonnegative-definite matrixes: D ~ W_p(V,n), where n=degrees and p=Ntraits
get_random_diffusivity_matrix = function(Ntraits, degrees=NULL, V=1){
	if(is.null(degrees)) degrees = Ntraits;
	if(degrees<Ntraits) stop(sprintf("ERROR: Degrees of freedom should be greater or equal to the number of traits (%d)",Ntraits))
	if(V<=0) stop(sprintf("ERROR: Scaling V must be positive (got V=%g)",V))
	X = matrix(stats::rnorm(n=degrees*Ntraits, mean=0, sd=V), ncol=Ntraits)
	D = t(X) %*% X
	return(D)
}

