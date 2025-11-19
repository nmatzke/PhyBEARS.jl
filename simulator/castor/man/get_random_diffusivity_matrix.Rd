\name{get_random_diffusivity_matrix}
\alias{get_random_diffusivity_matrix}
\title{
Create a random diffusivity matrix for a Brownian motion model.
}
\description{
Create a random diffusivity matrix for a Brownian motion model of multi-trait evolution. This may be useful for testing purposes. The diffusivity matrix is drawn from the Wishart distribution of symmetric, nonnegative-definite matrixes: 
\deqn{D = X^T \cdot X,\quad X[i,j]\sim N(0,V),\quad i=1,..,n, j=1,..,p,} 
where n is the degrees of freedom, p is the number of traits and V is a scalar scaling.
}
\usage{
get_random_diffusivity_matrix(Ntraits, degrees=NULL, V=1)
}
\arguments{
\item{Ntraits}{
The number of traits modelled. Equal to the number of rows and the number of columns of the returned matrix.
}
\item{degrees}{
Degrees of freedom for the Wishart distribution. Must be equal to or greater than \code{Ntraits}. Can also be \code{NULL}, which is the same as setting it equal to \code{Ntraits}.
}
\item{V}{
Positive number. A scalar scaling for the Wishart distribution.
}
}


%\details{
%}


\value{
A real-valued quadratic symmetric non-negative definite matrix of size Ntraits x Ntraits. Almost surely (in the probabilistic sense), this matrix will be positive definite.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_random_mk_transition_matrix}}, 
\code{\link{simulate_bm_model}}
}

\examples{
# generate a 5x5 diffusivity matrix
D = get_random_diffusivity_matrix(Ntraits=5)

# check that it is indeed positive definite
if(all(eigen(D)$values>0)){ 
  cat("Indeed positive definite\n");
}else{ 
  cat("Not positive definite\n"); 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BM model}
\keyword{random}
