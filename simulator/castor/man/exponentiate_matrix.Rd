\name{exponentiate_matrix}
\alias{exponentiate_matrix}
\title{
Exponentiate a matrix.
}
\description{
Calculate the exponential \eqn{\exp(T\cdot A)} of some quadratic real-valued matrix A for one or more scalar scaling factors T.
}
\usage{
exponentiate_matrix(A, scalings=1, max_absolute_error=1e-3, 
                    min_polynomials=1, max_polynomials=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{A}{
A real-valued quadratic matrix of size N x N.
}
\item{scalings}{
Vector of real-valued scalar scaling factors T, for each of which the exponential \eqn{\exp(T\cdot A)} should be calculated.
}
\item{max_absolute_error}{
Maximum allowed absolute error for the returned approximations. A smaller allowed error implies a greater computational cost as more matrix polynomials need to be included (see below). The returned approximations may have a greater error if the parameter \code{max_polynomials} is set too low.
}
\item{min_polynomials}{
Minimum number of polynomials to include in the approximations (see equation below), even if \code{max_absolute_error} may be satisfied with fewer polynomials. If you don't know how to choose this, in most cases the default is fine. Note that regardless of \code{min_polynomials} and \code{max_absolute_error}, the number of polynomials used will not exceed \code{max_polynomials}.
}
\item{max_polynomials}{
Maximum allowed number of polynomials to include in the approximations (see equation below). Meant to provide a safety limit for the amount of memory and the computation time required. For example, a value of 1000 means that up to 1000 matrices (powers of A) of size N x N may be computed and stored temporarily in memory. Note that if \code{max_polynomials} is too low, the requested accuracy (via \code{max_absolute_error}) may not be achieved. That said, for large trees more memory may be required to store the actual result rather than the intermediate polynomials, i.e. for purposes of saving RAM it doesn't make much sense to choose \code{max_polynomials} much smaller than the length of \code{scalings}.
}

}
\details{
Discrete character evolution Markov models often involve repeated exponentiations of the same transition matrix along each edge of the tree (i.e. with the scaling T being the edge length). Matrix exponentiation can become a serious computational bottleneck for larger trees or large matrices (i.e. spanning multiple discrete states). This function pre-calculates polynomials \eqn{A^p/p!} of the matrix, and then uses linear combinations of the same polynomials for each requested T: 
\deqn{
\exp(T\cdot A) = \sum_{p=0}^P T^p\frac{A^p}{p!} + ...
}
This function thus becomes very efficient when the number of scaling factors is large (e.g. >10,000).
The number of polynomials included is determined based on the specified \code{max_absolute_error}, and based on the largest (by magnitude) scaling factor requested. The function utilizes the balancing algorithm proposed by James et al (2014, Algorithm 3) and the scaling & squaring method (Moler and Van Loan, 2003) to improve the conditioning of the matrix prior to exponentiation.
}
\value{
A 3D numeric matrix of size N x N x S, where N is the number of rows & column of the input matrix A and S is the length of \code{scalings}. The [r,c,s]-th element of this matrix is the entry in the r-th row and c-th column of \eqn{\exp(scalings[s]\cdot A)}.
}

\author{Stilianos Louca}

\references{
R. James, J. Langou and B. R. Lowery (2014). On matrix balancing and eigenvector computation. arXiv:1401.5766

C. Moler and C. Van Loan (2003). Nineteen dubious ways to compute the exponential of a matrix, twenty-five years later. SIAM Review. 45:3-49.
}

%\seealso{
%\code{\link{ps_ouss}}
%}

\examples{
# create a random 5 x 5 matrix
A = get_random_mk_transition_matrix(Nstates=5, rate_model="ER")

# calculate exponentials exp(0.1*A) and exp(10*A)
exponentials = exponentiate_matrix(A, scalings=c(0.1,10))

# print 1st exponential: exp(0.1*A)
print(exponentials[,,1])

# print 2nd exponential: exp(10*A)
print(exponentials[,,2])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ algebra }
\keyword{ exponential }% __ONLY ONE__ keyword per line
