\name{get_stationary_distribution}
\alias{get_stationary_distribution}
\title{
Stationary distribution of Markov transition matrix.
}
\description{
Calculate the stationary probability distribution vector p for a transition matrix Q of a continuous-time Markov chain. That is, calculate \eqn{p\in[0,1]^n} such that \code{sum(p)==0} and \eqn{p^TQ=0}.
}
\usage{
get_stationary_distribution(Q)
}
\arguments{
\item{Q}{
A valid transition rate matrix of size Nstates x Nstates, i.e. a quadratic matrix in which every row sums up to zero.
}
}
\details{
A stationary distribution of a discrete-state continuous-time Markov chain is a probability distribution across states that remains constant over time, i.e. \eqn{p^TQ=0}. Note that in some cases (i.e. if Q is not irreducible), there may be multiple distinct stationary distributions. In that case,which one is returned by this function is unpredictable. Internally, p is estimated by stepwise minimization of the norm of \eqn{p^TQ}, starting with the vector p in which every entry equals 1/Nstates.
}
\value{
A numeric vector of size Nstates and with non-negative entries, satisfying the conditions \code{p\%*\%Q==0} and \code{sum(p)==1.0}.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{exponentiate_matrix}}
}

\examples{
# generate a random 5x5 Markov transition matrix
Q = get_random_mk_transition_matrix(Nstates=5, rate_model="ARD")

# calculate stationary probability distribution
p = get_stationary_distribution(Q)
print(p)

# test correctness (p*Q should be 0, apart from rounding errors)
cat(sprintf("max(abs(p*Q)) = \%g\n",max(abs(p \%*\% Q))))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{stationary distribution}
\keyword{Mk model}
