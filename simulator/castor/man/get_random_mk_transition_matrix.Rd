\name{get_random_mk_transition_matrix}
\alias{get_random_mk_transition_matrix}
\title{
Create a random transition matrix for an Mk model.
}
\description{
Create a random transition matrix for a fixed-rates continuous-time Markov model of discrete trait evolution ("Mk model"). This may be useful for testing purposes.
}
\usage{
get_random_mk_transition_matrix(Nstates, rate_model, min_rate=0, max_rate=1)
}
\arguments{
\item{Nstates}{
The number of distinct states represented in the transition matrix (number of rows & columns).
}
\item{rate_model}{
Rate model that the transition matrix must satisfy. Can be "ER" (all rates equal), "SYM" (transition rate i-->j is equal to transition rate j-->i), "ARD" (all rates can be different) or "SUEDE" (only stepwise transitions i-->i+1 and i-->i-1 allowed, all 'up' transitions are equal, all 'down' transitions are equal).
}
\item{min_rate}{
A non-negative number, specifying the minimum rate in off-diagonal entries of the transition matrix.
}
\item{max_rate}{
A non-negative number, specifying the maximum rate in off-diagonal entries of the transition matrix. Must not be smaller than \code{min_rate}.
}
}


%\details{
%}


\value{
A real-valued quadratic matrix of size Nstates x Nstates, representing a transition matrix for an Mk model. Each row will sum to 0. The [r,c]-th entry represents the transition rate r-->c. The number of unique off-diagonal rates will depend on the \code{rate_model} chosen.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{exponentiate_matrix}}, \code{\link{get_stationary_distribution}}
}

\examples{
# generate a 5x5 Markov transition rate matrix
Q = get_random_mk_transition_matrix(Nstates=5, rate_model="ARD")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Mk model}
\keyword{random}
