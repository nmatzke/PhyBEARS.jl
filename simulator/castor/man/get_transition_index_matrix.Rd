\name{get_transition_index_matrix}
\alias{get_transition_index_matrix}
\title{
Create an index matrix for a Markov transition model.
}
\description{
Create an index matrix encoding the parametric structure of the transition rates in a discrete-state continuous-time Markov model (e.g., Mk model of trait evolution). Such an index matrix is required by certain functions for mapping independent rate parameters to transition rates. For example, an index matrix may encode the information that each rate i-->j is equal to its reversed counterpart j-->i.
}
\usage{
get_transition_index_matrix(Nstates, rate_model)
}
\arguments{
\item{Nstates}{
Integer, the number of distinct states represented in the transition matrix (number of rows & columns).
}
\item{rate_model}{
Rate model that the transition matrix must satisfy. Can be "ER" (all rates equal), "SYM" (transition rate i-->j is equal to transition rate j-->i), "ARD" (all rates can be different) or "SUEDE" (only stepwise transitions i-->i+1 and i-->i-1 allowed, all 'up' transitions are equal, all 'down' transitions are equal).
}
}


\details{
The returned index matrix will include as many different positive integers as there are independent rate parameters in the requested rate model, plus potentially the value 0 (which has a special meaning, see below).
}


\value{
A named list with the following elements:
\item{index_matrix}{
Integer matrix of size Nstates x Nstates, with values between 0 and Nstates, assigning each entry in the transition matrix to an independent transition rate parameter. A value of 0 means that the corresponding rate is fixed to zero (if off-diagonal) or will be adjusted to ensure a valid Markov transition rate matrix (if on the diagonal).
}
\item{Nrates}{
Integer, the number of independent rate parameters in the model.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_random_mk_transition_matrix}}
}

%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Mk model}
\keyword{MuSSE}
\keyword{BiSSE}
\keyword{HiSSE}
