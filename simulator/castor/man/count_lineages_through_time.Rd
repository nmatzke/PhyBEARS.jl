\name{count_lineages_through_time}
\alias{count_lineages_through_time}
\title{
Count number of lineages through time (LTT).
}
\description{
Given a rooted timetree (i.e., a tree whose edge lengths represent time intervals), calculate the number of lineages represented in the tree at various time points, otherwise known as "lineages through time"" (LTT) curve. The root is interpreted as time 0, and the distance of any node or tip from the root is interpreted as time elapsed since the root. Optionally, the slopes and relative slopes of the LTT curve are also returned.
}
\usage{
count_lineages_through_time(  tree, 
                              Ntimes        = NULL, 
                              min_time      = NULL,
                              max_time      = NULL,
                              times         = NULL, 
                              include_slopes= FALSE,
                              ultrametric   = FALSE,
                              degree        = 1,
                              regular_grid  = TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo", where edge lengths represent time intervals (or similar).
}
\item{Ntimes}{
Integer, number of equidistant time points at which to count lineages. Can also be \code{NULL}, in which case \code{times} must be provided.
}
\item{min_time}{
Minimum time (distance from root) to consider. If \code{NULL}, this will be set to the minimum possible (i.e. 0). Only relevant if \code{times==NULL}.
}
\item{max_time}{
Maximum time (distance from root) to consider. If \code{NULL}, this will be set to the maximum possible. Only relevant if \code{times==NULL}.
}
\item{times}{
Integer vector, listing time points (in ascending order) at which to count lineages. Can also be \code{NULL}, in which case \code{Ntimes} must be provided.
}
\item{include_slopes}{
Logical, specifying whether the slope and the relative slope of the returned clades-per-time-point curve should also be returned.
}
\item{ultrametric}{
Logical, specifying whether the input tree is guaranteed to be ultrametric, even in the presence of some numerical inaccuracies causing some tips not have exactly the same distance from the root. If you know the tree is ultrametric, then this option helps the function choose a better time grid for the LTT.
}
\item{degree}{
Integer, specifying the "degree" of the LTT curve: LTT(t) will be the number of lineages in the tree at time t that have at least n descending tips in the tree. Typically \code{order=1}, which corresponds to the classical LTT curve.
}
\item{regular_grid}{
Logical, specifying whether the automatically generated time grid should be regular (equal distances between grid points). This option only matters if \code{times==NULL}. If \code{regular_grid==FALSE} and \code{times==NULL}, the time grid will be irregular, with grid point density being roughly proportional to the square root of the number of lineages at any particular time (i.e., the grid becomes finer towards the tips).
}
}


\details{
Given a sequence of time points between a tree's root and tips, this function essentially counts how many edges "cross" each time point (if \code{degree==1}). The slopes and relative slopes are calculated from this curve using finite differences.

Note that the classical LTT curve (\code{degree=1}) is non-decreasing over time, whereas higher-degree LTT's may be decreasing as well as increasing over time.
 
If \code{tree$edge.length} is missing, then every edge in the tree is assumed to be of length 1. The tree may include multifurcations as well as monofurcations (i.e. nodes with only one child). The tree need not be ultrametric, although in general this function only makes sense for dated trees (e.g., where edge lengths are time intervals or similar).

Either \code{Ntimes} or \code{times} must be non-\code{NULL}, but not both. If \code{times!=NULL}, then \code{min_time} and \code{max_time} must be \code{NULL}.
}


\value{
A list with the following elements:
\item{Ntimes}{
Integer, indicating the number of returned time points. Equal to the provided \code{Ntimes} if applicable.
}
\item{times}{
Numeric vector of size Ntimes, listing the time points at which the LTT was calculated. If \code{times} was provided as an argument to the function, then this will be the same as provided.
}
\item{lineages}{
Integer vector of size Ntimes, listing the number of lineages represented in the tree at each time point that have at least \code{degree} descending tips, i.e. the LTT curve.
}
\item{slopes}{
Numeric vector of size Ntimes, listing the slopes (finite-difference approximation of 1st derivative) of the LTT curve.
}
\item{relative_slopes}{
Numeric vector of size Ntimes, listing the relative slopes of the LTT curve, i.e. \code{slopes} divided by a sliding-window average of \code{lineages}.
}
}


\author{Stilianos Louca}

%\references{
%}

%\seealso{
%}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1), max_tips=1000)$tree

# calculate classical LTT curve
results = count_lineages_through_time(tree, Ntimes=100)

# plot classical LTT curve
plot(results$times, results$lineages, type="l", xlab="time", ylab="# clades")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{lineages through time}
