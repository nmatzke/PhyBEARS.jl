\name{get_trait_stats_over_time}
\alias{get_trait_stats_over_time}
\title{
Calculate mean & standard deviation of a numeric trait on a dated tree over time.
}
\description{
Given a rooted and dated phylogenetic tree, and a scalar numeric trait with known value on each node and tip of the tree, calculate the mean and the variance of the trait's states across the tree at discrete time points. For example, if the trait represents "body size", then this function calculates the mean body size of extant clades over time.
}
\usage{
get_trait_stats_over_time(tree,
                          states, 
                          Ntimes            = NULL,
                          times             = NULL,
                          include_quantiles = TRUE,
                          check_input       = TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo", where edge lengths represent time intervals (or similar).
}
\item{states}{
Numeric vector, specifying the trait's state at each tip and each node of the tree (in the order in which tips & nodes are indexed). May include \code{NA} or \code{NaN} if values are missing for some tips/nodes.
}
\item{Ntimes}{
Integer, number of equidistant time points for which to calculade clade counts. Can also be \code{NULL}, in which case \code{times} must be provided.
}
\item{times}{
Integer vector, listing time points (in ascending order) for which to calculate clade counts. Can also be \code{NULL}, in which case \code{Ntimes} must be provided.
}
\item{include_quantiles}{
Logical, specifying whether to include information on quantiles (e.g., median, CI95, CI50) of the trait over time, in addition to the means and standard deviations. This option increases computation time and memory needs for large trees, so if you only care about means and standard deviations you can set this to \code{FALSE}.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
If \code{tree$edge.length} is missing, then every edge in the tree is assumed to be of length 1. The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). The tree need not be ultrametric (e.g. may include extinct tips), although in general this function only makes sense if edge lengths correspond to time (or similar).

Either \code{Ntimes} or \code{times} must be non-\code{NULL}, but not both. \code{states} need not include names; if it does, then these are checked to be in the same order as in the tree (if \code{check_input==TRUE}).
}


\value{
A list with the following elements:
\item{Ntimes}{
Integer, indicating the number of returned time points. Equal to the provided \code{Ntimes} if applicable.
}
\item{times}{
Numeric vector of size Ntimes, listing the considered time points in increasing order. If \code{times} was provided as an argument to the function, then this will be the same as provided.
}
\item{clade_counts}{
Integer vector of size Ntimes, listing the number of tips or nodes considered at each time point.
}
\item{means}{
Numeric vector of size Ntimes, listing the arithmetic mean of trait states at each time point.
}
\item{stds}{
Numeric vector of size Ntimes, listing the standard deviation of trait states at each time point.
}
\item{medians}{
Numeric vector of size Ntimes, listing the median trait state at each time point. Only returned if \code{include_uantiles=TRUE}.
}
\item{CI50lower}{
Numeric vector of size Ntimes, listing the lower end of the equal-tailed 50\% range of trait states (i.e., the 25\% percentile) at each time point. Only returned if \code{include_uantiles=TRUE}.
}
\item{CI50upper}{
Numeric vector of size Ntimes, listing the upper end of the equal-tailed 50\% range of trait states (i.e., the 75\% percentile) at each time point. Only returned if \code{include_uantiles=TRUE}.
}
\item{CI95lower}{
Numeric vector of size Ntimes, listing the lower end of the equal-tailed 95\% range of trait states (i.e., the 2.5\% percentile) at each time point. Only returned if \code{include_uantiles=TRUE}.
}
\item{CI95upper}{
Numeric vector of size Ntimes, listing the upper end of the equal-tailed 95\% range of trait states (i.e., the 97.5\% percentile) at each time point. Only returned if \code{include_uantiles=TRUE}.
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

# simulate a numeric trait under Brownian-motion
trait = simulate_bm_model(tree, diffusivity=1)
states = c(trait$tip_states,trait$node_states)

# calculate trait stats over time
results = get_trait_stats_over_time(tree, states, Ntimes=100)

# plot trait stats over time (mean +/- std)
M = results$means
S = results$stds
matplot(x=results$times, 
        y=matrix(c(M-S,M+S),ncol=2,byrow=FALSE),
        main = "Simulated BM trait over time",
        lty = 1, col="black",
        type="l", xlab="time", ylab="mean +/- std")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{trait evolution}
