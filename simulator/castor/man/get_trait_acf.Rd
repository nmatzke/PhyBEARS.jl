\name{get_trait_acf}
\alias{get_trait_acf}
\title{
Phylogenetic autocorrelation function of a numeric trait.
}
\description{
Given a rooted phylogenetic tree and a numeric (typically continuous) trait with known value (state) on each tip, calculate the phylogenetic autocorrelation function (ACF) of the trait. The ACF is a function of phylogenetic distance x, where ACF(x) is the Pearson autocorrelation of the trait between two tips, provided that the tips have phylogenetic ("patristic") distance x. The function \code{get_trait_acf} also calculates the mean absolute difference and the mean relative difference of the trait between any two random tips at phylogenetic distance x (see details below).
}
\usage{
get_trait_acf(tree, 
              tip_states,
              Npairs            = 10000,
              Nbins             = NULL,
              min_phylodistance = 0,
              max_phylodistance = NULL,
              uniform_grid      = FALSE,
              phylodistance_grid= NULL)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
A numeric vector of size Ntips, specifying the value of the trait at each tip in the tree. Note that tip_states[i] (where i is an integer index) must correspond to the i-th tip in the tree.
}
\item{Npairs}{
Total number of random tip pairs to draw. A greater number of tip pairs will improve the accuracy of the estimated ACF within each distance bin. Tip pairs are drawn randomly with replacement. If \code{Npairs<=0}, then every tip pair is included exactly once.
}
\item{Nbins}{
Number of distance bins to consider within the range of phylogenetic distances encountered between tip pairs in the tree. A greater number of bins will increase the resolution of the ACF as a function of phylogenetic distance, but will decrease the number of tip pairs falling within each bin (which reduces the accuracy of the estimated ACF). If \code{NULL}, then \code{Nbins} is automatically and somewhat reasonably chosen based on the size of the input trees.
}
\item{min_phylodistance}{
Numeric, minimum phylogenetic distance to conssider. Only relevant if \code{phylodistance_grid} is \code{NULL}.
}
\item{max_phylodistance}{
Numeric, optional maximum phylogenetic distance to consider. If \code{NULL}, this is automatically set to the maximum phylodistance between any two tips.
}
\item{uniform_grid}{
Logical, specifying whether the phylodistance grid should be uniform, i.e., with equally sized phylodistance bins. If \code{FALSE}, then the grid is chosen non-uniformly (i.e., each bin has different size) such that each bin roughly contains the same number of tip pairs. Only relevant if \code{phylodistance_grid} is \code{NULL}. It is generally recommended to keep \code{uniform_grid=FALSE}, to avoid uneven estimation errors across bins.
}
\item{phylodistance_grid}{
Numeric vector, optional explicitly specified phylodistance bins (left boundaries thereof) on which to evaluate the ACF. Must contain non-negative numbers in strictly ascending order. Hence, the first bin will range from \code{phylodistance_grid[1]} to \code{phylodistance_grid[2]}, while the last bin will range from \code{tail(phylodistance_grid,1)} to \code{max_phylodistance}. Can be used as an alternative to \code{Nbins}. If non-\code{NULL}, then \code{Nbins}, \code{min_phylodistance} and \code{uniform_grid} are irrelevant. 
}
}
\details{
The phylogenetic autocorrelation function (ACF) of a trait can give insight into the evolutionary processes shaping its distribution across clades. An ACF that decays slowly with increasing phylogenetic distance indicates a strong phylogenetic conservatism of the trait, whereas a  rapidly decaying ACF indicates weak phylogenetic conservatism. Similarly, if the mean absolute difference in trait value between two random tips increases with phylogenetic distance, this indicates a phylogenetic autocorrelation of the trait (Zaneveld et al. 2014). Here, phylogenetic distance between tips refers to their patristic distance, i.e. the minimum cumulative edge length required to connect the two tips.

Since the phylogenetic distances between all possible tip pairs do not cover a continuoum (as there is only a finite number of tips), this function randomly draws tip pairs from the tree, maps them onto a finite set of equally-sized distance bins and then estimates the ACF for the centroid of each distance bin based on tip pairs in that bin. In practice, as a next step one would usually plot the estimated ACF (returned vector \code{autocorrelations}) over the centroids of the distance bins (returned vector \code{distances}).

Phylogenetic distance bins can be specified in two alternative ways: Either a set of bins (phylodistance grid) is automatically calculated based on the provided \code{Nbins}, \code{min_phylodistance}, \code{max_phylodistance} and \code{uniform_grid}, or a phylodistance grid is explicitly provided via \code{phylodistance_grid} and \code{max_phylodistance}.

The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). If \code{tree$edge.length} is missing, then every edge is assumed to have length 1. The input tree must be rooted at some node for technical reasons (see function \code{\link{root_at_node}}), but the choice of the root node does not influence the result.

This function assumes that each tip is assigned exactly one trait value. This might be problematic in situations where each tip covers a range of trait values, for example if tips are species and multiple individuals were sampled from each species. In that case, one might consider representing each individual as a separate tip in the tree, so that each tip has exactly one trait value.
}
\value{
A list with the following elements:
\item{phylodistances}{
Numeric vector of size Nbins, storing the center of each phylodistance bin in increasing order. This is equal to \code{0.5*(left_phylodistances+right_phylodistances)}. Typically, you will want to plot \code{autocorrelations} over \code{phylodistances}.
}
\item{left_phylodistances}{
Numeric vector of size Nbins, storing the left boundary of each phylodistance bin in increasing order.
}
\item{right_phylodistances}{
Numeric vector of size Nbins, storing the right boundary of each phylodistance bin in increasing order.
}
\item{autocorrelations}{
Numeric vector of size Nbins, storing the estimated Pearson autocorrelation of the trait for each distance bin.
}
\item{mean_abs_differences}{
Numeric vector of size Nbins, storing the mean absolute difference of the trait between tip pairs in each distance bin.
}
\item{mean_rel_differences}{
Numeric vector of size Nbins, storing the mean relative difference of the trait between tip pairs in each distance bin. The relative difference between two values \eqn{X} and \eqn{Y} is 0 if \eqn{X==Y}, and equal to 
\deqn{\frac{|X-Y|}{0.5\cdot (|X|+|Y|)}}
otherwise.
}
\item{Npairs_per_distance}{
Integer vector of size Nbins, storing the number of random tip pairs associated with each phylodistance bin.
}
}


\author{Stilianos Louca}

\references{
J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
}

\seealso{
\code{\link{consentrait_depth}}, 
\code{\link{geographic_acf}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=1000)$tree

# simulate continuous trait evolution on the tree
tip_states = simulate_bm_model(tree, diffusivity=1)$tip_states

# calculate autocorrelation function
ACF = get_trait_acf(tree, tip_states, Npairs=1e7, Nbins=10)

# plot ACF (autocorrelation vs phylogenetic distance)
plot(ACF$phylodistances, ACF$autocorrelations, type="l", xlab="distance", ylab="ACF")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{trait evolution}
