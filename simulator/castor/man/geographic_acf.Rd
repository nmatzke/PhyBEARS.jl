\name{geographic_acf}
\alias{geographic_acf}
\title{
Phylogenetic autocorrelation function of geographic locations.
}
\description{
Given a rooted phylogenetic tree and geographic coordinates (latitudes & longitudes) of each tip, calculate the phylogenetic autocorrelation function (ACF) of the geographic locations. The ACF is a function of phylogenetic distance x, i.e., ACF(x) is the autocorrelation between two tip locations conditioned on the tips having phylogenetic ("patristic") distance x.
}
\usage{
geographic_acf( trees,
                tip_latitudes,
                tip_longitudes,
                Npairs              = 10000,
                Nbins               = NULL,
                min_phylodistance   = 0,
                max_phylodistance   = NULL,
                uniform_grid        = FALSE,
                phylodistance_grid  = NULL)
}
\arguments{
\item{trees}{
Either a single rooted tree of class "phylo", or a list of multiple such trees.
}
\item{tip_latitudes}{
Either a numeric vector of size Ntips (if \code{trees} was a single tree), specifying the latitudes (decimal degrees) of the tree's tips, or a list of such numeric vectors (if \code{trees} contained multiple trees) specifying the latitudes of each tree's tips. Note that \code{tip_latitudes[k][i]} must correspond to the i-th tip in the k-th input tree, i.e. as listed in \code{trees[[k]]$tip.label}. By convention, positive latitudes correspond to the northern hemisphere.
}
\item{tip_longitudes}{
Similar to \code{tip_latitudes}, but listing the latitudes (decimal degrees) of each tip in each input tree. By convention, positive longitudes correspond to the hemisphere East of the prime meridian.
}
\item{Npairs}{
Maximum number of random tip pairs to draw from each tree. A greater number of tip pairs will improve the accuracy of the estimated ACF within each distance bin. Tip pairs are drawn randomly with replacement, if \code{Npairs} is lower than the number of tip pairs in a tree. If \code{Npairs=Inf}, then every tip pair of every tree is included exactly once (for small and moderately sized trees this is recommended).
}
\item{Nbins}{
Number of phylogenetic distance bins to consider. A greater number of bins will increase the resolution of the ACF as a function of phylogenetic distance, but will decrease the number of tip pairs falling within each bin (which reduces the accuracy of the estimated ACF). If \code{NULL}, then \code{Nbins} is automatically and somewhat reasonably chosen based on the size of the input trees.
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
The autocorrelation between random geographic locations is defined as the expectation of \eqn{<X,Y>}, where <> is the scalar product and \eqn{X} and \eqn{Y} are the unit vectors pointing towards the two random locations on the sphere. For comparison, for a spherical Brownian Motion model with constant diffusivity \eqn{D} and radius \eqn{r} the autocorrelation function is given by \eqn{ACF(t)=e^{-2Dt/r^2}} (see e.g. \code{\link{simulate_sbm}}). Note that this function assumes that Earth is a perfect sphere.

The phylogenetic autocorrelation function (ACF) of the geographic distribution of species can give insight into the dispersal processes shaping species distributions over global scales. An ACF that decays slowly with increasing phylogenetic distance indicates a strong phylogenetic conservatism of the location and thus slow dispersal, whereas a  rapidly decaying ACF indicates weak phylogenetic conservatism and thus fast dispersal. Similarly, if the mean distance between two random tips increases with phylogenetic distance, this indicates a phylogenetic autocorrelation of species locations. Here, phylogenetic distance between tips refers to their patristic distance, i.e. the minimum cumulative edge length required to connect the two tips.

Since the phylogenetic distances between all possible tip pairs do not cover a continuoum (as there is only a finite number of tips), this function randomly draws tip pairs from the tree, maps them onto a finite set of phylodistance bins and then estimates the ACF for the centroid of each bin based on tip pairs in that bin. In practice, as a next step one would usually plot the estimated ACF (returned vector \code{autocorrelations}) over the centroids of the phylodistance bins (returned vector \code{phylodistances}). When multiple trees are provided as input, then the ACF is first calculated separately for each tree, and then averaged across trees (weighted by the number of tip pairs included from each tree in each bin).

Phylogenetic distance bins can be specified in two alternative ways: Either a set of bins (phylodistance grid) is automatically calculated based on the provided \code{Nbins}, \code{min_phylodistance}, \code{max_phylodistance} and \code{uniform_grid}, or a phylodistance grid is explicitly provided via \code{phylodistance_grid} and \code{max_phylodistance}.

The trees may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). If edge lengths are missing from the trees, then every edge is assumed to have length 1. The input trees must be rooted at some node for technical reasons (see function \code{\link{root_at_node}}), but the choice of the root node does not influence the result.

This function assumes that each tip is assigned exactly one geographic location. This might be problematic in situations where each tip covers multiple geographic locations, for example if tips are species and multiple individuals were sampled from each species. In that case, one might consider representing each individual as a separate tip in the tree, so that each tip has exactly one geographic location.
}
\value{
A list with the following elements:
\item{success}{
Logical, indicating whether the calculation was successful. If \code{FALSE}, an additional element \code{error} (character) is returned that provides a brief description of the error that occurred; in that case all other return values may be undefined.
}
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
Numeric vector of size Nbins, storing the estimated geographic autocorrelation for each phylodistance bin.
}
\item{std_autocorrelations}{
Numeric vector of size Nbins, storing the standard deviation of geographic autocorrelations encountered in each phylodistance bin. Note that this is not the standard error of the estimated ACF; it is a measure for how different the geographic locations are between tip pairs within each phylodistance bin.
}
\item{mean_geodistances}{
Numeric vector of size Nbins, storing the mean geographic distance between tip pairs in each distance bin, in units of sphere radii. If you want geographic distances in km, you need to multiply these by Earth's mean radius in km (about 6371). If multiple input trees were provided, this is the average across all trees, weighted by the number of tip pairs included from each tree in each bin.
}
\item{std_geodistances}{
Numeric vector of size Nbins, storing the standard deviation of geographic distances between tip pairs in each distance bin, in units of sphere radii.
}
\item{Npairs_per_distance}{
Integer vector of size Nbins, storing the number of random tip pairs associated with each distance bin.
}
}


\author{Stilianos Louca}

%\references{
%J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
%}

\seealso{
\code{\link{consentrait_depth}}, 
\code{\link{get_trait_acf}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=1000)$tree

# simulate spherical Brownian Motion on the tree
simul = simulate_sbm(tree, radius=1, diffusivity=0.1)
tip_latitudes  = simul$tip_latitudes
tip_longitudes = simul$tip_longitudes

# calculate geographical autocorrelation function
ACF = geographic_acf(tree, tip_latitudes, tip_longitudes, Nbins=10)

# plot ACF (autocorrelation vs phylogenetic distance)
plot(ACF$phylodistances, ACF$autocorrelations, type="l", xlab="distance", ylab="ACF")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{trait evolution}
