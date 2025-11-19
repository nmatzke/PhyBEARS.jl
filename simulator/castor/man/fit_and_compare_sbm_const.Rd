\name{fit_and_compare_sbm_const}
\alias{fit_and_compare_sbm_const}
\title{
Fit and compare Spherical Brownian Motion models for diffusive geographic dispersal between two data sets.
}
\description{
Given two rooted phylogenetic trees and geographic coordinates of the trees' tips, fit a Spherical Brownian Motion (SBM) model of diffusive geographic dispersal with constant diffusivity to each tree and compare the fitted models. This function estimates the diffusivity (\eqn{D}) for each data set (i.e., each set of trees + tip-coordinates) via maximum-likelihood and assesses whether the log-difference between the two fitted diffusivities is statistically significant, under the null hypothesis that the two data sets exhibit the same diffusivity. Optionally, multiple trees can be used as input for each data set, under the assumption that dispersal occurred according to the same diffusivity in each tree of that dataset. For more details on how SBM is fitted to each data set see the function \code{\link{fit_sbm_const}}.
}
\usage{
fit_and_compare_sbm_const(  trees1, 
                            tip_latitudes1,
                            tip_longitudes1,
                            trees2,
                            tip_latitudes2,
                            tip_longitudes2,
                            radius,
                            planar_approximation    = FALSE,
                            only_basal_tip_pairs    = FALSE,
                            only_distant_tip_pairs  = FALSE,
                            min_MRCA_time           = 0,
                            max_MRCA_age            = Inf,
                            max_phylodistance       = Inf,
                            min_diffusivity         = NULL,
                            max_diffusivity         = NULL,
                            Nbootstraps             = 0,
                            Nsignificance           = 0,
                            SBM_PD_functor          = NULL,
                            verbose                 = FALSE,
                            verbose_prefix          = "")
}
\arguments{
\item{trees1}{
Either a single rooted tree or a list of rooted trees, of class "phylo", corresponding to the first data set on which an SBM model is to be fitted. Edge lengths are assumed to represent time intervals or a similarly interpretable phylogenetic distance.
}
\item{tip_latitudes1}{
Numeric vector listing the latitude (in decimal degrees) of each tip in each tree in the first data set. If \code{trees1} is a single tree, then \code{tip_latitudes1} must be a numeric vector of size Ntips, listing the latitudes for each tip in the tree.
If \code{trees1} is a list of Ntrees trees, then \code{tip_latitudes1} must be a list of length Ntrees, each element of which lists the latitudes for the corresponding tree (as a vector, similarly to the single-tree case).
}
\item{tip_longitudes1}{
Similar to \code{tip_latitudes1}, but listing longitudes (in decimal degrees) of each tip in each tree in the first data set.
}
\item{trees2}{
Either a single rooted tree or a list of rooted trees, of class "phylo", corresponding to the second data set on which an SBM model is to be fitted. Edge lengths are assumed to represent time intervals or a similarly interpretable phylogenetic distance.
}
\item{tip_latitudes2}{
Numeric vector listing the latitude (in decimal degrees) of each tip in each tree in the second data set, similarly to \code{tip_latitudes1}.
}
\item{tip_longitudes2}{
Numeric vector listing the longitude (in decimal degrees) of each tip in each tree in the second data set, similarly to \code{tip_longitudes1}.
}
\item{radius}{
Strictly positive numeric, specifying the radius of the sphere. For Earth, the mean radius is 6371 km.
}
\item{planar_approximation}{
Logical, specifying whether to estimate the diffusivity based on a planar approximation of the SBM model, i.e. by assuming that geographic distances between tips are as if tips are distributed on a 2D cartesian plane. This approximation is only accurate if geographical distances between tips are small compared to the sphere's radius.
}
\item{only_basal_tip_pairs}{
Logical, specifying whether to only compare immediate sister tips, i.e., tips connected through a single parental node.
}
\item{only_distant_tip_pairs}{
Logical, specifying whether to only compare tips at distinct geographic locations.
}
\item{min_MRCA_time}{
Numeric, specifying the minimum allowed time (distance from root) of the most recent common ancestor (MRCA) of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set \code{min_MRCA_time<=0} to disable this filter.
}
\item{max_MRCA_age}{
Numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set \code{max_MRCA_age=Inf} to disable this filter.
}
\item{max_phylodistance}{
Numeric, maximum allowed geodistance for an independent contrast to be included in the SBM fitting. Set \code{max_phylodistance=Inf} to disable this filter.
}
\item{min_diffusivity}{
Non-negative numeric, specifying the minimum possible diffusivity. If NULL, this is automatically chosen.
}
\item{max_diffusivity}{
Non-negative numeric, specifying the maximum possible diffusivity. If NULL, this is automatically chosen.
}
\item{Nbootstraps}{
Integer, specifying the number of parametric bootstraps to perform for calculating the confidence intervals of SBM diffusivities fitted to each data set. If <=0, no bootstrapping is performed.
}
\item{Nsignificance}{
Integer, specifying the number of simulations to perform for assessing the statistical significance of the linear difference and log-transformed difference between the diffusivities fitted to the two data sets, i.e. of \eqn{|D_1-D_2|} and of \eqn{|\log(D_1)-\log(D_2)|}. Set to 0 to not calculate statistical significances. See below for additional details.
}
\item{SBM_PD_functor}{
SBM probability density functor object. Used internally and for debugging purposes. Unless you know what you're doing, you should keep this \code{NULL}.
}
\item{verbose}{
Logical, specifying whether to print progress report messages to the screen.
}
\item{verbose_prefix}{
Character, specifying a prefix to include in front of progress report messages on each line. Only relevant if \code{verbose==TRUE}.
}
}


\details{
For details on the Spherical Brownian Motion model see \code{\link{fit_sbm_const}} and \code{\link{simulate_sbm}}. This function separately fits an SBM model with constant diffusivity to each of two data sets; internally, this function applies \code{fit_sbm_const} to each data set. 

If \code{Nsignificance>0}, the statistical significance of the linear difference (\eqn{|D_1-D_2|}) and log-transformed difference (\eqn{|\log(D_1)-\log(D_2)|}) of the two fitted diffusivities is assessed under the null hypothesis that both data sets were generated by the same common SBM model. The diffusivity of this common SBM model is estimated by fitting to both datasets at once, i.e. after merging the two datasets into a single dataset of trees and tip coordinates (see return variable \code{fit_common} below). For each of the \code{Nsignificance} random simulations of the common SBM model on the two tree sets, the diffusivities are again separately fitted on the two simulated sets and the resulting difference and log-difference is compared to those of the original data sets. The returned \code{lin_significance} (or \code{log_significance}) is the probability that the diffusivities would have a difference (or log-difference) larger than the observed one, if the two data sets had been generated under the common SBM model.

If \code{edge.length} is missing from one of the input trees, each edge in the tree is assumed to have length 1. Trees may include multifurcations as well as monofurcations, however multifurcations are internally expanded into bifurcations by adding dummy nodes.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether the fitting was successful for both data sets. If \code{FALSE}, then an additional return variable, \code{error}, will contain a description of the error; in that case all other return variables may be undefined.
}
\item{fit1}{
A named list containing the fitting results for the first data set, in the same format as returned by \code{\link{fit_sbm_const}}. In particular, the diffusivity fitted to the first data set will be stored in \code{fit1$diffusivity}.
}
\item{fit2}{
A named list containing the fitting results for the second data set, in the same format as returned by \code{\link{fit_sbm_const}}. In particular, the diffusivity fitted to the second data set will be stored in \code{fit2$diffusivity}.
}
\item{lin_difference}{
The absolute difference between the two diffusivities, i.e. \eqn{|D_1-D_2|}.
}
\item{log_difference}{
The absolute difference between the two log-transformed diffusivities, i.e. \eqn{|\log(D_1)-\log(D_2)|}.
}
\item{lin_significance}{
Numeric, statistical significance of the observed lin-difference under the null hypothesis that the two data sets were generated by a common SBM model. Only returned if \code{Nsignificance>0}.
}
\item{log_significance}{
Numeric, statistical significance of the observed log-difference under the null hypothesis that the two data sets were generated by a common SBM model. Only returned if \code{Nsignificance>0}.
}
\item{fit_common}{
A named list containing the fitting results for the two data sets combined, in the same format as returned by \code{\link{fit_sbm_const}}. The common diffusivity, \code{fit_common$diffusivity} is used for the random simulations when assessing the statistical significance of the lin-difference and log-difference of the separately fitted diffusivities. Only returned if \code{Nsignificance>0}.
}
}

\author{Stilianos Louca}

\references{
S. Louca (in review as of 2020). Phylogeographic estimation and simulation of global diffusive dispersal. Systematic Biology.
}

\seealso{
\code{\link{simulate_sbm}},
\code{\link{fit_sbm_const}},
\code{\link{fit_sbm_linear}},
\code{\link{fit_sbm_parametric}}
}

\examples{
\dontrun{
# simulate distinct SBM models on two random trees
radius = 6371   # Earth's radius
D1     = 1      # diffusivity on 1st tree
D2     = 3      # diffusivity on 2nd tree
tree1  = generate_random_tree(list(birth_rate_factor=1),max_tips=100)$tree
tree2  = generate_random_tree(list(birth_rate_factor=1),max_tips=100)$tree
sim1   = simulate_sbm(tree=tree1, radius=radius, diffusivity=D1)
sim2   = simulate_sbm(tree=tree2, radius=radius, diffusivity=D2)
tip_latitudes1  = sim1$tip_latitudes
tip_longitudes1 = sim1$tip_longitudes
tip_latitudes2  = sim2$tip_latitudes
tip_longitudes2 = sim2$tip_longitudes

# fit and compare SBM models between the two hypothetical data sets
fit = fit_and_compare_sbm_const(trees1          = tree1, 
                                tip_latitudes1  = tip_latitudes1, 
                                tip_longitudes1 = tip_longitudes1, 
                                trees2          = tree2,
                                tip_latitudes2  = tip_latitudes2, 
                                tip_longitudes2 = tip_longitudes2, 
                                radius          = radius,
                                Nbootstraps     = 0,
                                Nsignificance   = 100)

# print summary of results
cat(sprintf("Fitted D1 = \%g, D2 = \%g, significance of log-diff. = \%g\n",
            fit$fit1$diffusivity, fit$fit2$diffusivity, fit$log_significance))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SBM model}
