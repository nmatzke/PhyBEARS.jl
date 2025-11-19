\name{fit_sbm_const}
\alias{fit_sbm_const}
\title{
Fit a phylogeographic Spherical Brownian Motion model.
}
\description{
Given one or more rooted phylogenetic trees and geographic coordinates (latitudes & longitudes) for the tips of each tree, this function estimates the diffusivity of a Spherical Brownian Motion (SBM) model for the evolution of geographic location along lineages (Perrin 1928; Brillinger 2012). Estimation is done via maximum-likelihood and using independent contrasts between sister lineages.
}
\usage{
fit_sbm_const(trees, 
        tip_latitudes, 
        tip_longitudes, 
        radius,
        phylodistance_matrixes  = NULL,
        clade_states            = NULL,
        planar_approximation    = FALSE,
        only_basal_tip_pairs    = FALSE,
        only_distant_tip_pairs  = FALSE,
        min_MRCA_time           = 0,
        max_MRCA_age            = Inf,
        max_phylodistance       = Inf,
        no_state_transitions    = FALSE,
        only_state              = NULL,
        min_diffusivity         = NULL,
        max_diffusivity         = NULL,
        Nbootstraps             = 0, 
        NQQ                     = 0,
        SBM_PD_functor          = NULL,
        focal_diffusivities     = NULL)
}
\arguments{
\item{trees}{
Either a single rooted tree or a list of rooted trees, of class "phylo". The root of each tree is assumed to be the unique node with no incoming edge. Edge lengths are assumed to represent time intervals or a similarly interpretable phylogenetic distance. When multiple trees are provided, it is either assumed that their roots coincide in time (if \code{align_trees_at_root=TRUE}) or that each tree's youngest tip was sampled at present day (if \code{align_trees_at_root=FALSE}).
}
\item{tip_latitudes}{
Numeric vector of length Ntips, or a list of vectors, listing latitudes of tips in decimal degrees (from -90 to 90). If \code{trees} is a list of trees, then \code{tip_latitudes} should be a list of vectors of the same length as \code{trees}, listing tip latitudes for each of the input trees.
}
\item{tip_longitudes}{
Numeric vector of length Ntips, or a list of vectors, listing longitudes of tips in decimal degrees (from -180 to 180). If \code{trees} is a list of trees, then \code{tip_longitudes} should be a list of vectors of the same length as \code{trees}, listing tip longitudes for each of the input trees.
}
\item{radius}{
Strictly positive numeric, specifying the radius of the sphere. For Earth, the mean radius is 6371 km.
}
\item{phylodistance_matrixes}{
Numeric matrix, or a list of numeric matrixes, listing phylogenetic distances between tips for each tree. If \code{trees} is a list of trees, then \code{phylodistance_matrixes} should be a list of the same length as \code{trees}, whose n-th element should be a numeric matrix comprising as many rows and columns as there are tips in the n-th tree; the entry \code{phylodistance_matrixes[[n]][i,j]} is the phylogenetic distance between tips i and j in tree n. If \code{trees} is a single tree, then \code{phylodistance_matrixes} can be a single numeric matrix. If \code{NULL} (default), phylogenetic distances between tips are calculated based on the provided trees, otherwise phylogenetic distances are taken from \code{phylodistance_matrixes}; in the latter case the trees are only used for the topology (determining tip pairs for independent contrasts), but not for calculating phylogenetic distances.
}
\item{clade_states}{
Either NULL, or an integer vector of length Ntips+Nnodes, or a list of integer vectors, listing discrete states of every tip and node in the tree. If \code{trees} is a list of trees, then \code{clade_states} should be a list of vectors of the same length as \code{trees}, listing tip and node states for each of the input trees. For example, \code{clade_states[[2]][10]} specifies the state of the 10-th tip or node in the 2nd tree.
States may be, for example, geographic regions, sub-types, discrete traits etc, and can be used to restrict independent contrasts to tip pairs within the same state (see option \code{no_state_transitions}). 
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
\item{no_state_transitions}{
Logical, specifying whether to omit independent contrasts between tips whose shortest connecting paths include state transitions. If \code{TRUE}, only tips within the same state and with no transitions between them (as specified in \code{clade_states}) are compared. If \code{TRUE}, then \code{clade_states} must be provided.
}
\item{only_state}{
Optional integer, specifying the state in which tip pairs (and their connecting ancestral nodes) must be in order to be considered. If specified, then \code{clade_states} must be provided.
}
\item{min_diffusivity}{
Non-negative numeric, specifying the minimum possible diffusivity. If NULL, this is automatically chosen.
}
\item{max_diffusivity}{
Non-negative numeric, specifying the maximum possible diffusivity. If NULL, this is automatically chosen.
}
\item{Nbootstraps}{
Non-negative integer, specifying an optional number of parametric bootstraps to performs for estimating standard errors and confidence intervals.
}
\item{NQQ}{
Integer, optional number of simulations to perform for creating QQ plots of the theoretically expected distribution of geodistances vs. the empirical distribution of geodistances (across independent contrasts). The resolution of the returned QQ plot will be equal to the number of independent contrasts used for fitting. If <=0, no QQ plots will be calculated.
}
\item{SBM_PD_functor}{
SBM probability density functor object. Used internally for efficiency and for debugging purposes, and should be kept at its default value \code{NULL}.
}
\item{focal_diffusivities}{
Optional numeric vector, listing diffusivities of particular interest and for which the log-likelihoods should be returned. This may be used e.g. for diagnostic purposes, e.g. to see how "sharp" the likelihood peak is at the maximum-likelihood estimate.
}
}


\details{
For short expected transition distances this function uses the approximation formula by Ghosh et al. (2012). For longer expected transition distances the function uses a truncated approximation of the series representation of SBM transition densities (Perrin 1928). It is assumed that tips are sampled randomly without any biases for certain geographic regions. If you suspect strong geographic sampling biases, consider using the function \code{\link{fit_sbm_geobiased_const}}.

This function can use multiple trees to fit the diffusivity under the assumption that each tree is an independent realization of the same SBM process, i.e. all lineages in all trees dispersed with the same diffusivity.

If \code{edge.length} is missing from one of the input trees, each edge in the tree is assumed to have length 1. The tree may include multifurcations as well as monofurcations, however multifurcations are internally expanded into bifurcations by adding dummy nodes.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether the fitting was successful. If \code{FALSE}, then an additional return variable, \code{error}, will contain a description of the error; in that case all other return variables may be undefined.
}
\item{diffusivity}{
Numeric, the estimated diffusivity, in units distance^2/time. Distance units are the same as used for the \code{radius}, and time units are the same as the tree's edge lengths. For example, if the \code{radius} was specified in km and edge lengths are in Myr, then the estimated diffusivity will be in km^2/Myr.
}
\item{loglikelihood}{
Numeric, the log-likelihood of the data at the estimated diffusivity.
}
\item{Ncontrasts}{
Integer, number of independent contrasts (i.e., tip pairs) used to estimate the diffusivity. This is the number of independent data points used.
}
\item{phylodistances}{
Numeric vector of length \code{Ncontrasts}, listing the phylogenetic distances of the independent contrasts used in the fitting.
}
\item{geodistances}{
Numeric vector of length \code{Ncontrasts}, listing the geographical distances of the independent contrasts used in the fitting.
}
\item{focal_loglikelihoods}{
Numeric vector of the same length as \code{focal_diffusivities}, listing the log-likelihoods for the diffusivities provided in \code{focal_diffusivities}.
}
\item{standard_error}{
Numeric, estimated standard error of the estimated diffusivity, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50lower}{
Numeric, lower bound of the 50\% confidence interval for the estimated diffusivity (25-75\% percentile), based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50upper}{
Numeric, upper bound of the 50\% confidence interval for the estimated diffusivity, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95lower}{
Numeric, lower bound of the 95\% confidence interval for the estimated diffusivity (2.5-97.5\% percentile), based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95upper}{
Numeric, upper bound of the 95\% confidence interval for the estimated diffusivity, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{consistency}{
Numeric between 0 and 1, estimated consistency of the data with the fitted model. If \eqn{L} denotes the loglikelihood of new data generated by the fitted model (under the same model) and \eqn{M} denotes the expectation of \eqn{L}, then \code{consistency} is the probability that \eqn{|L-M|} will be greater or equal to \eqn{|X-M|}, where \eqn{X} is the loglikelihood of the original data under the fitted model.
Only returned if \code{Nbootstraps>0}. A low consistency (e.g., <0.05) indicates that the fitted model is a poor description of the data. See Lindholm et al. (2019) for background.
}
\item{QQplot}{
Numeric matrix of size Ncontrasts x 2, listing the computed QQ-plot. The first column lists quantiles of geodistances in the original dataset, the 2nd column lists quantiles of hypothetical geodistances simulated based on the fitted model.
}
\item{SBM_PD_functor}{
SBM probability density functor object. Used internally for efficiency and for debugging purposes.
}
}

\author{Stilianos Louca}

\references{
F. Perrin (1928). Etude mathematique du mouvement Brownien de rotation. 45:1-51.

D. R. Brillinger (2012). A particle migrating randomly on a sphere. in Selected Works of David Brillinger. Springer.

A. Ghosh, J. Samuel, S. Sinha (2012). A Gaussian for diffusion on the sphere. Europhysics Letters. 98:30003.

A. Lindholm, D. Zachariah, P. Stoica, T. B. Schoen (2019). Data consistency approach to model validation. IEEE Access. 7:59788-59796.

S. Louca (2021). Phylogeographic estimation and simulation of global diffusive dispersal. Systematic Biology. 70:340-359.
}

\seealso{
\code{\link{fit_sbm_geobiased_const}},
\code{\link{simulate_sbm}},
\code{\link{fit_sbm_parametric}},
\code{\link{fit_sbm_linear}},
\code{\link{fit_sbm_on_grid}}
}

\examples{
\dontrun{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=500)$tree

# simulate SBM on the tree
D = 1e4
simulation = simulate_sbm(tree, radius=6371, diffusivity=D)

# fit SBM on the tree
fit = fit_sbm_const(tree,simulation$tip_latitudes,simulation$tip_longitudes,radius=6371)
cat(sprintf('True D=\%g, fitted D=\%g\\n',D,fit$diffusivity))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SBM model}
\keyword{random}
\keyword{fitting}
\keyword{simulation}
