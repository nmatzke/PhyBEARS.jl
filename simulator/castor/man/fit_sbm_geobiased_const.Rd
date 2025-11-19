\name{fit_sbm_geobiased_const}
\alias{fit_sbm_geobiased_const}
\title{
Fit a phylogeographic Spherical Brownian Motion model with geographic sampling bias.
}
\description{
Given one or more rooted phylogenetic trees and geographic coordinates (latitudes & longitudes) for the tips of each tree, this function estimates the diffusivity of a Spherical Brownian Motion (SBM) model for the evolution of geographic location along lineages (Perrin 1928; Brillinger 2012), while correcting for geographic sampling biases. Estimation is done via maximum-likelihood and using independent contrasts between sister lineages, while correction for geographic sampling bias is done through an iterative simulation+fitting process until convergence.
}
\usage{
fit_sbm_geobiased_const(trees,
                        tip_latitudes,
                        tip_longitudes,
                        radius,
                        reference_latitudes    = NULL,
                        reference_longitudes   = NULL,
                        only_basal_tip_pairs   = FALSE,
                        only_distant_tip_pairs = FALSE,
                        min_MRCA_time          = 0,
                        max_MRCA_age           = Inf,
                        max_phylodistance      = Inf,
                        min_diffusivity        = NULL,
                        max_diffusivity        = NULL,
                        rarefaction            = 0.1,
                        Nsims                  = 100,
                        max_iterations         = 100,
                        Nbootstraps            = 0,
                        NQQ                    = 0,
                        Nthreads               = 1,
                        include_simulations    = FALSE,
                        SBM_PD_functor         = NULL,
                        verbose                = FALSE,
                        verbose_prefix         = "")
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
\item{reference_latitudes}{
Optional numeric vector, listing latitudes of reference coordinates based on which to calculate the geographic sampling density. If \code{NULL}, the geographic sampling density is estimated based on \code{tip_latitudes} and \code{tip_longitudes}.
}
\item{reference_longitudes}{
Optional numeric vector of the same length as \code{reference_latitudes}, listing latitudes of reference coordinates based on which to calculate the geographic sampling density. If \code{NULL}, the geographic sampling density is estimated based on \code{tip_latitudes} and \code{tip_longitudes}.
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
Non-negative numeric, specifying the minimum possible diffusivity. If \code{NULL}, this is automatically chosen.
}
\item{max_diffusivity}{
Non-negative numeric, specifying the maximum possible diffusivity. If \code{NULL}, this is automatically chosen.
}

\item{rarefaction}{
Numeric, between 00 and 1, specifying the fraction of extant lineages to sample from the simulated trees. Should be strictly smaller than 1, in order for geographic bias correction to have an effect. Note that regardless of \code{rarefaction}, the simulated trees will have the same size as the original trees.
}
\item{Nsims}{
Integer, number of SBM simulatons to perform per iteration for assessing the effects of geographic bias. Smaller trees require larger \code{Nsims} (due to higher stochasticity). This must be at least 2, although values of 100-1000 are recommended.
}
\item{max_iterations}{
Integer, maximum number of iterations (correction steps) to perform before giving up.
}

\item{Nbootstraps}{
Non-negative integer, specifying an optional number of parametric bootstraps to performs for estimating standard errors and confidence intervals.
}
\item{NQQ}{
Integer, optional number of simulations to perform for creating QQ plots of the theoretically expected distribution of geodistances vs. the empirical distribution of geodistances (across independent contrasts). The resolution of the returned QQ plot will be equal to the number of independent contrasts used for fitting. If <=0, no QQ plots will be calculated.
}
\item{Nthreads}{
Integer, number of parallel threads to use. Ignored on Windows machines.
}
\item{include_simulations}{
Logical, whether to include the trees and tip coordinates simulated under the final fitted SBM model, in the returned results. May be useful e.g. for checking model adequacy.
}
\item{SBM_PD_functor}{
SBM probability density functor object. Used internally for efficiency and for debugging purposes, and should be kept at its default value \code{NULL}.
}
\item{verbose}{
Logical, specifying whether to print progress reports and warnings to the screen.
}
\item{verbose_prefix}{
Character, specifying the line prefix for printing progress reports to the screen.
}
}


\details{
This function tries to estimate the true spherical diffusivity of an SBM model of geographic diffusive dispersal, while correcting for geographic sampling biases. This is done using an iterative refinement approach, by which trees and tip locations are repeatedly simulated under the current true diffusivity estimate and the diffusivity estimated from those simulated data are compared to the originally uncorrected diffusivity estimate. Trees are simulated according to a birth-death model with constant rates, fitted to the original input trees (a congruent birth-death model is chosen to match the requested \code{rarefaction}). Simulated trees are subsampled (rarefied) to match the original input tree sizes, with sampled lineages chosen randomly but in a geographically biased way that resembles the original geographic sampling density (e.g., as inferred from the \code{reference_latitudes} and \code{reference_longitudes}). Internally, this function repeatedly applies \code{\link{fit_sbm_const}} and \code{\link{simulate_sbm}}. If the true sampling fraction of the input trees is unknown, then it is advised to perform the analysis with a few alternative \code{rarefaction} values (e.g., 0.01 and 0.1) to verify the robustness of the estimates.

If \code{edge.length} is missing from one of the input trees, each edge in the tree is assumed to have length 1. The tree may include multifurcations as well as monofurcations, however multifurcations are internally expanded into bifurcations by adding dummy nodes.
}			

\value{
A list with the following elements:
\item{success}{
Logical, indicating whether the fitting was successful. If \code{FALSE}, then an additional return variable, \code{error}, will contain a description of the error; in that case all other return variables may be undefined.
}
\item{Nlat}{
Integer, number of latitude-tiles used for building a map of the geographic sampling biases.
}
\item{Nlon}{
Integer, number of longitude-tiles used for building a map of the geographic sampling biases.
}
\item{diffusivity}{
Numeric, the estimated true diffusivity, i.e. accounting for geographic sampling biases, in units distance^2/time. Distance units are the same as used for the \code{radius}, and time units are the same as the tree's edge lengths. For example, if the \code{radius} was specified in km and edge lengths are in Myr, then the estimated diffusivity will be in km^2/Myr.
}
\item{correction_factor}{
Numeric, estimated ratio between the true diffusivity and the original (uncorrected) diffusivity estimate.
}
\item{Niterations}{
Integer, the number of iterations performed until convergence.
}
\item{stopping_criterion}{
Character, a short description of the criterion by which the iteration was eventually halted.
}
\item{uncorrected_fit_diffusivity}{
Numeric, the originally estimated (uncorrected) diffusivity.
}
\item{last_sim_fit_diffusivity}{
Numeric, the mean uncorrected diffuvity estimated from the simulated data in the last iteration. Convergence means that \code{last_sim_fit_diffusivity} came close to \code{uncorrected_fit_diffusivity}.
}
\item{all_correction_factors}{
Numeric vector of length \code{Niterations}, listing the estimated correction factors in each iteration.
}
\item{all_diffusivity_estimates}{
Numeric vector of length \code{Niterations}, listing the mean uncorrected diffusivity estimated from the simulated data in each iteration.
}
\item{Ntrees}{
Integer, number of trees considered for the simulations. This might have smaller than \code{length(trees)}, if for some trees fitting a birth-death model was not possible.
}
\item{lambda}{
Numeric vector of length \code{Ntrees}, listing the birth rates used to simulate the trees.
}
\item{mu}{
Numeric vector of length \code{Ntrees}, listing the death rates used to simulate the trees.
}
\item{rarefaction}{
Numeric vector of length \code{Ntrees}, listing the rarefactions (sampling fractions) used to simulate the trees. These will typically be equal to the \code{rarefaction} provided by the function caller, but may differ for example if the congruence class did not include a birth-death model with the requested rarefaction.
}
\item{Ncontrasts}{
Integer, number of independent contrasts (i.e., tip pairs) used to estimate the diffusivity. This is the number of independent data points used.
}
\item{standard_error}{
Numeric, estimated standard error of the estimated true diffusivity, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50lower}{
Numeric, lower bound of the 50\% confidence interval for the estimated true diffusivity (25-75\% percentile), based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50upper}{
Numeric, upper bound of the 50\% confidence interval for the estimated true diffusivity, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95lower}{
Numeric, lower bound of the 95\% confidence interval for the estimated true diffusivity (2.5-97.5\% percentile), based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95upper}{
Numeric, upper bound of the 95\% confidence interval for the estimated true diffusivity, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{QQplot}{
Numeric matrix of size Ncontrasts x 2, listing the computed QQ-plot. The first column lists quantiles of geodistances in the original dataset, the 2nd column lists quantiles of hypothetical geodistances simulated based on the estimated true diffusivity.
}
\item{simulations}{
List, containing the trees and tip coordinates simulated under the final fitted SBM model, accounting for geographic biases. Each entry is itself a named list, containing the simulations corresponding to a specific input tree. In particular, \code{simulations[[t]]$sims[[r]]} is the r-th simulation performed that corresponds to the t-th input tree. Each simulation is again a named list, containing the elements \code{success} (logical), \code{tree} (of class \code{phylo}), \code{latitudes} (numeric vector) and \code{longitudes} (numeric vector). This data structure may be useful for testing the adequacy of the fitted SBM model; only use this if you know what you are doing. Only returned if \code{include_simulations} was \code{TRUE}.
}
\item{SBM_PD_functor}{
SBM probability density functor object. Used internally for efficiency and for debugging purposes. Most users can ignore this.
}
}

\author{Stilianos Louca}

\references{
F. Perrin (1928). Etude mathematique du mouvement Brownien de rotation. 45:1-51.

D. R. Brillinger (2012). A particle migrating randomly on a sphere. in Selected Works of David Brillinger. Springer.

A. Ghosh, J. Samuel, S. Sinha (2012). A Gaussian for diffusion on the sphere. Europhysics Letters. 98:30003.

S. Louca (2021). Phylogeographic estimation and simulation of global diffusive dispersal. Systematic Biology. 70:340-359.

S. Louca (in review as of 2021). The rates of global microbial dispersal.
}

\seealso{
\code{\link{simulate_sbm}},
\code{\link{fit_sbm_parametric}},
\code{\link{fit_sbm_linear}},
\code{\link{fit_sbm_on_grid}}
}

\examples{
\dontrun{
NFullTips   = 10000
diffusivity = 1
radius      = 6371

# generate tree and run SBM on it
cat(sprintf("Generating tree and simulating SBM (true D=\%g)..\n",diffusivity))
tree = castor::generate_tree_hbd_reverse(Ntips  = NFullTips,
                                         lambda = 5e-7,
                                         mu     = 2e-7,
                                         rho    = 1)$trees[[1]]
SBMsim = simulate_sbm(tree = tree, radius = radius, diffusivity = diffusivity)

# select subset of tips only found in certain geographic regions
min_abs_lat = 30
max_abs_lat = 80
min_lon     = 0
max_lon     = 90
keep_tips   = which((abs(SBMsim$tip_latitudes)<=max_abs_lat)
                    & (abs(SBMsim$tip_latitudes)>=min_abs_lat)
                    & (SBMsim$tip_longitudes<=max_lon)
                    & (SBMsim$tip_longitudes>=min_lon))
rarefaction     = castor::get_subtree_with_tips(tree, only_tips = keep_tips)
tree            = rarefaction$subtree
tip_latitudes   = SBMsim$tip_latitudes[rarefaction$new2old_tip]
tip_longitudes  = SBMsim$tip_longitudes[rarefaction$new2old_tip]
Ntips           = length(tree$tip.label)
rarefaction     = Ntips/NFullTips

# fit SBM while correcting for geographic sampling biases
fit = castor:::fit_sbm_geobiased_const(trees            = tree,
                                       tip_latitudes    = tip_latitudes,
                                       tip_longitudes   = tip_longitudes,
                                       radius           = radius,
                                       rarefaction      = Ntips/NFullTips,
                                       Nsims            = 10,
                                       Nthreads         = 4,
                                       verbose          = TRUE,
                                       verbose_prefix   = "  ")
if(!fit$success){
    cat(sprintf("ERROR: \%s\n",fit$error))
}else{
    cat(sprintf("Estimated true D = \%g\n",fit$diffusivity))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SBM model}
\keyword{random}
\keyword{fitting}
\keyword{simulation}
