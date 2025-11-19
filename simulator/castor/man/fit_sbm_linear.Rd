\name{fit_sbm_linear}
\alias{fit_sbm_linear}
\title{
Fit a phylogeographic Spherical Brownian Motion model with linearly varying diffusivity.
}
\description{
Given a rooted phylogenetic tree and geographic coordinates (latitudes & longitudes) for its tips, this function estimates the diffusivity of a Spherical Brownian Motion (SBM) model for the evolution of geographic location along lineages (Perrin 1928; Brillinger 2012), assuming that the diffusivity varies linearly over time. Estimation is done via maximum-likelihood and using independent contrasts between sister lineages. This function is designed to estimate the diffusivity over time, by fitting two parameters defining the diffusivity as a linear function of time. For fitting more general functional forms see \code{\link{fit_sbm_parametric}}.
}
\usage{
fit_sbm_linear(tree, 
              tip_latitudes,
              tip_longitudes,
              radius,
              clade_states          = NULL,
              planar_approximation  = FALSE,
              only_basal_tip_pairs  = FALSE,
              only_distant_tip_pairs= FALSE,
              min_MRCA_time         = 0,
              max_MRCA_age          = Inf,
              max_phylodistance     = Inf,
              no_state_transitions  = FALSE,
              only_state            = NULL,
              time1                 = 0,
              time2                 = NULL,
              Ntrials               = 1,
              Nthreads              = 1,
              Nbootstraps           = 0,
              Ntrials_per_bootstrap = NULL,
              Nsignificance         = 0,
              NQQ                   = 0,
              fit_control           = list(),
              SBM_PD_functor        = NULL,
              verbose               = FALSE,
              verbose_prefix        = "")
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge. Edge lengths are assumed to represent time intervals or a similarly interpretable phylogenetic distance.
}
\item{tip_latitudes}{
Numeric vector of length Ntips, listing latitudes of tips in decimal degrees (from -90 to 90). The order of entries must correspond to the order of tips in the tree (i.e., as listed in \code{tree$tip.label}).
}
\item{tip_longitudes}{
Numeric vector of length Ntips, listing longitudes of tips in decimal degrees (from -180 to 180). The order of entries must correspond to the order of tips in the tree (i.e., as listed in \code{tree$tip.label}).
}
\item{radius}{
Strictly positive numeric, specifying the radius of the sphere. For Earth, the mean radius is 6371 km.
}
\item{clade_states}{
Optional integer vector of length Ntips+Nnodes, listing discrete states of every tip and node in the tree. The order of entries must match the order of tips and nodes in the tree. States may be, for example, geographic regions, sub-types, discrete traits etc, and can be used to restrict independent contrasts to tip pairs within the same state (see option \code{no_state_transitions}).
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
Numeric, specifying the minimum allowed time (distance from root) of the most recent common ancestor (MRCA) of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at least this distance from the root. Set \code{min_MRCA_time=0} to disable this filter.
}
\item{max_MRCA_age}{
Numeric, specifying the maximum allowed age (distance from youngest tip) of the MRCA of sister tips considered in the fitting. In other words, an independent contrast is only considered if the two sister tips' MRCA has at most this age (time to present). Set \code{max_MRCA_age=Inf} to disable this filter.
}
\item{max_phylodistance}{
Numeric, maximum allowed geodistance for an independent contrast to be included in the SBM fitting. Set \code{max_phylodistance=Inf} to disable this filter.
}
\item{no_state_transitions}{
Logical, specifying whether to omit independent contrasts between tips whose shortest connecting paths include state transitions. If \code{TRUE}, only tips within the same state and with no transitions between them (as specified in \code{clade_states}) are compared.
}
\item{only_state}{
Optional integer, specifying the state in which tip pairs (and their connecting ancestral nodes) must be in order to be considered. If specified, then \code{clade_states} must be provided.
}
\item{time1}{
Optional numeric, specifying the first time point at which to estimate the diffusivity. By default this is set to root (i.e., time 0).
}
\item{time2}{
Optional numeric, specifying the first time point at which to estimate the diffusivity. By default this is set to the present day (i.e., the maximum distance of any tip from the root).
}
\item{Ntrials}{
Integer, specifying the number of independent fitting trials to perform, each starting from a random choice of model parameters. Increasing \code{Ntrials} reduces the risk of reaching a non-global local maximum in the fitting objective.
}
\item{Nthreads}{
Integer, specifying the number of parallel threads to use for performing multiple fitting trials simultaneously. This should generally not exceed the number of available CPUs on your machine. Parallel computing is not available on the Windows platform.
}
\item{Nbootstraps}{
Integer, specifying the number of parametric bootstraps to perform for estimating standard errors and confidence intervals of estimated model parameters. Set to 0 for no bootstrapping.
}
\item{Ntrials_per_bootstrap}{
Integer, specifying the number of fitting trials to perform for each bootstrap sampling. If \code{NULL}, this is set equal to \code{max(1,Ntrials)}. Decreasing \code{Ntrials_per_bootstrap} will reduce computation time, at the expense of potentially inflating the estimated confidence intervals; in some cases (e.g., for very large trees) this may be useful if fitting takes a long time and confidence intervals are very narrow anyway. Only relevant if \code{Nbootstraps>0}.
}
\item{Nsignificance}{
Integer, specifying the number of simulations to perform under a const-diffusivity model for assessing the statistical significance of the fitted slope. Set to 0 to not calculate the significance of the slope.
}
\item{NQQ}{
Integer, optional number of simulations to perform for creating QQ plots of the theoretically expected distribution of geodistances vs. the empirical distribution of geodistances (across independent contrasts). The resolution of the returned QQ plot will be equal to the number of independent contrasts used for fitting. If <=0, no QQ plots will be calculated.
}
\item{fit_control}{
 Named list containing options for the \code{nlminb} optimization routine, such as \code{iter.max}, \code{eval.max} or \code{rel.tol}. For a complete list of options and default values see the documentation of \code{nlminb} in the \code{stats} package.
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
This function is essentially a wrapper for the more general function \code{\link{fit_sbm_parametric}}, with the addition that it can estimate the statistical significance of the fitted linear slope. 

The statistical significance of the slope is the probability that a constant-diffusivity SBM model would generate data that would yield a fitted linear slope equal to or greater than the one fitted to the original data; the significance is estimated by simulating \code{Nsignificance} constant-diffusivity models and then fitting a linear-diffusivity model. The constant diffusivity assumed in these simulations is the maximum-likelihood diffusivity fitted internally using \code{\link{fit_sbm_const}}.

Note that estimation of diffusivity at older times is only possible if the timetree includes extinct tips or tips sampled at older times (e.g., as is often the case in viral phylogenies). If tips are only sampled once at present-day, i.e. the timetree is ultrametric, reliable diffusivity estimates can only be achieved near present times.

For short expected transition distances this function uses the approximation formula by Ghosh et al. (2012) to calculate the probability density of geographical transitions along edges. For longer expected transition distances the function uses a truncated approximation of the series representation of SBM transition densities (Perrin 1928).

If \code{edge.length} is missing from one of the input trees, each edge in the tree is assumed to have length 1. The tree may include multifurcations as well as monofurcations, however multifurcations are internally expanded into bifurcations by adding dummy nodes.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether the fitting was successful. If \code{FALSE}, then an additional return variable, \code{error}, will contain a description of the error; in that case all other return variables may be undefined.
}
\item{objective_value}{
The maximized fitting objective. Currently, only maximum-likelihood estimation is implemented, and hence this will always be the maximized log-likelihood.
}
\item{objective_name}{
The name of the objective that was maximized during fitting. Currently, only maximum-likelihood estimation is implemented, and hence this will always be ``loglikelihood''.
}
\item{times}{
Numeric vector of size 2, listing the two time points at which the diffusivity was estimated (\code{time1} and \code{time2}).
}
\item{diffusivities}{
Numeric vector of size 2, listing the fitted diffusivity at \code{time1} and \code{time2}. The fitted model assumes that the diffusivity varied linearly between those two time points.
}
\item{loglikelihood}{
The log-likelihood of the fitted linear model for the given data.
}
\item{NFP}{
Integer, number of fitted (i.e., non-fixed) model parameters. Will always be 2.
}
\item{Ncontrasts}{
Integer, number of independent contrasts used for fitting.
}
\item{AIC}{
The Akaike Information Criterion for the fitted model, defined as \eqn{2k-2\log(L)}, where \eqn{k} is the number of fitted parameters and \eqn{L} is the maximized likelihood.
}
\item{BIC}{
The Bayesian information criterion for the fitted model, defined as \eqn{\log(n)k-2\log(L)}, where \eqn{k} is the number of fitted parameters, \eqn{n} is the number of data points (number of independent contrasts), and \eqn{L} is the maximized likelihood.
}
\item{converged}{
Logical, specifying whether the maximum likelihood was reached after convergence of the optimization algorithm. Note that in some cases the maximum likelihood may have been achieved by an optimization path that did not yet converge (in which case it's advisable to increase \code{iter.max} and/or \code{eval.max}).
}
\item{Niterations}{
Integer, specifying the number of iterations performed during the optimization path that yielded the maximum likelihood.
}
\item{Nevaluations}{
Integer, specifying the number of likelihood evaluations performed during the optimization path that yielded the maximum likelihood.
}
\item{trial_start_objectives}{
Numeric vector of size \code{Ntrials}, listing the initial objective values (e.g., loglikelihoods) for each fitting trial, i.e. at the start parameter values.
}
\item{trial_objective_values}{
Numeric vector of size \code{Ntrials}, listing the final maximized objective values (e.g., loglikelihoods) for each fitting trial.
}
\item{trial_Nstart_attempts}{
Integer vector of size \code{Ntrials}, listing the number of start attempts for each fitting trial, until a starting point with valid likelihood was found.
}
\item{trial_Niterations}{
Integer vector of size \code{Ntrials}, listing the number of iterations needed for each fitting trial.
}
\item{trial_Nevaluations}{
Integer vector of size \code{Ntrials}, listing the number of likelihood evaluations needed for each fitting trial.
}
\item{standard_errors}{
Numeric vector of size 2, estimated standard error of the fitted diffusivity at the root and present, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50lower}{
Numeric vector of size 2, lower bound of the 50\% confidence interval (25-75\% percentile) for the fitted diffusivity at the root and present, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50upper}{
Numeric vector of size 2, upper bound of the 50\% confidence interval for the fitted diffusivity at the root and present, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95lower}{
Numeric vector of size 2, lower bound of the 95\% confidence interval (2.5-97.5\% percentile) for the fitted diffusivity at the root and present, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95upper}{
Numeric vector of size 2, upper bound of the 95\% confidence interval for the fitted diffusivity at the root and present, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{consistency}{
Numeric between 0 and 1, estimated consistency of the data with the fitted model. See the documentation of \code{\link{fit_sbm_const}} for an explanation. Only returned if \code{Nbootstraps>0}.
}
\item{significance}{
Numeric between 0 and 1, estimate statistical significance of the fitted linear slope. Only returned if \code{Nsignificance>0}.
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

S. Louca (2021). Phylogeographic estimation and simulation of global diffusive dispersal. Systematic Biology. 70:340-359.
}

\seealso{
\code{\link{simulate_sbm}},
\code{\link{fit_sbm_const}},
\code{\link{fit_sbm_parametric}},
\code{\link{fit_sbm_on_grid}}
}

\examples{
\dontrun{
# generate a random tree, keeping extinct lineages
tree_params = list(birth_rate_factor=1, death_rate_factor=0.95)
tree = generate_random_tree(tree_params,max_tips=1000,coalescent=FALSE)$tree

# calculate max distance of any tip from the root
max_time = get_tree_span(tree)$max_distance

# simulate time-dependent SBM on the tree
# we assume that diffusivity varies linearly with time
# in this example we measure distances in Earth radii
radius = 1
diffusivity_functor = function(times, params){
    return(params[1] + (times/max_time)*(params[2]-params[1]))
}
true_params = c(1, 2)
time_grid   = seq(0,max_time,length.out=2)
simulation  = simulate_sbm(tree,
                    radius = radius, 
                    diffusivity = diffusivity_functor(time_grid,true_params), 
                    time_grid = time_grid)

# fit time-independent SBM to get a rough estimate
fit_const = fit_sbm_const(tree,simulation$tip_latitudes,simulation$tip_longitudes,radius=radius)

# fit SBM model with linearly varying diffusivity
fit = fit_sbm_linear(tree,
            simulation$tip_latitudes,
            simulation$tip_longitudes,
            radius = radius,
            Ntrials = 10)
    
# compare fitted & true params
print(true_params)
print(fit$diffusivities)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SBM model}
\keyword{random}
\keyword{fitting}
\keyword{simulation}
