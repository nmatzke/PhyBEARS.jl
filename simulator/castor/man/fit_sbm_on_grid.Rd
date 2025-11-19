\name{fit_sbm_on_grid}
\alias{fit_sbm_on_grid}
\title{
Fit a phylogeographic Spherical Brownian Motion model with piecewise-linear diffusivity.
}
\description{
Given a rooted phylogenetic tree and geographic coordinates (latitudes & longitudes) for its tips, this function estimates the diffusivity of a Spherical Brownian Motion (SBM) model with time-dependent diffusivity for the evolution of geographic location along lineages (Perrin 1928; Brillinger 2012). Estimation is done via maximum-likelihood and using independent contrasts between sister lineages. This function is designed to estimate the diffusivity over time, approximated as a piecewise linear profile, by fitting the diffusivity on a discrete set of time points. The user thus provides a set of time points (\code{time_grid}), and \code{fit_sbm_on_grid} estimates the diffusivity on each time point, while assuming that the diffusivity varies linearly between time points.
}
\usage{
fit_sbm_on_grid(tree, 
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
              time_grid             = 0,
              guess_diffusivity     = NULL,
              min_diffusivity       = NULL,
              max_diffusivity       = Inf,
              Ntrials               = 1,
              Nthreads              = 1,
              Nbootstraps           = 0,
              Ntrials_per_bootstrap = NULL,
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
\item{time_grid}{
Numeric vector, specifying discrete time points (counted since the root) at which the \code{diffusivity} should be fitted; between these time points the diffusivity is assumed to vary linearly. This time grid should be fine enough to sufficiently capture the variation in the diffusivity over time, but must not be too big to avoid overfitting. If \code{NULL} or of size 1, then the diffusivity is assumed to be time-independent. Listed times must be strictly increasing, and should cover at least the full considered time interval (from 0 to the maximum distance of any tip from the root); otherwise, constant extrapolation is used to cover missing times. Note that time is measured in the same units as the tree's edge lengths.
}
\item{guess_diffusivity}{
Optional numeric vector, specifying a first guess for the diffusivity. Either of size 1 (the same first guess for all time points), or of the same length as \code{time_grid} (different first guess for each time point, \code{NA} are replaced with an automatically chosen first guess). If \code{NULL}, the first guess is chosen automatically.
}
\item{min_diffusivity}{
Optional numeric vector, specifying lower bounds for the fitted diffusivity. Either of size 1 (the same lower bound is assumed for all time points), or of the same length as \code{time_grid} (different lower bound for each time point, \code{NA} are replaced with an automatically chosen lower bound). If \code{NULL}, lower bounds are chosen automatically. 
}
\item{max_diffusivity}{
Optional numeric vector, specifying upper bounds for the fitted diffusivity. Either of size 1 (the same upper bound is assumed for all time points), or of the same length as \code{time_grid} (different upper bound for each time point, \code{NA} are replaced with infinity). If \code{NULL}, no upper bound is imposed.
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
Logical, specifying whether to print progress reports and warnings to the screen. Note that errors always cause a return of the function (see return values \code{success} and \code{error}).
}
\item{verbose_prefix}{
Character, specifying the line prefix for printing progress reports to the screen.
}
}


\details{
This function is designed to estimate the diffusivity profile over time, approximated by a piecewise linear function. Fitting is done by maximizing the likelihood of observing the given tip coordinates under the SBM model. Internally, this function uses \code{\link{fit_sbm_parametric}}.

It is generally advised to provide as much information to the function \code{fit_sbm_on_grid} as possible, including reasonable lower and upper bounds (\code{min_diffusivity} and \code{max_diffusivity}). It is important that the \code{time_grid} is sufficiently fine to capture the variation of the true diffusivity over time, since the likelihood is calculated under the assumption that the diffusivity varies linearly between grid points. However, depending on the size of the tree, the grid size must not be too large, since otherwise overfitting becomes very likely. The \code{time_grid} does not need to be uniform, i.e., you may want to use a finer grid in regions where there's more data (tips) available.

Note that estimation of diffusivity at older times is only possible if the timetree includes extinct tips or tips sampled at older times (e.g., as is often the case in viral phylogenies). If tips are only sampled once at present-day, i.e. the timetree is ultrametric, reliable diffusivity estimates can only be achieved near present times. If the tree is ultrametric, you should consider using \code{\link{fit_sbm_const}} instead.

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
\item{time_grid}{
Numeric vector, the time-grid on which the diffusivity was fitted.
}
\item{diffusivity}{
Numeric vector of size Ngrid (length of \code{time_grid}), listing the fitted diffusivities at the various time-grid points.
}
\item{loglikelihood}{
The log-likelihood of the fitted model for the given data.
}
\item{NFP}{
Integer, number of fitted (i.e., non-fixed) model parameters.
}
\item{Ncontrasts}{
Integer, number of independent contrasts used for fitting.
}
\item{phylodistances}{
Numeric vector of length Ncontrasts, listing phylogenetic (patristic) distances of the independent contrasts.
}
\item{geodistances}{
Numeric vector of length Ncontrasts, listing geographic (great circle) distances of the independent contrasts.
}
\item{child_times1}{
Numeric vector of length Ncontrasts, listing the times (distance from root) of the first tip in each independent contrast.
}
\item{child_times2}{
Numeric vector of length Ncontrasts, listing the times (distance from root) of the second tip in each independent contrast.
}
\item{MRCA_times}{
Numeric vector of length Ncontrasts, listing the times (distance from root) of the MRCA of the two tips in each independent contrast.
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
Numeric vector of size NP, estimated standard error of the parameters, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{medians}{
Numeric vector of size NP, median the estimated parameters across parametric bootstraps. Only returned if \code{Nbootstraps>0}.
}
\item{CI50lower}{
Numeric vector of size NP, lower bound of the 50\% confidence interval (25-75\% percentile) for the parameters, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50upper}{
Numeric vector of size NP, upper bound of the 50\% confidence interval for the parameters, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95lower}{
Numeric vector of size NP, lower bound of the 95\% confidence interval (2.5-97.5\% percentile) for the parameters, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95upper}{
Numeric vector of size NP, upper bound of the 95\% confidence interval for the parameters, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{consistency}{
Numeric between 0 and 1, estimated consistency of the data with the fitted model. See the documentation of \code{\link{fit_sbm_const}} for an explanation.
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
\code{\link{fit_sbm_linear}}
}

\examples{
\dontrun{
# generate a random tree, keeping extinct lineages
tree_params = list(birth_rate_factor=1, death_rate_factor=0.95)
tree = generate_random_tree(tree_params,max_tips=2000,coalescent=FALSE)$tree

# calculate max distance of any tip from the root
max_time = get_tree_span(tree)$max_distance

# simulate time-dependent SBM on the tree
# using a diffusivity that varies roughly exponentially with time
# In this example we measure distances in Earth radii
radius = 1
fine_time_grid = seq(from=0, to=max_time, length.out=10)
fine_D = 0.01 + 0.03*exp(-2*fine_time_grid/max_time)
simul = simulate_sbm(tree, 
                     radius     = radius, 
                     diffusivity= fine_D, 
                     time_grid  = fine_time_grid)

# fit time-dependent SBM on a time-grid of size 4
fit = fit_sbm_on_grid(tree,
            simul$tip_latitudes,
            simul$tip_longitudes,
            radius    = radius,
            time_grid = seq(from=0,to=max_time,length.out=4),
            Nthreads  = 3,  # use 3 CPUs
            Ntrials   = 30) # avoid local optima through multiple trials
    
# visually compare fitted & true params
plot(x      = fine_time_grid,
     y      = fine_D,
     type   = 'l',
     col    = 'black',
     xlab   = 'time',
     ylab   = 'D',
     ylim   = c(0,max(fine_D)))
lines(x     = fit$time_grid,
      y     = fit$diffusivity,
      type  = 'l',
      col   = 'blue')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SBM model}
\keyword{random}
\keyword{fitting}
\keyword{simulation}
