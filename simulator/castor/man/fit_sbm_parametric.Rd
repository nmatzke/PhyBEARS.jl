\name{fit_sbm_parametric}
\alias{fit_sbm_parametric}
\title{
Fit a time-dependent phylogeographic Spherical Brownian Motion model.
}
\description{
Given a rooted phylogenetic tree and geographic coordinates (latitudes & longitudes) for its tips, this function estimates the diffusivity of a Spherical Brownian Motion (SBM) model with time-dependent diffusivity for the evolution of geographic location along lineages (Perrin 1928; Brillinger 2012). Estimation is done via maximum-likelihood and using independent contrasts between sister lineages. This function is designed to estimate the diffusivity over time, by fitting a finite number of parameters defining the diffusivity as a function of time. The user thus provides the general functional form of the diffusivity that depends on time and NP parameters, and \code{fit_sbm_parametric} estimates each of the free parameters.
}
\usage{
fit_sbm_parametric(tree, 
              tip_latitudes,
              tip_longitudes,
              radius,
              param_values,
              param_guess,
              diffusivity,
              time_grid             = NULL,
              clade_states          = NULL,
              planar_approximation  = FALSE,
              only_basal_tip_pairs  = FALSE,
              only_distant_tip_pairs= FALSE,
              min_MRCA_time         = 0,
              max_MRCA_age          = Inf,
              max_phylodistance     = Inf,
              no_state_transitions  = FALSE,
              only_state            = NULL,
              param_min             = -Inf,
              param_max             = +Inf,
              param_scale           = NULL,
              Ntrials               = 1,
              max_start_attempts    = 1,
              Nthreads              = 1,
              Nbootstraps           = 0,
              Ntrials_per_bootstrap = NULL,
              NQQ                   = 0,
              fit_control           = list(),
              SBM_PD_functor        = NULL,
              focal_param_values    = NULL,
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
\item{param_values}{
Numeric vector of length NP, specifying fixed values for a some or all model parameters. For fitted (i.e., non-fixed) parameters, use \code{NaN} or \code{NA}. For example, the vector \code{c(1.5,NA,40)} specifies that the 1st and 3rd model parameters are fixed at the values 1.5 and 40, respectively, while the 2nd parameter is to be fitted. The length of this vector defines the total number of model parameters. If entries in this vector are named, the names are taken as parameter names. Names should be included if you'd like returned parameter vectors to have named entries, or if the \code{diffusivity} function queries parameter values by name (as opposed to numeric index).
}
\item{param_guess}{
Numeric vector of size NP, specifying a first guess for the value of each model parameter. For fixed parameters, guess values are ignored. Can be \code{NULL} only if all model parameters are fixed.
}
\item{diffusivity}{
Function specifying the diffusivity at any given time (time since the root) and for any given parameter values. This function must take exactly two arguments, the 1st one being a numeric vector (one or more times) and the 2nd one being a numeric vector of size NP (parameter values), and return a numeric vector of the same size as the 1st argument.
}
\item{time_grid}{
Numeric vector, specifying times (counted since the root) at which the \code{diffusivity} function should be evaluated. This time grid must be fine enough to capture the possible variation in the diffusivity over time, within the permissible parameter range. If of size 1, then the diffusivity is assumed to be time-independent.
Listed times must be strictly increasing, and should cover at least the full considered time interval (from 0 to the maximum distance of any tip from the root); otherwise, constant extrapolation is used to cover missing times. Can also be \code{NULL} or a vector of size 1, in which case the diffusivity is assumed to be time-independent. Note that time is measured in the same units as the tree's edge lengths.
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
\item{param_min}{
Optional numeric vector of size NP, specifying lower bounds for model parameters. If of size 1, the same lower bound is applied to all parameters. Use \code{-Inf} to omit a lower bound for a parameter. If \code{NULL}, no lower bounds are applied. For fixed parameters, lower bounds are ignored.
}
\item{param_max}{
Optional numeric vector of size NP, specifying upper bounds for model parameters. If of size 1, the same upper bound is applied to all parameters. Use \code{+Inf} to omit an upper bound for a parameter. If \code{NULL}, no upper bounds are applied. For fixed parameters, upper bounds are ignored.
}
\item{param_scale}{
Optional numeric vector of size NP, specifying typical scales for model parameters. If of size 1, the same scale is assumed for all parameters. If \code{NULL}, scales are determined automatically. For fixed parameters, scales are ignored. It is strongly advised to provide reasonable scales, as this facilitates the numeric optimization algorithm.
}
\item{Ntrials}{
Integer, specifying the number of independent fitting trials to perform, each starting from a random choice of model parameters. Increasing \code{Ntrials} reduces the risk of reaching a non-global local maximum in the fitting objective.
}
\item{max_start_attempts}{
Integer, specifying the number of times to attempt finding a valid start point (per trial) before giving up on that trial. Randomly choosen extreme start parameters may occasionally result in Inf/undefined likelihoods, so this option allows the algorithm to keep looking for valid starting points.
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
\item{focal_param_values}{
Optional numeric matrix having NP columns and an arbitrary number of rows, listing combinations of parameter values of particular interest and for which the log-likelihoods should be returned. This may be used e.g. for diagnostic purposes, e.g. to examine the shape of the likelihood function.
}
\item{verbose}{
Logical, specifying whether to print progress reports and warnings to the screen. Note that errors always cause a return of the function (see return values \code{success} and \code{error}).
}
\item{verbose_prefix}{
Character, specifying the line prefix for printing progress reports to the screen.
}
}


\details{
This function is designed to estimate a finite set of scalar parameters (\eqn{p_1,..,p_n\in\R}) that determine the diffusivity over time, by maximizing the likelihood of observing the given tip coordinates under the SBM model. For example, the investigator may assume that the diffusivity exponentially over time, i.e. can be described by \eqn{D(t)=A\cdot e^{-B t}} (where \eqn{A} and \eqn{B} are unknown coefficients and \eqn{t} is time since the root). In this case the model has 2 free parameters, \eqn{p_1=A} and \eqn{p_2=B}, each of which may be fitted to the tree.

It is generally advised to provide as much information to the function \code{fit_sbm_parametric} as possible, including reasonable lower and upper bounds (\code{param_min} and \code{param_max}), a reasonable parameter guess (\code{param_guess}) and reasonable parameter scales \code{param_scale}. If some model parameters can vary over multiple orders of magnitude, it is advised to transform them so that they vary across fewer orders of magnitude (e.g., via log-transformation). It is also important that the \code{time_grid} is sufficiently fine to capture the variation of the diffusivity over time, since the likelihood is calculated under the assumption that the diffusivity varies linearly between grid points.

Estimation of diffusivity at older times is only possible if the timetree includes extinct tips or tips sampled at older times (e.g., as is often the case in viral phylogenies). If tips are only sampled once at present-day, i.e. the timetree is ultrametric, reliable diffusivity estimates can only be achieved near present times. If the tree is ultrametric, you should consider using \code{\link{fit_sbm_const}} instead.

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
\item{param_fitted}{
Numeric vector of size NP (number of model parameters), listing all fitted or fixed model parameters in their standard order (see details above).
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
\item{guess_loglikelihood}{
The loglikelihood of the data for the initial parameter guess (\code{param_guess}).
}
\item{focal_loglikelihoods}{
A numeric vector of the same size as \code{nrow(focal_param_values)}, listing loglikelihoods for each of the focal parameter conbinations listed in \code{focal_loglikelihoods}.
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
\code{\link{fit_sbm_linear}}
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
                      radius      = radius, 
                      diffusivity = diffusivity_functor(time_grid,true_params), 
                      time_grid   = time_grid)

# fit time-independent SBM to get a rough estimate
fit_const = fit_sbm_const(tree,simulation$tip_latitudes,simulation$tip_longitudes,radius=radius)

# fit time-dependent SBM, i.e. fit the 2 parameters of the linear form
fit = fit_sbm_parametric(tree,
            simulation$tip_latitudes,
            simulation$tip_longitudes,
            radius = radius,
            param_values = c(NA,NA),
            param_guess = c(fit_const$diffusivity,fit_const$diffusivity),
            diffusivity = diffusivity_functor,
            time_grid = time_grid,
            Ntrials = 10)
    
# compare fitted & true params
print(true_params)
print(fit$param_fitted)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SBM model}
\keyword{random}
\keyword{fitting}
\keyword{simulation}
