\name{fit_hbd_psr_on_grid}
\alias{fit_hbd_psr_on_grid}
\title{
Fit pulled speciation rates of birth-death models on a time grid.
}
\description{
Given an ultrametric timetree, estimate the pulled speciation rate of homogenous birth-death (HBD) models that best explains the tree via maximum likelihood. Every HBD model is defined by some speciation and extinction rates (\eqn{\lambda} and \eqn{\mu}) over time, as well as the sampling fraction \eqn{\rho} (fraction of extant species sampled). ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction rates. For any given HBD model there exists an infinite number of alternative HBD models that predict the same deterministic lineages-through-time curve and yield the same likelihood for any given reconstructed timetree; these ``congruent'' models cannot be distinguished from one another solely based on the tree. 

Each congruence class is uniquely described by the ``pulled speciation rate'' (PSR), defined as the relative slope of the deterministic LTT over time, \eqn{PSR=-M^{-1}dM/d\tau} (where \eqn{\tau} is time before present). In other words, two HBD models are congruent if and only if they have the same PSR. This function is designed to estimate the generating congruence class for the tree, by fitting the PSR on a discrete time grid.
}
\usage{
fit_hbd_psr_on_grid(  tree, 
                      oldest_age            = NULL,
                      age0                  = 0,
                      age_grid              = NULL,
                      min_PSR               = 0,
                      max_PSR               = +Inf,
                      guess_PSR             = NULL,
                      fixed_PSR             = NULL,
                      splines_degree        = 1,
                      condition             = "auto",
                      relative_dt           = 1e-3,
                      Ntrials               = 1,
                      Nbootstraps           = 0,
                      Ntrials_per_bootstrap = NULL,
                      Nthreads              = 1,
                      max_model_runtime     = NULL,
                      fit_control           = list(),
                      verbose               = FALSE,
                      diagnostics           = FALSE,
                      verbose_prefix        = "")
}
\arguments{
\item{tree}{
A rooted ultrametric timetree of class "phylo", representing the time-calibrated phylogeny of a set of extant sampled species.
}
\item{oldest_age}{
Strictly positive numeric, specifying the oldest time before present (``age'') to consider when calculating the likelihood. If this is equal to or greater than the root age, then \code{oldest_age} is taken as the stem age, and the classical formula by Morlon et al. (2011) is used.  If \code{oldest_age} is less than the root age, the tree is split into multiple subtrees at that age by treating every edge crossing that age as the stem of a subtree, and each subtree is considered an independent realization of the HBD model stemming at that age. This can be useful for avoiding points in the tree close to the root, where estimation uncertainty is generally higher. If \code{oldest_age==NULL}, it is automatically set to the root age.
}
\item{age0}{
Non-negative numeric, specifying the youngest age (time before present) to consider for fitting. If \code{age0>0}, the tree essentially is trimmed at \code{age0}, omitting anything younger than \code{age0}, and the PSR is fitted to the trimmed tree while shifting time appropriately.
}
\item{age_grid}{
Numeric vector, listing ages in ascending order at which the PSR is allowed to vary independently. This grid must cover at least the age range from \code{age0} to \code{oldest_age}. If \code{NULL} or of length <=1 (regardless of value), then the PSR is assumed to be time-independent.
}
\item{min_PSR}{
Numeric vector of length Ngrid (=\code{max(1,length(age_grid))}), or a single numeric, specifying lower bounds for the fitted PSR at each point in the age grid. If a single numeric, the same lower bound applies at all ages. Note that the PSR is never negative.
}
\item{max_PSR}{
Numeric vector of length Ngrid, or a single numeric, specifying upper bounds for the fitted PSR at each point in the age grid. If a single numeric, the same upper bound applies at all ages. Use \code{+Inf} to omit upper bounds.
}
\item{guess_PSR}{
Initial guess for the PSR at each age-grid point. Either \code{NULL} (an initial guess will be computed automatically), or a single numeric (guessing the same PSR at all ages) or a numeric vector of size Ngrid specifying a separate guess at each age-grid point. To omit an initial guess for some but not all age-grid points, set their guess values to \code{NA}. Guess values are ignored for non-fitted (i.e., fixed) parameters.
}
\item{fixed_PSR}{
Optional fixed (i.e. non-fitted) PSR values on one or more age-grid points. Either \code{NULL} (PSR is not fixed anywhere), or a single numeric (PSR fixed to the same value at all grid points) or a numeric vector of size Ngrid (PSR fixed at one or more age-grid points, use \code{NA} for non-fixed values).
}
\item{splines_degree}{
Integer between 0 and 3 (inclusive), specifying the polynomial degree of the PSR between age-grid points. If 0, then the PSR is considered piecewise constant, if 1 then the PSR is considered piecewise linear, if 2 or 3 then the PSR is considered to be a spline of degree 2 or 3, respectively. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on. A degree of 0 is generally not recommended.
}
\item{condition}{
Character, either "crown", "stem", "auto", "stemN" or "crownN" (where N is an integer >=2), specifying on what to condition the likelihood. If "crown", the likelihood is conditioned on the survival of the two daughter lineages branching off at the root at that time. If "stem", the likelihood is conditioned on the survival of the stem lineage, with the process having started at \code{oldest_age}. Note that "crown" and "crownN"" really only make sense when \code{oldest_age} is equal to the root age, while "stem" is recommended if \code{oldest_age} differs from the root age. If "stem2", the condition is that the process yielded at least two sampled tips, and similarly for "stem3" etc. If "crown3", the condition is that a splitting occurred at the root age, both child clades survived, and in total yielded at least 3 sampled tips (and similarly for "crown4" etc). If "auto", the condition is chosen according to the recommendations mentioned earlier.
}
\item{relative_dt}{
Strictly positive numeric (unitless), specifying the maximum relative time step allowed for integration over time, when calculating the likelihood. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
}
\item{Ntrials}{
Integer, specifying the number of independent fitting trials to perform, each starting from a random choice of model parameters. Increasing \code{Ntrials} reduces the risk of reaching a non-global local maximum in the fitting objective.
}
\item{Nbootstraps}{
Integer, specifying an optional number of bootstrap samplings to perform, for estimating standard errors and confidence intervals of maximum-likelihood fitted parameters. If 0, no bootstrapping is performed. Typical values are 10-100. At each bootstrap sampling, a random timetree is generated under the birth-death model according to the fitted PSR, the parameters are estimated anew based on the generated tree, and subsequently compared to the original fitted parameters. Each bootstrap sampling will use roughly the same information and similar computational resources as the original maximum-likelihood fit (e.g., same number of trials, same optimization parameters, same initial guess, etc).
}
\item{Ntrials_per_bootstrap}{
Integer, specifying the number of fitting trials to perform for each bootstrap sampling. If \code{NULL}, this is set equal to \code{max(1,Ntrials)}. Decreasing \code{Ntrials_per_bootstrap} will reduce computation time, at the expense of potentially inflating the estimated confidence intervals; in some cases (e.g., for very large trees) this may be useful if fitting takes a long time and confidence intervals are very narrow anyway. Only relevant if \code{Nbootstraps>0}.
}
\item{Nthreads}{
Integer, specifying the number of parallel threads to use for performing multiple fitting trials simultaneously. This should generally not exceed the number of available CPUs on your machine. Parallel computing is not available on the Windows platform.
}
\item{max_model_runtime}{
Optional numeric, specifying the maximum number of seconds to allow for each evaluation of the likelihood function. Use this to abort fitting trials leading to parameter regions where the likelihood takes a long time to evaluate (these are often unlikely parameter regions).
}
\item{fit_control}{
 Named list containing options for the \code{nlminb} optimization routine, such as \code{iter.max}, \code{eval.max} or \code{rel.tol}. For a complete list of options and default values see the documentation of \code{nlminb} in the \code{stats} package.
}
\item{verbose}{
Logical, specifying whether to print progress reports and warnings to the screen. Note that errors always cause a return of the function (see return values \code{success} and \code{error}).
}
\item{diagnostics}{
Logical, specifying whether to print detailed information (such as model likelihoods) at every iteration of the fitting routine. For debugging purposes mainly.
}
\item{verbose_prefix}{
Character, specifying the line prefix for printing progress reports to the screen.
}
}


\details{
It is generally advised to provide as much information to the function \code{fit_hbd_psr_on_grid} as possible, including reasonable lower and upper bounds (\code{min_PSR} and \code{max_PSR}) and a reasonable parameter guess (\code{guess_PSR}). It is also important that the \code{age_grid} is sufficiently fine to capture the expected major variations of the PSR over time, but keep in mind the serious risk of overfitting when \code{age_grid} is too fine and/or the tree is too small.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether model fitting succeeded. If \code{FALSE}, the returned list will include an additional ``error'' element (character) providing a description of the error; in that case all other return variables may be undefined.
}
\item{objective_value}{
The maximized fitting objective. Currently, only maximum-likelihood estimation is implemented, and hence this will always be the maximized log-likelihood.
}
\item{objective_name}{
The name of the objective that was maximized during fitting. Currently, only maximum-likelihood estimation is implemented, and hence this will always be ``loglikelihood''.
}
\item{loglikelihood}{
The log-likelihood of the fitted model for the given timetree.
}
\item{fitted_PSR}{
Numeric vector of size Ngrid, listing fitted or fixed pulled speciation rates (PSR) at each age-grid point. Between grid points the fitted PSR should be interpreted as a piecewise polynomial function (natural spline) of degree \code{splines_degree}; to evaluate this function at arbitrary ages use the \code{castor} routine \code{\link{evaluate_spline}}.
}
\item{guess_PSR}{
Numeric vector of size Ngrid, specifying the initial guess for the PSR at each age-grid point.
}
\item{age_grid}{
The age-grid on which the PSR is defined. This will be the same as the provided \code{age_grid}, unless the latter was \code{NULL} or of length <=1.
}
\item{NFP}{
Integer, number of fitted (i.e., non-fixed) parameters. If none of the PSRs were fixed, this will be equal to Ngrid.
}
\item{AIC}{
The Akaike Information Criterion for the fitted model, defined as \eqn{2k-2\log(L)}, where \eqn{k} is the number of fitted parameters, and \eqn{L} is the maximized likelihood.
}
\item{BIC}{
The Bayesian information criterion for the fitted model, defined as \eqn{\log(n)k-2\log(L)}, where \eqn{k} is the number of fitted parameters, \eqn{n} is the number of data points (number of branching times), and \eqn{L} is the maximized likelihood.
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
\item{bootstrap_estimates}{
If \code{Nbootstraps>0}, this will be a numeric matrix of size \code{Nbootstraps} x Ngrid, listing the fitted PSR at each grid point and for each bootstrap.
}
\item{standard_errors}{
If \code{Nbootstraps>0}, this will be a numeric vector of size NGrid, listing bootstrap-estimated standard errors for the fitted PSR at each grid point.
}
\item{CI50lower}{
If \code{Nbootstraps>0}, this will be a numeric vector of size Ngrid, listing bootstrap-estimated lower bounds of the 50-percent confidence intervals for the fitted PSR at each grid point.
}
\item{CI50upper}{
Similar to \code{CI50lower}, listing upper bounds of 50-percentile confidence intervals.
}
\item{CI95lower}{
Similar to \code{CI50lower}, listing lower bounds of 95-percentile confidence intervals.
}
\item{CI95upper}{
Similar to \code{CI95lower}, listing upper bounds of 95-percentile confidence intervals.
}
}


\author{Stilianos Louca}

\references{
S. Louca et al. (2018). Bacterial diversification through geological time. Nature Ecology & Evolution. 2:1458-1467.

S. Louca and M. W. Pennell (2020). Extant timetrees are consistent with a myriad of diversification histories. Nature. 580:502-505.
}

\seealso{
\code{\link{simulate_deterministic_hbd}}

\code{\link{loglikelihood_hbd}}

\code{\link{fit_hbd_model_parametric}}

\code{\link{fit_hbd_model_on_grid}}

\code{\link{fit_hbd_pdr_parametric}}

\code{\link{fit_hbd_pdr_on_grid}}

\code{\link{fit_hbd_psr_on_best_grid_size}}

\code{\link{model_adequacy_hbd}}
}

\examples{
\dontrun{
# Generate a random tree with exponentially varying lambda & mu
Ntips     = 10000
rho       = 0.5 # sampling fraction
time_grid = seq(from=0, to=100, by=0.01)
lambdas   = 2*exp(0.1*time_grid)
mus       = 1.5*exp(0.09*time_grid)
sim       = generate_random_tree( parameters  = list(rarefaction=rho), 
                                  max_tips    = Ntips/rho,
                                  coalescent  = TRUE,
                                  added_rates_times     = time_grid,
                                  added_birth_rates_pc  = lambdas,
                                  added_death_rates_pc  = mus)
tree = sim$tree
root_age = castor::get_tree_span(tree)$max_distance
cat(sprintf("Tree has \%d tips, spans \%g Myr\n",length(tree$tip.label),root_age))

# Fit PSR on grid
oldest_age=root_age/2 # only consider recent times when fitting
Ngrid     = 10
age_grid  = seq(from=0,to=oldest_age,length.out=Ngrid)
fit = fit_hbd_psr_on_grid(tree, 
          oldest_age  = oldest_age,
          age_grid    = age_grid,
          min_PSR     = 0,
          max_PSR     = +100,
          condition   = "crown",
          Ntrials     = 10,
          Nthreads    = 4,
          max_model_runtime = 1) 	# limit model evaluation to 1 second
if(!fit$success){
  cat(sprintf("ERROR: Fitting failed: \%s\n",fit$error))
}else{
  cat(sprintf("Fitting succeeded:\nLoglikelihood=\%g\n",fit$loglikelihood))
  # plot fitted PSR
  plot( x     = fit$age_grid,
        y     = fit$fitted_PSR,
        main  = 'Fitted PSR',
        xlab  = 'age',
        ylab  = 'PSR',
        type  = 'b',
        xlim  = c(root_age,0)) 
        
 # plot deterministic LTT of fitted model
  plot( x     = fit$age_grid,
        y     = fit$fitted_LTT,
        main  = 'Fitted dLTT',
        xlab  = 'age',
        ylab  = 'lineages',
        type  = 'b',
        log   = 'y',
        xlim  = c(root_age,0)) 

  # get fitted PSR as a function of age
  PSR_fun = approxfun(x=fit$age_grid, y=fit$fitted_PSR)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{birth-death model}
\keyword{pulled speciation rate}
