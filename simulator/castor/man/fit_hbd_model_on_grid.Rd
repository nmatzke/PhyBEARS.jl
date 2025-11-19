\name{fit_hbd_model_on_grid}
\alias{fit_hbd_model_on_grid}
\title{
Fit a homogenous birth-death model on a discrete time grid.
}
\description{
Given an ultrametric timetree, fit a homogenous birth-death (HBD) model in which speciation and extinction rates (\eqn{\lambda} and \eqn{mu}) are defined on a fixed grid of discrete time points and assumed to vary polynomially between grid points. ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction rates (in the literature this is sometimes referred to simply as ``birth-death model''). Every HBD model is defined based on the values that \eqn{\lambda} and \eqn{\mu} take over time as well as the sampling fraction \eqn{\rho} (fraction of extant species sampled). This function estimates the values of \eqn{\lambda} and \eqn{\mu} at each grid point by maximizing the likelihood (Morlon et al. 2011) of the timetree under the resulting HBD model.
}
\usage{
fit_hbd_model_on_grid(tree, 
                      oldest_age        = NULL,
                      age0              = 0,
                      age_grid          = NULL,
                      min_lambda        = 0,
                      max_lambda        = +Inf,
                      min_mu            = 0,
                      max_mu            = +Inf,
                      min_rho0          = 1e-10,
                      max_rho0          = 1,
                      guess_lambda      = NULL,
                      guess_mu          = NULL,
                      guess_rho0        = 1,
                      fixed_lambda      = NULL,
                      fixed_mu          = NULL,
                      fixed_rho0        = NULL,
                      const_lambda      = FALSE,
                      const_mu          = FALSE,
                      splines_degree    = 1,
                      condition         = "auto",
                      relative_dt       = 1e-3,
                      Ntrials           = 1,
                      Nthreads          = 1,
                      max_model_runtime = NULL,
                      fit_control       = list())
}
\arguments{
\item{tree}{
A rooted ultrametric timetree of class "phylo", representing the time-calibrated reconstructed phylogeny of a set of extant sampled species.
}
\item{oldest_age}{
Strictly positive numeric, specifying the oldest time before present (``age'') to consider when calculating the likelihood. If this is equal to or greater than the root age, then \code{oldest_age} is taken as the stem age, and the classical formula by Morlon et al. (2011) is used.  If \code{oldest_age} is less than the root age, the tree is split into multiple subtrees at that age by treating every edge crossing that age as the stem of a subtree, and each subtree is considered an independent realization of the HBD model stemming at that age. This can be useful for avoiding points in the tree close to the root, where estimation uncertainty is generally higher. If \code{oldest_age==NULL}, it is automatically set to the root age.
}
\item{age0}{
Non-negative numeric, specifying the youngest age (time before present) to consider for fitting, and with respect to which \code{rho} is defined. If \code{age0>0}, then \code{rho0} refers to the sampling fraction at age \code{age0}, i.e. the fraction of lineages extant at \code{age0} that are included in the tree. See below for more details.
}
\item{age_grid}{
Numeric vector, listing ages in ascending order, on which \eqn{\lambda} and \eqn{\mu} are allowed to vary independently. This grid must cover at least the age range from \code{age0} to \code{oldest_age}. If \code{NULL} or of length <=1 (regardless of value), then \eqn{\lambda} and \eqn{\mu} are assumed to be time-independent.
}
\item{min_lambda}{
Numeric vector of length Ngrid (=\code{max(1,length(age_grid))}), or a single numeric, specifying lower bounds for the fitted \eqn{\lambda} at each point in the age grid. If a single numeric, the same lower bound applies at all ages.
}
\item{max_lambda}{
Numeric vector of length Ngrid, or a single numeric, specifying upper bounds for the fitted \eqn{\lambda} at each point in the age grid. If a single numeric, the same upper bound applies at all ages. Use \code{+Inf} to omit upper bounds.
}
\item{min_mu}{
Numeric vector of length Ngrid, or a single numeric, specifying lower bounds for the fitted \eqn{\mu} at each point in the age grid. If a single numeric, the same lower bound applies at all ages.
}
\item{max_mu}{
Numeric vector of length Ngrid, or a single numeric, specifying upper bounds for the fitted \eqn{\mu} at each point in the age grid. If a single numeric, the same upper bound applies at all ages. Use \code{+Inf} to omit upper bounds.
}
\item{min_rho0}{
Numeric, specifying a lower bound for the fitted sampling fraction \eqn{\rho} (fraction of extant species included in the tree).
}
\item{max_rho0}{
Numeric, specifying an upper bound for the fitted sampling fraction \eqn{\rho}.
}
\item{guess_lambda}{
Initial guess for \eqn{\lambda} at each age-grid point. Either \code{NULL} (an initial guess will be computed automatically), or a single numeric (guessing the same \eqn{\lambda} at all ages) or a numeric vector of size Ngrid specifying a separate guess for \eqn{\lambda} at each age-grid point. To omit an initial guess for some but not all age-grid points, set their guess values to \code{NA}. Guess values are ignored for non-fitted (i.e., fixed) parameters.
}
\item{guess_mu}{
Initial guess for \eqn{\mu} at each age-grid point. Either \code{NULL} (an initial guess will be computed automatically), or a single numeric (guessing the same \eqn{\mu} at all ages) or a numeric vector of size Ngrid specifying a separate guess for \eqn{\mu} at each age-grid point. To omit an initial guess for some but not all age-grid points, set their guess values to \code{NA}. Guess values are ignored for non-fitted (i.e., fixed) parameters.
}
\item{guess_rho0}{
Numeric, specifying an initial guess for the sampling fraction \eqn{\rho} at \code{age0}. Setting this to \code{NULL} or \code{NA} is the same as setting it to 1.
}
\item{fixed_lambda}{
Optional fixed (i.e. non-fitted) \eqn{\lambda} values on one or more age-grid points. Either \code{NULL} (\eqn{\lambda} is not fixed anywhere), or a single numeric (\eqn{\lambda} fixed to the same value at all grid points) or a numeric vector of size Ngrid (\eqn{\lambda} fixed on one or more age-grid points, use \code{NA} for non-fixed values).
}
\item{fixed_mu}{
Optional fixed (i.e. non-fitted) \eqn{\mu} values on one or more age-grid points. Either \code{NULL} (\eqn{\mu} is not fixed anywhere), or a single numeric (\eqn{\mu} fixed to the same value at all grid points) or a numeric vector of size Ngrid (\eqn{\mu} fixed on one or more age-grid points, use \code{NA} for non-fixed values).
}
\item{fixed_rho0}{
Numeric between 0 and 1, optionallly specifying a fixed value for the sampling fraction \eqn{\rho}. If \code{NULL} or \code{NA}, the sampling fraction \eqn{\rho} is estimated, however note that this may not always be meaningful (Stadler 2009, Stadler 2013).
}
\item{const_lambda}{
Logical, specifying whether \eqn{\lambda} should be assumed constant across the grid, i.e. time-independent. Setting \code{const_lambda=TRUE} reduces the number of free (i.e., independently fitted) parameters.
If \eqn{\lambda} is fixed on some grid points (i.e. via \code{fixed_lambda}), then only the non-fixed lambdas are assumed to be identical to one another.
}
\item{const_mu}{
Logical, specifying whether \eqn{\mu} should be assumed constant across the grid, i.e. time-independent. Setting \code{const_mu=TRUE} reduces the number of free (i.e., independently fitted) parameters. If \eqn{\mu} is fixed on some grid points (i.e. via \code{fixed_mu}), then only the non-fixed mus are assumed to be identical to one another.
}
\item{splines_degree}{
Integer between 0 and 3 (inclusive), specifying the polynomial degree of \eqn{\lambda} and \eqn{\mu} between age-grid points. If 0, then \eqn{\lambda} and \eqn{\mu} are considered piecewise constant, if 1 then \eqn{\lambda} and \eqn{\mu} are considered piecewise linear, if 2 or 3 then \eqn{\lambda} and \eqn{\mu} are considered to be splines of degree 2 or 3, respectively. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on. A degree of 0 is generally not recommended. The case \code{splines_degree=0} is also known as ``skyline'' model.
}
\item{condition}{
Character, either "crown", "stem", "auto", "stemN" or "crownN" (where N is an integer >=2), specifying on what to condition the likelihood. If "crown", the likelihood is conditioned on the survival of the two daughter lineages branching off at the root at that time. If "stem", the likelihood is conditioned on the survival of the stem lineage, with the process having started at \code{oldest_age}. Note that "crown" and "crownN"" really only make sense when \code{oldest_age} is equal to the root age, while "stem" is recommended if \code{oldest_age} differs from the root age. If "stem2", the condition is that the process yielded at least two sampled tips, and similarly for "stem3" etc. If "crown3", the condition is that a splitting occurred at the root age, both child clades survived, and in total yielded at least 3 sampled tips (and similarly for "crown4" etc). If "auto", the condition is chosen according to the recommendations mentioned earlier. "none" is generally not recommended.
}
\item{relative_dt}{
Strictly positive numeric (unitless), specifying the maximum relative time step allowed for integration over time, when calculating the likelihood. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
}
\item{Ntrials}{
Integer, specifying the number of independent fitting trials to perform, each starting from a random choice of model parameters. Increasing \code{Ntrials} reduces the risk of reaching a non-global local maximum in the fitting objective.
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
}


\details{
Warning: Unless well-justified constraints are imposed on either \eqn{\lambda} and/or \eqn{\mu} and \eqn{\rho}, it is generally impossible to reliably estimate \eqn{\lambda} and \eqn{\mu} from extant timetrees alone (Louca and Pennell, 2019). This routine (and any other software that claims to estimate \eqn{\lambda} and \eqn{\mu} solely from extant timetrees) should thus be used with great suspicion. If your only source of information is an extant timetree, and you have no a priori information on how \eqn{\lambda} or \eqn{\mu} might have looked like, you should consider using the more appropriate routines \code{\link{fit_hbd_pdr_on_grid}} and \code{\link{fit_hbd_psr_on_grid}} instead.

If \code{age0>0}, the input tree is essentially trimmed at \code{age0} (omitting anything younger than \code{age0}), and the various variables are fitted to this new (shorter) tree, with time shifted appropriately. For example, the fitted \code{rho0} is thus the sampling fraction at \code{age0}, i.e. the fraction of lineages extant at \code{age0} that are represented in the timetree.

It is generally advised to provide as much information to the function \code{fit_hbd_model_on_grid} as possible, including reasonable lower and upper bounds (\code{min_lambda}, \code{max_lambda}, \code{min_mu}, \code{max_mu}, \code{min_rho0} and \code{max_rho0}) and a reasonable parameter guess (\code{guess_lambda}, \code{guess_mu} and \code{guess_rho0}). It is also important that the \code{age_grid} is sufficiently fine to capture the expected major variations of \eqn{\lambda} and \eqn{\mu} over time, but keep in mind the serious risk of overfitting when \code{age_grid} is too fine and/or the tree is too small.
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
\item{fitted_lambda}{
Numeric vector of size Ngrid, listing fitted or fixed speciation rates \eqn{\lambda} at each age-grid point. Between grid points \eqn{\lambda} should be interpreted as a piecewise polynomial function (natural spline) of degree \code{splines_degree}; to evaluate this function at arbitrary ages use the \code{castor} routine \code{\link{evaluate_spline}}.
}
\item{fitted_mu}{
Numeric vector of size Ngrid, listing fitted or fixed extinction rates \eqn{\mu} at each age-grid point. Between grid points \eqn{\mu} should be interpreted as a piecewise polynomial function (natural spline) of degree \code{splines_degree}; to evaluate this function at arbitrary ages use the \code{castor} routine \code{\link{evaluate_spline}}.
}
\item{fitted_rho}{
Numeric, specifying the fitted or fixed sampling fraction \eqn{\rho}.
}
\item{guess_lambda}{
Numeric vector of size Ngrid, specifying the initial guess for \eqn{\lambda} at each age-grid point.
}
\item{guess_mu}{
Numeric vector of size Ngrid, specifying the initial guess for \eqn{\mu} at each age-grid point.
}
\item{guess_rho0}{
Numeric, specifying the initial guess for \eqn{\rho}.
}
\item{age_grid}{
The age-grid on which \eqn{\lambda} and \eqn{\mu} are defined. This will be the same as the provided \code{age_grid}, unless the latter was \code{NULL} or of length <=1.
}
\item{NFP}{
Integer, number of free (i.e., independently) fitted parameters. If none of the \eqn{\lambda}, \eqn{\mu} and \eqn{\rho} were fixed, and \code{const_lambda=FALSE} and \code{const_mu=FALSE}, then \code{NFP} will be equal to 2*Ngrid+1.
}
\item{AIC}{
The Akaike Information Criterion for the fitted model, defined as \eqn{2k-2\log(L)}, where \eqn{k} is the number of fitted parameters and \eqn{L} is the maximized likelihood.
}
\item{BIC}{
The Bayesian information criterion for the fitted model, defined as \eqn{\log(n)k-2\log(L)}, where \eqn{k} is the number of fitted parameters, \eqn{n} is the number of data points (number of branching times), and \eqn{L} is the maximized likelihood.
}
\item{condition}{
Character, specifying what conditioning was root for the likelihood (e.g. "crown" or "stem").
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
}


\author{Stilianos Louca}

\references{
T. Stadler (2009). On incomplete sampling under birth-death models and connections to the sampling-based coalescent. Journal of Theoretical Biology. 261:58-66.

T. Stadler (2013). How can we improve accuracy of macroevolutionary rate estimates? Systematic Biology. 62:321-329.

H. Morlon, T. L. Parsons, J. B. Plotkin (2011). Reconciling molecular phylogenies with the fossil record. Proceedings of the National Academy of Sciences. 108:16327-16332.

S. Louca et al. (2018). Bacterial diversification through geological time. Nature Ecology & Evolution. 2:1458-1467.

S. Louca and M. W. Pennell (2020). Extant timetrees are consistent with a myriad of diversification histories. Nature. 580:502-505.
}

\seealso{
\code{\link{simulate_deterministic_hbd}}

\code{\link{loglikelihood_hbd}}

\code{\link{fit_hbd_model_parametric}}

\code{\link{fit_hbd_pdr_on_grid}}

\code{\link{fit_hbd_pdr_parametric}}

\code{\link{fit_hbd_psr_on_grid}}
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


# Fit mu on grid
# Assume that lambda & rho are known
Ngrid     = 5
age_grid  = seq(from=0,to=root_age,length.out=Ngrid)
fit = fit_hbd_model_on_grid(tree, 
          age_grid    = age_grid,
          max_mu      = 100,
          fixed_lambda= approx(x=time_grid,y=lambdas,xout=sim$final_time-age_grid)$y,
          fixed_rho0   = rho,
          condition   = "crown",
          Ntrials     = 10,	# perform 10 fitting trials
          Nthreads    = 2,	# use two CPUs
          max_model_runtime = 1) 	# limit model evaluation to 1 second
if(!fit$success){
  cat(sprintf("ERROR: Fitting failed: \%s\n",fit$error))
}else{
  cat(sprintf("Fitting succeeded:\nLoglikelihood=\%g\n",fit$loglikelihood))
  
  # plot fitted & true mu
  plot( x     = fit$age_grid,
        y     = fit$fitted_mu,
        main  = 'Fitted & true mu',
        xlab  = 'age',
        ylab  = 'mu',
        type  = 'b',
        col   = 'red',
        xlim  = c(root_age,0))
  lines(x     = sim$final_time-time_grid,
        y     = mus,
        type  = 'l',
        col   = 'blue');
        
  # get fitted mu as a function of age
  mu_fun = approxfun(x=fit$age_grid, y=fit$fitted_mu)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{birth-death model}
\keyword{skyline model}
