\name{fit_hbd_model_parametric}
\alias{fit_hbd_model_parametric}
\title{
Fit a parametric homogenous birth-death model to a timetree.
}
\description{
Given an ultrametric timetree, fit a homogenous birth-death (HBD) model in which speciation and extinction rates (\eqn{\lambda} and \eqn{\mu}) are given as parameterized functions of time before present.  ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction rates (in the literature this is sometimes referred to simply as ``birth-death model''). Every HBD model is defined based on the values that \eqn{\lambda} and \eqn{\mu} take over time as well as the sampling fraction \eqn{\rho} (fraction of extant species sampled); in turn, \eqn{\lambda}, \eqn{\mu} and \eqn{\rho} can be parameterized by a finite set of parameters. This function estimates these parameters by maximizing the likelihood (Morlon et al. 2011) of the timetree under the resulting HBD model.
}
\usage{
fit_hbd_model_parametric( tree, 
                          param_values,
                          param_guess           = NULL,
                          param_min             = -Inf,
                          param_max             = +Inf,
                          param_scale           = NULL,
                          oldest_age            = NULL,
                          age0                  = 0,
                          lambda,
                          mu                    = 0,
                          rho0                  = 1,
                          age_grid              = NULL,
                          condition             = "auto",
                          relative_dt           = 1e-3,
                          Ntrials               = 1,
                          max_start_attempts    = 1,
                          Nthreads              = 1,
                          max_model_runtime     = NULL,
                          Nbootstraps           = 0,
                          Ntrials_per_bootstrap = NULL,
                          fit_algorithm         = "nlminb",
                          fit_control           = list(),
                          focal_param_values    = NULL,
                          verbose               = FALSE,
                          diagnostics           = FALSE,
                          verbose_prefix        = "")
}
\arguments{
\item{tree}{
A rooted ultrametric timetree of class "phylo", representing the time-calibrated reconstructed phylogeny of a set of extant sampled species.
}
\item{param_values}{
Numeric vector, specifying fixed values for a some or all model parameters. For fitted (i.e., non-fixed) parameters, use \code{NaN} or \code{NA}. For example, the vector \code{c(1.5,NA,40)} specifies that the 1st and 3rd model parameters are fixed at the values 1.5 and 40, respectively, while the 2nd parameter is to be fitted. The length of this vector defines the total number of model parameters. If entries in this vector are named, the names are taken as parameter names. Names should be included if you'd like returned parameter vectors to have named entries, or if the functions \code{lambda}, \code{mu} or \code{rho} query parameter values by name (as opposed to numeric index).
}
\item{param_guess}{
Numeric vector of size NP, specifying a first guess for the value of each model parameter. For fixed parameters, guess values are ignored. Can be \code{NULL} only if all model parameters are fixed.
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
\item{oldest_age}{
Strictly positive numeric, specifying the oldest time before present (``age'') to consider when calculating the likelihood. If this is equal to or greater than the root age, then \code{oldest_age} is taken as the stem age, and the classical formula by Morlon et al. (2011) is used.  If \code{oldest_age} is less than the root age, the tree is split into multiple subtrees at that age by treating every edge crossing that age as the stem of a subtree, and each subtree is considered an independent realization of the HBD model stemming at that age. This can be useful for avoiding points in the tree close to the root, where estimation uncertainty is generally higher. If \code{oldest_age==NULL}, it is automatically set to the root age.
}
\item{age0}{
Non-negative numeric, specifying the youngest age (time before present) to consider for fitting, and with respect to which \code{rho} is defined. If \code{age0>0}, then \code{rho0} refers to the sampling fraction at age \code{age0}, i.e. the fraction of lineages extant at \code{age0} that are included in the tree. See below for more details.
}
\item{lambda}{
Function specifying the speciation rate at any given age (time before present) and for any given parameter values. This function must take exactly two arguments, the 1st one being a numeric vector (one or more ages) and the 2nd one being a numeric vector of size NP (parameter values), and return a numeric vector of the same size as the 1st argument with strictly positive entries. Can also be a single number (i.e., lambda is fixed).
}
\item{mu}{
Function specifying the extinction rate at any given age and for any given parameter values. This function must take exactly two arguments, the 1st one being a numeric vector (one or more ages) and the 2nd one being a numeric vector of size NP (parameter values), and return a numeric vector of the same size as the 1st argument with non-negative entries. Can also be a single number (i.e., mu is fixed).
}
\item{rho0}{
Function specifying the sampling fraction (fraction of extant species sampled at \code{age0}) for any given parameter values. This function must take exactly one argument, a numeric vector of size NP (parameter values), and return a numeric between 0 (exclusive) and 1 (inclusive). Can also be a single number (i.e., rho0 is fixed).
}
\item{age_grid}{
Numeric vector, specifying ages at which the \code{lambda} and \code{mu} functionals should be evaluated. This age grid must be fine enough to capture the possible variation in \eqn{\lambda} and \eqn{\mu} over time, within the permissible parameter range. If of size 1, then lambda & mu are assumed to be time-independent.
Listed ages must be strictly increasing, and must cover at least the full considered age interval (from 0 to \code{oldest_age}). Can also be \code{NULL} or a vector of size 1, in which case the speciation rate and extinction rate is assumed to be time-independent.
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
\item{max_start_attempts}{
Integer, specifying the number of times to attempt finding a valid start point (per trial) before giving up on that trial. Randomly choosen extreme start parameters may occasionally result in Inf/undefined likelihoods, so this option allows the algorithm to keep looking for valid starting points.
}
\item{Nthreads}{
Integer, specifying the number of parallel threads to use for performing multiple fitting trials simultaneously. This should generally not exceed the number of available CPUs on your machine. Parallel computing is not available on the Windows platform.
}
\item{max_model_runtime}{
Optional numeric, specifying the maximum number of seconds to allow for each evaluation of the likelihood function. Use this to abort fitting trials leading to parameter regions where the likelihood takes a long time to evaluate (these are often unlikely parameter regions).
}
\item{Nbootstraps}{
Integer, specifying the number of parametric bootstraps to perform for estimating standard errors and confidence intervals of estimated model parameters. Set to 0 for no bootstrapping.
}
\item{Ntrials_per_bootstrap}{
Integer, specifying the number of fitting trials to perform for each bootstrap sampling. If \code{NULL}, this is set equal to \code{max(1,Ntrials)}. Decreasing \code{Ntrials_per_bootstrap} will reduce computation time, at the expense of potentially inflating the estimated confidence intervals; in some cases (e.g., for very large trees) this may be useful if fitting takes a long time and confidence intervals are very narrow anyway. Only relevant if \code{Nbootstraps>0}.
}
\item{fit_algorithm}{
Character, specifying which optimization algorithm to use. Either "nlminb" or "subplex" are allowed.
}
\item{fit_control}{
Named list containing options for the \code{nlminb} or \code{subplex} optimization routine, depending on the choice of \code{fit_algorithm}. For example, for "nlminb" commonly modified options are \code{iter.max}, \code{eval.max} or \code{rel.tol}. For a complete list of options and default values see the documentation of \code{nlminb} in the \code{stats} package or of \code{nloptr} in the \code{nloptr} package.
}
\item{focal_param_values}{
Optional numeric matrix having NP columns and an arbitrary number of rows, listing combinations of parameter values of particular interest and for which the log-likelihoods should be returned. This may be used for diagnostic purposes, e.g., to examine the shape of the likelihood function.
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
This function is designed to estimate a finite set of scalar parameters (\eqn{p_1,..,p_n\in\R}) that determine the speciation rate \eqn{\lambda}, the extinction rate \eqn{\mu} and the sampling fraction \eqn{\rho}, by maximizing the likelihood of observing a given timetree under the HBD model. For example, the investigator may assume that both \eqn{\lambda} and \eqn{\mu} vary exponentially over time, i.e. they can be described by \eqn{\lambda(t)=\lambda_o\cdot e^{-\alpha t}} and \eqn{\mu(t)=\mu_o\cdot e^{-\beta t}} (where \eqn{\lambda_o}, \eqn{\mu_o} are unknown present-day rates and \eqn{\alpha}, \eqn{\beta} are unknown factors, and \eqn{t} is time before present), and that the sampling fraction \eqn{\rho} is known. In this case the model has 4 free parameters, \eqn{p_1=\lambda_o}, \eqn{p_2=\mu_o}, \eqn{p_3=\alpha} and \eqn{p_4=\beta}, each of which may be fitted to the tree.

It is generally advised to provide as much information to the function \code{fit_hbd_model_parametric} as possible, including reasonable lower and upper bounds (\code{param_min} and \code{param_max}), a reasonable parameter guess (\code{param_guess}) and reasonable parameter scales \code{param_scale}. If some model parameters can vary over multiple orders of magnitude, it is advised to transform them so that they vary across fewer orders of magnitude (e.g., via log-transformation). It is also important that the \code{age_grid} is sufficiently fine to capture the variation of \code{lambda} and \code{mu} over time, since the likelihood is calculated under the assumption that both vary linearly between grid points.
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
\item{param_fitted}{
Numeric vector of size NP (number of model parameters), listing all fitted or fixed model parameters in their standard order (see details above). If \code{param_names} was provided, elements in \code{fitted_params} will be named.
}
\item{param_guess}{
Numeric vector of size NP, listing guessed or fixed values for all model parameters in their standard order. If \code{param_names} was provided, elements in \code{param_guess} will be named.
}
\item{loglikelihood}{
The log-likelihood of the fitted model for the given timetree.
}
\item{NFP}{
Integer, number of fitted (i.e., non-fixed) model parameters.
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
Numeric between 0 and 1, estimated consistency of the data with the fitted model (Lindholm et al. 2019). See the documentation of \code{\link{fit_hbds_model_on_grid}} for an explanation.
}
}


\author{Stilianos Louca}

\references{
H. Morlon, T. L. Parsons, J. B. Plotkin (2011). Reconciling molecular phylogenies with the fossil record. Proceedings of the National Academy of Sciences. 108:16327-16332.

S. Louca et al. (2018). Bacterial diversification through geological time. Nature Ecology & Evolution. 2:1458-1467.

A. Lindholm, D. Zachariah, P. Stoica, T. B. Schoen (2019). Data consistency approach to model validation. IEEE Access. 7:59788-59796.

S. Louca and M. W. Pennell (2020). Extant timetrees are consistent with a myriad of diversification histories. Nature. 580:502-505.
}

\seealso{
\code{\link{simulate_deterministic_hbd}}

\code{\link{loglikelihood_hbd}}

\code{\link{fit_hbd_model_on_grid}}

\code{\link{fit_hbd_pdr_on_grid}}

\code{\link{fit_hbd_pdr_parametric}}
}

\examples{
\dontrun{
# Generate a random tree with exponentially varying lambda & mu
Ntips     = 10000
rho       = 0.5 # sampling fraction
time_grid = seq(from=0, to=100, by=0.01)
lambdas   = 2*exp(0.1*time_grid)
mus       = 1.5*exp(0.09*time_grid)
tree      = generate_random_tree( parameters  = list(rarefaction=rho), 
                                  max_tips    = Ntips/rho,
                                  coalescent  = TRUE,
                                  added_rates_times     = time_grid,
                                  added_birth_rates_pc  = lambdas,
                                  added_death_rates_pc  = mus)$tree
root_age = castor::get_tree_span(tree)$max_distance
cat(sprintf("Tree has \%d tips, spans \%g Myr\n",length(tree$tip.label),root_age))

# Define a parametric HBD model, with exponentially varying lambda & mu
# Assume that the sampling fraction is known
# The model thus has 4 parameters: lambda0, mu0, alpha, beta
lambda_function = function(ages,params){
	return(params['lambda0']*exp(-params['alpha']*ages));
}
mu_function = function(ages,params){
	return(params['mu0']*exp(-params['beta']*ages));
}
rho_function = function(params){
	return(rho) # rho does not depend on any of the parameters
}

# Define an age grid on which lambda_function & mu_function shall be evaluated
# Should be sufficiently fine to capture the variation in lambda & mu
age_grid = seq(from=0,to=100,by=0.01)

# Perform fitting
# Lets suppose extinction rates are already known
cat(sprintf("Fitting model to tree..\n"))
fit = fit_hbd_model_parametric(	tree, 
                      param_values  = c(lambda0=NA, mu0=3, alpha=NA, beta=-0.09),
                      param_guess   = c(1,1,0,0),
                      param_min     = c(0,0,-1,-1),
                      param_max     = c(10,10,1,1),
                      param_scale   = 1, # all params are in the order of 1
                      lambda        = lambda_function,
                      mu            = mu_function,
                      rho0          = rho_function,
                      age_grid      = age_grid,
                      Ntrials       = 10,    # perform 10 fitting trials
                      Nthreads      = 2,     # use 2 CPUs
                      max_model_runtime = 1, # limit model evaluation to 1 second
                      fit_control       = list(rel.tol=1e-6))
if(!fit$success){
	cat(sprintf("ERROR: Fitting failed: \%s\n",fit$error))
}else{
	cat(sprintf("Fitting succeeded:\nLoglikelihood=\%g\n",fit$loglikelihood))
	print(fit)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{birth-death model}
