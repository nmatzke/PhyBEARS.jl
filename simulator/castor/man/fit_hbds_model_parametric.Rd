\name{fit_hbds_model_parametric}
\alias{fit_hbds_model_parametric}
\title{
Fit a parametric homogenous birth-death-sampling model to a timetree.
}
\description{
Given a timetree (potentially sampled through time and not necessarily ultrametric), fit a homogenous birth-death-sampling (HBDS) model in which speciation, extinction and lineage sampling occurs at some continuous (Poissonian) rates \eqn{\lambda}, \eqn{\mu} and \eqn{\psi}, which are given as parameterized functions of time before present. Sampled lineages are kept in the pool of extant lineages at some ``retention probability'' \eqn{\kappa}, which may also depend on time. In addition, this model can include concentrated sampling attempts (CSAs) at a finite set of discrete time points \eqn{t_1,..,t_m}. ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction/sampling rates. Every HBDS model is thus defined based on the values that \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} take over time, as well as the sampling probabilities \eqn{\rho_1,..,\rho_m} and retention probabilities \eqn{\kappa_1,..,\kappa_m} during the concentrated sampling attempts; each of these parameters, in turn, is assumed to be determined by a finite set of parameters. This function estimates these parameters by maximizing the corresponding likelihood of the timetree. Special cases of this model are sometimes known as ``birth-death-skyline plots'' in the literature (Stadler 2013). In epidemiology, these models are often used to describe the phylogenies of viral strains sampled over the course of the epidemic.
}
\usage{
fit_hbds_model_parametric(tree, 
                          param_values,
                          param_guess           = NULL,
                          param_min             = -Inf,
                          param_max             = +Inf,
                          param_scale           = NULL,
                          root_age              = NULL,
                          oldest_age            = NULL,
                          lambda                = 0,
                          mu                    = 0,
                          psi                   = 0,
                          kappa                 = 0,
                          age_grid              = NULL,
                          CSA_ages              = NULL,
                          CSA_probs             = NULL,
                          CSA_kappas            = 0,
                          condition             = "auto",
                          ODE_relative_dt       = 0.001,
                          ODE_relative_dy       = 1e-3,
                          CSA_age_epsilon       = NULL,
                          Ntrials               = 1,
                          max_start_attempts    = 1,
                          Nthreads              = 1,
                          max_model_runtime     = NULL,
                          Nbootstraps           = 0,
                          Ntrials_per_bootstrap = NULL,
                          fit_control           = list(),
                          focal_param_values    = NULL,
                          verbose               = FALSE,
                          diagnostics           = FALSE,
                          verbose_prefix        = "")
}
\arguments{
\item{tree}{
A timetree of class "phylo", representing the time-calibrated reconstructed phylogeny of a set of extant and/or extinct species. Tips of the tree are interpreted as terminally sampled lineages, while monofurcating nodes are interpreted as non-terminally sampled lineages, i.e., lineages sampled at some past time point and with subsequently sampled descendants.
}
\item{param_values}{
Numeric vector, specifying fixed values for a some or all model parameters. For fitted (i.e., non-fixed) parameters, use \code{NaN} or \code{NA}. For example, the vector \code{c(1.5,NA,40)} specifies that the 1st and 3rd model parameters are fixed at the values 1.5 and 40, respectively, while the 2nd parameter is to be fitted. The length of this vector defines the total number of model parameters. If entries in this vector are named, the names are taken as parameter names. Names should be included if the functions \code{lambda}, \code{mu}, \code{psi}, \code{kappa}, \code{CSA_psi} and \code{CSA_kappa} query parameter values by name (as opposed to numeric index).
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
\item{root_age}{
Positive numeric, specifying the age of the tree's root. Can be used to define a time offset, e.g. if the last tip was not actually sampled at the present. If \code{NULL}, this will be calculated from the tree and it will be assumed that the last tip was sampled at the present.
}
\item{oldest_age}{
Strictly positive numeric, specifying the oldest time before present (``age'') to consider when calculating the likelihood. If this is equal to or greater than the root age, then \code{oldest_age} is interpreted as the stem age.  If \code{oldest_age} is less than the root age, the tree is split into multiple subtrees at that age by treating every edge crossing that age as the stem of a subtree, and each subtree is considered an independent realization of the HBDS model stemming at that age. This can be useful for avoiding points in the tree close to the root, where estimation uncertainty is generally higher. If \code{oldest_age==NULL}, it is automatically set to the root age.
}
\item{lambda}{
Function specifying the speciation rate at any given age (time before present) and for any given parameter values. This function must take exactly two arguments, the 1st one being a numeric vector (one or more ages) and the 2nd one being a numeric vector of size NP (parameter values), and return a numeric vector of the same size as the 1st argument with strictly positive entries. Can also be a single numeric (i.e., lambda is fixed).
}
\item{mu}{
Function specifying the extinction rate at any given age and for any given parameter values. This function must take exactly two arguments, the 1st one being a numeric vector (one or more ages) and the 2nd one being a numeric vector of size NP (parameter values), and return a numeric vector of the same size as the 1st argument with non-negative entries. Can also be a single numeric (i.e., mu is fixed).
}
\item{psi}{
Function specifying the continuous (Poissonian) lineage sampling rate at any given age and for any given parameter values. This function must take exactly two arguments, the 1st one being a numeric vector (one or more ages) and the 2nd one being a numeric vector of size NP (parameter values), and return a numeric vector of the same size as the 1st argument with non-negative entries. Can also be a single numeric (i.e., psi is fixed).
}
\item{kappa}{
Function specifying the retention probability for continuously sampled lineages, at any given age and for any given parameter values. This function must take exactly two arguments, the 1st one being a numeric vector (one or more ages) and the 2nd one being a numeric vector of size NP (parameter values), and return a numeric vector of the same size as the 1st argument with non-negative entries. The retention probability is the probability of a sampled lineage remaining in the pool of extant lineages. Can also be a single numeric (i.e., kappa is fixed).
}
\item{age_grid}{
Numeric vector, specifying ages at which the \code{lambda}, \code{mu}, \code{psi} and \code{kappa} functionals should be evaluated. This age grid must be fine enough to capture the possible variation in \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} over time, within the permissible parameter range. Listed ages must be strictly increasing, and must cover at least the full considered age interval (from 0 to \code{oldest_age}). Can also be \code{NULL} or a vector of size 1, in which case \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} are assumed to be time-independent.
}
\item{CSA_ages}{
Optional numeric vector, listing ages (in ascending order) at which concentrated sampling attempts occurred. If \code{NULL}, it is assumed that no concentrated sampling attempts took place and that all tips were sampled according to the continuous sampling rate \code{psi}.
}
\item{CSA_probs}{
Function specifying the sampling probabilities during the various concentrated sampling attempts, depending on parameter values. Hence, for any choice of parameters, \code{CSA_probs} must return a numeric vector of the same size as \code{CSA_ages}. Can also be a single numeric (i.e., concentrated sampling probability is fixed).
}
\item{CSA_kappas}{
Function specifying the retention probabilities during the various concentrated sampling attempts, depending on parameter values. Hence, for any choice of parameters, \code{CSA_kappas} must return a numeric vector of the same size as \code{CSA_ages}. Can also be a single numeric (i.e., retention probability during concentrated samplings is fixed).
}
\item{condition}{
Character, either "crown", "stem", "none" or "auto", specifying on what to condition the likelihood. If "crown", the likelihood is conditioned on the survival of the two daughter lineages branching off at the root. If "stem", the likelihood is conditioned on the survival of the stem lineage. Note that "crown" really only makes sense when \code{oldest_age} is equal to the root age, while "stem" is recommended if \code{oldest_age} differs from the root age. "none" is usually not recommended. If "auto", the condition is chosen according to the above recommendations.
}
\item{ODE_relative_dt}{
Positive unitless number, specifying the default relative time step for the ordinary differential equation solvers. Typical values are 0.01-0.001.
}
\item{ODE_relative_dy}{
Positive unitless number, specifying the relative difference between subsequent simulated and interpolated values, in internally used ODE solvers. Typical values are 1e-2 to 1e-5. A smaller \code{ODE_relative_dy} increases interpolation accuracy, but also increases memory requirements and adds runtime (scaling with the tree's age span, not with Ntips).
}
\item{CSA_age_epsilon}{
Non-negative numeric, in units of time, specfying the age radius around a concentrated sampling attempt, within which to assume that sampling events were due to that concentrated sampling attempt. If \code{NULL}, this is chosen automatically based on the anticipated scale of numerical rounding errors. Only relevant if concentrated sampling attempts are included.
}
\item{Ntrials}{
Integer, specifying the number of independent fitting trials to perform, each starting from a random choice of model parameters. Increasing \code{Ntrials} reduces the risk of reaching a non-global local maximum in the fitting objective.
}
\item{max_start_attempts}{
Integer, specifying the number of times to attempt finding a valid start point (per trial) before giving up on that trial. Randomly chosen extreme start parameters may occasionally result in Inf/undefined likelihoods, so this option allows the algorithm to keep looking for valid starting points.
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
\item{fit_control}{
 Named list containing options for the \code{nlminb} optimization routine, such as \code{iter.max}, \code{eval.max} or \code{rel.tol}. For a complete list of options and default values see the documentation of \code{nlminb} in the \code{stats} package.
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
This function is designed to estimate a finite set of scalar parameters (\eqn{p_1,..,p_n\in\R}) that determine the speciation rate \eqn{\lambda}, the extinction rate \eqn{\mu}, the sampling rate \eqn{\psi}, the retention rate \eqn{\kappa}, the concentrated sampling probabilities \eqn{\rho_1,..,\rho_m} and the concentrated retention probabilities \eqn{\kappa_1,..,\kappa_m}, by maximizing the likelihood of observing a given timetree under the HBDS model. Note that the ages (times before present) of the concentrated sampling attempts are assumed to be known and are not fitted.

It is generally advised to provide as much information to the function \code{fit_hbds_model_parametric} as possible, including reasonable lower and upper bounds (\code{param_min} and \code{param_max}), a reasonable parameter guess (\code{param_guess}) and reasonable parameter scales \code{param_scale}. If some model parameters can vary over multiple orders of magnitude, it is advised to transform them so that they vary across fewer orders of magnitude (e.g., via log-transformation). It is also important that the \code{age_grid} is sufficiently fine to capture the variation of \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} over time, since the likelihood is calculated under the assumption that these functions vary linearly between grid points.

Note that in this function age always refers to time before present, i.e., present day age is 0 and age increases from tips to root. The functions \code{lambda}, \code{mu}, \code{psi} and \code{kappa} should be functions of age, not forward time. Similarly, concentrated sampling attempts (CSAs) are enumerated in order of increasing age, i.e., starting with the youngest CSA and moving towards older CSAs.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether model fitting succeeded. If \code{FALSE}, the returned list will include an additional ``\code{error}'' element (character) providing a description of the error; in that case all other return variables may be undefined.
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
\item{param_fitted}{
Numeric vector of size NP (number of model parameters), listing all fitted or fixed model parameters in their standard order (see details above). If \code{param_names} was provided, elements in \code{fitted_params} will be named.
}
\item{param_guess}{
Numeric vector of size NP, listing guessed or fixed values for all model parameters in their standard order. If \code{param_names} was provided, elements in \code{param_guess} will be named.
}
\item{guess_loglikelihood}{
The loglikelihood of the data for the initial parameter guess (\code{param_guess}).
}
\item{focal_loglikelihoods}{
A numeric vector of the same size as \code{nrow(focal_param_values)}, listing loglikelihoods for each of the focal parameter conbinations listed in \code{focal_loglikelihoods}.
}
\item{NFP}{
Integer, number of fitted (i.e., non-fixed) model parameters.
}
\item{Ndata}{
Number of data points used for fitting, i.e., the number of sampling and branching events that occurred between ages 0 and \code{oldest_age}.
}
\item{AIC}{
The Akaike Information Criterion for the fitted model, defined as \eqn{2k-2\log(L)}, where \eqn{k} is the number of fitted parameters and \eqn{L} is the maximized likelihood.
}
\item{BIC}{
The Bayesian information criterion for the fitted model, defined as \eqn{\log(n)k-2\log(L)}, where \eqn{k} is the number of fitted parameters, \eqn{n} is the number of data points (\code{Ndata}), and \eqn{L} is the maximized likelihood.
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
T. Stadler, D. Kuehnert, S. Bonhoeffer, A. J. Drummond (2013). Birth-death skyline plot reveals temporal changes of epidemic spread in HIV and hepatitis C virus (HCV). PNAS. 110:228-233.

A. Lindholm, D. Zachariah, P. Stoica, T. B. Schoen (2019). Data consistency approach to model validation. IEEE Access. 7:59788-59796.
}

\seealso{
\code{\link{generate_tree_hbds}},
\code{\link{simulate_deterministic_hbds}}
}

\examples{
\dontrun{
# Generate a random tree with exponentially varying lambda & mu and constant psi
# assume that all sampled lineages are removed from the pool (i.e. kappa=0)
time_grid = seq(from=0, to=100, by=0.01)
root_age  = 5
tree = generate_tree_hbds(max_time  = root_age,
                        time_grid   = time_grid,
                        lambda      = 2*exp(0.1*time_grid),
                        mu          = 0.1*exp(0.09*time_grid),
                        psi         = 0.1,
                        kappa       = 0)$tree
cat(sprintf("Tree has \%d tips\n",length(tree$tip.label)))


# Define a parametric HBDS model, with exponentially varying lambda & mu
# Assume that the sampling rate is constant but unknown
# The model thus has 5 parameters: lambda0, mu0, alpha, beta, psi
lambda_function = function(ages,params){
    return(params['lambda0']*exp(-params['alpha']*ages));
}
mu_function = function(ages,params){
    return(params['mu0']*exp(-params['beta']*ages));
}
psi_function = function(ages,params){
    return(rep(params['psi'],length(ages)))
}

# Define an age grid on which lambda_function & mu_function shall be evaluated
# Should be sufficiently fine to capture the variation in lambda & mu
age_grid = seq(from=0,to=root_age,by=0.01)

# Perform fitting
cat(sprintf("Fitting model to tree..\n"))
fit = fit_hbds_model_parametric(tree, 
                      root_age      = root_age,
                      param_values  = c(lambda0=NA, mu0=NA, alpha=NA, beta=NA, psi=NA),
                      param_guess   = c(1,1,0,0,0.5),
                      param_min     = c(0,0,-1,-1,0),
                      param_max     = c(10,10,1,1,10),
                      param_scale   = 1, # all params are in the order of 1
                      lambda        = lambda_function,
                      mu            = mu_function,
                      psi           = psi_function,
                      kappa         = 0,
                      age_grid      = age_grid,
                      Ntrials       = 4, # perform 4 fitting trials
                      Nthreads      = 2) # use 2 CPUs
if(!fit$success){
    cat(sprintf("ERROR: Fitting failed: \%s\n",fit$error))
}else{
    cat(sprintf("Fitting succeeded:\nLoglikelihood=\%g\n",fit$loglikelihood))
    # print fitted parameters
    print(fit$param_fitted)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{birth-death-sampling model}
