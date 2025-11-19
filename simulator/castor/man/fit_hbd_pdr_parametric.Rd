\name{fit_hbd_pdr_parametric}
\alias{fit_hbd_pdr_parametric}
\title{
Fit parameterized pulled diversification rates of birth-death models.
}
\description{
Given an ultrametric timetree, estimate the pulled diversification rate (PDR) of homogenous birth-death (HBD) models that best explains the tree via maximum likelihood, assuming that the PDR is given as a parameterized function of time before present. Every HBD model is defined by some speciation and extinction rates (\eqn{\lambda} and \eqn{\mu}) over time, as well as the sampling fraction \eqn{\rho} (fraction of extant species sampled). ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction rates. For any given HBD model there exists an infinite number of alternative HBD models that generate extant trees with the same probability distributions and yield the same likelihood for any given extant timetree; these ``congruent'' models cannot be distinguished from one another solely based on an extant timetree. 

Each congruence class is uniquely described by its PDR, defined as \eqn{PDR=\lambda-\mu+\lambda^{-1}d\lambda/d\tau} (where \eqn{\tau} is time before present) as well as the product \eqn{\rho\lambda_o} (where \eqn{\lambda_o} is the present-day speciation rate). That is, two HBD models are congruent if and only if they have the same PDR and the same product \eqn{\rho\lambda_o}. This function is designed to estimate the generating congruence class for the tree, by fitting a finite number of parameters defining the PDR and \eqn{\rho\lambda_o}.
}
\usage{
fit_hbd_pdr_parametric( tree, 
                        param_values,
                        param_guess        = NULL,
                        param_min          = -Inf,
                        param_max          = +Inf,
                        param_scale        = NULL,
                        oldest_age         = NULL,
                        age0               = 0,
                        PDR,
                        rholambda0,
                        age_grid           = NULL,
                        condition          = "auto",
                        relative_dt        = 1e-3,
                        Ntrials            = 1,
                        max_start_attempts = 1,
                        Nthreads           = 1,
                        max_model_runtime  = NULL,
                        fit_control        = list())
}
\arguments{
\item{tree}{
A rooted ultrametric timetree of class "phylo", representing the time-calibrated phylogeny of a set of extant sampled species.
}
\item{param_values}{
Numeric vector, specifying fixed values for a some or all model parameters. For fitted (i.e., non-fixed) parameters, use \code{NaN} or \code{NA}. For example, the vector \code{c(1.5,NA,40)} specifies that the 1st and 3rd model parameters are fixed at the values 1.5 and 40, respectively, while the 2nd parameter is to be fitted. The length of this vector defines the total number of model parameters. If entries in this vector are named, the names are taken as parameter names. Names should be included if you'd like returned parameter vectors to have named entries, or if the functions \code{PDR} or \code{rho} query parameter values by name (as opposed to numeric index).
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
Non-negative numeric, specifying the youngest age (time before present) to consider for fitting, and with respect to which \code{rholambda0} is defined. If \code{age0>0}, then \code{rholambda0} refers to the product of the sampling fraction at age \code{age0} and the speciation rate at age \code{age0}. See below for more details.
}
\item{PDR}{
Function specifying the pulled diversification rate at any given age (time before present) and for any given parameter values. This function must take exactly two arguments, the 1st one being a numeric vector (one or more ages) and the 2nd one being a numeric vector of size NP (parameter values), and return a numeric vector of the same size as the 1st argument. Can also be a single number (i.e., PDR is fixed).
}
\item{rholambda0}{
Function specifying the product \eqn{\rho\lambda_o} (sampling fraction times speciation rate at \code{age0}) for any given parameter values. This function must take exactly one argument, a numeric vector of size NP (parameter values), and return a strictly positive numeric. Can also be a single number (i.e., rholambda0 is fixed).
}
\item{age_grid}{
Numeric vector, specifying ages at which the \code{PDR} function should be evaluated. This age grid must be fine enough to capture the possible variation in the PDR over time, within the permissible parameter range. If of size 1, then the PDR is assumed to be time-independent.
Listed ages must be strictly increasing, and must cover at least the full considered age interval (from \code{age0} to \code{oldest_age}). Can also be \code{NULL} or a vector of size 1, in which case the PDR is assumed to be time-independent.
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
\item{max_start_attempts}{
Integer, specifying the number of times to attempt finding a valid start point (per trial) before giving up on that trial. Randomly choosen extreme start parameters may occasionally result in Inf/undefined likelihoods, so this option allows the algorithm to keep looking for valid starting points.
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
This function is designed to estimate a finite set of scalar parameters (\eqn{p_1,..,p_n\in\R}) that determine the PDR and the product \eqn{\rho\lambda_o} (sampling fraction times present-dat extinction rate), by maximizing the likelihood of observing a given timetree under the HBD model. For example, the investigator may assume that the PDR varies exponentially over time, i.e. can be described by \eqn{PDR(t)=A\cdot e^{-B t}} (where \eqn{A} and \eqn{B} are unknown coefficients and \eqn{t} is time before present), and that the product \eqn{\rho\lambda_o} is unknown. In this case the model has 3 free parameters, \eqn{p_1=A}, \eqn{p_2=B} and \eqn{p_3=\rho\lambda_o}, each of which may be fitted to the tree.

If \code{age0>0}, the input tree is essentially trimmed at \code{age0} (omitting anything younger than \code{age0}), and the PDR and \code{rholambda0} are fitted to this new (shorter) tree, with time shifted appropriately. The fitted \code{rholambda0} is thus the product of the sampling fraction at \code{age0} and the speciation rate at \code{age0}. Note that the sampling fraction at \code{age0} is simply the fraction of lineages extant at \code{age0} that are represented in the timetree. Most users will typically want to leave \code{age0=0}.

It is generally advised to provide as much information to the function \code{fit_hbd_pdr_parametric} as possible, including reasonable lower and upper bounds (\code{param_min} and \code{param_max}), a reasonable parameter guess (\code{param_guess}) and reasonable parameter scales \code{param_scale}. If some model parameters can vary over multiple orders of magnitude, it is advised to transform them so that they vary across fewer orders of magnitude (e.g., via log-transformation). It is also important that the \code{age_grid} is sufficiently fine to capture the variation of the PDR over time, since the likelihood is calculated under the assumption that both vary linearly between grid points.
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
Numeric vector of size NP, listing guessed or fixed values for all model parameters in their standard order.
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
}


\author{Stilianos Louca}

\references{
H. Morlon, T. L. Parsons, J. B. Plotkin (2011). Reconciling molecular phylogenies with the fossil record. Proceedings of the National Academy of Sciences. 108:16327-16332.

S. Louca et al. (2018). Bacterial diversification through geological time. Nature Ecology & Evolution. 2:1458-1467.

S. Louca and M. W. Pennell (2020). Extant timetrees are consistent with a myriad of diversification histories. Nature. 580:502-505.
}

\seealso{
\code{\link{simulate_deterministic_hbd}}

\code{\link{loglikelihood_hbd}}

\code{\link{fit_hbd_model_on_grid}}

\code{\link{fit_hbd_model_parametric}}

\code{\link{fit_hbd_pdr_on_grid}}

\code{\link{fit_hbd_psr_parametric}}

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
tree      = generate_random_tree( parameters  = list(rarefaction=rho), 
                                  max_tips    = Ntips/rho,
                                  coalescent  = TRUE,
                                  added_rates_times     = time_grid,
                                  added_birth_rates_pc  = lambdas,
                                  added_death_rates_pc  = mus)$tree
root_age = castor::get_tree_span(tree)$max_distance
cat(sprintf("Tree has \%d tips, spans \%g Myr\n",length(tree$tip.label),root_age))

# Define a parametric HBD congruence class, with exponentially varying PDR
# The model thus has 3 parameters
PDR_function = function(ages,params){
	return(params['A']*exp(-params['B']*ages));
}
rholambda0_function = function(params){
	return(params['rholambda0'])
}

# Define an age grid on which PDR_function shall be evaluated
# Should be sufficiently fine to capture the variation in the PDR
age_grid = seq(from=0,to=100,by=0.01)

# Perform fitting
cat(sprintf("Fitting class to tree..\n"))
fit = fit_hbd_pdr_parametric(	tree, 
                      param_values  = c(A=NA, B=NA, rholambda0=NA),
                      param_guess   = c(1,0,1),
                      param_min     = c(-10,-10,0),
                      param_max     = c(10,10,10),
                      param_scale   = 1, # all params are in the order of 1
                      PDR           = PDR_function,
                      rholambda0    = rholambda0_function,
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
