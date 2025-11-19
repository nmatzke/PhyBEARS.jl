\name{fit_tree_model}
\alias{fit_tree_model}
\title{
Fit a cladogenic model to an existing tree.
}
\description{
Fit the parameters of a tree generation model to an existing phylogenetic tree; branch lengths are assumed to be in time units. The fitted model is a stochastic cladogenic process in which speciations (births) and extinctions (deaths) are Poisson processes, as simulated by the function \code{\link{generate_random_tree}}. The birth and death rates of tips can each be constant or power-law functions of the number of extant tips. For example,
\deqn{
B = I + F\cdot N^E,
}
where \eqn{B} is the birth rate, \eqn{I} is the intercept, \eqn{F} is the power-law factor, \eqn{N} is the current number of extant tips and \eqn{E} is the power-law exponent. Each of the parameters I, F, E can be fixed or fitted. 

Fitting can be performed via maximum-likelihood estimation, based on the waiting times between subsequent speciation and/or extinction events represented in the tree. Alternatively, fitting can be performed using least-squares estimation, based on the number of lineages represented in the tree over time ("diversity-vs-time" curve, a.k.a. "lineages-through-time"" curve). Note that the birth and death rates are NOT per-capita rates, they are absolute rates of species appearance and disappearance per time.
}
\usage{
fit_tree_model( tree, 
                parameters          = list(),
                first_guess         = list(),
                min_age             = 0,
                max_age             = 0,
                age_centile         = NULL,
                Ntrials             = 1,
                Nthreads            = 1,
                coalescent          = FALSE,
                discovery_fraction  = NULL,
                fit_control         = list(),
                min_R2              = -Inf,
                min_wR2             = -Inf,
                grid_size           = 100,
                max_model_runtime   = NULL,
                objective           = 'LL')
}

\arguments{
\item{tree}{
A phylogenetic tree, in which branch lengths are assumed to be in time units. The tree may be a coalescent tree (i.e. only include extant clades) or a tree including extinct clades; the tree type influences what type of models can be fitted with each method.
}
\item{parameters}{
A named list specifying fixed and/or unknown birth-death model parameters, with one or more of the following elements:
\itemize{
\item\code{birth_rate_intercept}: Non-negative number. The intercept of the Poissonian rate at which new species (tips) are added. In units 1/time.
\item\code{birth_rate_factor}:
Non-negative number. The power-law factor of the Poissonian rate at which new species (tips) are added. In units 1/time.
\item\code{birth_rate_exponent}:
Numeric. The power-law exponent of the Poissonian rate at which new species (tips) are added. Unitless.
\item\code{death_rate_intercept}:
Non-negative number. The intercept of the Poissonian rate at which extant species (tips) go extinct. In units 1/time.
\item\code{death_rate_factor}:
Non-negative number. The power-law factor of the Poissonian rate at which extant species (tips) go extinct. In units 1/time.
\item\code{death_rate_exponent}:
Numeric. The power-law exponent of the Poissonian rate at which extant species (tips) go extinct. Unitless.
\item\code{resolution}:
Numeric. Resolution at which the tree was collapsed (i.e. every node of age smaller than this resolution replaced by a single tip). In units time. A resolution of 0 means the tree was not collapsed.
\item\code{rarefaction}:
Numeric. Species sampling fraction, i.e. fraction of extant species represented (as tips) in the tree. A rarefaction of 1, for example, implies that the tree is complete, i.e. includes all extant species. Rarefaction is assumed to have occurred after collapsing.
\item\code{extant_diversity}: The current total extant diversity, regardless of the rarefaction and resolution of the tree at hand. For example, if \code{resolution==0} and \code{rarefaction==0.5} and the tree has 1000 tips, then \code{extant_diversity} should be \code{2000}. If \code{resolution} is fixed at 0 and \code{rarefaction} is also fixed, this can be left \code{NULL} and will be inferred automatically by the function.
}
Each of the above elements can also be \code{NULL}, in which case the parameter is fitted. Elements can also be vectors of size 2 (specifying constraint intervals), in which case the parameters are fitted and constrained within the intervals specified. For example, to fit \code{death_rate_factor} while constraining it to the interval [1,2], set its value to \code{c(1,2)}.
}
\item{first_guess}{
A named list (with entries named as in \code{parameters}) specifying starting values for any of the fitted model parameters. Note that if \code{Ntrials>1}, then start values may be randomly modified in all but the first trial. For any parameters missing from \code{first_guess}, initial values are always randomly chosen. \code{first_guess} can also be \code{NULL}.
}
\item{min_age}{
Numeric. Minimum distance from the tree crown, for a node/tip to be considered in the fitting. If <=0 or \code{NULL}, this constraint is ignored. Use this option to omit most recent nodes.
}
\item{max_age}{
Numeric. Maximum distance from the tree crown, for a node/tip to be considered in the fitting. If <=0 or \code{NULL}, this constraint is ignored. Use this option to omit old nodes, e.g. with highly uncertain placements.
}
\item{age_centile}{
Numeric within 0 and 1. Fraction of youngest nodes/tips to consider for the fitting. This can be used as an alternative to max_age. E.g. if set to 0.6, then the 60\% youngest nodes/tips are considered. Either \code{age_centile} or \code{max_age} must be non-NULL, but not both.
}
\item{Ntrials}{
Integer. Number of fitting attempts to perform, each time using randomly varied start values for fitted parameters. The returned fitted parameter values will be taken from the trial with greatest achieved fit objective. A larger number of trials will decrease the chance of hitting a local non-global optimum during fitting.
}
\item{Nthreads}{
Number of threads to use for parallel execution of multiple fitting trials. On Windows, this option has no effect because Windows does not support forks.
}
\item{coalescent}{
Logical, specifying whether the input tree is a coalescent tree (and thus the coalescent version of the model should be fitted). Only available if \code{objective=='R2'}.
}
\item{discovery_fraction}{
Function handle, mapping age to the fraction of discovered lineages in a tree. That is, \code{discovery_fraction(tau)} is the probability that a lineage at age \code{tau}, that has an extant descendant today, will be represented (discovered) in the coalescent tree. In particular, \code{discovery_fraction(0)} equals the fraction of extant lineages represented in the tree. If this is provided, then \code{parameters$rarefaction} is fixed to 1, and \code{discovery_fraction} is applied after simulation. Only relevant if \code{coalescent==TRUE}. Experimental, so leave this \code{NULL} if you don't know what it means.
}
\item{fit_control}{
Named list containing options for the \code{stats::nlminb} optimization routine, such as \code{eval.max} (max number of evaluations), \code{iter.max} (max number of iterations) and \code{rel.tol} (relative tolerance for convergence).
}
\item{min_R2}{
Minimum coefficient of determination of the diversity curve (clade counts vs time) of the model when compared to the input tree, for a fitted model to be accepted. For example, if set to 0.5 then only fit trials achieving an R2 of at least 0.5 will be considered. Set this to \code{-Inf} to not filter fitted models based on the R2.
}
\item{min_wR2}{
Similar to \code{min_R2}, but applying to the weighted R2, where squared-error weights are proportional to the inverse squared diversities.
}
\item{grid_size}{
Integer. Number of equidistant time points to consider when calculating the R2 of a model's diversity-vs-time curve.
}
\item{max_model_runtime}{
Numeric. Maximum runtime (in seconds) allowed for each model evaluation during fitting. Use this to escape from badly parameterized models during fitting (this will likely cause the affected fitting trial to fail). If \code{NULL} or <=0, this option is ignored.
}
\item{objective}{
Character. Objective function to optimize during fitting. Can be either "LL" (log-likelihood of waiting times between speciation events and between extinction events), "R2" (coefficient of determination of diversity-vs-time curve), "wR2" (weighted R2, where weights of squared errors are proportional to the inverse diversities observed in the tree) or "lR2" (logarithmic R2, i.e. R2 calculated for the logarithm of the diversity-vs-time curve). Note that "wR2" will weight errors at lower diversities more strongly than "R2".
}
}


%\details{
%}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the fitting was successful.
}
\item{objective_value}{
Numeric. The achieved maximum value of the objective function (log-likelihood, R2 or weighted R2).
}
\item{parameters}{
A named list listing all model parameters (fixed and fitted).
}
\item{start_parameters}{
A named list listing the start values of all model parameters. In the case of multiple fitting trials, this will list the initial (non-randomized) guess.
}
\item{R2}{
Numeric. The achieved coefficient of determination of the fitted model, based on the diversity-vs-time curve.
}
\item{wR2}{
Numeric. The achieved weighted coefficient of determination of the fitted model, based on the diversity-vs-time curve. Weights of squared errors are proportional to the inverse squared diversities observed in the tree.
}
\item{lR2}{
Numeric. The achieved coefficient of determination of the fitted model on a log axis, i.e. based on the logarithm of the diversity-vs-time curve.
}
\item{Nspeciations}{
Integer. Number of speciation events (=nodes) considered during fitting. This only includes speciations visible in the tree.
}
\item{Nextinctions}{
Integer. Number of extinction events (=non-crown tips) considered during fitting. This only includes extinctions visible in the tree, i.e. tips whose distance from the root is lower than the maximum.
}
\item{grid_times}{
Numeric vector. Time points considered for the diversity-vs-time curve. Times will be constrained between \code{min_age} and \code{max_age} if these were specified.
}
\item{tree_diversities}{
Number of lineages represented in the tree through time, calculated for each of \code{grid_times}.
}
\item{model_diversities}{
Number of lineages through time as predicted by the model (in the deterministic limit), calculated for each of \code{grid_times}. If \code{coalescent==TRUE} then these are the number of lineages expected to be represented in the coalescent tree (this may be lower than the actual number of extant clades at any given time point, if the model includes extinctions).
}
\item{fitted_parameter_names}{
Character vector, listing the names of fitted (i.e. non-fixed) parameters.
}
\item{locally_fitted_parameters}{
Named list of numeric vectors, listing the fitted values for each parameter and for each fitting trial. For example, if \code{birth_rate_factor} was fitted, then \code{locally_fitted_parameters$birth_rate_factor} will be a numeric vector of size \code{Ntrials} (or less, if some trials failed or omitted), listing the locally-optimized values of the parameter for each considered fitting trial. Mainly useful for diagnostic purposes.
}
\item{objective}{
Character. The name of the objective function used for fitting ("LL", "R2" or "wR2").
}
\item{Ntips}{
The number of tips in the input tree.
}
\item{Nnodes}{
The number of nodes in the input tree.
}
\item{min_age}{
The minimum age of nodes/tips considered during fitting.
}
\item{max_age}{
The maximum age of nodes/tips considered during fitting.
}
\item{age_centile}{
Numeric or \code{NULL}, equal to the \code{age_centile} specified as input to the function.
}
}


\author{Stilianos Louca}

%\references{
%D. J. Aldous (2001). Stochastic models and descriptive statistics for phylogenetic trees, from Yule to today. Statistical Science. 16:23-34.
%
%M. Steel and A. McKenzie (2001). Properties of phylogenetic trees generated by Yule-type speciation models. Mathematical Biosciences. 170:91-112.
%}

\seealso{
\code{\link{generate_random_tree}},
\code{\link{simulate_diversification_model}}
\code{\link{reconstruct_past_diversification}}
}

\examples{
# Generate a tree using a simple speciation model
parameters = list(birth_rate_intercept  = 1, 
                  birth_rate_factor     = 0,
                  birth_rate_exponent   = 0,
                  death_rate_intercept  = 0,
                  death_rate_factor     = 0,
                  death_rate_exponent   = 0,
                  resolution            = 0,
                  rarefaction           = 1)
tree = generate_random_tree(parameters, max_tips=100)

# Fit model to the tree
fitting_parameters = parameters
fitting_parameters$birth_rate_intercept = NULL # fit only this parameter
fitting = fit_tree_model(tree,fitting_parameters)

# compare fitted to true value
T = parameters$birth_rate_intercept
F = fitting$parameters$birth_rate_intercept
cat(sprintf("birth_rate_intercept: true=\%g, fitted=\%g\n",T,F))
}
% Only 1 keyword per line
\keyword{fitting}
\keyword{tree model}
