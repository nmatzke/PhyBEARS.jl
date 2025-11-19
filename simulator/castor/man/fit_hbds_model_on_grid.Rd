\name{fit_hbds_model_on_grid}
\alias{fit_hbds_model_on_grid}
\title{
Fit a homogenous birth-death-sampling model on a discrete time grid.
}
\description{
Given a timetree (potentially sampled through time and not necessarily ultrametric), fit a homogenous birth-death-sampling (HBDS) model in which speciation, extinction and lineage sampling occurs at some continuous (Poissonian) rates \eqn{\lambda}, \eqn{\mu} and \eqn{\psi}, which are defined on a fixed grid of discrete time points and assumed to vary polynomially between grid points. Sampled lineages are kept in the pool of extant lineages at some ``retention probability'' \eqn{\kappa}, which may also depend on time. In addition, this model can include concentrated sampling attempts (CSAs) at a finite set of discrete time points \eqn{t_1,..,t_m}. ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction/sampling rates. Every HBDS model is thus defined based on the values that \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} take over time, as well as the sampling probabilities \eqn{\rho_1,..,\rho_m} and retention probabilities \eqn{\kappa_1,..,\kappa_m} during the concentrated sampling attempts. This function estimates the values of \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} on each grid point, as well as the \eqn{\rho_1,..,\rho_m} and \eqn{\kappa_1,..,\kappa_m}, by maximizing the corresponding likelihood of the timetree. Special cases of this model (when rates are piecewise constant through time) are sometimes known as ``birth-death-skyline plots'' in the literature (Stadler 2013). In epidemiology, these models are often used to describe the phylogenies of viral strains sampled over the course of the epidemic.
}
\usage{
fit_hbds_model_on_grid( tree, 
                        root_age                = NULL,
                        oldest_age              = NULL,
                        age_grid                = NULL,
                        CSA_ages                = NULL,
                        min_lambda              = 0,
                        max_lambda              = +Inf,
                        min_mu                  = 0,
                        max_mu                  = +Inf,
                        min_psi                 = 0,
                        max_psi                 = +Inf,
                        min_kappa               = 0,
                        max_kappa               = 1,
                        min_CSA_probs           = 0,
                        max_CSA_probs           = 1,
                        min_CSA_kappas          = 0,
                        max_CSA_kappas          = 1,
                        guess_lambda            = NULL,
                        guess_mu                = NULL,
                        guess_psi               = NULL,
                        guess_kappa             = NULL,
                        guess_CSA_probs         = NULL,
                        guess_CSA_kappas        = NULL,
                        fixed_lambda            = NULL,
                        fixed_mu                = NULL,
                        fixed_psi               = NULL,
                        fixed_kappa             = NULL,
                        fixed_CSA_probs         = NULL,
                        fixed_CSA_kappas        = NULL,
                        fixed_age_grid          = NULL,
                        const_lambda            = FALSE,
                        const_mu                = FALSE,
                        const_psi               = FALSE,
                        const_kappa             = FALSE,
                        const_CSA_probs         = FALSE,
                        const_CSA_kappas        = FALSE,
                        splines_degree          = 1,
                        condition               = "auto",
                        ODE_relative_dt         = 0.001,
                        ODE_relative_dy         = 1e-3,
                        CSA_age_epsilon         = NULL,
                        Ntrials                 = 1,
                        max_start_attempts      = 1,
                        Nthreads                = 1,
                        max_model_runtime       = NULL,
                        Nbootstraps             = 0,
                        Ntrials_per_bootstrap   = NULL,
                        fit_control             = list(),
                        focal_param_values      = NULL,
                        verbose                 = FALSE,
                        diagnostics             = FALSE,
                        verbose_prefix          = "")
}
\arguments{
\item{tree}{
A timetree of class "phylo", representing the time-calibrated reconstructed phylogeny of a set of extant and/or extinct species. Tips of the tree are interpreted as terminally sampled lineages, while monofurcating nodes are interpreted as non-terminally sampled lineages, i.e., lineages sampled at some past time point and with subsequently sampled descendants.
}
\item{root_age}{
Positive numeric, specifying the age of the tree's root. Can be used to define a time offset, e.g. if the last tip was not actually sampled at the present. If \code{NULL}, this will be calculated from the tree and it will be assumed that the last tip was sampled at the present.
}
\item{oldest_age}{
Strictly positive numeric, specifying the oldest time before present (``age'') to consider when calculating the likelihood. If this is equal to or greater than the root age, then \code{oldest_age} is interpreted as the stem age.  If \code{oldest_age} is less than the root age, the tree is split into multiple subtrees at that age by treating every edge crossing that age as the stem of a subtree, and each subtree is considered an independent realization of the HBDS model stemming at that age. This can be useful for avoiding points in the tree close to the root, where estimation uncertainty is generally higher. If \code{oldest_age==NULL}, it is automatically set to the root age.
}
\item{age_grid}{
Numeric vector, listing ages in ascending order, on which \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} are fitted and allowed to vary independently. This grid must cover at least the age range from the present (age 0) to \code{oldest_age}. If \code{NULL} or of length <=1 (regardless of value), then \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} are assumed to be time-independent.
}
\item{CSA_ages}{
Optional numeric vector, listing ages (in ascending order) at which concentrated sampling attempts (CSAs) occurred. If \code{NULL}, it is assumed that no concentrated sampling attempts took place and that all tips were sampled according to the continuous sampling rate \code{psi}.
}
\item{min_lambda}{
Numeric vector of length Ngrid (=\code{max(1,length(age_grid))}), or a single numeric, specifying lower bounds for the fitted speciation rate \eqn{\lambda} at each point in the age grid. If a single numeric, the same lower bound applies at all ages.
}
\item{max_lambda}{
Numeric vector of length Ngrid, or a single numeric, specifying upper bounds for the fitted speciation rate \eqn{\lambda} at each point in the age grid. If a single numeric, the same upper bound applies at all ages. Use \code{+Inf} to omit upper bounds.
}
\item{min_mu}{
Numeric vector of length Ngrid, or a single numeric, specifying lower bounds for the fitted extinction rate \eqn{\mu} at each point in the age grid. If a single numeric, the same lower bound applies at all ages.
}
\item{max_mu}{
Numeric vector of length Ngrid, or a single numeric, specifying upper bounds for the fitted extinction rate \eqn{\mu} at each point in the age grid. If a single numeric, the same upper bound applies at all ages. Use \code{+Inf} to omit upper bounds.
}
\item{min_psi}{
Numeric vector of length Ngrid, or a single numeric, specifying lower bounds for the fitted Poissonian sampling rate \eqn{\psi} at each point in the age grid. If a single numeric, the same lower bound applies at all ages.
}
\item{max_psi}{
Numeric vector of length Ngrid, or a single numeric, specifying upper bounds for the fitted Poissonian sampling rate \eqn{\psi} at each point in the age grid. If a single numeric, the same upper bound applies at all ages. Use \code{+Inf} to omit upper bounds.
}
\item{min_kappa}{
Numeric vector of length Ngrid, or a single numeric, specifying lower bounds for the fitted retention probability \eqn{\kappa} at each point in the age grid. If a single numeric, the same lower bound applies at all ages.
}
\item{max_kappa}{
Numeric vector of length Ngrid, or a single numeric, specifying upper bounds for the fitted retention probability \eqn{\kappa} at each point in the age grid. If a single numeric, the same upper bound applies at all ages. Use \code{+Inf} to omit upper bounds.
}
\item{min_CSA_probs}{
Numeric vector of length NCSA (=\code{length(CSA_ages)}), or a single numeric, specifying lower bounds for the fitted sampling probabilities \eqn{\rho_1},..,\eqn{\rho_m} at each concentrated sampling attempt. If a single numeric, the same lower bound applies at all CSAs. Note that, since \eqn{\rho_1}, \eqn{\rho_2}, ... are probabilities, \code{min_CSA_probs} should not be negative.
}
\item{max_CSA_probs}{
Numeric vector of length NCSA, or a single numeric, specifying upper bounds for the fitted sampling probabilities \eqn{\rho_1}, \eqn{\rho_2}, ... at each concentrated sampling attempt. If a single numeric, the same upper bound applies at all CSAs. Note that, since \eqn{\rho_1}, \eqn{\rho_2}, ... are probabilities, \code{max_CSA_probs} should not be greater than 1.
}
\item{min_CSA_kappas}{
Numeric vector of length NCSA, or a single numeric, specifying lower bounds for the fitted retention probabilities \eqn{\kappa_1}, \eqn{\kappa_2}, ... at each concentrated sampling attempt. If a single numeric, the same lower bound applies at all CSAs. Note that, since \eqn{\kappa_1}, \eqn{\kappa_2}, ... are probabilities, \code{min_CSA_kappas} should not be negative.
}
\item{max_CSA_kappas}{
Numeric vector of length NCSA, or a single numeric, specifying upper bounds for the fitted sampling probabilities \eqn{\kappa_1}, \eqn{\kappa_2}, ... at each concentrated sampling attempt. If a single numeric, the same upper bound applies at all CSAs. Note that, since \eqn{\kappa_1}, \eqn{\kappa_2}, .. are probabilities, \code{max_CSA_kappas} should not be greater than 1.
}
\item{guess_lambda}{
Initial guess for \eqn{\lambda} at each age-grid point. Either \code{NULL} (an initial guess will be computed automatically), or a single numeric (guessing the same \eqn{\lambda} at all ages) or a numeric vector of size Ngrid specifying a separate guess for \eqn{\lambda} at each age-grid point. To omit an initial guess for some but not all age-grid points, set their guess values to \code{NA}. Guess values are ignored for non-fitted (i.e., fixed) parameters.
}
\item{guess_mu}{
Initial guess for \eqn{\mu} at each age-grid point. Either \code{NULL} (an initial guess will be computed automatically), or a single numeric (guessing the same \eqn{\mu} at all ages) or a numeric vector of size Ngrid specifying a separate guess for \eqn{\mu} at each age-grid point. To omit an initial guess for some but not all age-grid points, set their guess values to \code{NA}. Guess values are ignored for non-fitted (i.e., fixed) parameters.
}
\item{guess_psi}{
Initial guess for \eqn{\psi} at each age-grid point. Either \code{NULL} (an initial guess will be computed automatically), or a single numeric (guessing the same \eqn{\psi} at all ages) or a numeric vector of size Ngrid specifying a separate guess for \eqn{\psi} at each age-grid point. To omit an initial guess for some but not all age-grid points, set their guess values to \code{NA}. Guess values are ignored for non-fitted (i.e., fixed) parameters.
}
\item{guess_kappa}{
Initial guess for \eqn{\kappa} at each age-grid point. Either \code{NULL} (an initial guess will be computed automatically), or a single numeric (guessing the same \eqn{\kappa} at all ages) or a numeric vector of size Ngrid specifying a separate guess for \eqn{\kappa} at each age-grid point. To omit an initial guess for some but not all age-grid points, set their guess values to \code{NA}. Guess values are ignored for non-fitted (i.e., fixed) parameters.
}
\item{guess_CSA_probs}{
Initial guess for the \eqn{\rho_1}, \eqn{\rho_2}, ... at each concentrated sampling attempt. Either \code{NULL} (an initial guess will be computed automatically), or a single numeric (guessing the same value at every CSA) or a numeric vector of size NCSA specifying a separate guess at each CSA. To omit an initial guess for some but not all CSAs, set their guess values to \code{NA}. Guess values are ignored for non-fitted (i.e., fixed) parameters.
}
\item{guess_CSA_kappas}{
Initial guess for the \eqn{\kappa_1}, \eqn{\kappa_2}, ... at each concentrated sampling attempt. Either \code{NULL} (an initial guess will be computed automatically), or a single numeric (guessing the same value at every CSA) or a numeric vector of size NCSA specifying a separate guess at each CSA. To omit an initial guess for some but not all CSAs, set their guess values to \code{NA}. Guess values are ignored for non-fitted (i.e., fixed) parameters.
}
\item{fixed_lambda}{
Optional fixed (i.e. non-fitted) \eqn{\lambda} values on one or more age-grid points. Either \code{NULL} (\eqn{\lambda} is not fixed anywhere), or a single numeric (\eqn{\lambda} fixed to the same value at all grid points) or a numeric vector of size Ngrid (if \code{fixed_age_grid=NULL}; \eqn{\lambda} fixed on one or more age-grid points, use \code{NA} for non-fixed values) or a numeric vector of the same size as \code{fixed_age_grid} (if \code{fixed_age_grid!=NULL}, in which case all entries in \code{fixed_lambda} must be finite numbers).
}
\item{fixed_mu}{
Optional fixed (i.e. non-fitted) \eqn{\mu} values on one or more age-grid points. Either \code{NULL} (\eqn{\mu} is not fixed anywhere), or a single numeric (\eqn{\mu} fixed to the same value at all grid points) or a numeric vector of size Ngrid (if \code{fixed_age_grid=NULL}; \eqn{\mu} fixed on one or more age-grid points, use \code{NA} for non-fixed values) or a numeric vector of the same size as \code{fixed_age_grid} (if \code{fixed_age_grid!=NULL}, in which case all entries in \code{fixed_mu} must be finite numbers).
}
\item{fixed_psi}{
Optional fixed (i.e. non-fitted) \eqn{\psi} values on one or more age-grid points. Either \code{NULL} (\eqn{\psi} is not fixed anywhere), or a single numeric (\eqn{\psi} fixed to the same value at all grid points) or a numeric vector of size Ngrid (if \code{fixed_age_grid=NULL}; \eqn{\psi} fixed on one or more age-grid points, use \code{NA} for non-fixed values) or a numeric vector of the same size as \code{fixed_age_grid} (if \code{fixed_age_grid!=NULL}, in which case all entries in \code{fixed_psi} must be finite numbers).
}
\item{fixed_kappa}{
Optional fixed (i.e. non-fitted) \eqn{\kappa} values on one or more age-grid points. Either \code{NULL} (\eqn{\kappa} is not fixed anywhere), or a single numeric (\eqn{\kappa} fixed to the same value at all grid points) or a numeric vector of size Ngrid (if \code{fixed_age_grid=NULL}; \eqn{\kappa} fixed on one or more age-grid points, use \code{NA} for non-fixed values) or a numeric vector of the same size as \code{fixed_age_grid} (if \code{fixed_age_grid!=NULL}, in which case all entries in \code{fixed_kappa} must be finite numbers).
}
\item{fixed_CSA_probs}{
Optional fixed (i.e. non-fitted) \eqn{\rho_1}, \eqn{\rho_2}, ... values on one or more age-grid points. Either \code{NULL} (none of the \eqn{\rho_1}, \eqn{\rho_2},... are fixed), or a single numeric (\eqn{\rho_1}, \eqn{\rho_2},... are fixed to the same value at all CSAs) or a numeric vector of size NCSA (one or more of the \eqn{\rho_1}, \eqn{\rho_2}, ... are fixed, use \code{NA} for non-fixed values).
}
\item{fixed_CSA_kappas}{
Optional fixed (i.e. non-fitted) \eqn{\kappa_1}, \eqn{\kappa_2}, ... values on one or more age-grid points. Either \code{NULL} (none of the \eqn{\kappa_1}, \eqn{\kappa_2},... are fixed), or a single numeric (\eqn{\kappa_1}, \eqn{\kappa_2},... are fixed to the same value at all CSAs) or a numeric vector of size NCSA (one or more of the \eqn{\kappa_1}, \eqn{\kappa_2}, ... are fixed, use \code{NA} for non-fixed values).
}
\item{fixed_age_grid}{
Optional numeric vector, specifying an age grid on which \code{fixed_lambda}, \code{fixed_mu}, \code{fixed_psi} and \code{fixed_kappa} (whichever is provided) are defined instead of on the \code{age_grid}. If \code{fixed_age_grid} is provided, then each of \code{fixed_lambda}, \code{fixed_mu}, \code{fixed_psi} and \code{fixed_kappa} must be defined (i.e. have a finite non-negative value) on every point in \code{fixed_age_grid}. Entries in \code{fixed_age_grid} must be in ascending order and must cover at least the ages 0 to \code{oldest_age}.

This option may be useful if you want to fit some parameters on a coarse grid, but want to specify (fix) some other parameters on a much finer grid. Also note that if \code{fixed_age_grid} is used, all parameters \code{lambda}, \code{mu}, \code{psi} and \code{kappa} are internally re-interpolated onto \code{fixed_age_grid} when evaluating the likelihood; hence, in general \code{fixed_age_grid} should be much finer than \code{age_grid}. In most situations you would probably want to keep the default \code{fixed_age_grid=NULL}.
}
\item{const_lambda}{
Logical, specifying whether \eqn{\lambda} should be assumed constant across the grid, i.e. time-independent. Setting \code{const_lambda=TRUE} reduces the number of free (i.e., independently fitted) parameters. If \eqn{\lambda} is fixed on some grid points (i.e. via \code{fixed_lambda}), then only the non-fixed lambdas are assumed to be identical to one another.
}
\item{const_mu}{
Logical, specifying whether \eqn{\mu} should be assumed constant across the grid, i.e. time-independent. Setting \code{const_mu=TRUE} reduces the number of free (i.e., independently fitted) parameters. If \eqn{\mu} is fixed on some grid points (i.e. via \code{fixed_mu}), then only the non-fixed mus are assumed to be identical to one another.
}
\item{const_psi}{
Logical, specifying whether \eqn{\psi} should be assumed constant across the grid, i.e. time-independent. Setting \code{const_psi=TRUE} reduces the number of free (i.e., independently fitted) parameters. If \eqn{\psi} is fixed on some grid points (i.e. via \code{fixed_psi}), then only the non-fixed psis are assumed to be identical to one another.
}
\item{const_kappa}{
Logical, specifying whether \eqn{\kappa} should be assumed constant across the grid, i.e. time-independent. Setting \code{const_kappa=TRUE} reduces the number of free (i.e., independently fitted) parameters. If \eqn{\kappa} is fixed on some grid points (i.e. via \code{fixed_kappa}), then only the non-fixed kappas are assumed to be identical to one another.
}
\item{const_CSA_probs}{
Logical, specifying whether the \eqn{\rho_1}, \eqn{\rho_2}, ... should be the same across all CSAs. Setting \code{const_CSA_probs=TRUE} reduces the number of free (i.e., independently fitted) parameters. If some of the \eqn{\rho_1}, \eqn{\rho_2}, ... are fixed (i.e. via \code{fixed_CSA_probs}), then only the non-fixed CSA_probs are assumed to be identical to one another.
}
\item{const_CSA_kappas}{
Logical, specifying whether the \eqn{\kappa_1}, \eqn{\kappa_2}, ... should be the same across all CSAs. Setting \code{const_CSA_kappas=TRUE} reduces the number of free (i.e., independently fitted) parameters. If some of the \eqn{\kappa_1}, \eqn{\kappa_2}, ... are fixed (i.e. via \code{fixed_CSA_kappas}), then only the non-fixed CSA_kappas are assumed to be identical to one another.
}

\item{splines_degree}{
Integer between 0 and 3 (inclusive), specifying the polynomial degree of \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} between age-grid points. If 0, then \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} are considered piecewise constant, if 1 they are considered piecewise linear, if 2 or 3 they are considered to be splines of degree 2 or 3, respectively. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on. A degree of 0 is generally not recommended. The case \code{splines_degree=0} is also known as ``skyline'' model.
}
\item{condition}{
Character, either "crown", "stem", "none" or "auto", specifying on what to condition the likelihood. If "crown", the likelihood is conditioned on the survival of the two daughter lineages branching off at the root. If "stem", the likelihood is conditioned on the survival of the stem lineage. Note that "crown" really only makes sense when \code{oldest_age} is equal to the root age, while "stem" is recommended if \code{oldest_age} differs from the root age. "none" is generally not recommended. If "auto", the condition is chosen according to the above recommendations.
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
Integer, specifying the number of times to attempt finding a valid start point (per trial) before giving up on that trial. Randomly choosen extreme start parameters may occasionally result in Inf/undefined likelihoods, so this option allows the algorithm to keep looking for valid starting points.
}
\item{Nthreads}{
Integer, specifying the number of parallel threads to use for performing multiple fitting trials simultaneously. This should generally not exceed the number of available CPUs on your machine. Parallel computing is not available on the Windows platform.
}
\item{max_model_runtime}{
Optional numeric, specifying the maximum number of seconds to allow for each evaluation of the likelihood function. Use this to abort fitting trials leading to parameter regions where the likelihood takes a long time to evaluate (these are often unlikely parameter regions).
}
\item{Nbootstraps}{
Integer, specifying the number of parametric bootstraps to perform for estimating standard errors and confidence intervals of estimated parameters. Set to 0 for no bootstrapping.
}
\item{Ntrials_per_bootstrap}{
Integer, specifying the number of fitting trials to perform for each bootstrap sampling. If \code{NULL}, this is set equal to \code{max(1,Ntrials)}. Decreasing \code{Ntrials_per_bootstrap} will reduce computation time, at the expense of potentially inflating the estimated confidence intervals; in some cases (e.g., for very large trees) this may be useful if fitting takes a long time and confidence intervals are very narrow anyway. Only relevant if \code{Nbootstraps>0}.
}
\item{fit_control}{
 Named list containing options for the \code{nlminb} optimization routine, such as \code{iter.max}, \code{eval.max} or \code{rel.tol}. For a complete list of options and default values see the documentation of \code{nlminb} in the \code{stats} package.
}
\item{focal_param_values}{
Optional list, listing combinations of parameter values of particular interest and for which the log-likelihoods should be returned. Every element of this list should itself be a named list, containing the elements \code{lambda}, \code{mu}, \code{psi} and \code{kappa} (each being a numeric vector of size NG) as well as the elements \code{CSA_probs} and \code{CSA_kappas} (each being a numeric vector of size NCSA). This may be used e.g. for diagnostic purposes, e.g. to examine the shape of the likelihood function.
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
Warning: In the absence of concentrated sampling attempts (NCSA=0), and without well-justified a priori constraints on either \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and/or \eqn{\kappa}, it is generally impossible to reliably estimate \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} from timetrees alone. This routine (and any other software that claims to estimate \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} solely from timetrees) should thus be treated with great suspicion. Many epidemiological models make the (often reasonable assumption) that \eqn{\kappa=0}; note that even in this case, one generally can't co-estimate \eqn{\lambda}, \eqn{\mu} and \eqn{\psi} from the timetree alone.

It is advised to provide as much information to the function \code{fit_hbds_model_on_grid} as possible, including reasonable lower and upper bounds (\code{min_lambda}, \code{max_lambda}, \code{min_mu}, \code{max_mu}, \code{min_psi}, \code{max_psi}, \code{min_kappa}, \code{max_kappa}) and reasonable parameter guesses. It is also important that the \code{age_grid} is sufficiently fine to capture the expected major variations of \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} over time, but keep in mind the serious risk of overfitting when \code{age_grid} is too fine and/or the tree is too small. The \code{age_grid} does not need to be uniform, i.e., you may want to use a finer grid in regions where there's more data (tips) available. If strong lower and upper bounds are not available and fitting takes a long time to run, consider using the option \code{max_model_runtime} to limit how much time the fitting allows for each evaluation of the likelihood.

Note that here "age" refers to time before present, i.e., age increases from tips to root and age 0 is present-day. CSAs are enumerated in the order of increasing age, i.e., from the present to the past. Similarly, the age grid specifies time points from the present towards the past.
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
\item{guess_loglikelihood}{
The log-likelihood of the guessed model for the given timetree.
}
\item{param_fitted}{
Named list, specifying the fixed and fitted model parameters. This list will contain the elements \code{lambda}, \code{mu}, \code{psi} and \code{kappa} (each being a numeric vector of size NG, listing \eqn{\lambda},\eqn{\mu}, \eqn{\psi} and \eqn{\kappa} at each age-grid point) as well as the elements \code{CSA_probs} and \code{CSA_kappas} (each being a numeric vector of size NCSA).
}
\item{param_guess}{
Named list, specifying the guessed model parameters. This list will contain the elements \code{lambda}, \code{mu}, \code{psi} and \code{kappa} (each being a numeric vector of size NG) as well as the elements \code{CSA_probs} and \code{CSA_kappas} (each being a numeric vector of size NCSA). Between grid points \eqn{\lambda} should be interpreted as a piecewise polynomial function (natural spline) of degree \code{splines_degree}; to evaluate this function at arbitrary ages use the \code{castor} routine \code{\link{evaluate_spline}}. The same also applies to \eqn{\mu}, \eqn{\psi} and \eqn{\kappa}.
}
\item{age_grid}{
Numeric vector of size NG, the age-grid on which \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} are defined. This will be the same as the provided \code{age_grid}, unless the latter was \code{NULL} or of length <=1.
}
\item{CSA_ages}{
Numeric vector of size NCSA, ting listhe ages at which concentrated sampling attempts occurred. This is the same as provided to the function.
}
\item{NFP}{
Integer, number of free (i.e., independently) fitted parameters. If none of the \eqn{\lambda}, \eqn{\mu} and \eqn{\rho} were fixed, and \code{const_lambda=FALSE} and \code{const_mu=FALSE}, then \code{NFP} will be equal to 2*Ngrid+1.
}
\item{Ndata}{
Integer, the number of data points (sampling and branching events) used for fitting.
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
\item{standard_errors}{
Named list specifying the standard errors of the parameters, based on parametric bootstrapping. This list will contain the elements \code{lambda}, \code{mu}, \code{psi} and \code{kappa} (each being a numeric vector of size NG) as well as the elements \code{CSA_probs} and \code{CSA_kappas} (each being a numeric vector of size NCSA). Only included if \code{Nbootstraps>0}. Note that the standard errors of non-fitted (i.e., fixed) parameters will be zero.
}
\item{CI50lower}{
Named list specifying the lower end of the 50\% confidence interval (i.e. the 25\% quantile) for each parameter, based on parametric bootstrapping. This list will contain the elements \code{lambda}, \code{mu}, \code{psi} and \code{kappa} (each being a numeric vector of size NG) as well as the elements \code{CSA_probs} and \code{CSA_kappas} (each being a numeric vector of size NCSA). Only included if \code{Nbootstraps>0}.
}
\item{CI50upper}{
Similar to \code{CI50lower}, but listing the upper end of the 50\% confidence interval (i.e. the 75\% quantile) for each parameter. For example, the confidence interval for \eqn{\lambda} at age \code{age_grid[1]} will be between \code{CI50lower$lambda[1]} and \code{CI50upper$lambda[1]}. Only included if \code{Nbootstraps>0}.
}
\item{CI95lower}{
Similar to \code{CI50lower}, but listing the lower end of the 95\% confidence interval (i.e. the 2.5\% quantile) for each parameter. Only included if \code{Nbootstraps>0}.
}
\item{CI95upper}{
Similar to \code{CI50upper}, but listing the upper end of the 95\% confidence interval (i.e. the 97.5\% quantile) for each parameter. Only included if \code{Nbootstraps>0}.
}
\item{consistency}{
Numeric between 0 and 1, estimated consistency of the data with the fitted model. If \eqn{L} denotes the loglikelihood of new data generated by the fitted model (under the same model) and \eqn{M} denotes the expectation of \eqn{L}, then \code{consistency} is the probability that \eqn{|L-M|} will be greater or equal to \eqn{|X-M|}, where \eqn{X} is the loglikelihood of the original data under the fitted model. Only returned if \code{Nbootstraps>0}. A low consistency (e.g., <0.05) indicates that the fitted model is a poor description of the data. See Lindholm et al. (2019) for background.
}
}


\author{Stilianos Louca}

\references{
T. Stadler, D. Kuehnert, S. Bonhoeffer, A. J. Drummond (2013). Birth-death skyline plot reveals temporal changes of epidemic spread in HIV and hepatitis C virus (HCV). PNAS. 110:228-233.

A. Lindholm, D. Zachariah, P. Stoica, T. B. Schoen (2019). Data consistency approach to model validation. IEEE Access. 7:59788-59796.
}

\seealso{
\code{\link{simulate_deterministic_hbds}}, 
\code{\link{fit_hbds_model_parametric}}
}

\examples{
\dontrun{
# define lambda & mu & psi as functions of time
# Assuming an exponentially varying lambda & mu, and a constant psi
time2lambda = function(times){ 2*exp(0.1*times) }
time2mu     = function(times){ 0.1*exp(0.09*times) }
time2psi    = function(times){ rep(0.2, times=length(times)) }

# define concentrated sampling attempts
CSA_times   = c(3,4)
CSA_probs   = c(0.1, 0.2)

# generate random tree based on lambda, mu & psi
# assume that all sampled lineages are removed from the pool (i.e. kappa=0)
time_grid = seq(from=0, to=100, by=0.01)
simul = generate_tree_hbds( max_time    = 5,
                            time_grid   = time_grid,
                            lambda      = time2lambda(time_grid),
                            mu          = time2mu(time_grid),
                            psi         = time2psi(time_grid),
                            kappa       = 0,
                            CSA_times   = CSA_times,
                            CSA_probs   = CSA_probs,
                            CSA_kappas  = 0)
tree     = simul$tree
root_age = simul$root_age
cat(sprintf("Tree has \%d tips\n",length(tree$tip.label)))

# Define an age grid on which lambda_function & mu_function shall be fitted
fit_age_grid = seq(from=0,to=root_age,length.out=3)

# Fit an HBDS model on a grid
# Assume that psi is known and that sampled lineages are removed from the pool
# Hence, we only fit lambda & mu & CSA_probs
cat(sprintf("Fitting model to tree..\n"))
fit = fit_hbds_model_on_grid(tree, 
                             root_age           = root_age,
                             age_grid           = fit_age_grid,
                             CSA_ages           = rev(simul$final_time - CSA_times),
                             fixed_psi          = time2psi(simul$final_time-fit_age_grid),
                             fixed_kappa        = 0,
                             fixed_CSA_kappas   = 0,
                             Ntrials            = 4,
                             Nthreads           = 4,
                             Nbootstraps        = 0,
                             verbose            = TRUE,
                             verbose_prefix     = "  ")
if(!fit$success){
    cat(sprintf("ERROR: Fitting failed: \%s\n",fit$error))
}else{
    # compare fitted lambda to true lambda
    plot(x=fit$age_grid, 
         y=fit$param_fitted$lambda, 
         type='l', 
         col='#000000', 
         xlim=c(root_age,0),
         xlab='age', 
         ylab='lambda')
    lines(x=simul$final_time-time_grid, 
          y=time2lambda(time_grid), 
          type='l', 
          col='#0000AA')
}


# compare true and fitted model in terms of their LTTs
LTT      = castor::count_lineages_through_time(tree, Ntimes=100, include_slopes=TRUE)
LTT$ages = root_age - LTT$times

cat(sprintf("Simulating deterministic HBDS (true model)..\n"))
age0 = 0.5 # reference age at which to equate LTTs
LTT0 = approx(x=LTT$ages, y=LTT$lineages, xout=age0)$y # tree LTT at age0
fsim = simulate_deterministic_hbds( age_grid        = fit$age_grid,
                                    lambda          = fit$param_fitted$lambda,
                                    mu              = fit$param_fitted$mu,
                                    psi             = fit$param_fitted$psi,
                                    kappa           = fit$param_fitted$kappa,
                                    CSA_ages        = fit$CSA_ages,
                                    CSA_probs       = fit$param_fitted$CSA_probs,
                                    CSA_kappas      = fit$param_fitted$CSA_kappas,
                                    requested_ages  = seq(0,root_age,length.out=200),
                                    age0            = age0,
                                    LTT0            = LTT0,
                                    splines_degree  = 1)
if(!fsim$success){
    cat(sprintf("ERROR: Could not simulate fitted model: \%s\n",fsim$error))
    stop()
}
plot(x=LTT$ages, y=LTT$lineages, type='l', col='#0000AA', lwd=2, xlim=c(root_age,0))
lines(x=fsim$ages, y=fsim$LTT, type='l', col='#000000', lwd=2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{birth-death-sampling model}
\keyword{skyline model}
