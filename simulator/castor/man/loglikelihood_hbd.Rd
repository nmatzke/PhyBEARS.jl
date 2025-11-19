\name{loglikelihood_hbd}
\alias{loglikelihood_hbd}
\title{
Galculate the log-likelihood of a homogenous birth-death model.
}
\description{
Given a rooted ultrametric timetree, and a homogenous birth-death (HBD) model, i.e., with speciation rate \eqn{\lambda}, extinction rate \eqn{\mu} and sampling fraction \eqn{\rho}, calculate the likelihood of the tree under the model. The speciation and extinction rates may be time-dependent. ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction rates (in the literature this is sometimes referred to simply as ``birth-death model''). Alternatively to \eqn{\lambda} and \eqn{\mu}, the likelihood may also be calculated based on the pulled diversification rate (PDR; Louca et al. 2018) and the product \eqn{\rho(0)\cdot\lambda(0)}, or based on the pulled speciation rate (PSR). In either case, the time-profiles of \eqn{\lambda}, \eqn{\mu}, the PDR or the PSR are specified as piecewise polynomially functions (splines), defined on a discrete grid of ages.
}
\usage{
loglikelihood_hbd(tree, 
                  oldest_age        = NULL,
                  age0              = 0,
                  rho0              = NULL,
                  rholambda0        = NULL,
                  age_grid          = NULL,
                  lambda            = NULL,
                  mu                = NULL,
                  PDR               = NULL,
                  PSR               = NULL,
                  splines_degree    = 1,
                  condition         = "auto",
                  max_model_runtime = -1,
                  relative_dt       = 1e-3)
}
\arguments{
\item{tree}{
A rooted ultrametric tree of class "phylo".
}
\item{oldest_age}{
Strictly positive numeric, specifying the oldest time before present (``age'') to consider when calculating the likelihood. If this is equal to or greater than the root age, then \code{oldest_age} is taken as the stem age, and the classical formula by Morlon et al. (2011) is used.  If \code{oldest_age} is less than the root age, the tree is split into multiple subtrees at that age by treating every edge crossing that age as the stem of a subtree, and each subtree is considered an independent realization of the HBD model stemming at that age. This can be useful for avoiding points in the tree close to the root, where estimation uncertainty is generally higher. If \code{oldest_age==NULL}, it is automatically set to the root age.
}
\item{age0}{
Non-negative numeric, specifying the youngest age (time before present) to consider for fitting, and with respect to which \code{rho} and \code{rholambda0} are defined. If \code{age0>0}, then \code{rho} refers to the sampling fraction at age \code{age0}, and \code{rholambda0} to the product between \code{rho} and the speciation rate at age \code{age0}. See below for more details.
}
\item{rho0}{
Numeric between 0 (exclusive) and 1 (inclusive), specifying the sampling fraction of the tree at \code{age0}, i.e. the fraction of lineages extant at \code{age0} that are included in the tree. Note that if \eqn{rho0<1}, lineages extant at \code{age0} are assumed to have been sampled randomly at equal probabilities. Can also be \code{NULL}, in which case \code{rholambda0} and \code{PDR} (see below) must be provided.
}
\item{rholambda0}{
Strictly positive numeric, specifying the product of the sampling fraction and the speciation rateat \code{age0}, units 1/time. Can be \code{NULL}, in which case \code{rarefaction}, \code{lambda} and \code{mu} must be provided.
}
\item{age_grid}{
Numeric vector, listing discrete ages (time before present) on which either \eqn{\lambda} and \eqn{\mu}, or the PDR, are specified. Listed ages must be strictly increasing, and must cover at least the full considered age interval (from \code{age0} to \code{oldest_age}). Can also be \code{NULL} or a vector of size 1, in which case the speciation rate, extinction rate and PDR are assumed to be time-independent.
}
\item{lambda}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing speciation rates (in units 1/time) at the ages listed in \code{age_grid}. Speciation rates should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). If \code{NULL}, then either \code{PDR} and \code{rholambda0}, or \code{PSR} alone, must be provided.
}
\item{mu}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing extinction rates (in units 1/time)at the ages listed in \code{age_grid}. Extinction rates should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). If \code{NULL}, then \code{PDR} and \code{rholambda0}, or \code{PSR} alone, must be provided.
}
\item{PDR}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing pulled diversification rates (in units 1/time) at the ages listed in \code{age_grid}. PDRs can be negative or positive, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). If \code{NULL}, then either \code{lambda} and \code{mu}, or \code{PSR} alone, must be provided.
}
\item{PSR}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing pulled speciation rates (in units 1/time) at the ages listed in \code{age_grid}. PSRs should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). If \code{NULL}, then either \code{lambda} and \code{mu}, or \code{PDR} and \code{rholambda0}, must be provided.
}
\item{splines_degree}{
Integer, either 0,1,2 or 3, specifying the polynomial degree of the provided \code{lambda}, \code{mu}, \code{PDR} and \code{PSR} (whichever applicable) between grid points in \code{age_grid}. For example, if \code{splines_degree==1}, then the provided \code{lambda}, \code{mu}, \code{PDR} and \code{PSR} are interpreted as piecewise-linear curves; if \code{splines_degree==2} they are interpreted as quadratic splines; if \code{splines_degree==3} they are interpreted as cubic splines. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on.
}
\item{condition}{
Character, either "crown", "stem", "auto" or "none" (the last one is only available if \code{lambda} and \code{mu} are given), specifying on what to condition the likelihood. If "crown", the likelihood is conditioned on the survival of the two daughter lineages branching off at the root. If "stem", the likelihood is conditioned on the survival of the stem lineage. Note that "crown" really only makes sense when \code{oldest_age} is equal to the root age, while "stem" is recommended if \code{oldest_age} differs from the root age. "none" is usually not recommended and is only available when \code{lambda} and \code{mu} are provided. If "auto", the condition is chosen according to the recommendations mentioned earlier.
}
\item{max_model_runtime}{
Numeric, maximum allowed runtime (in seconds) for evaluating the likelihood. If the likelihood calculation takes longer than this (appoximate) threshold, it halts and returns with an error. If negative (default), this option is ignored.
}
\item{relative_dt}{
Strictly positive numeric (unitless), specifying the maximum relative time step allowed for integration over time. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
}
}
\details{
If \code{age0>0}, the input tree is essentially trimmed at \code{age0} (omitting anything younger than \code{age0}), and the is likelihood calculated for the trimmed tree while shifting time appropriately. In that case, \code{rho0} is interpreted as the sampling fraction at \code{age0}, i.e. the fraction of lineages extant at \code{age0} that are repreented in the tree. Similarly, \code{rholambda0} is the product of the sampling fraction and \eqn{\lambda} at \code{age0}.

This function supports three alternative parameterizations of HBD models, either using the speciation and extinction rates and sampling fraction (\eqn{\lambda}, \eqn{\mu} and \eqn{\rho(\tau_o)} (for some arbitrary age \eqn{\tau_o}), or using the pulled diversification rate (PDR) and the product \eqn{\rho(\tau_o)\cdot\lambda(\tau_o} (sampling fraction times speciation rate at \eqn{\tau_o}), or using the pulled speciation rate (PSR). The latter two options should be interpreted as a parameterization of congruence classes, i.e. sets of models that have the same likelihood, rather than specific models, since multiple combinations of \eqn{\lambda}, \eqn{\mu} and \eqn{\rho(\tau_o)} can have identical PDRs, \eqn{\rho(\tau_o)\cdot\lambda(\tau_o)} and PSRs (Louca and Pennell, in review).

For large trees the asymptotic time complexity of this function is O(Nips). The tree may include monofurcations as well as multifurcations, and the likelihood formula accounts for those (i.e., as if monofurcations were omitted and multifurcations were expanded into bifurcations).
}

\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the calculation was successful. If \code{FALSE}, then the returned list includes an additional `\code{error}' element (character) containing a description of the error; all other return variables may be undefined.
}
\item{loglikelihood}{
Numeric. If \code{success==TRUE}, this will be the natural logarithm of the likelihood of the tree under the given model.
}
}


\author{Stilianos Louca}

\references{
H. Morlon, T. L. Parsons, J. B. Plotkin (2011). Reconciling molecular phylogenies with the fossil record. Proceedings of the National Academy of Sciences. 108:16327-16332.

S. Louca et al. (2018). Bacterial diversification through geological time. Nature Ecology & Evolution. 2:1458-1467.

S. Louca and M. W. Pennell (in review as of 2019)
}

\seealso{
\code{\link{simulate_deterministic_hbd}}

\code{\link{fit_hbd_model_parametric}}

\code{\link{fit_hbd_model_on_grid}}

\code{\link{fit_hbd_pdr_on_grid}}

\code{\link{fit_hbd_pdr_parametric}}
}

\examples{
# generate a random tree with constant rates
Ntips  = 100
params = list(birth_rate_factor=1, death_rate_factor=0.2, rarefaction=0.5)
tree   = generate_random_tree(params, max_tips=Ntips, coalescent=TRUE)$tree

# get the loglikelihood for an HBD model with the same parameters that generated the tree
# in particular, assuming time-independent speciation & extinction rates
LL = loglikelihood_hbd( tree, 
                        rho0      = params$rarefaction, 
                        age_grid  = NULL, # assume time-independent rates
                        lambda    = params$birth_rate_factor,
                        mu        = params$death_rate_factor)
if(LL$success){
  cat(sprintf("Loglikelihood for constant-rates model = \%g\n",LL$loglikelihood))
}

# get the likelihood for a model with exponentially decreasing (in forward time) lambda & mu
beta      = 0.01 # exponential decay rate of lambda over time
age_grid  = seq(from=0, to=100, by=0.1) # choose a sufficiently fine age grid
lambda    = 1*exp(beta*age_grid) # define lambda on the age grid
mu        = 0.2*lambda # assume similarly shaped but smaller mu
LL = loglikelihood_hbd( tree, 
                        rho0      = params$rarefaction, 
                        age_grid  = age_grid,
                        lambda    = lambda,
                        mu        = mu)
if(LL$success){
  cat(sprintf("Loglikelihood for exponential-rates model = \%g\n",LL$loglikelihood))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{homogenous birth-death model}
