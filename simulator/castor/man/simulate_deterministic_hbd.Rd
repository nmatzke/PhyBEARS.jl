\name{simulate_deterministic_hbd}
\alias{simulate_deterministic_hbd}
\title{
Simulate a deterministic homogenous birth-death model.
}
\description{
Given a homogenous birth-death (HBD) model, i.e., with speciation rate \eqn{\lambda}, extinction rate \eqn{\mu} and sampling fraction \eqn{\rho}, calculate various deterministic features of the model backwards in time, such as the total diversity over time. The speciation and extinction rates may be time-dependent. ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction rates (in the literature this is sometimes referred to simply as ``birth-death model''; Morlon et al. 2011). ``Deterministic'' refers to the fact that all calculated properties are completely determined by the model's parameters (i.e. non-random), as if an infinitely large tree was generated (aka. ``continuum limit'').

Alternatively to \eqn{\lambda}, one may provide the pulled diversification rate (PDR; Louca et al. 2018) and the speciation rate at some fixed age, \eqn{\lambda(\tau_o)}. Similarly, alternatively to \eqn{\mu}, one may provide the ratio of extinction over speciation rate, \eqn{\mu/\lambda}. In either case, the time-profiles of \eqn{\lambda}, \eqn{\mu}, \eqn{\mu/\lambda} or the PDR are specified as piecewise polynomial functions (splines), defined on a discrete grid of ages.
}
\usage{
simulate_deterministic_hbd(LTT0,
                           oldest_age,
                           age0           = 0,
                           rho0           = 1,
                           age_grid       = NULL,
                           lambda         = NULL,
                           mu             = NULL,
                           mu_over_lambda = NULL,
                           PDR            = NULL,
                           lambda0        = NULL,
                           splines_degree = 1,
                           relative_dt    = 1e-3,
                           allow_unreal   = FALSE)
}
\arguments{
\item{LTT0}{
The assumed number of sampled extant lineages at \code{age0}, defining the necessary initial condition for the simulation. If the HBD model is supposed to describe a specific timetree, then \code{LTT0} should correspond to the number of lineages in the tree ("lineages through time") at age \code{age0}.
}
\item{oldest_age}{
Strictly positive numeric, specifying the oldest time before present (``age'') to include in the simulation.
}
\item{age0}{
Non-negative numeric, specifying the age at which \code{LTT0}, \code{lambda0} and \code{rho} are given. Typically this will be 0, i.e., corresponding to the present.
}
\item{rho0}{
Numeric between 0 (exclusive) and 1 (inclusive), specifying the sampling fraction of the tree at \code{age0}, i.e. the fraction of lineages extant at \code{age0} that are included in the tree (aka. "rarefaction"). Note that if \code{rho0<1}, lineages extant at \code{age0} are assumed to have been sampled randomly at equal probabilities. Can also be \code{NULL}, in which case \code{rho0=1} is assumed.
}
\item{age_grid}{
Numeric vector, listing discrete ages (time before present) on which either \eqn{\lambda} and \eqn{\mu}, or the PDR and \eqn{\mu}, are specified. Listed ages must be strictly increasing, and must cover at least the full considered age interval (from \code{age0} to \code{oldest_age}). Can also be \code{NULL} or a vector of size 1, in which case the speciation rate, extinction rate and PDR are assumed to be time-independent.
}
\item{lambda}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing speciation rates (\eqn{\lambda}, in units 1/time) at the ages listed in \code{age_grid}. Speciation rates should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). If \code{NULL}, then \code{PDR} and \code{lambda0} must be provided.
}
\item{mu}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing extinction rates (\eqn{\mu}, in units 1/time) at the ages listed in \code{age_grid}. Extinction rates should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). Either \code{mu} or \code{mu_over_lambda} must be provided, but not both.
}
\item{mu_over_lambda}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing the ratio of extinction rates over speciation rates (\eqn{\mu/\lambda}) at the ages listed in \code{age_grid}. These ratios should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). Either \code{mu} or \code{mu_over_lambda} must be provided, but not both.
}
\item{PDR}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing pulled diversification rates (in units 1/time) at the ages listed in \code{age_grid}. PDRs can be negative or positive, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). If \code{NULL}, then \code{lambda} must be provided.
}
\item{lambda0}{
Non-negative numeric, specifying the speciation rate (in units 1/time) at \code{age0}. Either \code{lambda0} or \code{lambda} must be provided, but not both.
}
\item{splines_degree}{
Integer, either 0,1,2 or 3, specifying the polynomial degree of the provided \code{lambda}, \code{mu} and \code{PDR} between grid points in \code{age_grid}. For example, if \code{splines_degree==1}, then the provided \code{lambda}, \code{mu} and \code{PDR} are interpreted as piecewise-linear curves; if \code{splines_degree==2} they are interpreted as quadratic splines; if \code{splines_degree==3} they are interpreted as cubic splines. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on.
}
\item{relative_dt}{
Strictly positive numeric (unitless), specifying the maximum relative time step allowed for integration over time. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
}
\item{allow_unreal}{
Logical, specifying whether HBD models with unrealistic parameters (e.g., negative \eqn{\mu}) should be supported. This may be desired for example when examining model congruence classes with negative \eqn{\mu}.
}
}
\details{
This function supports the following alternative parameterizations of HBD models:
\itemize{
\item Using the speciation rate \eqn{\lambda} and extinction rate \eqn{\mu}.
\item Using the speciation rate \eqn{\lambda} and the ratio \eqn{\mu/\lambda}.
\item Using the pulled diversification rate (PDR), the extinction rate and the speciation rate given at some fixed \code{age0} (i.e. \code{lambda0}).
\item Using the PDR, the ratio \eqn{\mu/\lambda} and the speciation rate at some fixed \code{age0}.
}
The PDR is defined as \eqn{PDR = \lambda-\mu+\lambda^{-1}d\lambda/d\tau}, where \eqn{\tau} is age (time before present). To avoid ambiguities, only one of the above parameterizations is accepted at a time. The sampling fraction at \code{age0} (i.e., \code{rho0}) should always be provided; setting it to \code{NULL} is equivalent to setting it to 1.

Note that in the literature the sampling fraction usually refers to the fraction of lineages extant at present-day that have been sampled (included in the tree); this present-day sampling fraction is then used to parameterize birth-death cladogenic models. The sampling fraction can however be generalized to past times, by defining it as the probability that a lineage extant at any given time is included in the tree. The simulation function presented here allows for specifying this generalized sampling fraction at any age of choice, not just present-day.

The simulated LTT refers to a hypothetical tree sampled at age \code{age_grid[1]}, i.e. LTT(t) will be the number of lineages extant at age t that survived until age \code{age_grid[1]} and have been sampled, given that the fraction of sampled extant lineages at \code{age0} is \code{rho0}. Similarly, the returned Pextinct(t) (see below) is the probability that a lineage extant at age t would not survive until \code{age_grid[1]}. The same convention is used for the returned variables \code{Pmissing}, \code{shadow_diversity}, \code{PER}, \code{PSR}, \code{SER} and \code{PND}.
}

\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the calculation was successful. If \code{FALSE}, then the returned list includes an additional `\code{error}' element (character) providing a description of the error; all other return variables may be undefined.
}
\item{ages}{Numerical vector of size NG, listing discrete ages (time before present) on which all returned time-curves are specified. Listed ages will be in ascending order, will cover exactly the range \code{age_grid[1]} - \code{oldest_age}, may be irregularly spaced, and may be finer than the original provided \code{age_grid}. Note that \code{ages[1]} corresponds to the latest time point (closer to the tips), while \code{ages[NG]} corresponds to the oldest time point (\code{oldest_age}).}
\item{total_diversity}{Numerical vector of size NG, listing the predicted (deterministic) total diversity (number of extant species, denoted \eqn{N}) at the ages given in \code{ages[]}.}
\item{shadow_diversity}{Numerical vector of size NG, listing the predicted (deterministic) ``shadow diversity'' at the ages given in \code{ages[]}. The shadow diversity is defined as \eqn{N_s=N\cdot \rho(\tau_o)\lambda(\tau_o)/\lambda}, where \eqn{\tau_o} is \code{age0}.}
\item{Pmissing}{Numeric vector of size NG, listing the probability that a lineage, extant at a given age, will be absent from the tree either due to extinction or due to incomplete sampling.}
\item{Pextinct}{Numeric vector of size NG, listing the probability that a lineage, extant at a given age, will be fully extinct at present. Note that always \code{Pextinct<=Pmissing}.}
\item{LTT}{Numeric vector of size NG, listing the number of lineages represented in the tree at any given age, also known as ``lineages-through-time'' (LTT) curve. Note that \code{LTT} at \code{age0} will be equal to \code{LTT}, and that values in \code{LTT} will be non-increasing with age.}
\item{lambda}{Numeric vector of size NG, listing the speciation rate (in units 1/time) at the ages given in \code{ages[]}.}
\item{mu}{Numeric vector of size NG, listing the extinction rate (in units 1/time) at the ages given in \code{ages[]}.}
\item{diversification_rate}{Numeric vector of size NG, listing the net diversification rate (\eqn{\lambda-\mu}) at the ages given in \code{ages[]}.}
\item{PDR}{Numeric vector of size NG, listing the pulled diversification rate (PDR, in units 1/time) at the ages given in \code{ages[]}.}
\item{PND}{Numeric vector of size NG, listing the pulled normalized diversity (PND, in units 1/time) at the ages given in \code{ages[]}. The PND is defined as \eqn{PND=(N/N(\tau_o))\cdot\lambda(\tau_o)/\lambda}.}
\item{SER}{Numeric vector of size NG, listing the ``shadow extinction rate'' (SER, in units 1/time) at the ages given in \code{ages[]}. The SER is defined as \eqn{SER=\rho(\tau_o)\lambda(\tau_o)-PDR}.}
\item{PER}{Numeric vector of size NG, listing the ``pulled extinction rate'' (PER, in units 1/time) at the ages given in \code{ages[]}. The PER is defined as \eqn{SER=\lambda(\tau_o)-PDR} (Louca et al. 2018).}
\item{PSR}{Numeric vector of size NG, listing the ``pulled speciation rate'' (PSR, in units 1/time) at the ages given in \code{ages[]}. The PSR is defined as \eqn{PSR=\lambda\cdot(1-Pmissing)}.}
\item{rholambda0}{Non-negative numeric, specifying the product of the sampling fraction and the speciation rate at \code{age0}, \eqn{\rho\cdot\lambda(\tau_o)}.}
}


\author{Stilianos Louca}

\references{
H. Morlon, T. L. Parsons, J. B. Plotkin (2011). Reconciling molecular phylogenies with the fossil record. Proceedings of the National Academy of Sciences. 108:16327-16332.

S. Louca et al. (2018). Bacterial diversification through geological time. Nature Ecology & Evolution. 2:1458-1467.
}

\seealso{
\code{\link{loglikelihood_hbd}}
}

\examples{
# define an HBD model with exponentially decreasing speciation/extinction rates
Ntips     = 1000
beta      = 0.01 # exponential decay rate of lambda over time
oldest_age= 10
age_grid  = seq(from=0,to=oldest_age,by=0.1) # choose a sufficiently fine age grid
lambda    = 1*exp(beta*age_grid) # define lambda on the age grid
mu        = 0.2*lambda # assume similarly shaped but smaller mu

# simulate deterministic HBD model
simulation = simulate_deterministic_hbd(LTT0       = Ntips, 
                                        oldest_age = oldest_age, 
                                        rho0       = 0.5,
                                        age_grid   = age_grid,
                                        lambda     = lambda,
                                        mu         = mu)

# plot deterministic LTT
plot( x = simulation$ages, y = simulation$LTT, type='l',
      main='dLTT', xlab='age', ylab='lineages')
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{simulation}
\keyword{birth-death model}
