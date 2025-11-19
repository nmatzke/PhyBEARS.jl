\name{simulate_deterministic_hbds}
\alias{simulate_deterministic_hbds}
\title{
Simulate a deterministic homogenous birth-death-sampling model.
}
\description{
Given a homogenous birth-death-sampling (HBDS) model, i.e., with speciation rate \eqn{\lambda}, extinction rate \eqn{\mu}, continuous (Poissonian) sampling rate \eqn{\psi} and retention probability \eqn{\kappa}, calculate various deterministic features of the model backwards in time, such as the total population size and the LTT over time.
Continuously sampled lineages are kept in the pool of extant lineages at probability \eqn{\kappa}. The variables \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} may depend on time.
In addition, the model can include concentrated sampling attempts at a finite set of discrete time points \eqn{t_1,..,t_m}. ``Homogenous'' refers to the assumption that, at any given moment in time, all lineages exhibit the same speciation/extinction/sampling rates and retention proabbility. Every HBDS model is thus defined based on the values that \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} take over time, as well as the sampling probabilities \eqn{\psi_1,..,\psi_m} and retention probabilities \eqn{\kappa_1,..,\kappa_m} during the concentrated sampling attempts. Special cases of this model are sometimes known as ``birth-death-skyline plots'' in the literature (Stadler 2013). In epidemiology, these models are often used to describe the phylogenies of viral strains sampled over the course of the epidemic. A ``concentrated sampling attempt'' is a brief but intensified sampling period that lasted much less than the typical timescales of speciation/extinction. ``Deterministic'' refers to the fact that all calculated properties are completely determined by the model's parameters (i.e. non-random), as if an infinitely large tree was generated (aka. ``continuum limit'').
The time-profiles of \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} are specified as piecewise polynomial curves (splines), defined on a discrete grid of ages.
}
\usage{
simulate_deterministic_hbds(age_grid                        = NULL,
                            lambda                          = NULL,
                            mu                              = NULL,
                            psi                             = NULL,
                            kappa                           = NULL,
                            splines_degree                  = 1,
                            CSA_ages                        = NULL,
                            CSA_probs                       = NULL,
                            CSA_kappas                      = NULL,
                            requested_ages                  = NULL,
                            age0                            = 0,
                            N0                              = NULL,
                            LTT0                            = NULL,
                            ODE_relative_dt                 = 0.001,
                            ODE_relative_dy                 = 1e-4)
}
\arguments{
\item{age_grid}{
Numeric vector, listing discrete ages (time before present) on which either \eqn{\lambda} and \eqn{\mu}, or the PDR and \eqn{\mu}, are specified. Listed ages must be strictly increasing, and must cover at least the full considered age interval (from \code{age0} to \code{oldest_age}). Can also be \code{NULL} or a vector of size 1, in which case the speciation rate, extinction rate and PDR are assumed to be time-independent.
}
\item{lambda}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing speciation rates (\eqn{\lambda}, in units 1/time) at the ages listed in \code{age_grid}. Speciation rates should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}).
}
\item{mu}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing extinction rates (\eqn{\mu}, in units 1/time) at the ages listed in \code{age_grid}. Extinction rates should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}).
}
\item{psi}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing the continuous (Poissonian) sampling rate at the ages listed in \code{age_grid}. Sampling rates should be non-negative, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}).
}
\item{kappa}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing the retention probabilities following Poissonian sampling events, at the ages listed in \code{age_grid}. The listed values must be true probabilities, i.e. between 0 and 1, and are assumed to vary polynomially between grid points (see argument \code{splines_degree}). The retention probability is the probability that a continuously sampled lineage remains in the pool of extant lineages. Note that many epidemiological models assume kappa to be zero.
}
\item{splines_degree}{
Integer, either 0,1,2 or 3, specifying the polynomial degree of the provided \code{lambda}, \code{mu}, \code{psi} and \code{kappa} between grid points in \code{age_grid}. For example, if \code{splines_degree==1}, then the provided \code{lambda}, \code{mu}, \code{psi} and \code{kappa} are interpreted as piecewise-linear curves; if \code{splines_degree==2} they are interpreted as quadratic splines; if \code{splines_degree==3} they are interpreted as cubic splines. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on.
}
\item{CSA_ages}{
Optional numeric vector, listing the ages of concentrated sampling attempts, in ascending order. Concentrated sampling is performed in addition to any continuous (Poissonian) sampling specified by \code{psi}.
}
\item{CSA_probs}{
Optional numeric vector of the same size as \code{CSA_ages}, listing sampling probabilities at each concentrated sampling attempt. Note that in contrast to the sampling rates \code{psi}, the \code{CSA_probs} are interpreted as probabilities and must thus be between 0 and 1. \code{CSA_probs} must be provided if and only if \code{CSA_ages} is provided.
}
\item{CSA_kappas}{
Optional numeric vector of the same size as \code{CSA_ages}, listing retention probabilities at each concentrated sampling event, i.e. the probability at which a sampled lineage is kept in the pool of extant lineages. Note that the \code{CSA_kappas} are probabilities and must thus be between 0 and 1. \code{CSA_kappas} must be provided if and only if \code{CSA_ages} is provided.
}
\item{requested_ages}{
Optional numeric vector, listing ages (in ascending order) at which the various model variables are requested. If \code{NULL}, it will be set to \code{age_grid}.
}
\item{age0}{
Non-negative numeric, specifying the age at which \code{LTT0} and \code{pop_size0} are specified. Typically this will be 0, i.e., corresponding to the present.
}
\item{N0}{
Positive numeric, specifying the number of extant species (sampled or not) at \code{age0}. Used to determine the "scaling factor" for the returned population sizes and LTT. Either \code{pop_size0} or \code{LTT0} must be provided, but not both.
}
\item{LTT0}{
Positive numeric, specifying the number of lineages present in the tree at \code{age0}. Used to determine the "scaling factor" for the returned population sizes and LTT. Either \code{pop_size0} or \code{LTT0} must be provided, but not both.
}
\item{ODE_relative_dt}{
Positive unitless number, specifying the default relative time step for internally used ordinary differential equation solvers. Typical values are 0.01-0.001.
}
\item{ODE_relative_dy}{
Positive unitless number, specifying the relative difference between subsequent simulated and interpolated values, in internally used ODE solvers. Typical values are 1e-2 to 1e-5. A smaller \code{ODE_relative_dy} increases interpolation accuracy, but also increases memory requirements and adds runtime (scaling with the tree's age span, not with Ntips).
}
}
\details{
The simulated LTT refers to a hypothetical tree sampled at age 0, i.e. LTT(t) will be the number of lineages extant at age t that survived and were sampled until by the present day. Note that if a concentrated sampling attempt occurs at age \eqn{\tau}, then LTT(\eqn{\tau}) is the number of lineages in the tree right before the occurrence of the sampling attempt, i.e., in the limit where \eqn{\tau} is approached from above.

Note that in this function age always refers to time before present, i.e., present day age is 0, and age increases towards the root.
}

\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the calculation was successful. If \code{FALSE}, then the returned list includes an additional `\code{error}' element (character) providing a description of the error; all other return variables may be undefined.
}
\item{ages}{Numerical vector of size NG, listing discrete ages (time before present) on which all returned time-curves are specified. Will be equal to \code{requested_ages}, if the latter was provided.
}
\item{total_diversity}{Numerical vector of size NG, listing the predicted (deterministic) total diversity (number of extant species, denoted \eqn{N}) at the ages given in \code{ages[]}.
}
\item{LTT}{Numeric vector of size NG, listing the number of lineages represented in the tree at any given age, also known as ``lineages-through-time'' (LTT) curve. Note that \code{LTT} at \code{age0} will be equal to \code{LTT0} (if the latter was provided).
}
\item{nLTT}{Numeric vector of size NG, listing values of the normalized LTT at ages \code{ages[]}. The nLTT is calculated by dividing the LTT by its area-under-the-curve (AUC). The AUC is calculated by integrating the LTT over the time spanned by \code{ages} and using the trapezoid rule. Hence, the exact value of the AUC and of the nLTT depends on the span and resolution of \code{ages[]}. If you want the AUC to accurately refer to the entire area under the curve (i.e. across the full real axis), you should specify a sufficiently wide and sufficiently fine age grid (e.g., via \code{requested_ages}).
}
\item{Pmissing}{Numeric vector of size NG, listing the probability that a lineage, extant at a given age, will not be represented in the tree.}
\item{lambda}{
Numeric vector of size NG, listing the speciation rates at the ages given in \code{ages[]}.
}
\item{mu}{
Numeric vector of size NG, listing the extinctions rates at the ages given in \code{ages[]}.
}
\item{psi}{
Numeric vector of size NG, listing the Poissonian sampling rates at the ages given in \code{ages[]}.
}
\item{kappa}{
Numeric vector of size NG, listing the retention probabilities (for continuously sampled lineages) at the ages given in \code{ages[]}.
}
\item{PDR}{Numeric vector of size NG, listing the pulled diversification rate (PDR, in units 1/time) at the ages given in \code{ages[]}.}
\item{IPRP}{Numeric vector of size NG, listing the age-integrated pulled diversification rate at the ages given in \code{ages[]}, i.e. \eqn{IPDR(t)=\int_0^t PDR(s)ds}.}
\item{PSR}{Numeric vector of size NG, listing the ``pulled speciation rate'' (PSR, in units 1/time) at the ages given in \code{ages[]}. The PSR is defined as \eqn{PSR=\lambda\cdot(1-Pmissing)}.}
\item{PRP}{Numeric vector of size NG, listing the ``pulled retention probability'' (PRP) at the ages given in \code{ages[]}. The PRP is defined as \eqn{PRP=\kappa\cdot(1-Pmissing)}.}
\item{diversification_rate}{
Numeric vector of size NG, listing the net diversification rate (in units 1/time) at the ages given in \code{ages[]}.
}
\item{branching_density}{
Numeric vector of size NG, listing the deterministic branching density (PSR * nLTT, in units nodes/time) at the ages given in \code{ages[]}.
}
\item{sampling_density}{
Numeric vector of size NG, listing the deterministic sampling density (\eqn{\psi\cdot N/AUC}, in units tips/time, where AUC is the area-under-the-curve calculated for the LTT) at the ages given in \code{ages[]}.
}
\item{lambda_psi}{
Numeric vector of size NG, listing the product of the speciation rate and Poissonian sampling rate (in units 1/time^2) at the ages given in \code{ages[]}.
}
\item{kappa_psi}{
Numeric vector of size NG, listing the product of the continuous sampling rate and the continuous retention probability (in units 1/time) at the ages given in \code{ages[]}.
}
\item{Reff}{
Numeric vector of size NG, listing the effective reproduction ratio (\eqn{R_e=\lambda/(\mu+\psi(1-\kappa))}) at the ages given in \code{ages[]}.
}
\item{removal_rate}{
Numeric vector of size NG, listing the total removal rate (\eqn{\mu+\psi}), also known as ``become uninfectious rate'', at the ages given in \code{ages[]}.
}
\item{sampling_proportion}{
Numeric vector of size NG, listing the instantaneous sampling proportion (\eqn{\psi/(\mu+\psi)}) at the ages given in \code{ages[]}.
}
\item{CSA_pulled_probs}{
Numeric vector of size NG, listing the pulled concentrated sampling probabilities, \eqn{\tilde{\rho}_k=\rho_k/(1-E)}.
}
\item{CSA_psis}{
Numeric vector of size NG, listing the continuous (Poissonian) sampling rates during the concentrated sampling attempts, \eqn{\psi(t_1),..,\psi(t_m)}.
}
\item{CSA_PSRs}{
Numeric vector of size NG, listing the pulled speciation rates during the concentrated sampling attempts.
}
}


\author{Stilianos Louca}

\references{
T. Stadler, D. Kuehnert, S. Bonhoeffer, A. J. Drummond (2013). Birth-death skyline plot reveals temporal changes of epidemic spread in HIV and hepatitis C virus (HCV). PNAS. 110:228-233.
}

\seealso{
\code{\link{generate_tree_hbds}},
\code{\link{fit_hbds_model_parametric}},
\code{\link{simulate_deterministic_hbd}}
}

\examples{
# define an HBDS model with exponentially decreasing speciation/extinction rates
# and constant Poissonian sampling rate psi
oldest_age= 10
age_grid  = seq(from=0,to=oldest_age,by=0.1) # choose a sufficiently fine age grid
lambda    = 1*exp(0.01*age_grid) # define lambda on the age grid
mu        = 0.2*lambda # assume similarly shaped but smaller mu

# simulate deterministic HBD model
# scale LTT such that it is 100 at age 1
simulation = simulate_deterministic_hbds(age_grid   = age_grid,
                                         lambda     = lambda,
                                         mu         = mu,
                                         psi        = 0.1,
                                         age0       = 1,
                                         LTT0       = 100)
# plot deterministic LTT
plot( x = simulation$ages, y = simulation$LTT, type='l',
      main='dLTT', xlab='age', ylab='lineages', xlim=c(oldest_age,0))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{simulation}
\keyword{birth-death model}
