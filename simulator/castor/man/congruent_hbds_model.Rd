\name{congruent_hbds_model}
\alias{congruent_hbds_model}
\title{
Generate a congruent homogenous-birth-death-sampling model.
}
\description{
Given a homogenous birth-death-sampling (HBDS) model (or abstract congruence class), obtain the congruent model (or member of the congruence class) with a specific speciation rate \eqn{\lambda}, or extinction rate \eqn{\mu}, or sampling rate \eqn{\psi}, or effective reproduction ratio \eqn{R_e} or removal rate \eqn{\mu+\psi} (aka. "become uninfectious"" rate). All input and output time-profiles are specified as piecewise polynomial curves (splines), defined on a discrete grid of ages. This function allows exploration of a model's congruence class, by obtaining various members of the congruence class depending on the specified \eqn{\lambda}, \eqn{\mu}, \eqn{\psi}, \eqn{R_e} or removal rate. For more details on HBDS models and congruence classes see the documentation of \code{\link{simulate_deterministic_hbds}}.
}
\usage{
congruent_hbds_model(age_grid,
                     PSR,
                     PDR,
                     lambda_psi,
                     lambda             = NULL,
                     mu                 = NULL,
                     psi                = NULL,
                     Reff               = NULL,
                     removal_rate       = NULL,
                     lambda0            = NULL,
                     CSA_ages           = NULL,
                     CSA_pulled_probs   = NULL,
                     CSA_PSRs           = NULL,
                     splines_degree     = 1,
                     ODE_relative_dt    = 0.001,
                     ODE_relative_dy    = 1e-4)
}
\arguments{
\item{age_grid}{
Numeric vector, listing discrete ages (time before present) on which the various model variables (e.g., \eqn{\lambda}, \eqn{\mu} etc) are specified. Listed ages must be strictly increasing, and must include the present-day (i.e. age 0).
}
\item{PSR}{
Numeric vector, of the same size as \code{age_grid}, specifying the pulled speciation rate (PSR) (in units 1/time) at the ages listed in \code{age_grid}. The PSR is assumed to vary polynomially between grid points (see argument \code{splines_degree}). Can also be a single number, in which case \code{PSR} is assumed to be time-independent.
}
\item{PDR}{
Numeric vector, of the same size as \code{age_grid}, specifying the pulled diversification rate (PDR) (in units 1/time) at the ages listed in \code{age_grid}. The PDR is assumed to vary polynomially between grid points (see argument \code{splines_degree}). The PDR of a HBDS model is defined as \eqn{PDR=\lambda-\mu-\psi+(1/\lambda)d\lambda/dt} (where \eqn{t} denotes age). Can also be a single number, in which case \code{PDR} is assumed to be time-independent.
}
\item{lambda_psi}{
Numeric vector, of the same size as \code{age_grid}, specifying the product of speciation rate and sampling rate (\eqn{\lambda\psi}, in units 1/time^2) at the ages listed in \code{age_grid}. \eqn{\lambda\psi} is assumed to vary polynomially between grid points (see argument \code{splines_degree}). Can also be a single number, in which case \eqn{\lambda\psi} is assumed to be time-independent.
}
\item{lambda}{
Numeric vector, of the same size as \code{age_grid}, specifying the speciation rate (\eqn{\lambda}, in units 1/time) at the ages listed in \code{age_grid}. The speciation rate is assumed to vary polynomially between grid points (see argument \code{splines_degree}). Can also be a single number, in which case \eqn{\lambda} is assumed to be time-independent. By providing \eqn{\lambda}, one can select a specific model from the congruence class. Note that exactly one of \code{lambda}, \code{mu}, \code{psi}, \code{Reff} or \code{removal_rate} must be provided.
}
\item{mu}{
Numeric vector, of the same size as \code{age_grid}, specifying the extinction rate (\eqn{\mu}, in units 1/time) at the ages listed in \code{age_grid}. The extinction rate is assumed to vary polynomially between grid points (see argument \code{splines_degree}). Can also be a single number, in which case \eqn{\mu} is assumed to be time-independent. In an epidemiological context, \eqn{\mu} typically corresponds to the recovery rate plus the death rate of infected individuals. By providing \eqn{\mu} (together with \code{lambda0}, see below), one can select a specific model from the congruence class. Note that exactly one of \code{lambda}, \code{mu}, \code{psi}, \code{Reff} or \code{removal_rate} must be provided.
}
\item{psi}{
Numeric vector, of the same size as \code{age_grid}, specifying the (Poissonian) sampling rate (\eqn{\psi}, in units 1/time) at the ages listed in \code{age_grid}. The sampling rate is assumed to vary polynomially between grid points (see argument \code{splines_degree}). Can also be a single number, in which case \eqn{\psi} is assumed to be time-independent. By providing \eqn{\psi}, one can select a specific model from the congruence class. Note that exactly one of \code{lambda}, \code{mu}, \code{psi}, \code{Reff} or \code{removal_rate} must be provided.
}
\item{Reff}{
Numeric vector, of the same size as \code{age_grid}, specifying the effective reproduction ratio (\eqn{R_e}, unitless) at the ages listed in \code{age_grid}. The \eqn{R_e} is assumed to vary polynomially between grid points (see argument \code{splines_degree}). Can also be a single number, in which case \eqn{R_e} is assumed to be time-independent. By providing \eqn{R_e} (together with \code{lambda0}, see below), one can select a specific model from the congruence class. Note that exactly one of \code{lambda}, \code{mu}, \code{psi}, \code{Reff} or \code{removal_rate} must be provided.
}
\item{removal_rate}{
Numeric vector, of the same size as \code{age_grid}, specifying the removal rate (\eqn{\mu+\psi}, in units 1/time) at the ages listed in \code{age_grid}. IN an epidemiological context this is also known as "become uninfectious" rate. The removal rate is assumed to vary polynomially between grid points (see argument \code{splines_degree}). Can also be a single number, in which case the removal rate is assumed to be time-independent. By providing \eqn{\mu+\psi} (together with \code{lambda0}, see below), one can select a specific model from the congruence class. Note that exactly one of \code{lambda}, \code{mu}, \code{psi}, \code{Reff} or \code{removal_rate} must be provided.
}
\item{lambda0}{
Numeric, specifying the speciation rate at the present-day (i.e., at age 0). Must be provided if and only if one of \code{mu}, \code{Reff} or \code{removal_rate} is provided.
}
\item{CSA_ages}{
Optional numeric vector, listing the ages of concentrated sampling attempts, in ascending order. Concentrated sampling is assumed to occur in addition to any continuous (Poissonian) sampling specified by \code{psi}.
}
\item{CSA_pulled_probs}{
Optional numeric vector of the same size as \code{CSA_ages}, listing pulled sampling probabilities at each concentrated sampling attempt (CSA). Note that in contrast to the sampling rates \code{psi}, the \code{CSA_pulled_probs} are interpreted as probabilities and must thus be between 0 and 1. \code{CSA_pulled_probs} must be provided if and only if \code{CSA_ages} is provided.
}
\item{CSA_PSRs}{
Optional numeric vector of the same size as \code{CSA_ages}, specifying the pulled sampling rate (PSR) during each concentrated sampling attempt. While in principle the PSR is already provided by the argument \code{PSR}, the PSR may be non-continuous at CSAs, which makes a representation as piecewise polynomial function difficult; hence, you must explicitly provide the correct PSR at each CSA. \code{CSA_PSRs} must be provided if and only if \code{CSA_ages} is provided.
}
\item{splines_degree}{
Integer, either 0,1,2 or 3, specifying the polynomial degree of the provided time-dependent variables between grid points in \code{age_grid}. For example, if \code{splines_degree==1}, then the provided PDR, PSR and so on are interpreted as piecewise-linear curves; if \code{splines_degree==2} they are interpreted as quadratic splines; if \code{splines_degree==3} they are interpreted as cubic splines. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on.
}
\item{ODE_relative_dt}{
Positive unitless number, specifying the default relative time step for internally used ordinary differential equation solvers. Typical values are 0.01-0.001.
}
\item{ODE_relative_dy}{
Positive unitless number, specifying the relative difference between subsequent simulated and interpolated values, in internally used ODE solvers. Typical values are 1e-2 to 1e-5. A smaller \code{ODE_relative_dy} increases interpolation accuracy, but also increases memory requirements and adds runtime.
}
}
\details{
The PDR, PSR and the product \eqn{\lambda\psi} are variables that are invariant across the entire congruence class of an HBDS model, i.e. any two congruent models have the same PSR, PDR and product \eqn{\lambda\psi}. Reciprocally, any HBDS congruence class is fully determined by its PDR, PSR and \eqn{\lambda\psi}. This function thus allows "collapsing" a congruence class down to a single member (a specific HBDS model) by specifying one or more additional variables over time (such as \eqn{\lambda}, or \eqn{\psi}, or \eqn{\mu} and \eqn{\lambda_0}). Alternatively, this function may be used to obtain alternative models that are congruent to some reference model, for example to explore the robustness of specific epidemiological quantities of interest. The function returns a specific HBDS model in terms of the time profiles of various variables (such as \eqn{\lambda}, \eqn{\mu} and \eqn{\psi}).

In the current implementation it is assumed that any sampled lineages are immediately removed from the pool, that is, this function cannot accommodate models with a non-zero retention probability upon sampling. This is a common assumption in molecular epidemiology.
Note that in this function age always refers to time before present, i.e., present day age is 0, and age increases towards the root.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the calculation was successful. If \code{FALSE}, then the returned list includes an additional `\code{error}' element (character) providing a description of the error; all other return variables may be undefined.
}
\item{valid}{
Logical, indicating whether the returned model appears to be biologically valid (for example, does not have negative \eqn{\lambda}, \eqn{\mu} or \eqn{\psi}). In principle, a congruence class may include biologically invalid models, which might be returned depending on the input to \code{congruent_hbds_model}. Note that only biologically valid models can be subsequently simulated using \code{\link{simulate_deterministic_hbds}}.
}
\item{ages}{Numeric vector of size NG, specifying the discrete ages (time before present) on which all returned time-curves are specified. Will always be equal to \code{age_grid}.
}
\item{lambda}{Numeric vector of size NG, listing the speciation rates \eqn{\lambda} of the returned model at the ages given in \code{ages[]}.
}
\item{mu}{Numeric vector of size NG, listing the extinction rates \eqn{\mu} of the returned model at the ages given in \code{ages[]}.
}
\item{psi}{Numeric vector of size NG, listing the (Poissonian) sampling rates \eqn{\psi} of the returned model at the ages given in \code{ages[]}.
}
\item{lambda_psi}{Numeric vector of size NG, listing the product \eqn{\lambda\psi} at the ages given in \code{ages[]}.
}
\item{Reff}{Numeric vector of size NG, listing the effective reproduction ratio \eqn{R_e} of the returned model at the ages given in \code{ages[]}.
}
\item{removal_rate}{Numeric vector of size NG, listing the removal rate (\eqn{\mu+\psi}, aka. "become uninfectious" rate) of the returned model at the ages given in \code{ages[]}.
}
\item{Pmissing}{Numeric vector of size NG, listing the probability that a past lineage extant during \code{ages[]} will be missing from a tree generated by the model.
}
\item{CSA_probs}{Numeric vector of the same size as \code{CSA_ages}, listing the sampling probabilities at each of the CSAs.}
\item{CSA_Pmissings}{Numeric vector of the same size as \code{CSA_ages}, listing the probability that a past lineage extant during each of \code{CSA_ages[]} will be missing from a tree generated by the model.}
}


\author{Stilianos Louca}

\references{
T. Stadler, D. Kuehnert, S. Bonhoeffer, A. J. Drummond (2013). Birth-death skyline plot reveals temporal changes of epidemic spread in HIV and hepatitis C virus (HCV). PNAS. 110:228-233.

A. MacPherson, S. Louca, A. McLaughlin, J. B. Joy, M. W. Pennell (in review as of 2020). A general birth-death-sampling model for epidemiology and macroevolution. DOI:10.1101/2020.10.10.334383
}

\seealso{
\code{\link{generate_tree_hbds}},
\code{\link{fit_hbds_model_parametric}},
\code{\link{simulate_deterministic_hbds}}
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
sim = simulate_deterministic_hbds(age_grid  = age_grid,
                                  lambda    = lambda,
                                  mu        = mu,
                                  psi       = 0.1,
                                  age0      = 1,
                                  LTT0      = 100)
                                         
# calculate a congruent HBDS model with an alternative sampling rate
# use the previously simulated variables to define the congruence class
new_psi = 0.1*exp(-0.01*sim$ages) # consider a psi decreasing with age
congruent = congruent_hbds_model(age_grid   = sim$ages,
                                 PSR        = sim$PSR,
                                 PDR        = sim$PDR,
                                 lambda_psi = sim$lambda_psi,
                                 psi        = new_psi)
											 
# compare the deterministic LTT of the two models
# to confirm that the models are indeed congruent
if(!congruent$valid){
    cat("WARNING: Congruent model is not biologically valid\\n")
}else{
    # simulate the congruent model to get the LTT
    csim = simulate_deterministic_hbds(age_grid = congruent$ages,
                                       lambda   = congruent$lambda,
                                       mu       = congruent$mu,
                                       psi      = congruent$psi,
                                       age0     = 1,
                                       LTT0     = 100)

    # plot deterministic LTT of the original and congruent model
    plot( x = sim$ages, y = sim$LTT, type='l',
          main='dLTT', xlab='age', ylab='lineages', 
          xlim=c(oldest_age,0), col='red')
    lines(x= csim$ages, y=csim$LTT, 
          type='p', pch=21, col='blue')
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{simulation}
\keyword{birth-death model}
