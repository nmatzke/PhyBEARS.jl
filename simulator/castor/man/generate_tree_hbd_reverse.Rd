\name{generate_tree_hbd_reverse}
\alias{generate_tree_hbd_reverse}
\title{
Generate a tree from a birth-death model in reverse time.
}
\description{
Generate an ultrametric timetree (comprising only extant lineages) in reverse time (from present back to the root) based on the homogenous birth-death (HBD; Morlon et al., 2011) model, conditional on a specific number of extant species sampled and (optionally) conditional on the crown age or stem age. 

The probability distribution of such trees only depends on the congruence class of birth-death models (e.g., as specified by the pulled speciation rate) but not on the precise model within a congruence class (Louca and Pennell, 2019). Hence, in addition to allowing specification of speciation and extinction rates, this function can alternatively simulate trees simply based on some pulled speciation rate (PSR), or based on some pulled diversification rate (PDR) and the product \eqn{\rho\lambda_o} (present-day sampling fraction times present-day speciation rate). 

This function can be used to generate bootstrap samples after fitting an HBD model or HBD congruence class to a real timetree.
}
\usage{
generate_tree_hbd_reverse( Ntips,
                           stem_age         = NULL,
                           crown_age        = NULL,
                           age_grid         = NULL,
                           lambda           = NULL,
                           mu               = NULL,
                           rho              = NULL,
                           PSR              = NULL,
                           PDR              = NULL,
                           rholambda0       = NULL,
                           force_max_age    = Inf,
                           splines_degree   = 1,
                           relative_dt      = 1e-3,
                           Ntrees           = 1,
                           tip_basename     = "",
                           node_basename    = NULL,
                           edge_basename    = NULL)
}
\arguments{
\item{Ntips}{
Number of tips in the tree, i.e. number of extant species sampled at present day.
}
\item{stem_age}{
Numeric, optional stem age on which to condition the tree. If NULL or <=0, the tree is not conditioned on the stem age.
}
\item{crown_age}{
Numeric, optional crown age (aka. root age or MRCA age) on which to condition the tree. If NULL or <=0, the tree is not conditioned on the crown age. If both \code{stem_age} and \code{crown_age} are specified, only the crown age is used; in that case for consistency \code{crown_age} must not be greater than \code{stem_age}.
}
\item{age_grid}{
Numeric vector, listing discrete ages (time before present) on which the PSR is specified. Listed ages must be strictly increasing, and should cover at least the present day (age 0) as well as a sufficient duration into the past. If conditioning on the stem or crown age, that age must also be covered by \code{age_grid}. When not conditioning on crown nor stem age, and the generated tree ends up extending beyond the last time point in \code{age_grid}, the PSR will be extrapolated as a constant (with value equal to the last value in \code{PSR}) as necessary. \code{age_grid} also be \code{NULL} or a vector of size 1, in which case the PSR is assumed to be time-independent.
}
\item{lambda}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing speciation rates (\eqn{\lambda}, in units 1/time) at the ages listed in \code{age_grid}. Speciation rates must be non-negative, and are assumed to vary as a spline between grid points (see argument \code{splines_degree}). Can also be \code{NULL}, in which case either \code{PSR}, or \code{PDR} and \code{rholambda0}, must be provided.
}
\item{mu}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing extinction rates (\eqn{\mu}, in units 1/time) at the ages listed in \code{age_grid}. Extinction rates must be non-negative, and are assumed to vary as a spline between grid points (see argument \code{splines_degree}). Can also be \code{NULL}, in which case either \code{PSR}, or \code{PDR} and \code{rholambda0}, must be provided.
}
\item{rho}{
Numeric, sampling fraction at present day (fraction of extant species included in the tree). Can also be \code{NULL}, in which case either \code{PSR}, or \code{PDR} and \code{rholambda0}, must be provided.
}
\item{PSR}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing pulled speciation rates (\eqn{\lambda_p}, in units 1/time) at the ages listed in \code{age_grid}. The PSR must be non-negative (and strictly positive almost everywhere), and is assumed to vary as a spline between grid points (see argument \code{splines_degree}). Can also be \code{NULL}, in which case either \code{lambda} and \code{mu} and \code{rho}, or \code{PDR} and \code{rholambda0}, must be provided.
}
\item{PDR}{
Numeric vector, of the same size as \code{age_grid} (or size 1 if \code{age_grid==NULL}), listing pulled diversification rates (\eqn{r_p}, in units 1/time) at the ages listed in \code{age_grid}. The PDR is assumed to vary polynomially between grid points (see argument \code{splines_degree}). Can also be \code{NULL}, in which case either \code{lambda} and \code{mu} and \code{rho}, or \code{PSR}, must be provided.
}
\item{rholambda0}{
Strictly positive numeric, specifying the product \eqn{\rho\lambda_o} (present-day species sampling fraction times present-day speciation rate). Can also be \code{NULL}, in which case \code{PSR} must be provided.
}
\item{force_max_age}{
Numeric, specifying an optional maximum allowed age for the tree's root. If the tree ends up expanding past that age, all remaining lineages are forced to coalesce at that age. This is not statistically consistent with the provided HBD model (in fact it corresponds to a modified HBD model with a spike in the PSR at that time). This argument merely provides a way to prevent excessively large trees if the PSR is close to zero at older ages and when not conditioning on the stem nor crown age, while still keeping the original statistical properties at younger ages. To disable this feature set \code{force_max_age} to \code{Inf}.
}
\item{splines_degree}{
Integer, either 0,1,2 or 3, specifying the polynomial degree of the provided rates \code{PSR}, \code{PDR}, \code{lambda}, \code{mu} and \code{rho} between grid points in \code{age_grid}. For example, if \code{splines_degree==1}, then the provided rates are interpreted as piecewise-linear curves; if \code{splines_degree==2} the rates are interpreted as quadratic splines; if \code{splines_degree==3} the rates are interpreted as cubic splines. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on. If your \code{age_grid} is fine enough, then \code{splines_degree=1} is usually sufficient.
}
\item{relative_dt}{
Strictly positive numeric (unitless), specifying the maximum relative time step allowed for integration over time. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.
}
\item{Ntrees}{
Integer, number of trees to generate. The computation time per tree is lower if you generate multiple trees at once.
}
\item{tip_basename}{
Character. Prefix to be used for tip labels (e.g. "tip."). If empty (""), then tip labels will be integers "1", "2" and so on.
}
\item{node_basename}{
Character. Prefix to be used for node labels (e.g. "node."). If \code{NULL}, no node labels will be included in the tree.
}
\item{edge_basename}{
Character. Prefix to be used for edge labels (e.g. "edge."). Edge labels (if included) are stored in the character vector \code{edge.label}. If \code{NULL}, no edge labels will be included in the tree.
}
}
\details{
This function requires that the BD model, or the BD congruence class (Louca and Pennell, 2019), is specified using one of the following sets of arguments:
\itemize{
\item Using the speciation rate \eqn{\lambda}, the extinctin rate \eqn{\mu}, and the present-day sampling fraction \eqn{\rho}.
\item Using the pulled diversification rate (PDR) and the product \eqn{\rho\lambda(0)}. The PDR is defined as \eqn{r_p=\lambda-\mu+\frac{1}{\lambda}\frac{d\lambda}{d\tau}}, where \eqn{\tau} is age (time before present), \eqn{\lambda(\tau)} is the speciation rate at age \eqn{\tau} and \eqn{\mu(\tau)} is the extinction rate.
\item Using the pulled speciation rate (PSR). The PSR (\eqn{\lambda_p}) is defined as \eqn{\lambda_p(\tau) = \lambda(\tau)\cdot\Phi(\tau)}, where and \eqn{\Phi(\tau)} is the probability that a lineage extant at age \eqn{\tau} will survive until the present and be represented in the tree. 
}
Concurrently using/combining more than one the above parameterization methods is not supported.

Either the PSR, or the PDR and \code{rholambda0}, provide sufficient information to fully describe the probability distribution of the tree (Louca and Pennell, 2019).
For example, the probability distribution of generated trees only depends on the PSR, and not on the specific speciation rate \eqn{\lambda} or extinction rate \eqn{\mu} (various combinations of \eqn{\lambda} and \eqn{\mu} can yield the same PSR; Louca and Pennell, 2019). To calculate the PSR and PDR for any arbitrary \eqn{\lambda}, \eqn{\mu} and \eqn{\rho} you can use the function \code{\link{simulate_deterministic_hbd}}.

When not conditioning on the crown age, the age of the root of the generated tree will be stochastic (i.e., non-fixed). This function then assumes a uniform prior distribution (in a sufficiently large time interval) for the origin of the forward HBD process that would have generated the tree, based on a generalization of the EBDP algorithm provided by (Stadler, 2011). When conditioning on stem or crown age, this function is based on the algorithm proposed by Hoehna (2013, Eq. 8).

Note that HBD trees can also be generated using the function \code{\link{generate_random_tree}}. That function, however, generates trees in forward time, and hence when conditioning on the final number of tips the total duration of the simulation is unpredictable; consequently, speciation and extinction rates cannot be specified as functions of "age" (time before present). The function presented here provides a means to generate trees with a fixed number of tips, while specifying \eqn{\lambda}, \eqn{\mu}, \eqn{\lambda_p} or \eqn{r_p} as functions of age (time before present).
}

\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the simulation was successful. If \code{FALSE}, then the returned list includes an additional `\code{error}' element (character) providing a description of the error; all other return variables may be undefined.
}
\item{trees}{A list of length \code{Ntrees}, listing the generated trees. Each tree will be an ultrametric timetree of class "phylo".}
}


\author{Stilianos Louca}

\references{
H. Morlon, T. L. Parsons, J. B. Plotkin (2011). Reconciling molecular phylogenies with the fossil record. Proceedings of the National Academy of Sciences. 108:16327-16332.

T. Stadler (2011). Simulating trees with a fixed number of extant species. Systematic Biology. 60:676-684.

S. Hoehna (2013). Fast simulation of reconstructed phylogenies under global time-dependent birth-death processes. Bioinformatics. 29:1367-1374.

S. Louca and M. W. Pennell (in review as of 2019). Phylogenies of extant species are consistent with an infinite array of diversification histories.
}

\seealso{
\code{\link{loglikelihood_hbd}},
\code{\link{simulate_deterministic_hbd}},
\code{\link{generate_random_tree}}
}

\examples{
# EXAMPLE 1: Generate trees based on some speciation and extinction rate
# In this example we assume an exponentially decreasing speciation rate
#   and a temporary mass extinction event

# define parameters
age_grid = seq(0,100,length.out=1000)
lambda   = 0.1 + exp(-0.5*age_grid)
mu       = 0.05 + exp(-(age_grid-5)^2)
rho      = 0.5 # species sampling fraction at present-day

# generate a tree with 100 tips and no specific crown or stem age
sim = generate_tree_hbd_reverse(Ntips       = 100, 
                                age_grid    = age_grid, 
                                lambda      = lambda,
                                mu          = mu,
                                rho         = rho)
if(!sim$success){
    cat(sprintf("Tree generation failed: \%s\n",sim$error))
}else{
    cat(sprintf("Tree generation succeeded\n"))
    tree = sim$trees[[1]]
}


########################
# EXAMPLE 2: Generate trees based on the pulled speciation rate
# Here we condition the tree on some fixed crown (MRCA) age

# specify the PSR on a sufficiently fine and wide age grid
age_grid  = seq(0,1000,length.out=10000)
PSR       = 0.1+exp(-0.1*age_grid) # exponentially decreasing PSR

# generate a tree with 100 tips and MRCA age 10
sim = generate_tree_hbd_reverse(Ntips       = 100, 
                                age_grid    = age_grid, 
                                PSR         = PSR, 
                                crown_age   = 10)
if(!sim$success){
    cat(sprintf("Tree generation failed: \%s\n",sim$error))
}else{
    cat(sprintf("Tree generation succeeded\n"))
    tree = sim$trees[[1]]
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{simulation}
\keyword{birth-death model}
