\name{generate_tree_hbds}
\alias{generate_tree_hbds}
\title{
Generate a tree from a birth-death-sampling model in forward time.
}
\description{
Generate a random timetree according to a homogenous birth-death-sampling model with arbitrary time-varying speciation/extinction/sampling rates. Lineages split (speciate) or die (go extinct) at Poissonian rates and independently of each other. Lineages are sampled continuously (i.e., at Poissonian rates) in time and/or during concentrated sampling attempts (i.e., at specific time points). Sampled lineages are assumed to continue in the pool of extant lineages at some given "retention probability". The final tree can be restricted to sampled lineages only, but may optionally include extant (non-sampled) as well as extinct lineages. Speciation, extinction and sampling rates as well as retention probabilities may depend on time. This function may be used to simulate trees commonly encountered in viral epidemiology, where sampled patients are assumed to exit the pool of infectious individuals.
}
\usage{
generate_tree_hbds( max_sampled_tips             = NULL, 
                    max_sampled_nodes            = NULL, 
                    max_extant_tips              = NULL,
                    max_extinct_tips             = NULL,
                    max_tips                     = NULL,
                    max_time                     = NULL,
                    include_extant               = FALSE,
                    include_extinct              = FALSE,
                    as_generations               = FALSE,
                    time_grid                    = NULL,
                    lambda                       = NULL,
                    mu                           = NULL,
                    psi                          = NULL,
                    kappa                        = NULL,
                    splines_degree               = 1,
                    CSA_times                    = NULL,
                    CSA_probs                    = NULL,
                    CSA_kappas                   = NULL,
                    no_full_extinction           = FALSE,
                    max_runtime                  = NULL,
                    tip_basename                 = "",
                    node_basename                = NULL,
                    edge_basename                = NULL,
                    include_birth_times          = FALSE,
                    include_death_times          = FALSE)
}
\arguments{
\item{max_sampled_tips}{
Integer, maximum number of sampled tips. The simulation is halted once this number is reached. If \code{NULL} or <=0, this halting criterion is ignored.
}
\item{max_sampled_nodes}{
Integer, maximum number of sampled nodes, i.e., of lineages that were sampled but kept in the pool of extant lineages. The simulation is halted once this number is reached. If \code{NULL} or <=0, this halting criterion is ignored.
}
\item{max_extant_tips}{
Integer, maximum number of extant tips. The simulation is halted once the number of concurrently extant tips reaches this threshold. If \code{NULL} or <=0, this halting criterion is ignored.
}
\item{max_extinct_tips}{
Integer, maximum number of extant tips. The simulation is halted once this number is reached. If \code{NULL} or <=0, this halting criterion is ignored.
}
\item{max_tips}{
Integer, maximum number of tips (extant+extinct+sampled). The simulation is halted once this number is reached. If \code{NULL} or <=0, this halting criterion is ignored.
}
\item{max_time}{
Numeric, maximum duration of the simulation. If \code{NULL} or <=0, this halting criterion is ignored.
}
\item{include_extant}{
Logical, specifying whether to include extant tips (i.e., neither extinct nor sampled) in the final tree.
}
\item{include_extinct}{
Logical, specifying whether to include extant tips (i.e., neither extant nor sampled) in the final tree.
}
\item{as_generations}{
Logical, specifying whether edge lengths should correspond to generations. If \code{FALSE}, then edge lengths correspond to time. If \code{TRUE}, then the time between two subsequent events (speciation, extinction, sampling) is counted as "one generation".
}
\item{time_grid}{
Numeric vector, specifying time points (in ascending order) on which the rates \code{lambda}, \code{mu} and \code{psi} are provided. Rates are interpolated polynomially between time grid points as needed (according to splines_degree).
The time grid should generally cover the maximum possible simulation time, otherwise it will be polynomially extrapolated as needed.
}
\item{lambda}{
Numeric vector, of the same size as \code{time_grid} (or size 1 if \code{time_grid==NULL}), listing per-lineage speciation (birth) rates (\eqn{\lambda}, in units 1/time) at the times listed in \code{time_grid}. Speciation rates must be non-negative, and are assumed to vary as a spline between grid points (see argument \code{splines_degree}). Can also be a single numeric, in which case \eqn{\lambda} is assumed to be constant over time.
}
\item{mu}{
Numeric vector, of the same size as \code{time_grid} (or size 1 if \code{time_grid==NULL}), listing per-lineage extinction (death) rates (\eqn{\mu}, in units 1/time) at the times listed in \code{time_grid}. Extinction rates must be non-negative, and are assumed to vary as a spline between grid points (see argument \code{splines_degree}). Can also be a single numeric, in which case \eqn{\mu} is assumed to be constant over time. If omitted, the extinction rate is assumed to be zero.
}
\item{psi}{
Numeric vector, of the same size as \code{time_grid} (or size 1 if \code{time_grid==NULL}), listing per-lineage sampling rates (\eqn{\psi}, in units 1/time) at the times listed in \code{time_grid}. Sampling rates must be non-negative, and are assumed to vary as a spline between grid points (see argument \code{splines_degree}). Can also be a single numeric, in which case \eqn{\psi} is assumed to be constant over time. If omitted, the continuous sampling rate is assumed to be zero.
}
\item{kappa}{
Numeric vector, of the same size as \code{time_grid} (or size 1 if \code{time_grid==NULL}), listing retention probabilities (\eqn{\kappa}, unitless) of continuously (Poissonian) sampled lineages at the times listed in \code{time_grid}. Retention probabilities must be true probabilities (i.e., between 0 and 1), and are assumed to vary as a spline between grid points (see argument \code{splines_degree}). Can also be a single numeric, in which case \eqn{\kappa} is assumed to be constant over time. If omitted, the retention probability is assumed to be zero (a common assumption in epidemiology).
}
\item{splines_degree}{
Integer, either 0,1,2 or 3, specifying the polynomial degree of the provided \code{lambda}, \code{mu} and \code{psi} between grid points in \code{age_grid}. For example, if \code{splines_degree==1}, then the provided \code{lambda}, \code{mu} and \code{psi} are interpreted as piecewise-linear curves; if \code{splines_degree==2} the \code{lambda}, \code{mu} and \code{psi} are interpreted as quadratic splines; if \code{splines_degree==3} the \code{lambda}, \code{mu} and \code{psi} is interpreted as cubic splines. If your \code{age_grid} is fine enough, then \code{splines_degree=1} is usually sufficient.
}
\item{CSA_times}{
Optional numeric vector, listing times of concentrated sampling attempts, in ascending order. Concentrated sampling is performed in addition to any continuous (Poissonian) sampling specified by \code{psi}.
}
\item{CSA_probs}{
Optional numeric vector of the same size as \code{CSA_times}, listing sampling probabilities at each concentrated sampling time. Note that in contrast to the sampling rates \code{psi}, the \code{CSA_probs} are interpreted as probabilities and must thus be between 0 and 1. \code{CSA_probs} must be provided if and only if \code{CSA_times} is provided.
}
\item{CSA_kappas}{
Optional numeric vector of the same size as \code{CSA_times}, listing sampling retention probabilities at each concentrated sampling time, i.e. the probability at which a sampled lineage is kept in the pool of extant lineages. Note that the \code{CSA_kappas} are probabilities and must thus be between 0 and 1. \code{CSA_kappas} must be provided if and only if \code{CSA_times} is provided.
}
\item{no_full_extinction}{
Logical, specifying whether to prevent complete extinction of the tree. Full extinction is prevented by temporarily disabling extinctions whenever the number of extant tips is 1. Note that, strictly speaking, the trees generated do not exactly follow the proper probability distribution when \code{no_full_extinction} is \code{TRUE}.
}
\item{max_runtime}{
Numeric, optional maximum computation time (in seconds) to allow for the simulation. Use this to avoid occasional explosions of runtimes, for example due to very large generated trees. Aborted simulations will return with the flag \code{success=FALSE} (i.e., no tree is returned at all).
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
\item{include_birth_times}{
Logical. If \code{TRUE}, then the times of speciation events (in order of occurrence) will also be returned.
}
\item{include_death_times}{
Logical. If \code{TRUE}, then the times of extinction events (in order of occurrence) will also be returned.
}
}
\details{
The simulation proceeds in forward time, starting with a single root. Speciation/extinction and continuous (Poissonian) sampling events are drawn at exponentially distributed time steps, according to the rates specified by \code{lambda}, \code{mu} and \code{psi}. Sampling also occurs at the optional \code{CSA_times}. Only extant lineages are sampled at any time point, and sampled lineages are removed from the pool of extant lineages at probability 1-\code{kappa}.

The simulation halts as soon as one of the halting criteria are met, as specified by the options \code{max_sampled_tips}, \code{max_sampled_nodes}, \code{max_extant_tips}, \code{max_extinct_tips}, \code{max_tips} and \code{max_time}, or if no extant tips remain, whichever occurs first. Note that in some scenarios (e.g., if extinction rates are very high) the simulation may halt too early and the generated tree may only contain a single tip (i.e., the root lineage); in that case, the simulation will return an error (see return value \code{success}).

The function returns a single generated tree, as well as supporting information such as which tips are extant, extinct or sampled.
}

\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the simulation was successful. If \code{FALSE}, then the returned list includes an additional `\code{error}' element (character) providing a description of the error; all other return variables may be undefined.
}
\item{tree}{The generated timetree, of class "phylo". Note that this tree need not be ultrametric, for example if sampling occurs at multiple time points.
}
\item{root_time}{
Numeric, giving the time at which the tree's root was first split during the simulation. Note that this may be greater than 0, i.e., if the tips of the final tree do not coalesce all the way back to the simulation's start.
}
\item{final_time}{
Numeric, giving the final time at the end of the simulation.
}
\item{root_age}{
Numeric, giving the age (time before present) at the tree's root. This is equal to \code{final_time}-\code{root_time}.
}
\item{Nbirths}{
Integer, the total number of speciation (birth) events that occured during the simulation.
}
\item{Ndeaths}{
Integer, the total number of extinction (death) events that occured during the simulation.
}
\item{Nsamplings}{
Integer, the total number of sampling events that occured during the simulation.
}
\item{Nretentions}{
Integer, the total number of sampling events that occured during the simulation and for which lineages were kept in the pool of extant lineages.
}
\item{sampled_clades}{
Integer vector, specifying indices (from 1 to Ntips+Nnodes) of sampled tips and nodes in the final tree (regardless of whether their lineages were subsequently retained or removed from the pool).
}
\item{retained_clades}{
Integer vector, specifying indices (from 1 to Ntips+Nnodes) of sampled tips and nodes in the final tree that were retained, i.e., not removed from the pool following sampling.
}
\item{extant_tips}{
Integer vector, specifying indices (from 1 to Ntips) of extant (non-sampled and non-extinct) tips in the final tree. Will be empty if \code{include_extant==FALSE}.
}
\item{extinct_tips}{
Integer vector, specifying indices (from 1 to Ntips) of extinct (non-sampled and non-extant) tips in the final tree. Will be empty if \code{include_extinct==FALSE}.
}
}


\author{Stilianos Louca}

\references{
T. Stadler (2010). Sampling-through-time in birth--death trees. Journal of Theoretical Biology. 267:396-404.

T. Stadler et al. (2013). Birth--death skyline plot reveals temporal changes of epidemic spread in HIV and hepatitis C virus (HCV). PNAS. 110:228-233.
}

\seealso{
\code{\link{generate_tree_hbd_reverse}},
\code{\link{generate_gene_tree_msc}},
\code{\link{generate_random_tree}},
\code{\link{fit_hbds_model_parametric}},
\code{\link{simulate_deterministic_hbds}}
}

\examples{
# define time grid on which lambda, mu and psi will be specified
time_grid = seq(0,100,length.out=1000)

# specify the time-dependent extinction rate mu on the time-grid
mu_grid = 0.5*time_grid/(10+time_grid)

# define additional concentrated sampling attempts
CSA_times  = c(5,7,9)
CSA_probs  = c(0.5, 0.5, 0.5)
CSA_kappas = c(0.2, 0.1, 0.1)

# generate tree with a constant speciation & sampling rate,
# time-variable extinction rate and additional discrete sampling points
# assuming that all continuously sampled lineages are removed from the pool
simul = generate_tree_hbds( max_time        = 10,
                            include_extant  = FALSE,
                            include_extinct = FALSE,
                            time_grid       = time_grid,
                            lambda          = 1,
                            mu              = mu_grid,
                            psi             = 0.1,
                            kappa           = 0,
                            CSA_times       = CSA_times,
                            CSA_probs       = CSA_probs,
                            CSA_kappas      = CSA_kappas);
if(!simul$success){
    cat(sprintf("ERROR: Could not simulate tree: \%s\n",simul$error))
}else{
    # simulation succeeded. print some basic info about the generated tree
    tree = simul$tree
    cat(sprintf("Generated tree has \%d tips\n",length(tree$tip.label)))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree model}
\keyword{simulation}
\keyword{birth-death-sampling model}
