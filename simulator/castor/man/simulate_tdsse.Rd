\name{simulate_tdsse}
\alias{simulate_tdsse}
\title{
Simulate a time-dependent Discrete-State Speciation and Extinction (tdSSE) model.
}
\description{
Simulate a random phylogenetic tree in forward time based on a Poissonian speciation/extinction (birth/death) process, whereby birth and death rates are determined by a co-evolving discrete trait. New species are added (born) by splitting of a randomly chosen extant tip. The discrete trait, whose values determine birth/death rates, can evolve in two modes: (A) Anagenetically, i.e. according to a discrete-space continuous-time Markov process along each edge, with fixed or time-dependent transition rates between states, and/or (B) cladogenetically, i.e. according to fixed or time-dependent transition probabilities between states at each speciation event. This model class includes the Multiple State Speciation and Extinction (MuSSE) model described by FitzJohn et al. (2009), as well as the Cladogenetic SSE (ClaSSE) model described by Goldberg and Igis (2012). Optionally, the model can be turned into a Hidden State Speciation and Extinction model (Beaulieu and O'meara, 2016), by replacing the simulated tip/node states with "proxy" states, thus hiding the original states actually influencing speciation/extinction rates. This function is similar to \code{\link{simulate_dsse}}, the main difference being that state-specific speciation/extinction rates as well as state transition rates can be time-dependent.
}
\usage{
simulate_tdsse(Nstates,
               NPstates             = NULL,
               proxy_map            = NULL,
               time_grid            = NULL,
               parameters           = list(),
               splines_degree       = 1,
               start_state          = NULL,
               max_tips             = NULL, 
               max_time             = NULL,
               max_events           = NULL,
               sampling_fractions   = NULL,
               reveal_fractions     = NULL,
               coalescent           = TRUE,
               as_generations       = FALSE,
               no_full_extinction   = TRUE,
               Nsplits              = 2, 
               tip_basename         = "", 
               node_basename        = NULL,
               include_birth_times  = FALSE,
               include_death_times  = FALSE,
               include_labels       = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Nstates}{
Integer, specifying the number of possible discrete states a tip can have, influencing speciation/extinction rates. For example, if \code{Nstates==2} then this corresponds to the common Binary State Speciation and Extinction (BiSSE) model (Maddison et al., 2007). In the case of a HiSSE model, \code{Nstates} refers to the total number of diversification rate categories, as described by Beaulieu and O'meara (2016).
}
\item{NPstates}{
Integer, optionally specifying a number of "proxy-states" that are observed instead of the underlying speciation/extinction-modulating states. To simulate a HiSSE model, this should be smaller than \code{Nstates}. Each state corresponds to a different proxy-state, as defined using the variable \code{proxy_map} (see below). For BiSSE/MuSSE with no hidden states, \code{NPstates} can be set to either \code{NULL} or equal to \code{Nstates}, and proxy-states are equivalent to states.
}
\item{proxy_map}{
Integer vector of size \code{Nstates} and with values in 1,..,\code{NPstates}, specifying the correspondence between states (i.e. diversification-rate categories) and (observed) proxy-states, in a HiSSE model. Specifically, \code{proxy_map[s]} indicates which proxy-state the state s is represented by. Each proxy-state can represent multiple states (i.e. proxies are ambiguous), but each state must be represented by exactly one proxy-state. For non-HiSSE models, set this to \code{NULL}. See below for more details.
}
\item{time_grid}{
Numeric vector listing discrete times in ascending order, used to define the time-dependent rates of the model. The time grid should generally cover the maximum possible simulation time, otherwise it will be polynomially extrapolated (according to \code{splines_degree}).
}
\item{parameters}{
A named list specifying the time-dependent model parameters, including optional anagenetic and/or cladogenetic transition rates between states, as well as the mandatory state-dependent birth/death rates (see details below).
}
\item{splines_degree}{
Integer, either 0, 1, 2 or 3, specifying the polynomial degree of time-dependent model parameters (birth_rates, death_rates, transition_rates) between time-grid points. For example, \code{splines_degree=1} means that rates are to be considered linear between adjacent grid points.
}
\item{start_state}{
Integer within 1,..,\code{Nstates}, specifying the initial state, i.e. of the first lineage created. If left unspecified, this is chosen randomly and uniformly among all possible states.
}
\item{max_tips}{
Maximum number of tips in the generated tree, prior to any subsampling. If \code{coalescent=TRUE}, this refers to the number of extant tips, prior to subsampling. Otherwise, it refers to the number of extinct + extant tips, prior to subsampling. If \code{NULL} or <=0, the number of tips is not limited, so you should use \code{max_time} and/or \code{max_time_eq} and/or \code{max_events} to stop the simulation.
}
\item{max_time}{
Numeric, maximum duration of the simulation. If \code{NULL} or <=0, this constraint is ignored.
}
\item{max_events}{
Integer, maximum number of speciation/extinction/transition events before halting the simulation. If \code{NULL}, this constraint is ignored.
}
\item{sampling_fractions}{
A single number, or a numeric vector of size \code{NPstates}, listing tip sub-sampling fractions, depending on proxy-state. \code{sampling_fractions[p]} is the probability of including a tip in the final tree, if its proxy-state is p. If \code{NULL}, all tips (or all extant tips, if \code{coalescent==TRUE}) are included in the tree. If a single number, all tips are included with the same probability, i.e. regardless of their proxy-state.
}
\item{reveal_fractions}{
Numeric vector of size \code{NPstates}, listing reveal fractions of tip proxy-states, depending on proxy state. \code{reveal_fractions[p]} is the probability of knowing a tip's proxy-state, if its proxy state is p. Can also be NULL, in which case all tip proxy states will be known.
}
\item{coalescent}{
Logical, specifying whether only the coalescent tree (i.e. the tree spanning the extant tips) should be returned. If \code{coalescent==FALSE} and the death rate is non-zero, then the tree may include non-extant tips (i.e. tips whose distance from the root is less than the total time of evolution). In that case, the tree will not be ultrametric.
}
\item{as_generations}{
Logical, specifying whether edge lengths should correspond to generations. If FALSE, then edge lengths correspond to time.
}
\item{no_full_extinction}{
Logical, specifying whether to prevent complete extinction of the tree. Full extinction is prevented by temporarily disabling extinctions whenever the number of extant tips is 1. if \code{no_full_extinction==FALSE} and death rates are non-zero, the tree may go extinct during the simulation; if \code{coalescent==TRUE}, then the returned tree would be empty, hence the function will return unsuccessfully (i.e. \code{success} will be \code{FALSE}). By default \code{no_full_extinction} is \code{TRUE}, however in some special cases it may be desirable to allow full extinctions to ensure that the generated trees are statistically distributed exactly according to the underlying cladogenetic model.
}
\item{Nsplits}{
Integer greater than 1. Number of child-tips to generate at each diversification event. If set to 2, the generated tree will be bifurcating. If >2, the tree will be multifurcating.
}
\item{tip_basename}{
Character. Prefix to be used for tip labels (e.g. "tip."). If empty (""), then tip labels will be integers "1", "2" and so on.
}
\item{node_basename}{
Character. Prefix to be used for node labels (e.g. "node."). If \code{NULL}, no node labels will be included in the tree.
}
\item{include_birth_times}{
Logical. If \code{TRUE}, then the times of speciation events (in order of occurrence) will also be returned.
}
\item{include_death_times}{
Logical. If \code{TRUE}, then the times of extinction events (in order of occurrence) will also be returned.
}
\item{include_labels}{
Logical, specifying whether to include tip-labels and node-labels (if available) as names in the returned state vectors (e.g. \code{tip_states} and \code{node_states}). In any case, returned states are always listed in the same order as tips and nodes in the tree. Setting this to \code{FALSE} may increase computational efficiency for situations where labels are not required.
}
}


\details{
The function \code{simulate_tdsse} can be used to simulate a diversification + discrete-trait evolutionary process, in which birth/death (speciation/extinction) rates at each tip are determined by a tip's current "state". Lineages can transition between states anagenetically along each edge (according to some Markov transition rates) and/or cladogenetically at each speciation event (according to some transition probabilities). The speciation and extinction rates, as well as the transition rates, may be specified as time-dependent variables, defined as piecewise polynomial functions (natural splines) on a temporal grid.

In the following, Ngrid refers to the length of the vector \code{time_grid}.
The argument \code{parameters} should be a named list including one or more of the following elements:
\itemize{
\item\code{birth_rates}: Numeric 2D matrix of size Nstates x Ngrid, listing the per-capita birth rate (speciation rate) at each state and at each time-grid point. Can also be a single number (same birth rate for all states and at all times).
\item\code{death_rates}: Numeric 2D matrix of size Nstates x Ngrid, listing the per-capita death rate (extinction rate) at each state and at each time-grid point. Can also be a single number (same death rate for all states and at all times) or NULL (no deaths).
\item\code{transition_matrix_A}: Either a 3D numeric array of size Nstates x Nstates x Ngrid, or a 2D numeric matrix of size Nstates x Nstates, listing anagenetic transition rates between states along an edge. If a 3D array, then \code{transition_matrix_A[r,c,t]} is the infinitesimal rate for transitioning from state \code{r} to state \code{c} at time \code{time_grid[t]}. If a 2D matrix, \code{transition_matrix_A[r,c]} is the time-independent infintesimal rate for transitioning from state \code{r} to state \code{c}. At each time point (i.e., a fixed \code{t}), non-diagonal entries in \code{transition_matrix_A[,,t]} must be non-negative, diagonal entries must be non-positive, and the sum of each row must be zero.
\item\code{transition_matrix_C}: Either a 3D numeric array of size Nstates x Nstates x Ngrid, or a 2D numeric matrix of size Nstates x Nstates, listing cladogenetic transition probabilities between states during a speciation event, seperately for each child. If a 3D array, then \code{transition_matrix_C[r,c,t]} is the probability that a child emerging at time \code{time_grid[t]} will have state \code{c}, conditional upon the occurrence of a speciation event, given that the parent had state \code{r}, and independently of all other children. If a 2D matrix, then \code{transition_matrix_C[r,c]} is the (time-independent) probability that a child will have state \code{c}, conditional upon the occurrence of a speciation event, given that the parent had state \code{r}, and independently of all other children. Entries must be non-negative, and for any fixed \code{t} the sum of each row in \code{transition_matrix[,,t]} must be one.
}

If \code{max_time==NULL} and \code{max_events==NULL}, then the returned tree will always contain \code{max_tips} tips. If at any moment during the simulation the tree only includes a single extant tip, and if \code{no_full_extinction=TRUE} the death rate is temporarily set to zero to prevent the complete extinction of the tree. If \code{max_tips==NULL}, then the simulation is ran as long as specified by \code{max_time} and/or \code{max_events}. If neither \code{max_time}, \code{max_tips} nor \code{max_events} is \code{NULL}, then the simulation halts as soon as the time reaches \code{max_time}, or the number of tips (extant tips if \code{coalescent} is \code{TRUE}) reaches \code{max_tips}, or the number of speciation/extinction/transition events reaches \code{max_events} whichever occurs first. If \code{max_tips!=NULL} and \code{Nsplits>2}, then the last diversification even may generate fewer than \code{Nsplits} children, in order to keep the total number of tips within the specified limit. Note that this code generates trees in forward time, and halts as soon as one of the halting conditions is met; the halting condition chosen affects the precise distribution from which the generated trees are drawn (Stadler 2011).

For additional information on simulating HiSSE models see the related function \code{\link{simulate_dsse}}.

The parameter \code{transition_matrix_C} can be used to define ClaSSE models (Goldberg and Igic, 2012) or BiSSE-ness models (Magnuson-Ford and Otto, 2012), although care must be taken to properly define the transition probabilities. Here, cladogenetic transitions occur at probabilities that are defined conditionally upon a speciation event, whereas in other software they may be defined as probability rates.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the simulation was successful. If \code{FALSE}, an additional element \code{error} (of type character) is included containing an explanation of the error; in that case the value of any of the other elements is undetermined.
}
\item{tree}{
A rooted bifurcating (if \code{Nsplits==2}) or multifurcating (if \code{Nsplits>2}) tree of class "phylo", generated according to the specified birth/death model. 

If \code{coalescent==TRUE} or if all death rates are zero, and only if \code{as_generations==FALSE}, then the tree will be ultrametric. If \code{as_generations==TRUE} and \code{coalescent==FALSE}, all edges will have unit length.
}
\item{root_time}{
Numeric, giving the time at which the tree's root was first split during the simulation. 
Note that if \code{coalescent==TRUE}, this may be later than the first speciation event during the simulation.
}
\item{final_time}{
Numeric, giving the final time at the end of the simulation. If \code{coalescent==TRUE}, then this may be greater than the total time span of the tree (since the root of the coalescent tree need not correspond to the first speciation event).
}
\item{Nbirths}{
Numeric vector of size Nstates, listing the total number of birth events (speciations) that occurred at each state. The sum of all entries in \code{Nbirths} may be lower than the total number of tips in the tree if death rates were non-zero and \code{coalescent==TRUE}, or if \code{Nsplits>2}.
}
\item{Ndeaths}{
Numeric vector of size Nstates, listing the total number of death events (extinctions) that occurred at each state.
}
\item{Ntransitions_A}{
2D numeric matrix of size Nstates x Nstates, listing the total number of anagenetic transition events that occurred between each pair of states. For example, \code{Ntransitions_A[1,2]} is the number of anagenetic transitions (i.e., within a species) that occured from state 1 to state 2.
}
\item{Ntransitions_C}{
2D numeric matrix of size Nstates x Nstates, listing the total number of cladogenetic transition events that occurred between each pair of states. For example, \code{Ntransitions_C[1,2]} is the number of cladogenetic transitions (i.e., from a parent to a child) that occured from state 1 to state 2 during some speciation event. Note that each speciation event will have caused \code{Nsplits} transitions, and that the emergence of a child with the same state as the parent is counted as a transition between the same state (diagonal entries in \code{Ntransitions_C}).
}
\item{tip_states}{
Integer vector of size Ntips and with values in 1,..,Nstates, listing the state of each tip in the tree.
}
\item{node_states}{
Integer vector of size Nnodes and with values in 1,..,Nstates, listing the state of each node in the tree.
}
\item{tip_proxy_states}{
Integer vector of size Ntips and with values in 1,..,NPstates, listing the proxy state of each tip in the tree. Only included in the case of HiSSE models.
}
\item{node_proxy_states}{
Integer vector of size Nnodes and with values in 1,..,NPstates, listing the proxy state of each node in the tree. Only included in the case of HiSSE models.
}
\item{start_state}{
Integer, specifying the state of the first lineage (either provided during the function call, or generated randomly).
}
\item{birth_times}{
Numeric vector, listing the times of speciation events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root.
}
\item{death_times}{
Numeric vector, listing the times of extinction events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root.
}
}

\author{Stilianos Louca}

\references{
W. P. Maddison, P. E. Midford, S. P. Otto (2007). Estimating a binary character's effect on speciation and extinction. Systematic Biology. 56:701-710.

R. G. FitzJohn, W. P. Maddison, S. P. Otto (2009). Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology. 58:595-611

R. G. FitzJohn (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution. 3:1084-1092

E. E. Goldberg, B. Igic (2012). Tempo and mode in plant breeding system evolution. Evolution. 66:3701-3709.

K. Magnuson-Ford, S. P. Otto (2012). Linking the investigations of character evolution and species diversification. The American Naturalist. 180:225-245.

J. M. Beaulieu and B. C. O'Meara (2016). Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic Biology. 65:583-601.

T. Stadler (2011). Simulating trees with a fixed number of extant species. Systematic Biology. 60:676-684.
}

\seealso{
\code{\link{simulate_dsse}},
\code{\link{simulate_musse}},
\code{\link{fit_musse}}
}

\examples{\dontrun{
# prepare params for time-dependent BiSSE model
# include time-dependent speciation & extinction rates
# as well as time-dependent anagenetic transition rates
Nstates          = 2
reveal_fractions = c(1,0.5)
rarefaction      = 0.5 # species sampling fraction

time2lambda1 = function(times) rep(1,times=length(times))
time2lambda2 = function(times) rep(2,times=length(times))
time2mu1     = function(times) 0.5 + 2.5*exp(-((times-8)**2)/2)
time2mu2     = function(times) 1 + 2*exp(-((times-12)**2)/2)
time_grid    = seq(from=0, to=100, length.out=1000)

time2Q12    = function(times) 1*exp(0.1*times)
time2Q21    = function(times) 2*exp(-0.1*times)
QA          = array(0, dim=c(Nstates,Nstates,length(time_grid)))
QA[1,2,]    = time2Q12(time_grid)
QA[2,1,]    = time2Q21(time_grid)
QA[1,1,]    = -QA[1,2,]
QA[2,2,]    = -QA[2,1,]

parameters = list()
parameters$birth_rates = rbind(time2lambda1(time_grid), time2lambda2(time_grid))
parameters$death_rates = rbind(time2mu1(time_grid), time2mu2(time_grid))
parameters$transition_matrix_A = QA

# simulate time-dependent BiSSE model
cat(sprintf("Simulating tMuSSE model..\n"))
sim = castor::simulate_tdsse(Nstates            = Nstates,
                            time_grid           = time_grid,
                            parameters          = parameters, 
                            splines_degree      = 1,
                            max_tips            = 10000/rarefaction,
                            sampling_fractions  = rarefaction,
                            reveal_fractions    = reveal_fractions,
                            coalescent          = TRUE,
                            no_full_extinction  = TRUE)
if(!sim$success){
    cat(sprintf("ERROR: \%s\n",sim$error))
}else{
    # print some summary info about the generated tree
    tree        = sim$tree
    Ntips       = length(tree$tip.label)
    root_age    = get_tree_span(tree)$max_distance
    root_time   = sim$final_time - root_age
    tip_states  = sim$tip_states
    Nknown_tips = sum(!is.na(tip_states))
    cat(sprintf("Note: Simulated tree has root_age = \%g\n",root_age))
    cat(sprintf("Note: \%d tips have known state\n", Nknown_tips));
}
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{tree model}
\keyword{simulation}
\keyword{MuSSE}
\keyword{BiSSE}
