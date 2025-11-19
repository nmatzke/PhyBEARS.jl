\name{simulate_dsse}
\alias{simulate_dsse}
\alias{simulate_musse}

\title{
Simulate a Discrete-State Speciation and Extinction (dSSE) model.
}
\description{
Simulate a random phylogenetic tree in forward time based on a Poissonian speciation/extinction (birth/death) process, with optional Poissonian sampling over time, whereby birth/death/sampling rates are determined by a co-evolving discrete trait. New species are added (born) by splitting of a randomly chosen extant tip. The discrete trait, whose values determine birth/death/sampling rates over time, can evolve in two modes: (A) Anagenetically, i.e. according to a discrete-space continuous-time Markov process along each edge, with fixed transition rates between states, and/or (B) cladogenetically, i.e. according to fixed transition probabilities between states at each speciation event. Poissonian lineage sampling is assumed to lead to a removal of lineages from the pool of extant tips (as is common in epidemiology). 

This model class includes the Multiple State Speciation and Extinction (MuSSE) model described by FitzJohn et al. (2009), as well as the Cladogenetic SSE (ClaSSE) model described by Goldberg and Igis (2012). Optionally, the model can be turned into a Hidden State Speciation and Extinction model (Beaulieu and O'meara, 2016), by replacing the simulated tip/node states with "proxy" states, thus hiding the original states actually influencing speciation/extinction rates.
}
\usage{
simulate_dsse( Nstates,
               NPstates                 = NULL,
               proxy_map                = NULL,
               parameters               = list(),
               start_state              = NULL,
               max_tips                 = NULL, 
               max_extant_tips          = NULL,
               max_Psampled_tips        = NULL,
               max_time                 = NULL,
               max_time_eq              = NULL,
               max_events               = NULL,
               sampling_fractions       = NULL,
               reveal_fractions         = NULL,
               sampling_rates           = NULL,
               coalescent               = TRUE,
               as_generations           = FALSE,
               no_full_extinction       = TRUE,
               tip_basename             = "", 
               node_basename            = NULL,
               include_event_times      = FALSE,
               include_rates            = FALSE,
               include_labels           = TRUE)
               
simulate_musse(Nstates,  NPstates = NULL, proxy_map = NULL, 
            parameters = list(), start_state = NULL, 
            max_tips = NULL, max_extant_tips = NULL, max_Psampled_tips = NULL,
            max_time = NULL, max_time_eq = NULL, max_events = NULL, 
            sampling_fractions = NULL, reveal_fractions = NULL, sampling_rates = NULL,
            coalescent = TRUE, as_generations = FALSE, no_full_extinction = TRUE, 
            tip_basename = "", node_basename = NULL, 
            include_event_times = FALSE, include_rates = FALSE, include_labels = TRUE)
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
Integer vector of size \code{Nstates} and with values in 1,..\code{NPstates}, specifying the correspondence between states (i.e. diversification-rate categories) and (observed) proxy-states, in a HiSSE model. Specifically, \code{proxy_map[s]} indicates which proxy-state the state s is represented by. Each proxy-state can represent multiple states (i.e. proxies are ambiguous), but each state must be represented by exactly one proxy-state. For non-HiSSE models, set this to \code{NULL}. See below for more details.
}
\item{parameters}{
A named list specifying the dSSE model parameters, such as the anagenetic and/or cladogenetic transition rates between states and the state-dependent birth/death rates (see details below).
}
\item{start_state}{
Integer within 1,..,\code{Nstates}, specifying the initial state, i.e. of the first lineage created. If left unspecified, this is chosen randomly and uniformly among all possible states.
}
\item{max_tips}{
Integer, maximum number of tips (extant + Poissonian-sampled if \code{coalescent==TRUE}, or extant+extinct+Poissonian-sampled if \code{coalescent==FALSE}) in the generated tree, shortly before any present-day sampling. If \code{NULL} or <=0, the number of tips is not limited, so you should use another stopping criterion such as \code{max_time} and/or \code{max_time_eq} and/or \code{max_events} to stop the simulation.
}
\item{max_extant_tips}{
Integer, maximum number of extant tips in the generated tree, shortly before to any present-day sampling. If \code{NULL} or <=0, this constraint is ignored.
}
\item{max_Psampled_tips}{
Integer, maximum number of Poissonian-sampled tips in the generated tree. If \code{NULL} or <=0, this constraint is ignored.
}
\item{max_time}{
Numeric, maximum duration of the simulation. If \code{NULL} or <=0, this constraint is ignored.
}
\item{max_time_eq}{
Numeric, maximum duration of the simulation, counting from the first point at which speciation/extinction equilibrium is reached, i.e. when (birth rate - death rate) changed sign for the first time. If \code{NULL} or <0, this constraint is ignored.
}
\item{max_events}{
Integer, maximum number of speciation/extinction/transition events before halting the simulation. If \code{NULL}, this constraint is ignored.
}
\item{sampling_fractions}{
A single number, or a numeric vector of size \code{NPstates}, listing the sampling fractions for extant tips at the end of the simulation (i.e., at "present-day")", depending on proxy-state. \code{sampling_fractions[p]} is the probability of including an extant tip in the final tree, if its proxy-state is p. If a single number, all extant tips are sampled with the same probability, i.e. regardless of their proxy-state. If \code{NULL}, this is the same as setting \code{sampling_fractions} to 1, i.e., all extant tips are sampled at the end of the simulation.
}
\item{reveal_fractions}{
Numeric vector of size \code{NPstates}, listing reveal fractions of tip proxy-states, depending on proxy state. \code{reveal_fractions[p]} is the probability of knowing a tip's proxy-state, if its proxy state is p. Can also be NULL, in which case all tip proxy states will be known.
}
\item{sampling_rates}{
Numeric vector of size \code{NPstates}, listing Poissonian sampling rates of lineages over time, depending on proxy state. Hence, \code{sampling_rates[p]} is the sampling rate of a lineage if its proxy state is p. Can also be a single numeric, thus applying the same sampling rate to all lineages regardless of proxy state. Can also be \code{NULL}, in which case Poissonian sampling is not included.
}
\item{coalescent}{
Logical, specifying whether only the coalescent tree (i.e. the tree spanning the sampled tips) should be returned. If \code{coalescent==FALSE} and the death rate is non-zero, then the tree may include extinct tips.
}
\item{as_generations}{
Logical, specifying whether edge lengths should correspond to generations. If FALSE, then edge lengths correspond to time.
}
\item{no_full_extinction}{
Logical, specifying whether to prevent complete extinction of the tree. Full extinction is prevented by temporarily disabling extinctions and Poissonian samplings whenever the number of extant tips is 1. if \code{no_full_extinction==FALSE} and death rates and/or Poissonian sampling rates are non-zero, the tree may go extinct during the simulation; if \code{coalescent==TRUE}, then the returned could end up empty, hence the function will return unsuccessfully (i.e. \code{success} will be \code{FALSE}). By default \code{no_full_extinction} is \code{TRUE}, however in some special cases it may be desirable to allow full extinctions to ensure that the generated trees are statistically distributed exactly according to the underlying cladogenetic model.
}
\item{tip_basename}{
Character. Prefix to be used for tip labels (e.g. "tip."). If empty (""), then tip labels will be integers "1", "2" and so on.
}
\item{node_basename}{
Character. Prefix to be used for node labels (e.g. "node."). If \code{NULL}, no node labels will be included in the tree.
}
\item{include_event_times}{
Logical. If \code{TRUE}, then the times of speciation and extinction events (each in order of occurrence) will also be returned.
}
\item{include_rates}{
Logical. If \code{TRUE}, then the per-capita birth & death rates of all tips and nodes will also be returned.
}
\item{include_labels}{
Logical, specifying whether to include tip-labels and node-labels (if available) as names in the returned state vectors (e.g. \code{tip_states} and \code{node_states}). In any case, returned states are always listed in the same order as tips and nodes in the tree. Setting this to \code{FALSE} may increase computational efficiency for situations where labels are not required.
}
}


\details{
The function \code{simulate_dsse} can be used to simulate a diversification + discrete-trait evolutionary process, in which birth/death (speciation/extinction) and Poissonian sampling rates at each tip are determined by a tip's current "state". Lineages can transition between states anagenetically along each edge (according to fixed Markov transition rates) and/or cladogenetically at each speciation event (according to fixed transition probabilities). In addition to Poissonian sampling through time (commonly included in epidemiological models), extant tips can also be sampled at the end of the simulation (i.e. at "present-day") according to some state-specific \code{sampling_fractions} (common in macroevolution).

The function \code{simulate_musse} is a simplified variant meant to simulate MuSSE/HiSSE models in the absence of cladogenetic state transitions, and is included mainly for backward-compatibility reasons. The input arguments for \code{simulate_musse} are identical to \code{simulate_dsse}, with the exception that the \code{parameters} argument must include slightly different elements (explained below). Note that the standard MuSSE/HiSSE models published by FitzJohn et al. (2009) and Beaulieu and O'meara (2016) did not include Poissonian sampling through time, i.e. sampling of extant lineages was only done once at present-day.

For \code{simulate_dsse}, the argument \code{parameters} should be a named list including one or more of the following elements:
\itemize{
\item\code{birth_rates}: Numeric vector of size Nstates, listing the per-capita birth rate (speciation rate) at each state. Can also be a single number (all states have the same birth rate).
\item\code{death_rates}: Numeric vector of size Nstates, listing the per-capita death rate (extinction rate) at each state. Can also be a single number (all states have the same death rate).
\item\code{transition_matrix_A}: 2D numeric matrix of size Nstates x Nstates, listing anagenetic transition rates between states along an edge. Hence, \code{transition_matrix_A[r,c]} is the probability rate for transitioning from state \code{r} to state \code{c}. Non-diagonal entries must be non-negative, diagonal entries must be non-positive, and the sum of each row must be zero.
\item\code{transition_matrix_C}: 2D numeric matrix of size Nstates x Nstates, listing cladogenetic transition probabilities between states during a speciation event, seperately for each child. Hence, \code{transition_matrix_C[r,c]} is the probability that a child will have state \code{c}, conditional upon the occurrence of a speciation event, given that the parent had state \code{r}, and independently of all other children. Entries must be non-negative, and the sum of each row must be one.
}

For \code{simulate_musse}, the argument \code{parameters} should be a named list including one or more of the following elements:
\itemize{
\item\code{birth_rates}: Same as for \code{simulate_dsse}.
\item\code{death_rates}: Same as for \code{simulate_dsse}.
\item\code{transition_matrix}: 2D numeric matrix of size Nstates x Nstates, listing anagenetic transition rates between states. This is equivalent to \code{transition_matrix_A} in \code{simulate_dsse}.
}

Note that this code generates trees in forward time, and halts as soon as one of the enabled halting conditions is met; the halting conditions chosen affects the precise probability distribution from which the generated trees are drawn (Stadler 2011).
If at any moment during the simulation the tree only includes a single extant tip, and if \code{no_full_extinction=TRUE}, the death and sampling rate are temporarily set to zero to prevent the complete extinction of the tree. The tree will be ultrametric if \code{coalescent==TRUE} (or death rates were zero) and Poissonian sampling was not included.

HiSSE models (Beaulieu and O'meara, 2016) are closely related to BiSSE/MuSSE models, the main difference being the fact that the actual diversification-modulating states are not directly observed. Hence, this function is also able to simulate HiSSE models, with appropriate choice of the input variables \code{Nstates}, \code{NPstates} and \code{proxy_map}. For example, \code{Nstates=4}, \code{NPstates=2} and \code{proxy_map=c(1,2,1,2)} specifies that states 1 and 3 are represented by proxy-state 1, and states 2 and 4 are represented by proxy-state 2. This is the original case described by Beaulieu and O'meara (2016); in their terminology, there would be 2 "hidden"" states ("0" and "1") and 2 "observed" (proxy) states ("A" and "B"), and the 4 diversification rate categories (\code{Nstates=4}) would be called "0A", "1A", "0B" and "1B", respectively. The somewhat different terminology used here allows for easier generalization to an arbitrary number of diversification-modulating states and an arbitrary number of proxy states. For example, if there are 6 diversification modulating states, represented by 3 proxy-states as 1->A, 2->A, 3->B, 4->C, 5->C, 6->C, then one would set \code{Nstates=6}, \code{NPstates=3} and \code{proxy_map=c(1,1,2,3,3,3)}.

The parameter \code{transition_matrix_C} can be used to define ClaSSE models (Goldberg and Igic, 2012) or BiSSE-ness models (Magnuson-Ford and Otto, 2012), although care must be taken to properly define the transition probabilities. Here, cladogenetic transitions occur at probabilities that are defined conditionally upon a speciation event, whereas in other software they may be defined as probability rates.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the simulation was successful. If \code{FALSE}, an additional element \code{error} (of type character) is included containing an explanation of the error; in that case the value of any of the other elements is undetermined.
}
\item{tree}{
A rooted bifurcating tree of class "phylo", generated according to the specified birth/death model. 

If \code{coalescent==TRUE} or if all death rates are zero, and only if \code{as_generations==FALSE} and in the absence of Poissonian sampling, then the tree will be ultrametric. If \code{as_generations==TRUE} and \code{coalescent==FALSE}, all edges will have unit length.
}
\item{root_time}{
Numeric, giving the time at which the tree's root was first split during the simulation. 
Note that if \code{coalescent==TRUE}, this may be later than the first speciation event during the simulation.
}
\item{final_time}{
Numeric, giving the final time at the end of the simulation. If \code{coalescent==TRUE}, then this may be greater than the total time span of the tree (since the root of the coalescent tree need not correspond to the first speciation event).
}
\item{equilibrium_time}{
Numeric, giving the first time where the sign of (death rate - birth rate) changed from the beginning of the simulation, i.e. when speciation/extinction equilibrium was reached. May be infinite if the simulation stopped before reaching this point.
}
\item{Nbirths}{
Integer vector of size Nstates, listing the total number of birth events (speciations) that occurred at each state. The sum of all entries in \code{Nbirths} may be lower than the total number of tips in the tree if death rates were non-zero and \code{coalescent==TRUE}.
}
\item{Ndeaths}{
Integer vector of size Nstates, listing the total number of death events (extinctions) that occurred at each state.
}
\item{NPsamplings}{
Integer vector of size Nstates, listing the total number of Poissonian sampling events that occurred at each state.
}
\item{Ntransitions_A}{
2D numeric matrix of size Nstates x Nstates, listing the total number of anagenetic transition events that occurred between each pair of states. For example, \code{Ntransitions_A[1,2]} is the number of anagenetic transitions (i.e., within a species) that occured from state 1 to state 2.
}
\item{Ntransitions_C}{
2D numeric matrix of size Nstates x Nstates, listing the total number of cladogenetic transition events that occurred between each pair of states. For example, \code{Ntransitions_C[1,2]} is the number of cladogenetic transitions (i.e., from a parent to a child) that occured from state 1 to state 2 during some speciation event. Note that each speciation event will have caused 2 transitions (one per child), and that the emergence of a child with the same state as the parent is counted as a transition between the same state (diagonal entries in \code{Ntransitions_C}).
}
\item{NnonsampledExtant}{
Integer, specifying the number of extant tips not sampled at the end, i.e., omitted from the tree.
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
\item{extant_tips}{
Integer vector, listing the indices of any extant tips in the tree.
}
\item{extinct_tips}{
Integer vector, listing the indices of any extinct tips in the tree. Note that if \code{coalescent==TRUE}, this vector will be empty.
}
\item{Psampled_tips}{
Integer vector, listing the indices of any Poissonian-sampled tips in the tree.
}
\item{birth_times}{
Numeric vector, listing the times of speciation events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root. Only returned if \code{include_event_times==TRUE}.
}
\item{death_times}{
Numeric vector, listing the times of extinction events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root. Only returned if \code{include_event_times==TRUE}.
}
\item{Psampling_times}{
Numeric vector, listing the times of Poissonian sampling events during tree growth, in order of occurrence. Only returned if \code{include_event_times==TRUE}.
}
\item{clade_birth_rates}{
Numeric vector of size Ntips+Nnodes, listing the per-capita birth rate of each tip and node in the tree. Only included if \code{include_rates==TRUE}.
}
\item{clade_death_rates}{
Numeric vector of size Ntips+Nnodes, listing the per-capita death rate of each tip and node in the tree. Only included if \code{include_rates==TRUE}.
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

S. Louca and M. W. Pennell (2020). A general and efficient algorithm for the likelihood of diversification and discrete-trait evolutionary models. Systematic Biology. 69:545-556.
}

\seealso{
\code{\link{simulate_tdsse}}, 
\code{\link{fit_musse}}
}

\examples{
# Simulate a tree under a classical BiSSE model
# I.e., anagenetic transitions between two states, no Poissonian sampling through time.
A = get_random_mk_transition_matrix(Nstates=2, rate_model="ER", max_rate=0.1)
parameters = list(birth_rates         = c(1,1.5),
                  death_rates         = 0.5,
                  transition_matrix_A = A)
simulation = simulate_dsse( Nstates         = 2, 
                            parameters      = parameters, 
                            max_extant_tips = 1000, 
                            include_rates   = TRUE)
tree       = simulation$tree
Ntips      = length(tree$tip.label)

# plot distribution of per-capita birth rates of tips
rates = simulation$clade_birth_rates[1:Ntips]
barplot(table(rates)/length(rates), 
        xlab="rate", 
        main="Distribution of pc birth rates across tips (BiSSE model)")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{tree model}
\keyword{simulation}
\keyword{MuSSE}
\keyword{BiSSE}
