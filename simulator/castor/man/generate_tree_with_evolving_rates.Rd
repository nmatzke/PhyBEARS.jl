\name{generate_tree_with_evolving_rates}
\alias{generate_tree_with_evolving_rates}
\title{
Generate a random tree with evolving speciation/extinction rates.
}
\description{
Generate a random phylogenetic tree via simulation of a Poissonian speciation/extinction (birth/death) process. New species are added (born) by splitting of a randomly chosen extant tip. Per-capita birth and death rates (aka. speciation and extinction rates) evolve under some stochastic process (e.g. Brownian motion) along each edge. Thus, the probability rate of a tip splitting or going extinct depends on the tip, with closely related tips having more similar per-capita birth and death rates.
}
\usage{
generate_tree_with_evolving_rates(parameters           = list(),
                                  rate_model           = 'BM',
                                  max_tips             = NULL, 
                                  max_time             = NULL,
                                  max_time_eq          = NULL,
                                  coalescent           = TRUE,
                                  as_generations       = FALSE,
                                  tip_basename         = "", 
                                  node_basename        = NULL,
                                  include_event_times  = FALSE,
                                  include_rates        = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{parameters}{
A named list specifying the model parameters for the evolving birth/death rates. The precise entries expected depend on the chosen \code{rate_model} (see details below).
}

\item{rate_model}{
Character, specifying the model for the evolving per-capita birth/death rates. Must be one of the following: 'BM' (Brownian motion constrained to a finite interval via reflection), 'Mk' (discrete-state continuous-time Markov chain with fixed transition rates).
}

\item{max_tips}{
Maximum number of tips of the tree to be generated. If \code{coalescent=TRUE}, this refers to the number of extant tips. Otherwise, it refers to the number of extinct + extant tips. If \code{NULL} or <=0, the number of tips is unlimited (so be careful).
}
\item{max_time}{
Maximum duration of the simulation. If \code{NULL} or <=0, this constraint is ignored.
}
\item{max_time_eq}{
Maximum duration of the simulation, counting from the first point at which speciation/extinction equilibrium is reached, i.e. when (birth rate - death rate) changed sign for the first time. If \code{NULL} or <0, this constraint is ignored.
}
\item{coalescent}{
Logical, specifying whether only the coalescent tree (i.e. the tree spanning the extant tips) should be returned. If \code{coalescent==FALSE} and the death rate is non-zero, then the tree may include non-extant tips (i.e. tips whose distance from the root is less than the total time of evolution). In that case, the tree will not be ultrametric.
}
\item{as_generations}{
Logical, specifying whether edge lengths should correspond to generations. If FALSE, then edge lengths correspond to time.
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
Logical. If \code{TRUE}, then the bper-capita birth & death rates of all tips and nodes will also be returned.
}
}


\details{
If \code{max_time==NULL}, then the returned tree will always contain \code{max_tips} tips. In particular, if at any moment during the simulation the tree only includes a single extant tip, the death rate is temporarily set to zero to prevent the complete extinction of the tree. If \code{max_tips==NULL}, then the simulation is ran as long as specified by \code{max_time}. If neither \code{max_time} nor \code{max_tips} is \code{NULL}, then the simulation halts as soon as the time exceeds \code{max_time} or the number of tips (extant tips if \code{coalescent} is \code{TRUE}) exceeds \code{max_tips}.

If \code{rate_model=='BM'}, then per-capita birth rates (speciation rates) and per-capita death rates (extinction rates) evolve according to Brownian Motion, constrained to a finite interval via reflection. Note that speciation and extinction rates are only updated at branching points, i.e. during speciation events, while waiting times until speciation/extinction are based on rates at the previous branching point. The argument \code{parameters} should be a named list including one or more of the following elements:
\itemize{
\item\code{birth_rate_diffusivity}: Non-negative number. Diffusivity constant for the Brownian motion model of the evolving per-capita birth rate. In units 1/time^3. See \code{\link{simulate_bm_model}} for an explanation of the diffusivity parameter.
\item\code{min_birth_rate_pc}:
Non-negative number. The minimum allowed per-capita birth rate of a clade. In units 1/time. By default this is 0.
\item\code{max_birth_rate_pc}:
Non-negative number. The maximum allowed per-capita birth rate of a clade. In units 1/time. By default this is 1.
\item\code{death_rate_diffusivity}: Non-negative number. Diffusivity constant for the Brownian motion model of the evolving per-capita death rate. In units 1/time^3. See \code{\link{simulate_bm_model}} for an explanation of the diffusivity parameter.
\item\code{min_death_rate_pc}:
Non-negative number. The minimum allowed per-capita death rate of a clade. In units 1/time. By default this is 0.
\item\code{max_death_rate_pc}:
Non-negative number. The maximum allowed per-capita death rate of a clade. In units 1/time. By default this is 1.
\item\code{root_birth_rate_pc}:
Non-negative number, between \code{min_birth_rate_pc} and \code{max_birth_rate_pc}, specifying the initial per-capita birth rate of the root. If left unspecified, this will be chosen randomly and uniformly within the allowed interval.
\item\code{root_death_rate_pc}:
Non-negative number, between \code{min_death_rate_pc} and \code{max_death_rate_pc}, specifying the initial per-capita death rate of the root. If left unspecified, this will be chosen randomly and uniformly within the allowed interval.
\item\code{rarefaction}:
Numeric between 0 and 1. Rarefaction to be applied at the end of the simulation (fraction of random tips kept in the tree). 
Note that if \code{coalescent==FALSE}, rarefaction may remove both extant as well as extinct clades. Set \code{rarefaction=1} to not perform any rarefaction.
}

If \code{rate_model=='Mk'}, then speciation/extinction rates are determined by a tip's current "state", which evolves according to a continuous-time discrete-state Markov chain (Mk model) with constant transition rates. The argument \code{parameters} should be a named list including one or more of the following elements:
\itemize{
\item\code{Nstates}: Number of possible discrete states a tip can have. For example, if \code{Nstates} then this corresponds to the common Binary State Speciation and Extinction (BiSSE) model (Maddison et al., 2007). By default this is 1.
\item\code{state_birth_rates}: Numeric vector of size Nstates, listing the per-capita birth rate (speciation rate) at each state. Can also be a single number (all states have the same birth rate).
\item\code{state_death_rates}: Numeric vector of size Nstates, listing the per-capita death rate (extinction rate) at each state. Can also be a single number (all states have the same death rate).
\item\code{transition_matrix}: 2D numeric matrix of size Nstates x Nstates. Transition rate matrix for the Markov chain model of birth/death rate evolution.
\item\code{start_state}:
Integer within 1,..,\code{Nstates}, specifying the initial state of the first created lineage. If left unspecified, this is chosen randomly and uniformly among all possible states.
\item\code{rarefaction}: Same as when \code{rate_model=='BM'}.
}
Note: The option \code{rate_model=='Mk'} is deprecated and included for backward compatibility purposes only. To generate a tree with Markov transitions between states (known as Multiple State Speciation and Extinction model), use the command \command{simulate_dsse} instead.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the simulation was successful. If \code{FALSE}, an additional element \code{error} (of type character) is included containing an explanation of the error; in that case the value of any of the other elements is undetermined.
}
\item{tree}{
A rooted bifurcating tree of class "phylo", generated according to the specified birth/death model. 

If \code{coalescent==TRUE} or if all death rates are zero, and only if \code{as_generations==FALSE}, then the tree will be ultrametric. If \code{as_generations==TRUE} and \code{coalescent==FALSE}, all edges will have unit length.
}
\item{root_time}{
Numeric, giving the time at which the tree's root was first split during the simulation. 
Note that if \code{coalescent==TRUE}, this may be later than the first speciation event during the simulation.
}
\item{final_time}{
Numeric, giving the final time at the end of the simulation. If \code{coalescent==TRUE}, then this may be greater than the total time span of the tree (since the root of the coalescent tree need not correspond to the first speciation event).
}
\item{equilibrium_time}{
Numeric, giving the first time where the sign of (death rate - birth rate) changed from the beginning of the simulation, i.e. when speciation/extinction equilibrium was reached. May be infinite if the simulation stoped before reaching this point.
}
\item{Nbirths}{
Total number of birth events (speciations) that occurred during tree growth. This may be lower than the total number of tips in the tree if death rates were non-zero and \code{coalescent==TRUE}.
}
\item{Ndeaths}{
Total number of deaths (extinctions) that occurred during tree growth.
}
\item{birth_times}{
Numeric vector, listing the times of speciation events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root. Only returned if \code{include_event_times==TRUE}.
}
\item{death_times}{
Numeric vector, listing the times of extinction events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root. Only returned if \code{include_event_times==TRUE}.
}
\item{birth_rates_pc}{
Numeric vector of length Ntips+Nnodes, listing the per-capita birth rate of each tip and node in the tree. The length of an edge in the tree was thus drawn from an exponential distribution with rate equal to the per-capita birth rate of the child tip or node.
}
\item{death_rates_pc}{
Numeric vector of length Ntips+Nnodes, listing the per-capita death rate of each tip and node in the tree.
}
\item{states}{
Integer vector of size Ntips+Nnodes, listing the discrete state of each tip and node in the tree. Only included if \code{rate_model=="Mk"}.
}
\item{start_state}{
Integer, specifying the initial state of the first created lineage (either provided during the function call, or generated randomly). Only included if \code{rate_model=="Mk"}.
}
\item{root_birth_rate_pc}{
Numeric, specifying the initial per-capita birth rate of the root (either provided during the function call, or generated randomly). Only included if \code{rate_model=="BM"}.
}
\item{root_death_rate_pc}{
Numeric, specifying the initial per-capita death rate of the root (either provided during the function call, or generated randomly). Only included if \code{rate_model=="BM"}.
}
}

\author{Stilianos Louca}

\references{
D. J. Aldous (2001). Stochastic models and descriptive statistics for phylogenetic trees, from Yule to today. Statistical Science. 16:23-34.

W. P. Maddison, P. E. Midford, S. P. Otto (2007). Estimating a binary character's effect on speciation and extinction. Systematic Biology. 56:701-710.
}

\seealso{
\code{\link{simulate_dsse}}
}

\examples{
# Example 1
# Generate tree, with rates evolving under Brownian motion
parameters = list(birth_rate_diffusivity  = 1,
                  min_birth_rate_pc       = 1,
                  max_birth_rate_pc       = 2,
                  death_rate_diffusivity  = 0.5,
                  min_death_rate_pc       = 0,
                  max_death_rate_pc       = 1)
simulation = generate_tree_with_evolving_rates(parameters,
                                               rate_model='BM',
                                               max_tips=1000,
                                               include_rates=TRUE)
tree  = simulation$tree
Ntips = length(tree$tip.label)

# plot per-capita birth & death rates of tips
plot( x=simulation$birth_rates_pc[1:Ntips], 
      y=simulation$death_rates_pc[1:Ntips], 
      type='p', 
      xlab="pc birth rate", 
      ylab="pc death rate", 
      main="Per-capita birth & death rates across tips (BM model)",
      las=1)

      
######################
# Example 2
# Generate tree, with rates evolving under a binary-state model
Q = get_random_mk_transition_matrix(Nstates=2, rate_model="ER", max_rate=0.1)
parameters = list(Nstates = 2,
                  state_birth_rates = c(1,1.5),
                  state_death_rates = 0.5,
                  transition_matrix = Q)
simulation = generate_tree_with_evolving_rates(parameters,
                                               rate_model='Mk',
                                               max_tips=1000,
                                               include_rates=TRUE)
tree  = simulation$tree
Ntips = length(tree$tip.label)

# plot distribution of per-capita birth rates of tips
rates = simulation$birth_rates_pc[1:Ntips]
barplot(table(rates)/length(rates), 
        xlab="rate", 
        main="Distribution of pc birth rates across tips (Mk model)")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{tree model}
\keyword{simulation}
\keyword{BM model}
\keyword{Mk model}
