\name{simulate_diversification_model}
\alias{simulate_diversification_model}
\title{
Simulate a deterministic uniform speciation/extinction model.
}
\description{
Simulate a speciation/extinction cladogenic model for diversity over time, in the derministic limit. Speciation (birth) and extinction (death) rates can each be constant or power-law functions of the number of extant species. For example,
\deqn{
B = I + F\cdot N^E,
}
where \eqn{B} is the birth rate, \eqn{I} is the intercept, \eqn{F} is the power-law factor, \eqn{N} is the current number of extant species and \eqn{E} is the power-law exponent. Optionally, the model can account for incomplete taxon sampling (rarefaction of tips) and for the effects of collapsing a tree at a non-zero resolution (i.e. clustering closely related tips into a single tip).
}
\usage{
simulate_diversification_model( times,
                                parameters            = list(),
                                added_rates_times     = NULL,
                                added_birth_rates_pc  = NULL,
                                added_death_rates_pc  = NULL,
                                added_periodic        = FALSE,
                                start_time            = NULL,
                                final_time            = NULL,
                                start_diversity	      = 1,
                                final_diversity       = NULL,
                                reverse               = FALSE,
                                include_coalescent    = FALSE,
                                include_event_rates   = FALSE,
                                include_Nevents       = FALSE,
                                max_runtime           = NULL)
}
\arguments{
\item{times}{
Numeric vector, listing the times for which to calculate diversities, as predicted by the model. Values must be in ascending order.
}
\item{parameters}{
A named list specifying the birth-death model parameters, with one or more of the following entries:
\itemize{
\item\code{birth_rate_intercept}: Non-negative number. The intercept of the Poissonian rate at which new species (tips) are added. In units 1/time.
\item\code{birth_rate_factor}:
Non-negative number. The power-law factor of the Poissonian rate at which new species (tips) are added. In units 1/time.
\item\code{birth_rate_exponent}:
Numeric. The power-law exponent of the Poissonian rate at which new species (tips) are added. Unitless.
\item\code{death_rate_intercept}:
Non-negative number. The intercept of the Poissonian rate at which extant species (tips) go extinct. In units 1/time.
\item\code{death_rate_factor}:
Non-negative number. The power-law factor of the Poissonian rate at which extant species (tips) go extinct. In units 1/time.
\item\code{death_rate_exponent}:
Numeric. The power-law exponent of the Poissonian rate at which extant species (tips) go extinct. Unitless.
\item\code{resolution}: Non-negative number. Time resolution at which the final tree is assumed to be collapsed. Units are time units. E.g. if this is 10, then all nodes of age 10 or less, are assumed to be collapsed into (represented by) a single tip. This can be used to model OTU trees, obtained after clustering strains by some similarity (=age) threshold. Set to 0 to disable collapsing. If left unspecified, this is set to 0.
\item\code{rarefaction}: Numeric between 0 and 1, specifying the fraction of tips kept in the final tree after random subsampling. Rarefaction is assumed to occur after collapsing at the specified resolution (if applicable). This can be used to model incomplete taxon sampling. If left unspecified, this is set to 1.
}
}
\item{added_rates_times}{
Numeric vector, listing time points (in ascending order) for a custom per-capita birth and/or death rates time series (see \code{added_birth_rates_pc} and \code{added_death_rates_pc} below). Can also be \code{NULL}, in which case the custom time series are ignored.
}
\item{added_birth_rates_pc}{
Numeric vector of the same size as \code{added_rates_times}, listing per-capita birth rates to be added to the power law part. Added rates are interpolated linearly between time points in \code{added_rates_times}. Can also be \code{NULL}, in which case this option is ignored and birth rates are purely described by the power law.
}
\item{added_death_rates_pc}{
Numeric vector of the same size as \code{added_rates_times}, listing per-capita death rates to be added to the power law part. Added rates are interpolated linearly between time points in \code{added_rates_times}. Can also be \code{NULL}, in which case this option is ignored and death rates are purely described by the power law.
}
\item{added_periodic}{
Logical, indicating whether \code{added_birth_rates_pc} and \code{added_death_rates_pc} should be extended periodically if needed (i.e. if not defined for the entire simulation time). If \code{FALSE}, added birth & death rates are extended with zeros.
}
\item{start_time}{
Numeric. Start time of the tree (<=\code{times[1]}). Can also be \code{NULL}, in which case it is set to the first value in \code{times}.
}
\item{final_time}{
Numeric. Final (ending) time of the tree (>=\code{max(times)}). Can also be \code{NULL}, in which case it is set to the last value in \code{times}.
}
\item{start_diversity}{
Numeric. Total diversity at \code{start_time}. Only relevant if \code{reverse==FALSE}.
}
\item{final_diversity}{
Numeric. Total diversity at \code{final_time}, i.e. the final diversity of the tree (total extant species at age 0). Only relevant if \code{reverse==TRUE}.
}
\item{reverse}{
Logical. If \code{TRUE}, then the tree model is simulated in backward time direction. In that case, \code{final_diversity} is interpreted as the known diversity at the last time point, and all diversities at previous time points are calculated based on the model. If \code{FALSE}, then the model is simulated in forward-time, with initial diversity given by \code{start_diversity}.
}
\item{include_coalescent}{
Logical, specifying whether the diversity corresponding to a coalescent tree (i.e. the tree spanning only extant tips) should also be calculated. If \code{coalescent==TRUE} and the death rate is non-zero, then the coalescent diversities will generally be lower than the total diversities.
}
\item{include_event_rates}{
Logical. If \code{TRUE}, then the birth (speciation) and death (extinction) rates (for each time point) are included as returned values. This comes at a moderate computational overhead.
}
\item{include_Nevents}{
Logical. If \code{TRUE}, then the cumulative birth (speciation) and death (extinction) events (for each time point) are included as returned values. This comes at a moderate computational overhead.
}
\item{max_runtime}{
Numeric. Maximum runtime (in seconds) allowed for the simulation. If this time is surpassed, the simulation aborts.
}
}


\details{
The simulation is deterministic, meaning that diversification is modeled using ordinary differential equations, not as a stochastic process. The simulation essentially computes the deterministic diversity over time, not an actual tree. For stochastic cladogenic simulations yielding a random tree, see \code{\link{generate_random_tree}} and \code{\link{simulate_dsse}}.

In the special case where per-capita birth and death rates are constant (i.e. \eqn{I=0} and \eqn{E=1} for birth and death rates), this function uses an explicit analytical solution to the underlying differential equations, and is thus much faster than in the general case.

If \code{rarefaction<1} and \code{resolution>0}, collapsing of closely related tips (at the resolution specified) is assumed to take place prior to rarefaction (i.e., subsampling applies to the already collapsed tips).
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the simulation was successful. If the simulation aborted due to runtime constraints (option \code{max_runtime}), \code{success} will be \code{FALSE}.
}
\item{total_diversities}{
Numeric vector of the same size as \code{times}, listing the total diversity (extant at each the time) for each time point in \code{times}.
}
\item{coalescent_diversities}{
Numeric vector of the same size as \code{times}, listing the coalescent diversity (i.e. as seen in the coalescent tree spanning only extant species) for each time point in \code{times}. Only included if \code{include_coalescent==TRUE}.
}
\item{birth_rates}{
Numeric vector of the same size as \code{times}, listing the speciation (birth) rate at each time point. Only included if \code{include_event_rates==TRUE}.
}
\item{death_rates}{
Numeric vector of the same size as \code{times}, listing the extinction (death) rate at each time point. Only included if \code{include_event_rates==TRUE}.
}
\item{Nbirths}{
Numeric vector of the same size as \code{times}, listing the cumulative number of speciation (birth) events up to each time point. Only included if \code{include_Nevents==TRUE}.
}
\item{Ndeaths}{
Numeric vector of the same size as \code{times}, listing the cumulative number of extinction (death) events up to each time point. Only included if \code{include_Nevents==TRUE}.
}
}

\author{Stilianos Louca}


\seealso{
\code{\link{generate_random_tree}},
\code{\link{count_lineages_through_time}}
}

\examples{
# Generate a tree
max_time = 100
parameters = list(birth_rate_intercept  = 10, 
                  birth_rate_factor     = 0,
                  birth_rate_exponent   = 0,
                  death_rate_intercept  = 0,
                  death_rate_factor     = 0,
                  death_rate_exponent   = 0,
                  resolution            = 20,
                  rarefaction           = 0.5)
generator = generate_random_tree(parameters,max_time=max_time)
tree = generator$tree
final_total_diversity = length(tree$tip.label)+generator$Nrarefied+generator$Ncollapsed

# Calculate diversity-vs-time curve for the tree
times = seq(from=0,to=0.99*max_time,length.out=50)
tree_diversities = count_lineages_through_time(tree, times=times)$lineages

# simulate diversity curve based on deterministic model
simulation = simulate_diversification_model(times,
                                            parameters,
                                            reverse=TRUE,
                                            final_diversity=final_total_diversity,
                                            include_coalescent=TRUE)
model_diversities = simulation$coalescent_diversities

# compare diversities in the tree to the simulated ones
plot(tree_diversities,model_diversities,xlab="tree diversities",ylab="simulated diversities")
abline(a=0,b=1,col="#A0A0A0") # show diagonal for reference
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{tree model}
\keyword{simulation}