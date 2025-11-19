\name{generate_random_tree}
\alias{generate_random_tree}
\title{
Generate a tree using a Poissonian speciation/extinction model.
}
\description{
Generate a random timetree via simulation of a Poissonian speciation/extinction (birth/death) process. New species are added (born) by splitting of a randomly chosen extant tip. The tree-wide birth and death rates of tips can each be constant or power-law functions of the number of extant tips. For example,
\deqn{
B = I + F\cdot N^E,
}
where \eqn{B} is the tree-wide birth rate (species generation rate), \eqn{I} is the intercept, \eqn{F} is the power-law factor, \eqn{N} is the current number of extant tips and \eqn{E} is the power-law exponent. Optionally, the per-capita (tip-specific) birth and death rates can be extended by adding a custom time series provided by the user.
}
\usage{
generate_random_tree(parameters           = list(),
                     max_tips             = NULL, 
                     max_extant_tips      = NULL,
                     max_time             = NULL,
                     max_time_eq          = NULL,
                     coalescent           = TRUE,
                     as_generations       = FALSE,
                     no_full_extinction   = TRUE,
                     Nsplits              = 2,
                     added_rates_times    = NULL,
                     added_birth_rates_pc = NULL,
                     added_death_rates_pc = NULL,
                     added_periodic       = FALSE,
                     tip_basename         = "", 
                     node_basename        = NULL,
                     edge_basename        = NULL,
                     include_birth_times  = FALSE,
                     include_death_times  = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{parameters}{
A named list specifying the birth-death model parameters, with one or more of the following entries:

\code{birth_rate_intercept}: Non-negative number. The intercept of the Poissonian rate at which new species (tips) are added. In units 1/time. By default this is 0.

\code{birth_rate_factor}:
Non-negative number. The power-law factor of the Poissonian rate at which new species (tips) are added. In units 1/time. By default this is 0.

\code{birth_rate_exponent}:
Numeric. The power-law exponent of the Poissonian rate at which new species (tips) are added. Unitless. By default this is 1.

\code{death_rate_intercept}:
Non-negative number. The intercept of the Poissonian rate at which extant species (tips) go extinct. In units 1/time. By default this is 0.

\code{death_rate_factor}:
Non-negative number. The power-law factor of the Poissonian rate at which extant species (tips) go extinct. In units 1/time. By default this is 0.

\code{death_rate_exponent}:
Numeric. The power-law exponent of the Poissonian rate at which extant species (tips) go extinct. Unitless. By default this is 1.

\code{resolution}:
Non-negative numeric, specifying the resolution (in time units) at which to collapse the final tree by combining closely related tips. Any node whose age is smaller than this threshold, will be represented by a single tip. Set \code{resolution=0} to not collapse tips (default).

\code{rarefaction}:
Numeric between 0 and 1. Rarefaction to be applied to the final tree (fraction of random tips kept in the tree). Note that if \code{coalescent==FALSE}, rarefaction may remove both extant as well as extinct clades. Set \code{rarefaction=1} to not perform any rarefaction (default).
}
\item{max_tips}{
Integer, maximum number of tips of the tree to be generated. If \code{coalescent=TRUE}, this refers to the number of extant tips. Otherwise, it refers to the number of extinct + extant tips. If \code{NULL} or <=0, this halting condition is ignored.
}
\item{max_extant_tips}{
Integer, maximum number of extant lineages allowed at any moment during the simulation. If this number is reached, the simulation is halted. If \code{NULL} or <=0, this halting condition is ignored.
}
\item{max_time}{
Numeric, maximum duration of the simulation. If \code{NULL} or <=0, this constraint is ignored.
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
\item{no_full_extinction}{
Logical, specifying whether to prevent complete extinction of the tree. Full extinction is prevented by temporarily disabling extinctions whenever the number of extant tips is 1. Note that, strictly speaking, the trees generated do not exactly follow the proper probability distribution when \code{no_full_extinction} is \code{TRUE}.
}
\item{Nsplits}{
Integer greater than 1. Number of child-tips to generate at each diversification event. If set to 2, the generated tree will be bifurcating. If >2, the tree will be multifurcating.
}
\item{added_rates_times}{
Numeric vector, listing time points (in ascending order) for the custom per-capita birth and/or death rates time series (see \code{added_birth_rates_pc} and \code{added_death_rates_pc} below). Can also be \code{NULL}, in which case the custom time series are ignored.
}
\item{added_birth_rates_pc}{
Numeric vector of the same size as \code{added_rates_times}, listing per-capita birth rates to be added to the power law part. Can also be \code{NULL}, in which case this option is ignored and birth rates are purely described by the power law.
}
\item{added_death_rates_pc}{
Numeric vector of the same size as \code{added_rates_times}, listing per-capita death rates to be added to the power law part. Can also be \code{NULL}, in which case this option is ignored and death rates are purely described by the power law.
}
\item{added_periodic}{
Logical, indicating whether \code{added_birth_rates_pc} and \code{added_death_rates_pc} should be extended periodically if needed (i.e. if not defined for the entire simulation time). If \code{FALSE}, added birth & death rates are extended with zeros.
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
If \code{max_time==NULL}, then the returned tree will always contain \code{max_tips} tips. In particular, if at any moment during the simulation the tree only includes a single extant tip, the death rate is temporarily set to zero to prevent the complete extinction of the tree. If \code{max_tips==NULL}, then the simulation is ran as long as specified by \code{max_time}. If neither \code{max_time} nor \code{max_tips} is \code{NULL}, then the simulation halts as soon as the time exceeds \code{max_time} or the number of tips (extant tips if \code{coalescent} is \code{TRUE}) exceeds \code{max_tips}. If \code{max_tips!=NULL} and \code{Nsplits>2}, then the last diversification even may generate fewer than \code{Nsplits} children, in order to keep the total number of tips within the specified limit.

If \code{rarefaction<1} and \code{resolution>0}, collapsing of closely related tips (at the resolution specified) takes place prior to rarefaction (i.e., subsampling applies to the already collapsed tips).

Both the per-capita birth and death rates can be made into completely arbitrary functions of time, by setting all power-law coefficients to zero and providing custom time series \code{added_birth_rates_pc} and \code{added_death_rates_pc}.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the tree was successfully generated. If \code{FALSE}, the only other value returned is \code{error}.
}
\item{tree}{
A rooted bifurcating (if \code{Nsplits==2}) or multifurcating (if \code{Nsplits>2}) tree of class "phylo", generated according to the specified birth/death model. If \code{coalescent==TRUE} or if all death rates are zero, and only if \code{as_generations==FALSE}, then the tree will be ultrametric. If \code{as_generations==TRUE} and \code{coalescent==FALSE}, all edges will have unit length.
}
\item{root_time}{
Numeric, giving the time at which the tree's root was first split during the simulation. Note that if \code{coalescent==TRUE}, this may be later than the first speciation event during the simulation.
}
\item{final_time}{
Numeric, giving the final time at the end of the simulation. Note that if \code{coalescent==TRUE}, then this may be greater than the total time span of the tree (since the root of the coalescent tree need not correspond to the first speciation event).
}
\item{root_age}{
Numeric, giving the age (time before present) at the tree's root. This is equal to \code{final_time}-\code{root_time}.
}
\item{equilibrium_time}{
Numeric, giving the first time where the sign of (death rate - birth rate) changed from the beginning of the simulation, i.e. when speciation/extinction equilibrium was reached. May be infinite if the simulation stoped before reaching this point.
}
\item{extant_tips}{
Integer vector, listing indices of extant tips in the tree. If \code{coalescent==TRUE}, all tips will be extant.
}
\item{Nbirths}{
Total number of birth events (speciations) that occurred during tree growth. This may be lower than the total number of tips in the tree if death rates were non-zero and \code{coalescent==TRUE}, or if \code{Nsplits>2}.
}
\item{Ndeaths}{
Total number of deaths (extinctions) that occurred during tree growth.
}
\item{Ncollapsed}{
Number of tips removed from the tree while collapsing at the resolution specified.
}
\item{Nrarefied}{
Number of tips removed from the tree due to rarefaction.
}
\item{birth_times}{
Numeric vector, listing the times of speciation events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root.
}
\item{death_times}{
Numeric vector, listing the times of extinction events during tree growth, in order of occurrence. Note that if \code{coalescent==TRUE}, then \code{speciation_times} may be greater than the phylogenetic distance to the coalescent root.
}
\item{error}{
Character, containing an explanation of ther error that occurred. Only included if \code{success==FALSE}.
}
}

\author{Stilianos Louca}

\references{
D. J. Aldous (2001). Stochastic models and descriptive statistics for phylogenetic trees, from Yule to today. Statistical Science. 16:23-34.

M. Steel and A. McKenzie (2001). Properties of phylogenetic trees generated by Yule-type speciation models. Mathematical Biosciences. 170:91-112.
}

%\seealso{
%}

\examples{
# Simple speciation model
parameters = list(birth_rate_intercept=1)
tree = generate_random_tree(parameters,max_tips=100)$tree

# Exponential growth rate model
parameters = list(birth_rate_factor=1)
tree = generate_random_tree(parameters,max_tips=100)$tree
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{tree model}
\keyword{simulation}
