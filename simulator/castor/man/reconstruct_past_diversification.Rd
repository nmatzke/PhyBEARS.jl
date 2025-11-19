\name{reconstruct_past_diversification}
\alias{reconstruct_past_diversification}
\title{
Reconstruct past diversification dynamics from a diversity time series.
}
\description{
Given a time series of past diversities (coalescent or not), this function estimates instantaneous birth (speciation) and death (extinction) rates that would lead to the observed diversity time series. The function is based on a deterministic model (or the continuum limit of a stochastic cladogenic model), in which instantaneous birth and death rates lead to a predictable growth of a tree (one new species per birth event). The reconstruction is non-parametric, i.e. does not rely on fitting a parameterized model. The reconstruction is only accurate in the deterministic limit, i.e. for high diversities where the stochastic nature of the cladogenic process diminishes. Of particular importance is the case where the time series is coalescent, i.e. represents the diversity (lineages-through-time) that would be represented in a coalescent tree with extinctions. 

Note: This function is included for legacy reasons mainly. In most cases users should instead use the functions \code{\link{fit_hbd_model_on_grid}} and \code{\link{fit_hbd_model_parametric}} to fit birth-death models, or the functions \code{\link{fit_hbd_pdr_on_grid}}, \code{\link{fit_hbd_pdr_parametric}} and \code{\link{fit_hbd_psr_on_grid}} to fit BD model congruence classes (aka. ``pulled variables'') to a tree.
}
\usage{
reconstruct_past_diversification( times,
                                  diversities,
                                  birth_rates_pc            = NULL,
                                  rarefaction               = NULL,
                                  discovery_fractions       = NULL,
                                  discovery_fraction_slopes = NULL,
                                  max_age                   = NULL,
                                  coalescent                = FALSE,
                                  smoothing_span            = 0,
                                  smoothing_order           = 1)
}
\arguments{
\item{times}{
Numeric vector, listing the times at which diversities are given. Values must be in ascending order.
}
\item{diversities}{
Numeric vector of the same size as \code{times}, listing diversities (coalescent or not) at each time point.
}
\item{birth_rates_pc}{
Numeric vector of the same size as \code{times}, listing known or assumed per-capita birth rates (speciation rates). Can also be of size 1, in which case the same per-capita birth rate is assumed throughout. Alternatively if \code{coalescent==TRUE}, then this vector can also be empty, in which case a constant per-capita birth rate is assumed and estimated from the slope of the coalescent diversities at the last time point. The last alternative is not available when \code{coalescent==FALSE}.
}
\item{rarefaction}{
Numeric between 0 and 1.  Optional rarefaction fraction assumed for the diversities at the very end. Set to 1 to assume no rarefaction was performed.
}
\item{discovery_fractions}{
Numeric array of size Ntimes, listing the fractions of extant lineages represented in the tree over time. Hence, \code{discovery_fraction[t]} is the probability that a lineage at time \code{times[t]} with extant representatives will be represented in the tree. Can be used as an alternative to \code{rarefaction}, for example if discovery of extant species is non-random or phylogenetically biased. Experimental, so leave this \code{NULL} if you don't know what it means.
}
\item{discovery_fraction_slopes}{
Numeric array of size Ntimes, listing the 1st derivative of \code{discovery_fractions} (w.r.t. time) over time. If \code{NULL}, this will be estimated from \code{discovery_fractions} via basic finite differences if needed. Experimental, so leave this \code{NULL} if you don't know what it means.
}
\item{max_age}{
Numeric. Optional maximum distance from the end time to be considered. If \code{NULL} or <=0 or \code{Inf}, all provided time points are considered.
}
\item{coalescent}{
Logical, indicating whether the provided diversities are from a coalescent tree (only including clades with extant representatives) or total diversities (extant species at each time point).
}
\item{smoothing_span}{
Non-negative integer. Optional sliding window size (number of time points) for smoothening the diversities time series via Savitzky-Golay-filter. If <=2, no smoothing is done. Smoothening the time series can reduce the effects of noise on the reconstructed diversity dynamics.
}
\item{smoothing_order}{
Integer between 1 and 4. Polynomial order of the Savitzky-Golay smoothing filter to be applied. Only relevant if \code{smoothing_span>2}. A value of 1 or 2 is typically recommended.
}
}


\details{
This function can be used to fit a birth-death model to a coalescent diversity time series \eqn{N_c(\tau)} at various ages \eqn{\tau}, also known as ``lineages-through-time'' curve. The reconstruction of the total diversity \eqn{N(\tau)} is based on the following formulas:
\deqn{
E(\tau)=1+\frac{\nu(\tau)}{\beta(\tau)},\\
}
\deqn{
N(\tau)=\frac{N_c}{1-E(\tau)},
}
\deqn{
\nu(\tau)=\frac{1}{N_c(\tau)}\frac{dN_c(\tau)}{d\tau}
}
where \eqn{E(\tau)} is the probability that a clade of size 1 at age \eqn{\tau} went extinct by the end of the time series and \eqn{\beta} is the per-capita birth rate. If the per-capita birth rate is not explicitly provided for each time point (see argument \code{birth_rate_pc}), the function assumes that the per-capita birth rate (speciation rate) is constant at all times. If \code{birth_rates_pc==NULL} and \code{coalescent==TRUE}, the constant speciation rate is estimated as
\deqn{
\beta = -\frac{\nu(0)}{\rho},
}
where \eqn{\rho} is the fraction of species kept after rarefaction (see argument \code{rarefaction}). 

Assuming a constant speciation rate may or may not result in accurate estimates of past total diversities and other quantities. If a time-varying speciation rate is suspected but not known, additional information on past diversification dynamics may be obtained using modified (``pulled'') quantities that partly resemble the classical extinction rate, diversification rate and total diversity. Such quantities are the ``pulled diversification rate'':
\deqn{
\eta(\tau) = \delta(\tau) - \beta(\tau) + \frac{1}{\beta(\tau)}\frac{d\beta}{d\tau},
}
the ``pulled extinction rate'':
\deqn{
\delta_p(\tau) = \delta(\tau) + (\beta_o-\beta(\tau)) - \frac{1}{\beta(\tau)}\frac{d\beta}{d\tau},
}
and the ``pulled total diversity'':
\deqn{
N_p(\tau) = N(\tau)\cdot\frac{\beta_o}{\beta(\tau)},
}
where \eqn{\beta_o} is the provided or estimated (if not provided) speciation rate at the last time point. The advantage of these quantities is that they can be estimated from the coalescent diversities (lineages-through-time) without any assumptions on how \eqn{\beta} and \eqn{\delta} varied over time. The disadvantage is that they differ from their ``non-pulled'' quantities (\eqn{\beta-\delta}, \eqn{\delta} and \eqn{N}), in cases where \eqn{\beta} varied over time.
}


\value{
A named list with the following elements:
\item{success}{
Logical, specifying whether the reconstruction was successful. If \code{FALSE}, the remaining elements may not be defined.
}
\item{Ntimes}{
Integer. Number of time points for which reconstruction is returned.
}
\item{total_diversities}{
Numeric vector of the same size as \code{times}, listing the total diversity at each time point (number of extant lineages at each time point). If \code{coalescent==FALSE}, then these are the same as the \code{diversities} passed to the function.
}
\item{coalescent_diversities}{
Numeric vector of the same size as \code{times}, listing the coalescent diversities at each time point (number of species with at least one extant descendant at the last time point). If \code{coalescent==TRUE}, then these are the same as the \code{diversities} passed to the function.
}
\item{birth_rates}{
Numeric vector of the same size as \code{times}, listing the estimated birth rates (speciation events per time unit).
}
\item{death_rates}{
Numeric vector of the same size as \code{times}, listing the estimated death rates (extinction events per time unit).
}
\item{Psurvival}{
Numeric vector of the same size as \code{times}, listing the estimated fraction of lineages at each time point that eventually survive. \code{Psurvival[i]} is the probability that a clade of size 1 at time \code{times[i]} will be extant by the end of the time series. May be \code{NULL} in some cases.
}
\item{Pdiscovery}{
Numeric vector of the same size as \code{times}, listing the estimated fraction of lineages at each time point that are eventually discovered, provided that they survive. \code{Pdiscovery[i]} is the probability that a clade of size 1 at time \code{times[i]} that is extant by the end of the time series, will be discovered. May be \code{NULL} in some cases.
}
\item{Prepresentation}{
Numeric vector of the same size as \code{times}, listing the estimated fraction of lineages at each time point that eventually survive and are discovered. \code{Prepresentation[i]} is the probability that a clade of size 1 at time \code{times[i]} will be extant by the end of the time series and visible in the coalescent tree after rarefaction. Note that Prepresentation = Psurvival * Pdiscovery. May be \code{NULL} in some cases.
}
\item{total_births}{
Numeric, giving the estimated total number of birth events that occurred between times \code{T-max_age} and \code{T}, where \code{T} is the last time point of the time series.
}
\item{total_deaths}{
Numeric, giving the estimated total number of death events that occurred between times \code{T-max_age} and \code{T}, where \code{T} is the last time point of the time series.
}
\item{last_birth_rate_pc}{
The provided or estimated (if not provided) speciation rate at the last time point. This corresponds to the birth rate divided by the estimated true diversity (prior to rarefaction) at the last time point.
}
\item{last_death_rate_pc}{
The estimated extinction rate at the last time point. This corresponds to the death rate divided by the estimated true diversity (prior to rarefaction) at the last time point.
}
\item{pulled_diversification_rates}{
Numeric vector of the same size as \code{times}, listing the estimated pulled diversification rates.
}
\item{pulled_extinction_rates}{
Numeric vector of the same size as \code{times}, listing the estimated pulled extinction rates.
}
\item{pulled_total_diversities}{
Numeric vector of the same size as \code{times}, listing the estimated pulled total diversities.
}
}

\author{Stilianos Louca}

\references{
Louca et al (2018). Bacterial diversification through geological time. Nature Ecology & Evolution. 2:1458-1467.
}

\seealso{
\code{\link{generate_random_tree}},
\code{\link{fit_tree_model}},
\code{\link{count_lineages_through_time}},
\code{\link{fit_hbd_model_parametric}},
\code{\link{fit_hbd_model_on_grid}}
}

\examples{
#####################################################
# EXAMPLE 1

# Generate a coalescent tree
params = list(birth_rate_intercept  = 0, 
              birth_rate_factor     = 1,
              birth_rate_exponent   = 1,
              death_rate_intercept  = 0,
              death_rate_factor     = 0.05,
              death_rate_exponent   = 1.3,
              rarefaction           = 1)
simulation = generate_random_tree(params,max_time_eq=1,coalescent=TRUE)
tree = simulation$tree
time_span = simulation$final_time - simulation$root_time
cat(sprintf("Generated tree has \%d tips, spans \%g time units\n",length(tree$tip.label),time_span))

# Calculate diversity time series from the tree
counter = count_lineages_through_time(tree, times=seq(0,0.99*time_span,length.out=100))

# print coalescent diversities
print(counter$lineages)

# reconstruct diversification dynamics based on diversity time series
results = reconstruct_past_diversification( counter$times,
                                            counter$lineages,
                                            coalescent      = TRUE,
                                            smoothing_span  = 3,
                                            smoothing_order = 1)
                                            
# print reconstructed total diversities
print(results$total_diversities)
                                                  
# plot coalescent and reconstructed true diversities
matplot(x     = counter$times, 
        y     = matrix(c(counter$lineages,results$total_diversities), ncol=2, byrow=FALSE),
        type  = "b", 
        xlab  = "time", 
        ylab  = "# clades",
        lty   = c(1,2), pch = c(1,0), col = c("red","blue"))
legend( "topleft", 
        legend  = c("coalescent (simulated)","true (reconstructed)"), 
        col     = c("red","blue"), lty = c(1,2), pch = c(1,0));
        
        
        
#####################################################
# EXAMPLE 2

# Generate a non-coalescent tree
params = list(birth_rate_intercept  = 0, 
              birth_rate_factor     = 1,
              birth_rate_exponent   = 1,
              death_rate_intercept  = 0,
              death_rate_factor     = 0.05,
              death_rate_exponent   = 1.3,
              rarefaction           = 1)
simulation = generate_random_tree(params,max_time_eq=1,coalescent=FALSE)
tree = simulation$tree
time_span = simulation$final_time - simulation$root_time
cat(sprintf("Generated tree has \%d tips, spans \%g time units\n",length(tree$tip.label),time_span))

# Calculate diversity time series from the tree
counter = count_lineages_through_time(tree, times=seq(0,0.99*time_span,length.out=100))

# print true diversities
print(counter$lineages)

# reconstruct diversification dynamics based on diversity time series
results = reconstruct_past_diversification( counter$times,
                                            counter$lineages,
                                            birth_rates_pc  = params$birth_rate_factor,
                                            coalescent      = FALSE,
                                            smoothing_span  = 3,
                                            smoothing_order = 1)
                                            
# print coalescent diversities
print(results$coalescent_diversities)
                                                  
# plot coalescent and reconstructed true diversities
matplot(x     = counter$times, 
        y     = matrix(c(results$coalescent_diversities,counter$lineages), ncol=2, byrow=FALSE),
        type  = "b", 
        xlab  = "time", 
        ylab  = "# clades",
        lty   = c(1,2), pch = c(1,0), col = c("red","blue"))
legend( "topleft", 
        legend  = c("coalescent (reconstructed)","true (simulated)"), 
        col     = c("red","blue"), lty = c(1,2), pch = c(1,0));
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{tree model}
\keyword{fitting}
