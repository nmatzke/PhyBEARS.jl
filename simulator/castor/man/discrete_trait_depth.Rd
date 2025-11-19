\name{discrete_trait_depth}
\alias{discrete_trait_depth}
\title{
Calculate phylogenetic depth of a discrete trait.
}
\description{
Given a rooted phylogenetic tree and the state of a discrete trait at each tip, calculate the mean phylogenetic depth at which the trait is conserved across clades, using a modification of the consenTRAIT metric introduced by Martiny et al (2013). This is the mean depth of clades that are "maximally uniform" in the trait (see below for details).
}
\usage{
discrete_trait_depth(tree, 
                     tip_states, 
                     min_fraction         = 0.9, 
                     count_singletons     = TRUE,
                     singleton_resolution = 0,
                     weighted             = FALSE, 
                     Npermutations        = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
A vector of size Ntips specifying the state at each tip. Note that \code{tip_states[i]} (where i is an integer index) must correspond to the i-th tip in the tree. This vector may be of any base data type, although character or integer are the most typical types.
}
\item{min_fraction}{
Minimum fraction of tips in a clade that must have the dominant state, for the clade to be considered "uniform" in the trait.
}
\item{count_singletons}{
Logical, specifying whether to consider singleton clades in the statistics (e.g., tips not part of a larger uniform clade). The phylogenetic depth of singletons is taken to be half the length of their incoming edge, as proposed by Martiny et al (2013). If FALSE, singletons are ignored. If you suspect a high risk of false positives in the detection of a trait, it may be worth setting \code{count_singletons} to \code{FALSE} to avoid skewing the distribution of conservation depths towards shallower depths due to false positives.
}
\item{singleton_resolution}{
Numeric, specifying the phylogenetic resolution at which to resolve singletons. A clade will be considered a singleton if the distance of the clade's root to all descending tips is below this threshold.
}
\item{weighted}{
Whether to weight uniform clades by their number of tips in the dominant state. If \code{FALSE}, each uniform clades is weighted equally.
}
\item{Npermutations}{
Number of random permutations for estimating the statistical significance of the mean trait depth. If zero (default), the statistical significance is not calculated.
}
}
\details{
The depth of a clade is defined as the average distance of its tips to the clade's root. The "dominant" state of a clade is defined as the most frequent state among all of the clade's tips. A clade is considered "uniform" in the trait if the frequency of its dominant state is equal to or greater than \code{min_fraction}. The clade is "maximally uniform" if it is uniform and not descending from another uniform clade. The mean depth of the trait is defined as the average phylogenetic depth of all considered maximal uniform clades (whether a maximally uniform clade is considered in this statistic depends on \code{count_singletons} and \code{singleton_resolution}). A greater mean depth means that the trait tends to be conserved in deeper-rooting clades.

This function implements a modification of the "consenTRAIT" metric proposed by Martiny et al. (2013) for measuring the mean phylogenetic depth at which a binary trait is conserved across clades. Note that the original consenTRAIT metric by Martiny et al. (2013) does not treat the two states of a binary trait ("presence" and "absence") equally, whereas the function \code{discrete_trait_depth} does. If you want the original consenTRAIT metric for a binary trait, see the function \code{\link{consentrait_depth}}.

The statistical significance of the calculated mean depth, i.e. the probability of encountering such a mean dept or higher by chance, is estimated based on a null model in which each tip is re-assigned a state by randomly reshuffling the original \code{tip_states}. A low P value indicates that the trait exhibits a phylogenetic signal, whereas a high P value means that there is insufficient evidence in the data to suggest a phylogenetic signal (i.e., the trait's phylogenetic conservatism is indistinguishable from the null model of zero conservatism).

The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). If \code{tree$edge.length} is missing, then every edge is assumed to have length 1.
}
\value{
A list with the following elements:
\item{unique_states}{
Vector of the same type as \code{tip_states} and of length \code{Nstates}, listing the unique possible states of the trait.
}
\item{mean_depth}{
Numeric, specifying the mean phylogenetic depth of the trait, i.e., the mean depth of considered maximally uniform clades.
}
\item{var_depth}{
Numeric, specifying the variance of phylogenetic depths of considered maximally uniform clades.
}
\item{min_depth}{
Numeric, specifying the minimum phylogenetic depth of considered maximally uniform clades.
}
\item{max_depth}{
Numeric, specifying the maximum phylogenetic depth of considered maximally uniform clades.
}
\item{Nmax_uniform}{
Number of considered maximal uniform clades.
}
\item{mean_depth_per_state}{
Numeric vector of size Nstates. Mean depth of considered maximally uniform clades, separately for each state and in the same order as \code{unique_states}. Hence, \code{mean_depth_per_state[s]} lists the mean depth of considered maximally uniform clades whose dominant state is \code{unique_states[s]}.
}
\item{var_depth_per_state}{
Numeric vector of size Nstates. Variance of depths of considered maximally uniform clades, separately for each state and in the same order as \code{unique_states}
}
\item{min_depth_per_state}{
Numeric vector of size Nstates. Minimum phylogenetic depth of considered maximally uniform clades, separately for each state and in the same order as \code{unique_states}
}
\item{max_depth_per_state}{
Numeric vector of size Nstates. Maximum phylogenetic depth of considered maximally uniform clades, separately for each state and in the same order as \code{unique_states}
}
\item{Nmax_uniform_per_state}{
Integer vector of size Nstates. Number of considered maximally uniform clades, seperately for each state and in the same order as \code{unique_states}
}
\item{P}{
Statistical significance (P-value) of mean_depth, under a null model of random tip states (see details above). This is the probability that, under the null model, the \code{mean_depth} would be at least as high as observed in the data.
}
\item{mean_random_depth}{
Mean random mean_depth, under the null model of random tip states (see details above).
}
}


\author{Stilianos Louca}

\references{
A. C. Martiny, K. Treseder and G. Pusch (2013). Phylogenetic trait conservatism of functional traits in microorganisms. ISME Journal. 7:830-838.
}

\seealso{
\code{\link{get_trait_acf}}, 
\code{\link{consentrait_depth}}
}

\examples{\dontrun{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=1000)$tree

# simulate discrete trait evolution on the tree
# consider a trait with 3 discrete states
Q = get_random_mk_transition_matrix(Nstates=3, rate_model="ARD", max_rate=0.1)
tip_states = simulate_mk_model(tree, Q)$tip_states

# calculate phylogenetic conservatism of trait
results = discrete_trait_depth(tree, tip_states, count_singletons=FALSE, weighted=TRUE)
cat(sprintf("Mean depth = \%g, std = \%g\n",results$mean_depth,sqrt(results$var_depth)))
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{trait evolution}
