\name{consentrait_depth}
\alias{consentrait_depth}
\title{
Calculate phylogenetic depth of a binary trait using the consenTRAIT metric.
}
\description{
Given a rooted phylogenetic tree and presences/absences of a binary trait for each tip, calculate the mean phylogenetic depth at which the trait is conserved across clades, in terms of the consenTRAIT metric introduced by Martiny et al (2013). This is the mean depth of clades that are positive in the trait (i.e. in which a sufficient fraction of tips exhibits the trait).
}
\usage{
consentrait_depth(tree, 
                  tip_states, 
                  min_fraction        = 0.9, 
                  count_singletons    = TRUE,
                  singleton_resolution= 0,
                  weighted            = FALSE, 
                  Npermutations       = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
A numeric vector of size Ntips indicating absence (value <=0) or presence (value >0) of a particular trait at each tip of the tree. Note that tip_states[i] (where i is an integer index) must correspond to the i-th tip in the tree.
}
\item{min_fraction}{
Minimum fraction of tips in a clade exhibiting the trait, for the clade to be considered "positive" in the trait. In the original paper by Martiny et al (2013), this was 0.9.
}
\item{count_singletons}{
Logical, specifying whether to include singletons in the statistics (tips positive in the trait, but not part of a larger positive clade). The phylogenetic depth of singletons is taken to be half the length of their incoming edge, as proposed by Martiny et al (2013). If FALSE, singletons are ignored.
}
\item{singleton_resolution}{
Numeric, specifying the phylogenetic resolution at which to resolve singletons. Any clade found to be positive in a trait will be considered a singleton if the distance of the clade's root to all descending tips is below this threshold.
}
\item{weighted}{
Whether to weight positive clades by their number of positive tips. If FALSE, each positive clades is weighted equally, as proposed by Martiny et al (2013).
}
\item{Npermutations}{
Number of random permutations for estimating the statistical significance of the mean trait depth. If zero (default), the statistical significance is not calculated.
}
}
\details{
This function calculates the "consenTRAIT" metric (or variants thereof) proposed by Martiny et al. (2013) for measuring the mean phylogenetic depth at which a binary trait (e.g. presence/absence of a particular metabolic function) is conserved across clades. A greater mean depth means that the trait tends to be conserved in deeper-rooting clades. In their original paper, Martiny et al. proposed to consider a trait as conserved in a clade (i.e. marking a clade as "positive" in the trait) if at least 90\% of the clade's tips exhibit the trait (i.e. are "positive" in the trait). This fraction can be controlled using the \code{min_fraction} parameter. The depth of a clade is taken as the average distance of its tips to the clade's root. 

Note that the consenTRAIT metric does not treat "presence" and "absence" equally, i.e., if one were to reverse all presences and absences then the consenTRAIT metric will generally have a different value. This is because the focus is on the presence of the trait (e.g., presence of a metabolic function, or presence of a morphological feature).

The default parameters of this function reflect the original choices made by Martiny et al. (2013), however in some cases it may be sensible to adjust them. For example, if you suspect a high risk of false positives in the detection of a trait, it may be worth setting \code{count_singletons} to \code{FALSE} to avoid skewing the distribution of conservation depths towards shallower depths due to false positives.

The statistical significance of the calculated mean depth, i.e. the probability of encountering such a mean dept or higher by chance, is estimated based on a null model in which each tip is re-assigned a presence or absence of the trait by randomly reshuffling the original \code{tip_states}.

The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). If \code{tree$edge.length} is missing, then every edge is assumed to have length 1.
}
\value{
A list with the following elements:
\item{mean_depth}{
Mean phylogenetic depth of clades that are positive in the trait.
}
\item{var_depth}{
Variance of phylogenetic depths of clades that are positive in the trait.
}
\item{min_depth}{
Minimum phylogenetic depth of clades that are positive in the trait.
}
\item{max_depth}{
Maximum phylogenetic depth of clades that are positive in the trait.
}
\item{Npositives}{
Number of clades that are positive in the trait.
}
\item{P}{
Statistical significance (P-value) of mean_depth, under a null model of random trait presences/absences (see details above). This is the probability that, under the null model, the \code{mean_depth} would be at least as high as observed in the data.
}
\item{mean_random_depth}{
Mean random mean_depth, under a null model of random trait presences/absences (see details above).
}
\item{positive_clades}{
Integer vector, listing indices of tips and nodes (from 1 to Ntips+Nnodes) that were found to be positive in the trait and counted towards the statistic.
}
\item{positives_per_clade}{
Integer vector of size Ntips+Nnodes, listing the number of descending tips per clade (tip or node) that were positive in the trait.
}
\item{mean_depth_per_clade}{
Numeric vector of size Ntips+Nnodes, listing the mean phylogenetic depth of each clade (tip or node), i.e. the average distance to all its descending tips.
}
}


\author{Stilianos Louca}

\references{
A. C. Martiny, K. Treseder and G. Pusch (2013). Phylogenetic trait conservatism of functional traits in microorganisms. ISME Journal. 7:830-838.
}

\seealso{
\code{\link{get_trait_acf}}, 
\code{\link{discrete_trait_depth}}
}

\examples{\dontrun{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=1000)$tree

# simulate binary trait evolution on the tree
Q = get_random_mk_transition_matrix(Nstates=2, rate_model="ARD", max_rate=0.1)
tip_states = simulate_mk_model(tree, Q)$tip_states

# change states from 1/2 to 0/1 (presence/absence)
tip_states = tip_states - 1

# calculate phylogenetic conservatism of trait
results = consentrait_depth(tree, tip_states, count_singletons=FALSE, weighted=TRUE)
cat(sprintf("Mean depth = \%g, std = \%g\n",results$mean_depth,sqrt(results$var_depth)))
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{trait evolution}
