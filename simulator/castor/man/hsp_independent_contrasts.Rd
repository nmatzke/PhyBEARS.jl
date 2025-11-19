\name{hsp_independent_contrasts}
\alias{hsp_independent_contrasts}
\title{
Hidden state prediction via phylogenetic independent contrasts.
}
\description{
Reconstruct ancestral states of a continuous (numeric) trait for nodes and predict unknown (hidden) states for tips on a tree using phylogenetic independent contrasts.
}
\usage{
hsp_independent_contrasts(tree, tip_states, weighted=TRUE, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}

\item{tip_states}{
A numeric vector of size Ntips, specifying the state of each tip in the tree. \code{tip_states} can include \code{NA} to indicate an unknown tip state that is to be predicted.
}
\item{weighted}{
Logical, specifying whether to weight transition costs by the inverted edge lengths during ancestral state reconstruction. This corresponds to the "weighted squared-change parsimony" reconstruction by Maddison (1991) for a Brownian motion model of trait evolution.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
Any \code{NA} entries in \code{tip_states} are interpreted as unknown (hidden) states to be estimated. Prior to ancestral state reconstruction, the tree is temporarily prunned, keeping only tips with known state. The function then uses a postorder traversal algorithm to calculate the intermediate "X" variables (a state estimate for each node) introduced by Felsenstein (1985) in his phylogenetic independent contrasts method. Note that these are only local estimates, i.e. for each node the estimate is only based on the tip states in the subtree descending from that node (see discussion in Garland and Ives, 2000). The states of tips with hidden state are set to those of the most recent ancestor with reconstructed state, as described by Zaneveld and Thurber (2014).

This function has asymptotic time complexity O(Nedges).
If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. This is the same as setting \code{weighted=FALSE}. The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). 

Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. The vector \code{tip_states} need not include item names; if it does, however, they are checked for consistency (if \code{check_input==TRUE}).

This function is meant for reconstructing ancestral states in all nodes of a tree as well as predicting the states of tips with an a priory unknown state. If the state of all tips is known and only ancestral state reconstruction is needed, consider using the function \code{\link{asr_independent_contrasts}} for improved efficiency.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether HSP was successful. If \code{FALSE}, some return values may be \code{NULL}.
}
\item{states}{
A numeric vector of size Ntips+Nnodes, listing the reconstructed state of each tip and node. The entries in this vector will be in the order in which tips and nodes are indexed in \code{tree$edge}.
}
\item{total_sum_of_squared_changes}{
The total sum of squared changes in tree, minimized by the (optionally weighted) squared-change parsimony algorithm. This is equation 7 in (Maddison, 1991). Note that for the root, phylogenetic independent contrasts is equivalent to Maddison's squared-change parsimony.
}
}


\author{Stilianos Louca}

\references{
J. Felsenstein (1985). Phylogenies and the comparative method. The American Naturalist. 125:1-15.

T. Jr. Garland and A. R. Ives (2000). Using the past to predict the present: Confidence intervals for regression equations in phylogenetic comparative methods. The American Naturalist. 155:346-364.

W. P. Maddison (1991). Squared-change parsimony reconstructions of ancestral states for continuous-valued characters on a phylogenetic tree. Systematic Zoology. 40:304-314.

J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
}

\seealso{
\code{\link{asr_squared_change_parsimony}}
\code{\link{hsp_max_parsimony}}, 
\code{\link{hsp_mk_model}}, 
}

\examples{
# generate random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# simulate a continuous trait
tip_states = simulate_ou_model(tree, stationary_mean=0, spread=1, decay_rate=0.001)$tip_states

# print tip states
print(as.vector(tip_states))

# set half of the tips to unknown state
tip_states[sample.int(Ntips,size=as.integer(Ntips/2),replace=FALSE)] = NA

# reconstruct all tip states via weighted PIC
estimated_states = hsp_independent_contrasts(tree, tip_states, weighted=TRUE)$states

# print estimated tip states
print(estimated_states[1:Ntips])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum parsimony}
\keyword{ancestral state reconstruction}
\keyword{hidden state prediction}
\keyword{BM model}
