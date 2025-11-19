\name{asr_squared_change_parsimony}
\alias{asr_squared_change_parsimony}
\title{
Squared-change parsimony ancestral state reconstruction.
}
\description{
Reconstruct ancestral states for a continuous (numeric) trait using squared-change maximum parsimony (Maddison, 1991). Transition costs can optionally be weighted by the inverse edge lengths ("weighted squared-change parsimony" by Maddison).
}
\usage{
asr_squared_change_parsimony(tree, tip_states,  weighted=TRUE, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
A numeric vector of size Ntips, specifying the known state of each tip in the tree.
}
\item{weighted}{
Logical, specifying whether to weight transition costs by the inverted edge lengths. This corresponds to the "weighted squared-change parsimony" reconstruction by Maddison (1991) for a Brownian motion model of trait evolution.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
The function traverses the tree in postorder (tips-->root) to calculate the quadratic parameters described by Maddison (1991) and obtain the globally parsimonious squared-change parsimony state for the root. The function then reroots at each node, updates all affected quadratic parameters in the tree and calculates the node's globally parsimonious squared-change parsimony state.
The function has asymptotic time complexity O(Nedges).

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. This is the same as setting \code{weighted=FALSE}. The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). Edges with length 0 will be adjusted internally to some tiny length if needed (if \code{weighted==TRUE}).

Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. The vector \code{tip_states} need not include item names; if it does, however, they are checked for consistency (if \code{check_input==TRUE}).

If \code{weighted==FALSE}, then this function yields the same ancestral state reconstructions as

\code{ape::ace(tip_states, tree, type="continuous", method="ML", model="BM", CI=FALSE)} 

in the \code{ape} package (v. 0.5-64), assuming the tree as unit edge lengths. If \code{weighted==TRUE}, then this function yields the same ancestral state reconstructions as the maximum likelihood estimates under a Brownian motion model, as implemented by the \code{Rphylopars} package (v. 0.2.10):

\code{Rphylopars::anc.recon(tip_states, tree, vars=FALSE, CI=FALSE).}

}


\value{
A list with the following elements:
\item{ancestral_states}{
A numeric vector of size Nnodes, listing the reconstructed state of each node. The entries in this vector will be in the order in which nodes are indexed in the tree.
}
\item{total_sum_of_squared_changes}{
The total sum of squared changes, minimized by the (optionally weighted) squared-change parsimony algorithm. This is equation 7 in (Maddison, 1991).
}
}


\author{Stilianos Louca}

\references{
W. P. Maddison (1991). Squared-change parsimony reconstructions of ancestral states for continuous-valued characters on a phylogenetic tree. Systematic Zoology. 40:304-314.
}

\seealso{
\code{\link{asr_independent_contrasts}}
\code{\link{asr_max_parsimony}},
\code{\link{asr_mk_model}}
}

\examples{
# generate random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# simulate a continuous trait
tip_states = simulate_ou_model(tree, stationary_mean=0, spread=1, decay_rate=0.001)$tip_states

# reconstruct node states based on simulated tip states
node_states = asr_squared_change_parsimony(tree, tip_states, weighted=TRUE)$ancestral_states
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum parsimony}
\keyword{ancestral state reconstruction}
\keyword{BM model}
