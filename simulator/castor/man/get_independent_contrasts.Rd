\name{get_independent_contrasts}
\alias{get_independent_contrasts}
\title{
Phylogenetic independent contrasts for continuous traits.
}
\description{
Calculate phylogenetic independent contrasts (PICs) for one or more continuous traits on a phylogenetic tree, as described by Felsenstein (1985). The trait states are assumed to be known for all tips of the tree. PICs are commonly used to calculate correlations between multiple traits, while accounting for shared evolutionary history at the tips. This function also returns an estimate for the state of the root or, equivalently, the phylogenetically weighted mean of the tip states (Garland et al., 1999).
}
\usage{
get_independent_contrasts(tree, 
                          tip_states,
                          scaled = TRUE, 
                          only_bifurcations = FALSE,
                          check_input = TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
A numeric vector of size Ntips, or a 2D numeric matrix of size Ntips x Ntraits, specifying the numeric state of each trait at each tip in the tree.
}
\item{scaled}{
Logical, specifying whether to divide (standardize) PICs by the square root of their expected variance, as recommended by Felsenstein (1985).
}
\item{only_bifurcations}{
Logical, specifying whether to only calculate PICs for bifurcating nodes. If \code{FALSE}, then multifurcations are temporarily expanded to bifurcations, and an additional PIC is calculated for each created bifurcation. If \code{TRUE}, then multifurcations are not expanded and PICs will not be calculated for them.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
If the tree is bifurcating, then one PIC is returned for each node. If multifurcations are present and \code{only_bifurcations==FALSE}, these are internally expanded to bifurcations and an additional PIC is returned for each such bifurcation. PICs are never returned for monofurcating nodes. Hence, in general the number of returned PICs is the number of bifurcations in the tree, potentially after multifurcations have been expanded to bifurcations (if \code{only_bifurcations==FALSE}).

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. The tree may include multifurcations (i.e. nodes with more than 2 children) as well as monofurcations (i.e. nodes with only one child). Edges with length 0 will be adjusted internally to some tiny length (chosen to be much smaller than the smallest non-zero length).

Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. The vector \code{tip_states} need not include item names; if it does, however, they are checked for consistency (if \code{check_input==TRUE}).

The function has asymptotic time complexity O(Nedges x Ntraits). It is more efficient to calculate PICs of multiple traits with the same function call, than to calculate PICs for each trait separately. For a single trait, this function is equivalent to the function \code{ape::pic}, with the difference that it can handle multifurcating trees.
}


\value{
A list with the following elements:
\item{PICs}{
A numeric vector (if \code{tip_states} is a vector) or a numeric matrix (if \code{tip_states} is a matrix), listing the phylogenetic independent contrasts for each trait and for each bifurcating node (potentially after multifurcations have been expanded). If a matrix, then \code{PICs[:,T]} will list the PICs for the T-th trait.
Note that the order of elements in this vector (or rows, if \code{PICs} is a matrix) is not necesssarily the order of nodes in the tree, and that \code{PICs} may contain fewer or more elements (or rows) than there were nodes in the input tree.
}
\item{distances}{
Numeric vector of the same size as \code{PICs}. The ``evolutionary distances'' (or time) corresponding to the PICs under a Brownian motion model of trait evolution. These roughly correspond to the cumulative edge lengths between sister nodes from which PICs were calculated; hence their units are the same as those of edge lengths. They do not take into account the actual trait values. See Felsenstein (1985) for details.
}
\item{nodes}{
Integer vector of the same size as \code{PICs}, listing the node indices for which PICs are returned. If \code{only_bifurcations==FALSE}, then this vector may contain \code{NA}s, corresponding to temporary nodes created during expansion of multifurcations.

If \code{only_bifurcations==TRUE}, then this vector will only list nodes that were bifurcating in the input tree. In that case, \code{PICs[1]} will correspond to the node with name \code{tree$node.label[nodes[1]]}, whereas \code{PICs[2]} will correspond to the node with name \code{tree$node.label[nodes[2]]}, and so on.
}
\item{root_state}{
Numeric vector of size Ntraits, listing the globally estimated state for the root or, equivalently, the phylogenetically weighted mean of the tip states.
}
\item{root_standard_error}{
Numeric vector of size Ntraits, listing the phylogenetically estimated standard errors of the root state under a Brownian motion model. The standard errors have the same units as the traits and depend both on the tree topology as well as the tip states. Calculated according to the procedure described by Garland et al. (1999, page 377).
}
\item{root_CI95}{
Numeric vector of size Ntraits, listing the radius (half width) of the 95\% confidence interval of the root state. Calculated according to the procedure described by Garland et al. (1999, page 377). Note that in contrast to the CI95 returned by the \code{ace} function in the \code{ape} package (v. 0.5-64), \code{root_CI95} has the same units as the traits and depends both on the tree topology as well as the tip states.
}
}


\author{Stilianos Louca}

\references{
J. Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.

T. Garland Jr., P. E. Midford, A. R. Ives (1999). An introduction to phylogenetically based statistical methods, with a new method for confidence intervals on ancestral values. American Zoologist. 39:374-388. 
}

\seealso{
\code{\link{asr_independent_contrasts}}
}

\examples{
# generate random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# simulate a continuous trait
tip_states = simulate_bm_model(tree, diffusivity=0.1, include_nodes=FALSE)$tip_states;

# calculate PICs
results = get_independent_contrasts(tree, tip_states, scaled=TRUE, only_bifurcations=TRUE)

# assign PICs to the bifurcating nodes in the input tree
PIC_per_node = rep(NA, tree$Nnode)
valids = which(!is.na(results$nodes))
PIC_per_node[results$nodes[valids]] = results$PICs[valids]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BM model}
