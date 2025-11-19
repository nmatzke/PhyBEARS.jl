\name{date_tree_red}
\alias{date_tree_red}
\title{
Date a tree based on relative evolutionary divergences.
}
\description{
Given a rooted phylogenetic tree and a single node ('anchor') of known age (distance from the present), rescale all edge lengths so that the tree becomes ultrametric and edge lengths correspond to time intervals. The function is based on relative evolutionary divergences (RED), which measure the relative position of each node between the root and its descending tips (Parks et al. 2018). If no anchor node is provided, the root is simply assumed to have age 1. This function provides a heuristic quick-and-dirty way to date a phylogenetic tree.
}
\usage{
date_tree_red(tree, anchor_node = NULL, anchor_age = 1)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{anchor_node}{
Integer, ranging between 1 and Nnodes. Index of the node to be used as dating anchor. If \code{NULL}, the tree's root is used as anchor.
}
\item{anchor_age}{
Positive numeric. Age of the anchor node.
}
}


\details{
The RED of a node measures its relative placement between the root and the node's descending tips (Parks et al. 2018). The root's RED is set to 0. Traversing from root to tips (preorder traversal), for each node the RED is set to \eqn{P+(a/(a+b))\cdot(1-P)}, where \eqn{P} is the RED of the node's parent, \eqn{a} is the edge length connecting the node to its parent, and \eqn{b} is the average distance from the node to its descending tips. The RED of all tips is set to 1.

For each edge, the RED difference between child & parent is used to set the new length of that edge, multiplied by some common scaling factor to translate RED units into time units. The scaling factor is chosen such that the new distance of the anchor node from its descending tips equals \code{anchor_age}. All tips will have age 0. The topology of the dated tree, as well as tip/node/edge indices, remain unchanged.

This function provides a heuristic approach to making a tree ultrametric, and has not been derived from a specific evolutionary model. In particular, its statistical properties are unknown to the author.

The time complexity of this function is O(Nedges).
The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). If \code{tree$edge.length} is \code{NULL}, then all edges in the input tree are assumed to have length 1.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether the dating was successful. If \code{FALSE}, all other return values (except for \code{error}) may be undefined.
}
\item{tree}{
A new rooted tree of class "phylo", representing the dated tree.
}
\item{REDs}{
Numeric vector of size Nnodes, listing the RED of each node in the input tree.
}
\item{error}{
Character, listing any error message if \code{success==FALSE}.
}
}

\author{Stilianos Louca}

\references{
D. H. Parks, M. Chuvochina et al. (2018). A proposal for a standardized bacterial taxonomy based on genome phylogeny. bioRxiv 256800. DOI:10.1101/256800
}

\seealso{
\code{\link{congruent_divergence_times}}
}

\examples{
# generate a random non-ultrametric tree
params = list(birth_rate_intercept=1, death_rate_intercept=0.8)
tree = generate_random_tree(params, max_time=1000, coalescent=FALSE)$tree

# make ultrametric, by setting the root to 2 million years
dated_tree = date_tree_red(tree, anchor_age=2e6)
}
\keyword{dating}
\keyword{ultrametric}
\keyword{relative evolutionary divergence}
