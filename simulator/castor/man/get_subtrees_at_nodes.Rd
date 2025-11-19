\name{get_subtrees_at_nodes}
\alias{get_subtrees_at_nodes}
\title{
Extract subtrees descending from specific nodes.
}
\description{
Given a tree and a list of focal nodes, extract the subtrees descending from those focal nodes, with the focal nodes becoming the roots of the extracted subtrees.
}
\usage{
get_subtrees_at_nodes(tree, nodes)
}
\arguments{
\item{tree}{
A tree of class "phylo".
}
\item{nodes}{
Character vector or integer vector specifying the names or indices, respectively, of the focal nodes at which to extract the subtrees. If an integer vector, entries must be between 1 and \code{tree$Nnode}. If a character vector, each entry must be a valid entry in \code{tree$node.label}.
}
}


\details{
The input tree need not be rooted, however "descendance" from a focal node is inferred based on the direction of edges in \code{tree$edge}. The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child).
}


\value{
A list with the following elements:
\item{subtrees}{
List of the same length as \code{nodes}, with each element being a new tree of class "phylo", containing the subtrees descending from the focal nodes. Each subtree will be rooted at the corresponding focal node.
}
\item{new2old_tip}{
List of the same length as \code{nodes}, with the n-th element being an integer vector with values in 1,..,Ntips, mapping tip indices of the n-th subtree to tip indices in the original tree. In particular, \code{tree$tip.label[new2old_tip[[n]]]} will be equal to \code{subtrees[[n]]$tip.label}.
}
\item{new2old_node}{
List of the same length as \code{nodes}, with the n-th element being an integer vector with values in 1,..,Nnodes, mapping node indices of the n-th subtree to node indices in the original tree. 

For example, \code{new2old_node[[2]][1]} is the index that the 1st node of the 2nd subtree had within the original tree. In particular, \code{tree$node.label[new2old_node[[n]]]} will be equal to \code{subtrees[[n]]$node.label} (if node labels are available).
}
\item{new2old_edge}{
List of the same length as \code{nodes}, with the n-th element being an integer vector with values in 1,..,Nedges, mapping edge indices of the n-th subtree to edge indices in the original tree. In particular, \code{tree$edge.length[new2old_edge[[n]]]} will be equal to \code{subtrees[[n]]$edge.length} (if edge lengths are available).
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_subtree_at_node}},

\code{\link{get_subtree_with_tips}}
}

\examples{
# generate a random tree
Ntips = 1000
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# extract subtrees descending from random nodes
nodes = sample.int(tree$Nnode,size=10)
subtrees = get_subtrees_at_nodes(tree, nodes)$subtrees

# print summaries of extracted subtrees
for(n in length(nodes)){
  cat(sprintf("Subtree at \%d-th node has \%d tips\n",nodes[n],length(subtrees[[n]]$tip.label)))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{pruning}
\keyword{subtree}
