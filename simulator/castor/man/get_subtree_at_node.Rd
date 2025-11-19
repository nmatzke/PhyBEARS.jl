\name{get_subtree_at_node}
\alias{get_subtree_at_node}
\title{
Extract a subtree descending from a specific node.
}
\description{
Given a tree and a focal node, extract the subtree descending from the focal node and place the focal node as the root of the extracted subtree.
}
\usage{
get_subtree_at_node(tree, node)
}
\arguments{
\item{tree}{
A tree of class "phylo".
}
\item{node}{
Character or integer specifying the name or index, respectively, of the focal node at which to extract the subtree. If an integer, it must be between 1 and \code{tree$Nnode}. If a character, it must be a valid entry in \code{tree$node.label}.
}
}


\details{
The input tree need not be rooted, however "descendance" from the focal node is inferred based on the direction of edges in \code{tree$edge}. The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child).
}


\value{
A list with the following elements:
\item{subtree}{
A new tree of class "phylo", containing the subtree descending from the focal node. This tree will be rooted, with the new root being the focal node.
}
\item{new2old_tip}{
Integer vector of length Ntips_kept (=number of tips in the extracted subtree) with values in 1,..,Ntips, mapping tip indices of the subtree to tip indices in the original tree. In particular, \code{tree$tip.label[new2old_tip]} will be equal to \code{subtree$tip.label}.
}
\item{new2old_node}{
Integer vector of length Nnodes_kept (=number of nodes in the extracted subtree) with values in 1,..,Nnodes, mapping node indices of the subtree to node indices in the original tree. 

For example, \code{new2old_node[1]} is the index that the first node of the subtree had within the original tree. In particular, \code{tree$node.label[new2old_node]} will be equal to \code{subtree$node.label} (if node labels are available).
}
\item{new2old_edge}{
Integer vector of length Nedges_kept (=number of edges in the extracted subtree), with values in 1,..,Nedges, mapping edge indices of the subtree to edge indices in the original tree. In particular, \code{tree$edge.length[new2old_edge]} will be equal to \code{subtree$edge.length} (if edge lengths are available).
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_subtree_with_tips}}
}

\examples{
# generate a random tree
Ntips = 1000
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# extract subtree descending from a random node
node = sample.int(tree$Nnode,size=1)
subtree = get_subtree_at_node(tree, node)$subtree

# print summary of subtree
cat(sprintf("Subtree at \%d-th node has \%d tips\n",node,length(subtree$tip.label)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{pruning}
\keyword{subtree}
