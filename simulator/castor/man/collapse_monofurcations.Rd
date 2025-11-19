\name{collapse_monofurcations}
\alias{collapse_monofurcations}
\title{
Remove monofurcations from a tree.
}
\description{
Eliminate monofurcations (nodes with only a single child) from a phylogenetic tree, by connecting their incoming and outgoing edge.
}
\usage{
collapse_monofurcations(tree, force_keep_root=TRUE, as_edge_counts=FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{force_keep_root}{
Logical, indicating whether the root node should always be kept (i.e., even if it only has a single child).
}
\item{as_edge_counts}{
Logical, indicating whether all edges should be assumed to have length 1. If \code{TRUE}, the outcome is the same as if the tree had no edges.
}
}


\details{
All tips in the input tree retain their original indices, however the returned tree may include fewer nodes and edges. Edge and node indices may change.

If \code{tree$edge.length} is missing, then all edges in the input tree are assumed to have length 1.
}


\value{
A list with the following elements:
\item{tree}{
A new tree of class "phylo", containing only bifurcations (and multifurcations, if these existed in the input tree). The number of nodes in this tree, Nnodes_new, may be lower than of the input tree.
}
\item{new2old_node}{
Integer vector of length Nnodes_new, mapping node indices in the new tree to node indices in the old tree.
}
\item{Nnodes_removed}{
Integer. Number of nodes (monofurcations) removed from the tree.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{multifurcations_to_bifurcations}}
}

\examples{
# generate a random tree
Ntips = 1000
tree = generate_random_tree(list(birth_rate_intercept=1), max_tips=Ntips)$tree

# prune the tree to generate random monofurcations
random_tips = sample.int(n=Ntips, size=0.5 * Ntips, replace=FALSE)
tree = get_subtree_with_tips(tree, only_tips=random_tips, collapse_monofurcations=FALSE)$subtree

# collapse monofurcations
new_tree = collapse_monofurcations(tree)$tree

# print summary of old and new tree
cat(sprintf("Old tree has \%d nodes\n",tree$Nnode))
cat(sprintf("New tree has \%d nodes\n",new_tree$Nnode))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multifurcations}
