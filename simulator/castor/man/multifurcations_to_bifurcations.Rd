\name{multifurcations_to_bifurcations}
\alias{multifurcations_to_bifurcations}
\title{
Expand multifurcations to bifurcations.
}
\description{
Eliminate multifurcations from a phylogenetic tree, by replacing each multifurcation with multiple bifurcations.
}
\usage{
multifurcations_to_bifurcations(tree, dummy_edge_length=0, 
                                new_node_basename="node.", 
                                new_node_start_index=NULL)
}
\arguments{
\item{tree}{
A tree of class "phylo".
}
\item{dummy_edge_length}{
Non-negative numeric. Length to be used for new (dummy) edges when breaking multifurcations into bifurcations. Typically this will be 0, but can also be a positive number if zero edge lengths are not desired in the returned tree.
}
\item{new_node_basename}{
Character. Name prefix to be used for added nodes (e.g. "node." or "new.node."). Only relevant if the input tree included node labels.
}
\item{new_node_start_index}{
Integer. First index for naming added nodes. Can also be \code{NULL}, in which case this is set to Nnodes+1, where Nnodes is the number of nodes in the input tree.
}
}


\details{
For each multifurcating node (i.e. with more than 2 children), all children but one will be placed on new bifurcating nodes, connected to the original node through one or more dummy edges.

The input tree need not be rooted, however descendance from each node is inferred based on the direction of edges in \code{tree$edge}. The input tree may include multifurcations (i.e. nodes with more than 2 children) as well as monofurcations (i.e. nodes with only one child). Monofurcations are kept in the returned tree. 

All tips and nodes in the input tree retain their original indices, however the returned tree may include additional nodes and edges. Edge indices may change.

If \code{tree$edge.length} is missing, then all edges in the input tree are assumed to have length 1. The returned tree will include \code{edge.length}, with all new edges having length equal to \code{dummy_edge_length}.
}


\value{
A list with the following elements:
\item{tree}{
A new tree of class "phylo", containing only bifurcations (and monofurcations, if these existed in the input tree). 
}
\item{old2new_edge}{
Integer vector of length Nedges, mapping edge indices in the old tree to edge indices in the new tree.
}
\item{Nnodes_added}{
Integer. Number of nodes added to the new tree.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{collapse_monofurcations}}
}

\examples{
# generate a random multifurcating tree
Ntips = 1000
tree = generate_random_tree(list(birth_rate_intercept=1), Ntips, Nsplits=5)$tree

# expand multifurcations to bifurcations
new_tree = multifurcations_to_bifurcations(tree)$tree

# print summary of old and new tree
cat(sprintf("Old tree has \%d nodes\n",tree$Nnode))
cat(sprintf("New tree has \%d nodes\n",new_tree$Nnode))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multifurcations}
