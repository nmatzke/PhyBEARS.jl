\name{merge_nodes_to_multifurcations}
\alias{merge_nodes_to_multifurcations}
\title{
Merge specific nodes into multifurcations.
}
\description{
Given a rooted tree, merge one or more nodes ``upwards'' into their parent nodes, thus effectively generating multifurcations. Multiple generations of nodes (i.e., successive branching points) can be merged into a single "absorbing ancestor".
}
\usage{
merge_nodes_to_multifurcations( tree, 
                                nodes_to_merge,
                                merge_with_parents  = FALSE,
                                keep_ancestral_ages = FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{nodes_to_merge}{
Integer vector or character vector, listing nodes in the tree that should be merged with their parents (if \code{merge_with_parents=TRUE}) or with their children (if \code{merge_with_parents=FALSE}). If an integer vector, it must contain values in 1,..,Nnodes. If a character vector, it must list node labels, and the tree itself must also include node labels.
}
\item{merge_with_parents}{
Logical, specifying whether the nodes listed in \code{nodes_to_merge} should be merged with their parents. If \code{FALSE}, the specified nodes will be merged with their children (whenever these are not tips).
}
\item{keep_ancestral_ages}{
Logical, specifying whether the generated multifurcations should have the same age as the absorbing ancestor. If \code{FALSE}, then the age of a multifurcation will be the average of the absorbing ancestor's age and the ages of its merged child nodes (but constrained from below by the ages of non-merged descendants to avoid negative edge lengths). If \code{TRUE}, then the ages of multifurcations will be biased towards the root, since their age will be that of the absorbing ancestor.
}
}


\details{
All tips in the input tree are kept and retain their original indices, however the returned tree will include fewer nodes and edges. Edge and node indices may change. When a node is merged into its parent, the incoming edge is lost, and the parent's age remains unchanged. 

Nodes are merged in order from root to tips. Hence, if a node B is merged into ("absorbed by") its parent node A, and child node C is merged into node B, then effectively C ends up merged into node A (node A is the "absorbing ancestor").

If \code{tree$edge.length} is missing, then all edges in the input tree are assumed to have length 1.
}


\value{
A list with the following elements:
\item{tree}{
A new tree of class "phylo". The number of nodes in this tree, Nnodes_new, will generally be lower than of the input tree.
}
\item{new2old_node}{
Integer vector of length Nnodes_new, mapping node indices in the new tree to node indices in the old tree. Note that nodes merged with their parents are not represented in this list.
}
\item{old2new_node}{
Integer vector of length Nnodes, mapping node indices in the old tree to node indices in the new tree. Nodes merged with their parents (and thus missing from the new tree) will have value 0.
}
\item{Nnodes_removed}{
Integer. Number of nodes removed from the tree, due to being merged into their parents.
}
\item{Nedges_removed}{
Integer. Number of edges removed from the tree.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{multifurcations_to_bifurcations}}, 
\code{\link{collapse_monofurcations}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1), max_tips=Ntips)$tree

# merge a few nodes with their parents,
# thus obtaining a multifurcating tree
nodes_to_merge = c(1,3,4)
new_tree = merge_nodes_to_multifurcations(tree, nodes_to_merge)$tree

# print summary of old and new tree
cat(sprintf("Old tree has \%d nodes\n",tree$Nnode))
cat(sprintf("New tree has \%d nodes\n",new_tree$Nnode))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multifurcations}
