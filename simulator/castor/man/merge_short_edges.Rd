\name{merge_short_edges}
\alias{merge_short_edges}
\title{
Eliminate short edges in a tree by merging nodes into multifurcations.
}
\description{
Given a rooted phylogenetic tree and an edge length threshold, merge nodes/tips into multifurcations when their incoming edges are shorter than the threshold. 
}
\usage{
merge_short_edges(tree, 
                  edge_length_epsilon	= 0,
                  force_keep_tips     = TRUE,
                  new_tip_prefix      = "ex.node.tip.")
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{edge_length_epsilon}{
Non-negative numeric, specifying the maximum edge length for an edge to be considered ``short'' and thus to be eliminated. Typically 0 or some small positive number.
}
\item{force_keep_tips}{
Logical. If \code{TRUE}, then tips are always kept, even if their incoming edges are shorter than \code{edge_length_epsilon}. If \code{FALSE}, then tips with short incoming edges are removed from the tree; in that case some nodes may become tips.
}
\item{new_tip_prefix}{
Character or \code{NULL}, specifying the prefix to use for new tip labels stemming from nodes. Only relevant if \code{force_keep_tips==FALSE}. If \code{NULL}, then labels of tips stemming from nodes will be the node labels from the original tree (in this case the original tree should include node labels).
}
}


\details{
The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). Whenever a short edge is eliminated, the edges originating from its child are elongated according to the short edge's length. The corresponding grand-children become children of the short edge's parent. Short edges are eliminated in a depth-first-search manner, i.e. traversing from the root to the tips. 

Note that existing monofurcations are retained. If \code{force_keep_tips==FALSE}, then new monofurcations may also be introduced due to tips being removed.

This function is conceptually similar to the function \code{ape::di2multi}.
}


\value{
A list with the following elements:
\item{tree}{
A new rooted tree of class "phylo", containing the (potentially multifurcating) tree.
}
\item{new2old_clade}{
Integer vector of length equal to the number of tips+nodes in the new tree, with values in 1,..,Ntips+Nnodes, mapping tip/node indices of the new tree to tip/node indices in the original tree.
}
\item{new2old_edge}{
Integer vector of length equal to the number of edges in the new tree, with values in 1,..,Nedges, mapping edge indices of the new tree to edge indices in the original tree.
}
\item{Nedges_removed}{
Integer. Number of edges that have been eliminated.
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
Ntips = 100
tree = generate_random_tree(list(birth_rate_factor=1),max_tips=Ntips)$tree

# set some edge lengths to zero
tree$edge.length[sample.int(n=Ntips, size=10, replace=FALSE)] = 0

# print number of edges
cat(sprintf("Original tree has \%d edges\n",nrow(tree$edge)))

# eliminate any edges of length zero
merged = merge_short_edges(tree, edge_length_epsilon=0)$tree

# print number of edges
cat(sprintf("New tree has \%d edges\n",nrow(merged$edge)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multifurcations}
