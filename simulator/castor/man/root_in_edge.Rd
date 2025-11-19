\name{root_in_edge}
\alias{root_in_edge}
\title{
Root or re-root a tree in the middle of an edge.
}
\description{
Given a tree (rooted or unrooted), this function places the new root on some specified edge, effectively adding one more node, one more edge and changing the direction of edges as required.
}
\usage{
root_in_edge( tree, 
              root_edge, 
              location                = 0.5,
              new_root_name           = NULL, 
              collapse_monofurcations = TRUE)
}
\arguments{
\item{tree}{
A tree object of class "phylo". Can be unrooted or rooted.
}
\item{root_edge}{
Integer, index of the edge into which the new root is to be placed. Must be between 1 and Nedges.
}
\item{location}{
Numeric, between 0 and 1, specifying the relative location along the \code{root_edge} at which to place the root (relative to the edge length, measured from the upstream node). For example, \code{location=0.5} means the root is placed in the middle of the edge, while \code{location=0.1} means that it will be place closer to the upstream node (i.e., closer to \code{tree$edge[root_edge,1]}).
}
\item{new_root_name}{
Character, optional, specifying the node name to use for the new root. Only used if \code{tree$node.label} is not \code{NULL}.
}
\item{collapse_monofurcations}{
Logical, specifying whether monofurcations in the rerooted tree (e.g. stemming from the old root) should be collapsed by connecting incoming edges with outgoing edges.
}
}

\details{
The input tree may include an arbitrary number of incoming and outgoing edges per node (but only one edge per tip), and the direction of these edges can be arbitrary. Of course, the undirected graph defined by all edges must still be a valid tree.

The number of tips in the rerooted tree remains unchanged, the number of nodes is increased by 1.
Node indices may be modified. Tip indices always remain unchanged.

The asymptotic time complexity of this function is O(Nedges).
}

\value{
A tree object of class "phylo", representing the (re-)rooted phylogenetic tree. The element \code{tip.label} is the same as for the input tree, but all other elements may have changed.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{root_via_outgroup}},
\code{\link{root_at_node}},
\code{\link{root_at_midpoint}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# reroot the tree inside some arbitrary edge
focal_edge = 120
tree = root_in_edge(tree, focal_edge)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Rerooting }
