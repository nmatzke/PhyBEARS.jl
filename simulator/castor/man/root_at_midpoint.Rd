\name{root_at_midpoint}
\alias{root_at_midpoint}
\title{
Root or re-root a tree at the midpoint node.
}
\description{
Given a tree (rooted or unrooted), this function changes the direction of edges (\code{tree$edge}) such that the new root satisfies a "midpoint"" criterion. The number of tips and the number of nodes remain unchanged. The root can either be placed on one of the existing nodes (this node will be the one whose maximum distance to any tip is minimized) or in the middle of one of the existing edges (chosen to be in the middle of the longest path between any two tips).
}
\usage{
root_at_midpoint( tree, 
                  split_edge      = FALSE,
                  update_indices  = TRUE,
                  as_edge_counts  = FALSE,
                  is_rooted       = FALSE)
}
\arguments{
\item{tree}{
A tree object of class "phylo". Can be unrooted or rooted (but see option \code{is_rooted}).
}
\item{split_edge}{
Logical, specifying whether to place the new root in the middle of an edge (in the middle of the longest path of any two tips), thereby creating a new node. If \code{FALSE}, then the root will be placed on one of the existing nodes; note that the resulting tree may no longer be bifurcating at the root.
}
\item{update_indices}{
Logical, specifying whether to update the node indices such that the new root is the first node in the list, as is common convention. This will modify \code{tree$node.label} (if it exists) and also the node indices listed in \code{tree$edge}. Note that this option is only relevant if \code{split_edge=FALSE}; if \code{split_edge=TRUE} then \code{update_indices} will always be assumed \code{TRUE}.
}
\item{as_edge_counts}{
Logical, specifying whether phylogenetic distances should be measured as cumulative edge counts. This is the same if all edges had length 1.
}
\item{is_rooted}{
Logical, specifying whether the input tree can be assumed to be rooted. If you are not certain that the tree is rooted, set this to \code{FALSE}.
}
}

\details{
The input tree may include an arbitrary number of incoming and outgoing edges per node (but only one edge per tip), and the direction of these edges can be arbitrary. Of course, the undirected graph defined by all edges must still be a valid tree. Only set \code{is_rooted=TRUE} if you are sure that the input tree is rooted.

If \code{update_indices==FALSE} and \code{split_edge=FALSE}, then node indices remain unchanged. If \code{update_indices==TRUE} (default) or \code{split_edge=TRUE}, then node indices are modified such that the new root is the first node (i.e. with index Ntips+1 in \code{edge} and with index 1 in \code{node.label}), as is common convention. Setting \code{update_indices=FALSE} (when \code{split_edge=FALSE}) reduces the computation required for rerooting. Tip indices always remain unchanged.

The asymptotic time complexity of this function is O(Nedges).
}

\value{
A tree object of class "phylo", with the \code{edge} element modified such that the maximum distance of the root to any tip is minimized. The elements \code{tip.label}, \code{edge.length} and \code{root.edge} (if they exist) are the same as for the input tree. If \code{update_indices==FALSE}, then the element \code{node.label} will also remain the same.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{root_via_outgroup}},
\code{\link{root_at_node}},
\code{\link{root_in_edge}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# reroot the tree at its midpoint node
tree = root_at_midpoint(tree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Rerooting }
