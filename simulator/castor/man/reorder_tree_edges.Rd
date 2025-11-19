\name{reorder_tree_edges}
\alias{reorder_tree_edges}
\title{
Reorder tree edges in preorder or postorder.
}
\description{
Given a rooted tree, this function reorders the rows in \code{tree$edge} so that they are listed in preorder (root-->tips) or postorder (tips-->root) traversal.
}
\usage{
reorder_tree_edges(tree, root_to_tips=TRUE, 
                   depth_first_search=TRUE,
                   index_only=FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{root_to_tips}{
Logical, specifying whether to sort edges in preorder traversal (root-->tips), rather than in postorder traversal (tips-->roots).
}
\item{depth_first_search}{
Logical, specifying whether the traversal (or the reversed traversal, if \code{root_to_tips} is \code{FALSE}) should be in depth-first-search format rather than breadth-first-search format.
}
\item{index_only}{
Whether the function should only return a vector listing the reordered row indices of the edge matrix, rather than a modified tree.
}
}

\details{
This function does not change the tree structure, nor does it affect tip/node indices and names. It merely changes the order in which edges are listed in the matrix \code{tree$edge}, so that edges are listed in preorder or postorder traversal. Preorder traversal guarantees that each edge is listed before any of its descending edges. Likewise, postorder guarantees that each edge is listed after any of its descending edges.

With options \code{root_to_tips=TRUE} and \code{depth_first_search=TRUE}, this function is analogous to the function \code{reorder} in the \code{ape} package with option \code{order="cladewise"}.

The tree can include multifurcations (nodes with more than 2 children) as well as monofurcations (nodes with 1 child).
This function has asymptotic time complexity O(Nedges).
}

\value{
If \code{index_only==FALSE}, a tree object of class "phylo", with the rows in \code{edge} reordered such that they are listed in direction root-->tips (if \code{root_to_tips==TRUE}) or tips-->root. The vector \code{tree$edge.length} will also be updated in correspondence. Tip and node indices and names remain unchanged.

If \code{index_only=TRUE}, an integer vector (X) of size Nedges, listing the reordered row indices of \code{tree$edge}, i.e. such that \code{tree$edge[X,]} would be the reordered edge matrix.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_tree_traversal_root_to_tips}}
}

\examples{
\dontrun{
# generate a random tree
tree = generate_random_tree(list(birth_rate_factor=1), max_tips=100)$tree

# get new tree with reordered edges
postorder_tree = reorder_tree_edges(tree, root_to_tips=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{traversal}
