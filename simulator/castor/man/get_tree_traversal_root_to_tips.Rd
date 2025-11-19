\name{get_tree_traversal_root_to_tips}
\alias{get_tree_traversal_root_to_tips}
\title{
Traverse tree from root to tips.
}
\description{
Create data structures for traversing a tree from root to tips, and for efficient retrieval of a node's outgoing edges and children.
}
\usage{
get_tree_traversal_root_to_tips(tree, include_tips)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{include_tips}{
Include tips in the tarversal queue. If FALSE, then only nodes are included in the queue.
}
}
\details{
Many dynamic programming algorithms for phylogenetics involve traversing the tree in a certain direction (root to tips or tips to root), and efficient (O(1) complexity) access to a node's direct children can significantly speed up those algorithms. This function is meant to provide data structures that allow traversing the tree's nodes (and optionally tips) in such an order that each node is traversed prior to its descendants (root-->tips) or such that each node is traversed after its descendants (tips-->root). This function is mainly meant for use in other algorithms, and is probably of little relevance to the average user.

The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child).

The asymptotic time and memory complexity of this function is O(Ntips), where Ntips is the number of tips in the tree.
}
\value{
A list with the following elements:
  \item{queue}{
 An integer vector of size Nnodes (if \code{include_tips} was \code{FALSE}) or of size Nnodes+Ntips (if \code{include_tips} was \code{TRUE}), listing indices of nodes (and optionally tips) in the order root-->tips described above. In particular, \code{queue[1]} will be the index of the tree's root (typically Ntips+1).
}
  \item{edges}{
	An integer vector of size Nedges (\code{=nrow(tree$edge)}), listing indices of edges (corresponding to \code{tree$edge}) such that outgoing edges of the same node are listed in consequtive order.
}
	\item{node2first_edge}{
	An integer vector of size Nnodes listing the location of the first outgoing edge of each node in \code{edges}. That is, \code{edges[node2first_edge[n]]} points to the first outgoing edge of node n in \code{tree$edge}.
	}
	\item{node2last_edge}{
	An integer vector of size Nnodes listing the location of the last outgoing edge of each node in \code{edges}. That is, \code{edges[node2last_edge[n]]} points to the last outgoing edge of node n in \code{tree$edge}. The total number of outgoing edges of a node is thus given by \code{1+node2last_edge[n]-node2first_edge[n]}.
	}
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{reorder_tree_edges}}
}

\examples{
\dontrun{
# generate a random tree
tree = generate_random_tree(list(birth_rate_factor=1), max_tips=100)$tree

# get tree traversal
traversal = get_tree_traversal_root_to_tips(tree, include_tips=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{traversal}
