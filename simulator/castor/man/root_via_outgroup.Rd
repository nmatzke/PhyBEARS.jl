\name{root_via_outgroup}
\alias{root_via_outgroup}
\title{
Root or re-root a tree based on an outgroup tip.
}
\description{
Given a tree (rooted or unrooted) and a specific tip (``outgroup''), this function changes the direction of edges (\code{tree$edge}) such that the outgroup's parent node becomes the root. The number of tips and the number of nodes remain unchanged.
}
\usage{
root_via_outgroup(tree, outgroup, update_indices=TRUE)
}
\arguments{
\item{tree}{
A tree object of class "phylo". Can be unrooted or rooted.
}
\item{outgroup}{
Character or integer specifying the name or index, respectively, of the outgroup tip. If an integer, it must be between 1 and Ntips. If a character, it must be a valid entry in \code{tree$tip.label}.
}
\item{update_indices}{
Logical, specifying whether to update the node indices such that the new root is the first node in the list (as is common convention). This will modify \code{tree$node.label} (if it exists) and also the node indices listed in \code{tree$edge}.
}
}

\details{
The input tree may include an arbitrary number of incoming and outgoing edges per node (but only one edge per tip), and the direction of these edges can be arbitrary. Of course, the undirected graph defined by all edges must still be a valid tree.
The asymptotic time complexity of this function is O(Nedges).

If \code{update_indices==FALSE}, then node indices remain unchanged.
If \code{update_indices==TRUE} (default), then node indices are modified such that the new root is the first node (i.e. with index Ntips+1 in \code{edge} and with index 1 in \code{node.label}). This is common convention, but it may be undesirable in some cases. Setting \code{update_indices=FALSE} also reduces the computation required for rerooting. Tip indices always remain unchanged.
}

\value{
A tree object of class "phylo", with the \code{edge} element modified such that the outgroup tip's parent node is root. The elements \code{tip.label}, \code{edge.length} and \code{root.edge} (if they exist) are the same as for the input tree. If \code{update_indices==FALSE}, then the element \code{node.label} will also remain the same.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{root_at_node}},
\code{\link{root_at_midpoint}},
\code{\link{root_in_edge}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# reroot the tree using the 1st tip as outgroup
outgroup = 1
tree = root_via_outgroup(tree, outgroup, update_indices=FALSE)

# find new root index
cat(sprintf("New root is \%d\n",find_root(tree)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Rerooting }
