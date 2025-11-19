\name{find_root}
\alias{find_root}
\title{
Find the root of a tree.
}
\description{
Find the root of a phylogenetic tree. The root is defined as the unique node with no parent.
}
\usage{
find_root(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A tree of class "phylo". If the tree is not rooted, the function will return \code{NA}.
}
}


\details{
By convention, the root of a "phylo" tree is typically the first node (i.e. with index Ntips+1), however this is not always guaranteed. This function finds the root of a tree by searching for the node with no parent. If no such node exists, \code{NA} is returned. If multiple such nodes exist, \code{NA} is returned. If any node has more than 1 parent, \code{NA} is returned. Hence, this function can be used to test if a tree is rooted purely based on the edge structure, assuming that the tree is connected (i.e. not a forest).

The asymptotic time complexity of this function is O(Nedges), where Nedges is the number of edges in the tree.
}

\value{
Index of the tree's root, as listed in \code{tree$edge}. An integer ranging from Ntips+1 to Ntips+Nnodes, where Ntips and Nnodes is the number of tips and nodes in the tree, respectively. By convention, the root will typically be Ntips+1 but this is not guaranteed.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{find_root_of_monophyletic_tips}},
\code{\link{root_at_node}},
\code{\link{root_at_midpoint}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# reroot the tree at the 20-th node
new_root_node = 20
tree = root_at_node(tree, new_root_node, update_indices=FALSE)

# find new root index and compare with expectation
cat(sprintf("New root is \%d, expected at \%d\n",find_root(tree),new_root_node+Ntips))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }
