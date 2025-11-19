\name{trim_tree_at_height}
\alias{trim_tree_at_height}
\title{
Trim a rooted tree down to a specific height.
}
\description{
Given a rooted phylogenetic tree and a maximum allowed distance from the root (``height''), remove tips and nodes and shorten the remaining terminal edges so that the tree's height does not exceed the specified threshold. This corresponds to drawing the tree in rectangular layout and trimming everything beyond a specific phylogenetic distance from the root. Tips or nodes at the end of trimmed edges are kept, and the affected edges are shortened. 
}
\usage{
trim_tree_at_height(tree, height = Inf, by_edge_count = FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{height}{
Numeric, specifying the phylogenetic distance from the root at which to trim.
}
\item{by_edge_count}{
Logical. Instead of considering edge lengths, consider edge counts as phylogenetic distance. This is the same as if all edges had length equal to 1.
}
}


\details{
The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child).

Tip labels and uncollapsed node labels of the collapsed tree are inheritted from the original tree. Labels of tips that used to be nodes (i.e. of which all descendants have been removed) will be the node labels from the original tree. If the input tree has no node names, it is advised to first add node names to avoid \code{NA} in the resulting tip names.
}


\value{
A list with the following elements:
\item{tree}{
A new rooted tree of class "phylo", representing the trimmed tree.
}
\item{Nedges_trimmed}{
Integer. Number of edges trimmed (shortened).
}
\item{Nedges_removed}{
Integer. Number of edges removed.
}
\item{new2old_clade}{
Integer vector of length equal to the number of tips+nodes in the trimmed tree, with values in 1,..,Ntips+Nnodes, mapping tip/node indices of the trimmed tree to tip/node indices in the original tree. In particular,

\code{c(tree$tip.label,tree$node.label)[new2old_clade]}

will be equal to:

\code{c(trimmed_tree$tip.label,trimmed_tree$node.label)}.
}
\item{new2old_edge}{
Integer vector of length equal to the number of edges in the trimmed tree, with values in 1,..,Nedges, mapping edge indices of the trimmed tree to edge indices in the original tree. In particular, \code{tree$edge.length[new2old_edge]} will be equal to \code{trimmed_tree$edge.length} (if edge lengths are available).
}
\item{new_edges_trimmed}{
Integer vector, listing edge indices in the trimmed tree that we originally longer edges and have been trimmed. In other words, these are the edges that "crossed" the trimming height.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{split_tree_at_height}}
}

\examples{
# generate a random tree, include node names
tree = generate_random_tree(list(birth_rate_intercept=1),
                            max_time=1000,
                            node_basename="node.")$tree

# print number of tips
cat(sprintf("Simulated tree has \%d tips\n",length(tree$tip.label)))

# trim tree at height 500
trimmed = trim_tree_at_height(tree, height=500)$tree

# print number of tips in trimmed tree
cat(sprintf("Trimmed tree has \%d tips\n",length(trimmed$tip.label)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{pruning}
