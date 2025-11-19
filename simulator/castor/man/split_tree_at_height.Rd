\name{split_tree_at_height}
\alias{split_tree_at_height}
\title{
Split a tree into subtrees at a specific height.
}
\description{
Given a rooted phylogenetic tree and a specific distance from the root (``height''), split the tree into subtrees at the specific height. This corresponds to drawing the tree in rectangular layout and trimming everything below the specified phylogenetic distance from the root: What is obtained is a set of separated subtrees. The tips of the original tree are spread across those subtrees.
}
\usage{
split_tree_at_height(tree, height = 0, by_edge_count = FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{height}{
Numeric, specifying the phylogenetic distance from the root at which to split the tree. If <=0, the original tree is returned as the sole subtree.
}
\item{by_edge_count}{
Logical. Instead of considering edge lengths, consider edge counts as phylogenetic distance. This is the same as if all edges had length equal to 1.
}
}


\details{
This function can be used to generate multiple smaller trees from one large tree, with each subtree having a time span equal to or lower than a certain threshold. The input tree may include multifurcations (i.e. nodes with more than 2 children) as well as monofurcations (i.e. nodes with only one child).

Note that while edges are cut exactly at the specified distance from the root, the cut point does not become the root node of the obtained subtree; rather, the first node encountered after the cut will become the subtree's root. The length of the remaining edge segment leading into this node will be used as \code{root.edge} in the returned subtree.
}


\value{
A list with the following elements:
\item{Nsubtrees}{
Integer, the number of subtrees obtained.
}
\item{subtrees}{
A list of length Nsubtrees, each element of which is a named list containing the following elements:
\itemize{
\item \code{tree}: A rooted tree of class "phylo", representing a subtree obtained from the original tree.
\item \code{new2old_clade}: An integer vector of length NStips+NSnodes (where NStips is the number of tips and NSnodes the number of nodes of the subtree), mapping subtree tip and node indices (i.e., 1,..,NStips+NSnodes) to tip and node indices in the original tree.
\item \code{new2old_edge}: Integer vector of length NSedges (=number of edges in the subtree), mapping subtree edge indices (i.e., 1,..,NSedges) to edge indices in the original tree.
}
}
\item{clade2subtree}{
Integer vector of length Ntips+Nnodes and containing values from 1 to Nsubtrees, mapping tip and node indices of the original tree to their assigned subtree.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{trim_tree_at_height}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),
                            max_tips=100)$tree

# split tree halfway towards the root
root_age = get_tree_span(tree)$max_distance
splitting = split_tree_at_height(tree, height=0.5*root_age)

# print number of subtrees obtained
cat(sprintf("Obtained \%d subtrees\n",splitting$Nsubtrees))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{pruning}
