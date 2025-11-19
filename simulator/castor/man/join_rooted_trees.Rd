\name{join_rooted_trees}
\alias{join_rooted_trees}
\title{
Join two rooted trees.
}
\description{
Given two rooted phylogenetic trees, place one tree (tree2) onto an edge of the other tree (tree1), so that tree2 becomes a monophyletic group of the final joined tree. As a special case, this function can join two trees at their roots, i.e. so that both are disjoint monophyletic clades of the final tree, splitting at the new root.
}
\usage{
join_rooted_trees(  tree1, 
                    tree2,
                    target_edge1,
                    target_edge_length1,
                    root_edge_length2)
}
\arguments{
\item{tree1}{
A rooted tree of class "phylo".
}
\item{tree2}{
A rooted tree of class "phylo". This tree will become a monophyletic subclade of the final joined tree.
}
\item{target_edge1}{
Integer, edge index in \code{tree1} onto which \code{tree2} is to be joined. If <=0, then this refers to the hypothetical edge leading into the root of \code{tree1}, in which case both trees will become disjoint monophyletic subclades of the final joined tree.
}
\item{target_edge_length1}{
Numeric, length of the edge segment in \code{tree1} from the joining-point to the next child node, i.e. how far from the child of \code{target_edge1} should the joining occur. If \code{target_edge1<=0}, then \code{target_edge_length1} is the distance of the root of \code{tree1} from the final joined tree's root.
}
\item{root_edge_length2}{
Numeric, length of the edge leading into the root of \code{tree2}, i.e. the distance from the joining point to the root of \code{tree2}.
}
}


\details{
The input trees may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). If any of the input trees does not have edge lengths (i.e., \code{edge.length} is \code{NULL}), then its edge lengths are assumed to all be 1.

The tips of the two input trees will become the tips of the final joined tree. The nodes of the two input trees will become nodes of the final joined tree, however one additional node will be added at the joining point. Tip labels and node labels (if available) of the joined tree are inheritted from the two input trees.
}


\value{
A list with the following elements:
\item{tree}{
A new rooted tree of class "phylo", representing the joined tree.
}
\item{clade1_to_clade}{
Integer vector of length Ntips1+Nnodes1, mapping tip/node indices of the input \code{tree1} to tip/node indices in the final joined tree.
}
\item{clade2_to_clade}{
Integer vector of length Ntips2+Nnodes2, mapping tip/node indices of the input \code{tree2} to tip/node indices in the final joined tree.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{split_tree_at_height}}
}

\examples{
# generate two random trees, include tip & node names
tree1 = generate_random_tree(list(birth_rate_intercept=1),
                             max_tips=10,
                             tip_basename="tip1.",
                             node_basename="node1.")$tree
tree2 = generate_random_tree(list(birth_rate_intercept=1),
                             max_tips=5,
                             tip_basename="tip2.",
                             node_basename="node2.")$tree

# join trees at their roots
# each subtree's root should have distance 1 from the new root
joined_tree = join_rooted_trees(tree1, 
                                tree2,
                                target_edge1=0,
                                target_edge_length1=1,
                                root_edge_length2=1)$tree
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{joining trees}
