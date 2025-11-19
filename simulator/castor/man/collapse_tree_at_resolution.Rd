\name{collapse_tree_at_resolution}
\alias{collapse_tree_at_resolution}
\title{
Collapse nodes of a tree at a phylogenetic resolution.
}
\description{
Given a rooted tree and a phylogenetic resolution threshold, collapse all nodes whose distance to all descending tips does not exceed the threshold (or whose sum of descending edge lengths does not exceed the threshold), into new tips. This function can be used to obtain a "coarser" version of the tree, or to cluster closely related tips into a single tip.
}
\usage{
collapse_tree_at_resolution(tree, 
                            resolution             = 0, 
                            by_edge_count          = FALSE,
                            shorten                = TRUE,
                            rename_collapsed_nodes = FALSE,
                            criterion              = 'max_tip_depth')
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{resolution}{
Numeric, specifying the phylogenetic resolution at which to collapse the tree. This is the maximum distance a descending tip can have from a node, such that the node is collapsed into a new tip. If set to 0 (default), then only nodes whose descending tips are identical to the node will be collapsed.
}
\item{by_edge_count}{
Logical. Instead of considering edge lengths, consider edge counts as phylogenetic distance between nodes and tips. This is the same as if all edges had length equal to 1.
}
\item{shorten}{
Logical, indicating whether collapsed nodes should be turned into tips at the same location (thus potentially shortening the tree). If \code{FALSE}, then the incoming edge of each collapsed node is extended by some length L, where L is the distance of the node to its farthest descending tip (thus maintaining the height of the tree).
}
\item{rename_collapsed_nodes}{
Logical, indicating whether collapsed nodes should be renamed using a representative tip name (the farthest descending tip). See details below.
}
\item{criterion}{
Character, specifying the criterion to use for collapsing (i.e. how to interpret \code{resolution}). '\code{max_tip_depth}': Collapse nodes based on their maximum distance to any descending tip. '\code{sum_tip_paths}': Collapse nodes based on the sum of descending edges (each edge counted once). '\code{max_tip_pair_dist}': Collapse nodes based on the maximum distance between any pair of descending tips.
}
}


\details{
The tree is traversed from root to tips and nodes are collapsed into tips as soon as the criterion equals or falls below the resolution threshold.

The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child).
Tip labels and uncollapsed node labels of the collapsed tree are inheritted from the original tree. If \code{rename_collapsed_nodes==FALSE}, then labels of collapsed nodes will be the node labels from the original tree (in this case the original tree should include node labels). If \code{rename_collapsed_nodes==TRUE}, each collapsed node is given the label of its farthest descending tip. If \code{shorten==TRUE}, then edge lengths are the same as in the original tree. If \code{shorten==FALSE}, then edges leading into collapsed nodes may be longer than before.
}


\value{
A list with the following elements:
\item{tree}{
A new rooted tree of class "phylo", containing the collapsed tree.
}
\item{root_shift}{
Numeric, indicating the phylogenetic distance between the old and the new root. Will always be non-negative.
}
\item{collapsed_nodes}{
Integer vector, listing indices of collapsed nodes in the original tree (subset of 1,..,Nnodes).
}
\item{farthest_tips}{
Integer vector of the same length as \code{collapsed_nodes}, listing indices of the farthest tips for each collapsed node. Hence, \code{farthest_tips[n]} will be the index of a tip in the original tree that descended from node \code{collapsed_nodes[n]} and had the greatest distance from that node among all descending tips.
}
\item{new2old_clade}{
Integer vector of length equal to the number of tips+nodes in the collapsed tree, with values in 1,..,Ntips+Nnodes, mapping tip/node indices of the collapsed tree to tip/node indices in the original tree.
}
\item{new2old_edge}{
Integer vector of length equal to the number of edges in the collapsed tree, with values in 1,..,Nedges, mapping edge indices of the collapsed tree to edge indices in the original tree.
}
}

\author{Stilianos Louca}

%\references{
%}

%\seealso{
%}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=1000)$tree

# print number of nodes
cat(sprintf("Simulated tree has \%d nodes\n",tree$Nnode))

# collapse any nodes with tip-distances < 20
collapsed = collapse_tree_at_resolution(tree, resolution=20)$tree

# print number of nodes
cat(sprintf("Collapsed tree has \%d nodes\n",collapsed$Nnode))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{pruning}
\keyword{subtree}
