\name{shift_clade_times}
\alias{shift_clade_times}
\title{
Shift the time of specific nodes & tips.
}
\description{
Given a rooted tree, shift the times (distance from root) of specific tips & nodes.
}
\usage{
shift_clade_times(  tree,
                    clades_to_shift,
                    time_shifts,
                    shift_descendants       = FALSE,
                    negative_edge_lengths   = "error")
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{clades_to_shift}{
Integer or character vector, listing the tips and/or nodes whose time is to be shifted. If an integer vector, values must correspond to indices and must be in the range 1,..,Ntips+Nnodes. If a character vector, values must correspond to tip and/or node labels in the tree; if node labels are listed, the tree must contain node labels (attribute \code{node.label}).
}
\item{time_shifts}{
Numeric vector of the same length as \code{clades_to_shift}, specifying the time shifts to apply to every tip/node listed in \code{clades_to_shift}. Values can be negative (shift towards the root) or positive (shift away from the root).
}
\item{shift_descendants}{
Logical, specifying whether to shift the entire descending subclade when shifting a node. If \code{FALSE}, the descending tips & nodes retain their original time (unless negative edges are created, see option \code{negative_edge_lengths}).
}
\item{negative_edge_lengths}{
Character, specifying whether and how to fix negative edge lengths resulting from excessive shifting. Must be either "\code{error}", "\code{allow}" (allow and don't fix negative edge lengths), "\code{move_all_descendants}" (move all descendants forward as needed, to make the edge non-negative), "\code{move_all_ancestors}" (move all ancestors backward as needed, to make the edge non-negative), "\code{move_child}" (only move children to younger ages as needed, traversing the tree root->tips) or "\code{move_parent}" (only move parents to older ages as needed, traversing the tree tips->root). Note that "\code{move_child}" could result in tips moving, if an ancestral node is shifted too much towards younger ages. Similarly, "\code{move_parent}" could result in the root moving towards an older age if some descendant was shifted too far towards the root.
}
}


\details{
The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). The input tree does not need to be ultrametric, but edge lengths are interpreted as time. If edge lengths are missing from the tree, it is assumed that each edge has length 1.

All tips, nodes and edges are kept and indexed as in the input tree; the only thing that changes are the edgen lengths.

Note that excessive shifting can result in negative edge lengths, which can be corrected in a variety of alternative ways (see option \code{negative_edge_lengths}). However, to avoid changing the overall span of the tree (root age and tip times) in an effort to fix negative edge lengths, you should generally not shift a clade beyond the boundaries of the tree (i.e., resulting in a negative time or a time beyond its descending tips).
}


\value{
A list with the following elements:
\item{success}{
Logical, specifying whether the operation was successful. If \code{FALSE}, an additional variable \code{error} is returned, briefly specifying the error, but all other return variables may be undefined.
}
\item{tree}{
A new rooted tree of class "phylo", representing the tree with shifted clade times.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_all_distances_to_root}},
\code{\link{trim_tree_at_height}},
\code{\link{get_tree_span}}
}

\examples{
# generate a random tree, include node names
tree = generate_random_tree(list(birth_rate_intercept=1),
                            max_tips=20,
                            node_basename="node.")$tree

# shift a few nodes backward in time,
# changing as few of the remaining node timings as possible
clades_to_shift = c("node.2",   "node.5",   "node.6")
time_shifts     = c(-0.5,       -0.2,       -0.3)
new_tree        = shift_clade_times(tree, 
                                    clades_to_shift, 
                                    time_shifts,
                                    shift_descendants=FALSE,
                                    negative_edge_lengths="move_parent")$tree
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{phylogenetic distance}
\keyword{dating}
