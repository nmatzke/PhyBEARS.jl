\name{get_subtree_with_tips}
\alias{get_subtree_with_tips}
\title{
Extract a subtree spanning a specific subset of tips.
}
\description{
Given a rooted tree and a subset of tips, extract the subtree containing only those tips. The root of the tree is kept.
}
\usage{
get_subtree_with_tips(tree, 
                      only_tips               = NULL, 
                      omit_tips               = NULL, 
                      collapse_monofurcations = TRUE, 
                      force_keep_root         = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{only_tips}{
Either a character vector listing tip names to keep, or an integer vector listing tip indices to keep (between 1 and Ntips). Can also be \code{NULL}. Tips listed in \code{only_tips} not found in the tree will be silently ignored.
}
\item{omit_tips}{
Either a character vector listing tip names to omit, or an integer vector listing tip indices to omit (between 1 and Ntips). Can also be \code{NULL}. Tips listed in \code{omit_tips} not found in the tree will be silently ignored.
}
\item{collapse_monofurcations}{
A logical specifying whether nodes with a single outgoing edge remaining should be collapsed (removed). Incoming and outgoing edge of such nodes will be concatenated into a single edge, connecting the parent (or earlier) and child (or later) of the node. In that case, the returned tree will have edge lengths that reflect the concatenated edges.
}
\item{force_keep_root}{
Logical, specifying whether to keep the root even if \code{collapse_monofurcations==TRUE} and the root of the subtree is left with a single child. If \code{FALSE}, and \code{collapse_monofurcations==TRUE}, the root may be removed and one of its descendants may become root.
}
}


\details{
If both \code{only_tips} and \code{omit_tips} are \code{NULL}, then all tips are kept and the tree remains unchanged. If both \code{only_tips} and \code{omit_tips} are non-\code{NULL}, then only tips listed in \code{only_tips} and not listed in \code{omit_tips} will be kept. If \code{only_tips} and/or \code{omit_tips} is a character vector listing tip names, then \code{tree$tip.label} must exist.

If the input tree does not include \code{edge.length}, each edge in the input tree is assumed to have length 1. The root of the tree (which is always kept) is assumed to be the unique node with no incoming edge. The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child).

The asymptotic time complexity of this function is O(Nnodes+Ntips), where Ntips is the number of tips and Nnodes the number of nodes in the input tree.

When \code{only_tips==NULL}, \code{omit_tips!=NULL}, \code{collapse_monofurcations==TRUE} and \code{force_keep_root==FALSE}, this function is analogous to the function \code{drop.tip} in the \code{ape} package with option \code{trim_internal=TRUE} (v. 0.5-64).
}


\value{
A list with the following elements:
\item{subtree}{
A new tree of class "phylo", containing only the tips specified by \code{tips_to_keep} and the nodes & edges connecting those tips to the root. The returned tree will include \code{edge.lengh} as a member variable, listing the lengths of the remaining (possibly concatenated) edges. 
}
\item{root_shift}{
Numeric, indicating the phylogenetic distance between the old and the new root. Will always be non-negative.
}
\item{new2old_tip}{
Integer vector of length Ntips_kept (=number of tips in the extracted subtree) with values in 1,..,Ntips, mapping tip indices of the subtree to tip indices in the original tree. In particular, \code{tree$tip.label[new2old_tip]} will be equal to \code{subtree$tip.label}.
}
\item{new2old_node}{
Integer vector of length Nnodes_kept (=number of nodes in the extracted subtree) with values in 1,..,Nnodes, mapping node indices of the subtree to node indices in the original tree. 

For example, \code{new2old_node[1]} is the index that the first node of the subtree had within the original tree. In particular, \code{tree$node.label[new2old_node]} will be equal to \code{subtree$node.label} (if node labels are available).
}
\item{old2new_tip}{
Integer vector of length Ntips, with values in 1,..,Ntips_kept, mapping tip indices of the original tree to tip indices in the subtree (a value of 0 is used whenever a tip is absent in the subtree). This is essentially the inverse of the mapping \code{new2old_tip}.
}
\item{old2new_node}{
Integer vector of length Nnodes, with values in 1,..,Nnodes_kept, mapping node indices of the original tree to node indices in the subtree (a value of 0 is used whenever a node is absent in the subtree). This is essentially the inverse of the mapping \code{new2old_node}.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_subtree_at_node}}
}

\examples{
# generate a random tree
Ntips = 1000
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# choose a random subset of tips
tip_subset = sample.int(Ntips, size=as.integer(Ntips/10), replace=FALSE)

# extract subtree spanning the chosen tip subset
subtree = get_subtree_with_tips(tree, only_tips=tip_subset)$subtree

# print summary of subtree
cat(sprintf("Subtree has \%d tips and \%d nodes\n",length(subtree$tip.label),subtree$Nnode))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{pruning}
\keyword{subtree}
