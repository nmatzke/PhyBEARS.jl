\name{find_root_of_monophyletic_tips}
\alias{find_root_of_monophyletic_tips}
\title{
Find the node or tip that, as root, would make a set of target tips monophyletic.
}
\description{
Given a tree (rooted or unrooted) and a specific set of target tips, this function finds the tip or node that, if turned into root, would make a set of target tips a monophyletic group that either descends from a single child of the new root (if \code{as_MRCA==FALSE}) or whose MRCA is the new root (if \code{as_MRCA==TRUE}).
}
\usage{
find_root_of_monophyletic_tips(tree, monophyletic_tips, as_MRCA=TRUE, is_rooted=FALSE)
}
\arguments{
\item{tree}{
A tree object of class "phylo". Can be unrooted or rooted.
}
\item{monophyletic_tips}{
Character or integer vector, specifying the names or indices, respectively, of the target tips that should be turned monophyletic. If an integer vector, its elements must be between 1 and Ntips. If a character vector, its elements must be elements in \code{tree$tip.label}.
}
\item{as_MRCA}{
Logical, specifying whether the new root should become the MRCA of the target tips. If \code{FALSE}, the new root is chosen such that the MRCA of the target tips is the child of the new root.
}
\item{is_rooted}{
Logical, specifying whether the input tree can be assumed to be rooted. If you are sure that the input tree is rooted, set this to \code{TRUE} for computational efficiency, otherwise to be on the safe side set this to \code{FALSE}.
}
}

\details{
The input tree may include an arbitrary number of incoming and outgoing edges per node (but only one edge per tip), and the direction of these edges can be arbitrary. Of course, the undirected graph defined by all edges must still be a valid tree (i.e. a connected acyclic graph). Note that this function does not change the tree, it just determines which tip or node should be made root for the target tips to be a monophyletic group.


The asymptotic time complexity of this function is O(Nedges).
}

\value{
A single integer between 1 and (Ntips+Nnodes), specifying the index of the tip or node that, if made root, would make the target tips monophyletic. If this was not possible, \code{NA} is returned.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{find_root}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# pick a random node and find all descending tips
MRCA = sample.int(tree$Nnode,size=1)
monophyletic_tips = get_subtree_at_node(tree, MRCA)$new2old_tip

# change root of tree (change edge directions)
tree = root_at_node(tree, new_root_node=10, update_indices=FALSE)

# determine root that would make target tips monophyletic
new_root = find_root_of_monophyletic_tips(tree, monophyletic_tips, as_MRCA=TRUE, is_rooted=FALSE)

# compare expectation with result
cat(sprintf("MRCA = \%d, new root node=\%d\n",MRCA,new_root-Ntips))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Rerooting }
