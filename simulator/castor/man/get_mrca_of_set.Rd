\name{get_mrca_of_set}
\alias{get_mrca_of_set}
\title{
Most recent common ancestor of a set of tips/nodes.
}
\description{
Given a rooted phylogenetic tree and a set of tips and/or nodes ("descendants"), calculate the most recent common ancestor (MRCA) of those descendants.
}
\usage{
get_mrca_of_set(tree, descendants)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{descendants}{
An integer vector or character vector, specifying the tips/nodes for which to find the MRCA. If an integer vector, it must list indices of tips (from 1 to Ntips) and/or nodes (from Ntips+1 to Ntips+Nnodes), where Ntips and Nnodes is the number of tips and nodes in the tree, respectively. If a character vector, it must list tip and/or node names. In this case \code{tree} must include \code{tip.label}, as well as \code{node.label} if nodes are included in \code{descendants}.
}
}
\details{
The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child).  Duplicate entries in \code{descendants} are ignored.
}
\value{
An integer in 1,..,(Ntips+Nnodes), representing the MRCA using the same index as in \code{tree$edge}. If the MRCA is a tip, then this index will be in 1,..,Ntips. If the MRCA is a node, then this index will be in (Ntips+1),..,(Ntips+Nnodes).
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_pairwise_mrcas}}, \code{\link{get_tips_for_mrcas}}
}

\examples{
# generate a random tree
Ntips = 1000
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# pick 3 random tips or nodes
descendants = sample.int(n=(Ntips+tree$Nnode), size=3, replace=FALSE)

# calculate MRCA of picked descendants
mrca = get_mrca_of_set(tree, descendants)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MRCA}
