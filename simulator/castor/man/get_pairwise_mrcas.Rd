\name{get_pairwise_mrcas}
\alias{get_pairwise_mrcas}
\title{
Get most recent common ancestors of tip/node pairs.
}
\description{
Given a rooted phylogenetic tree and one or more pairs of tips and/or nodes, for each pair of tips/nodes find the most recent common ancestor (MRCA). If one clade is descendant of the other clade, the latter will be returned as MRCA.
}
\usage{
get_pairwise_mrcas(tree, A, B, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{A}{
An integer vector or character vector of size Npairs, specifying the first of the two members of each pair of tips/nodes for which to find the MRCA. If an integer vector, it must list indices of tips (from 1 to Ntips) and/or nodes (from Ntips+1 to Ntips+Nnodes). If a character vector, it must list tip and/or node names.
}
\item{B}{
An integer vector or character vector of size Npairs, specifying the second of the two members of each pair of tips/nodes for which to find the MRCA. If an integer vector, it must list indices of tips (from 1 to Ntips) and/or nodes (from Ntips+1 to Ntips+Nnodes). If a character vector, it must list tip and/or node names.
}
\item{check_input}{
Logical, whether to perform basic validations of the input data. If you know for certain that your input is valid, you can set this to \code{FALSE} to reduce computation time.
}
}
\details{
The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). If \code{tree$edge.length} is missing, then each edge is assumed to be of length 1. Note that in some cases the MRCA of two tips may be a tip, namely when both tips are the same.

If \code{A} and/or \code{B} is a character vector, then \code{tree$tip.label} must exist. If node names are included in \code{A} and/or \code{B}, then \code{tree$node.label} must also exist.

The asymptotic average time complexity of this function is O(Nedges), where Nedges is the number of edges in the tree.
}
\value{
An integer vector of size Npairs with values in 1,..,Ntips (tips) and/or in (Ntips+1),..,(Ntips+Nnodes) (nodes), with the i-th element being the index of the MRCA of tips/nodes \code{A[i]} and \code{B[i]}. 
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_mrca_of_set}}, \code{\link{get_tips_for_mrcas}}
}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# pick 3 random pairs of tips or nodes
Npairs = 3
A = sample.int(n=(Ntips+tree$Nnode), size=Npairs, replace=FALSE)
B = sample.int(n=(Ntips+tree$Nnode), size=Npairs, replace=FALSE)

# calculate MRCAs
MRCAs = get_pairwise_mrcas(tree, A, B)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MRCA}
