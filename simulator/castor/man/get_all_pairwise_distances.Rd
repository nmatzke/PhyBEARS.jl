\name{get_all_pairwise_distances}
\alias{get_all_pairwise_distances}
\title{
Get distances between all pairs of tips and/or nodes.
}
\description{
Calculate phylogenetic ("patristic") distances between all pairs of tips or nodes in the tree, or among a subset of tips/nodes requested.
}
\usage{
get_all_pairwise_distances( tree,
                            only_clades     = NULL, 
                            as_edge_counts  = FALSE,
                            check_input     = TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{only_clades}{
Optional integer vector or character vector, listing tips and/or nodes to which to restrict pairwise distance calculations. If an integer vector, it must list indices of tips (from 1 to Ntips) and/or nodes (from Ntips+1 to Ntips+Nnodes). If a character vector, it must list tip and/or node names. 

For example, if \code{only_clades=c('apple','lemon','pear')}, then only the distance between `apple' and `lemon', between `apple' and 'pear', and between `lemon' and `pear' are calculated. If \code{only_clades==NULL}, then this is equivalent to \code{only_clades=c(1:(Ntips+Nnodes))}.
}
\item{check_input}{
Logical, whether to perform basic validations of the input data. If you know for certain that your input is valid, you can set this to \code{FALSE} to reduce computation time.
}
\item{as_edge_counts}{
Logical, specifying whether distances should be calculated in terms of edge counts, rather than cumulative edge lengths. This is the same as if all edges had length 1.
}
}
\details{
The "patristic distance" between two tips and/or nodes is the shortest cumulative branch length that must be traversed along the tree in order to reach one tip/node from the other.This function returns a square distance matrix, containing the patristic distance between all possible pairs of tips/nodes in the tree (or among the ones provided in \code{only_clades}).

If \code{tree$edge.length} is missing, then each edge is assumed to be of length 1; this is the same as setting \code{as_edge_counts=TRUE}. The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). The input tree must be rooted at some node for technical reasons (see function \code{\link{root_at_node}}), but the choice of the root node does not influence the result. If \code{only_clades} is a character vector, then \code{tree$tip.label} must exist. If node names are included in \code{only_clades}, then \code{tree$node.label} must also exist.

The asymptotic average time complexity of this function for a balanced binary tree is O(NC*NC*Nanc + Ntips), where NC is the number of tips/nodes considered (e.g., the length of \code{only_clades}) and Nanc is the average number of ancestors per tip.
}
\value{
A 2D numeric matrix of size NC x NC, where NC is the number of tips/nodes considered, and with the entry in row r and column c listing the distance between the r-th and the c-th clade considered (e.g., between clades \code{only_clades[r]} and \code{only_clades[c]}). Note that if \code{only_clades} was specified, then the rows and columns in the returned distance matrix correspond to the entries in \code{only_clades} (i.e., in the same order). If \code{only_clades} was \code{NULL}, then the rows and columns in the returned distance matrix correspond to tips (1,..,Ntips) and nodes (Ntips+1,..,Ntips+Nnodes)
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_all_distances_to_root}},
\code{\link{get_pairwise_distances}}
}

\examples{
# generate a random tree
Ntips = 100
tree  = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# calculate distances between all internal nodes
only_clades = c((Ntips+1):(Ntips+tree$Nnode))
distances = get_all_pairwise_distances(tree, only_clades)

# reroot at some other node
tree = root_at_node(tree, new_root_node=20, update_indices=FALSE)
new_distances = get_all_pairwise_distances(tree, only_clades)

# verify that distances remained unchanged
plot(distances,new_distances,type='p')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{phylogenetic distance}
