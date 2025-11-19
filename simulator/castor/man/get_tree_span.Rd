\name{get_tree_span}
\alias{get_tree_span}
\title{
Get min and max distance of any tip to the root.
}
\description{
Given a rooted phylogenetic tree, calculate the minimum and maximum phylogenetic distance (cumulative branch length) of any tip from the root.
}
\usage{
get_tree_span(tree, as_edge_count=FALSE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{as_edge_count}{
Logical, specifying whether distances should be counted in number of edges, rather than cumulative edge length. This is the same as if all edges had length 1.
}
}
\details{
If \code{tree$edge.length} is missing, then every edge in the tree is assumed to be of length 1. The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child).
The asymptotic average time complexity of this function is O(Nedges), where Nedges is the number of edges in the tree.
}

\value{
A named list with the following elements:
\item{min_distance}{
Minimum phylogenetic distance that any of the tips has to the root.
}
\item{max_distance}{
Maximum phylogenetic distance that any of the tips has to the root.
}
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_pairwise_distances}}
}

\examples{
# generate a random tree
Ntips   = 1000
params  = list(birth_rate_intercept=1, death_rate_intercept=0.5)
tree    = generate_random_tree(params, max_tips=Ntips, coalescent=FALSE)$tree

# calculate min & max tip distances from root
tree_span = get_tree_span(tree)
cat(sprintf("Tip min dist = \%g, max dist = \%g\n",
            tree_span$min_distance,
            tree_span$max_distance))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{phylogenetic distance}
