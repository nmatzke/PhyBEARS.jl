\name{get_all_node_depths}
\alias{get_all_node_depths}
\title{
Get the phylogenetic depth of each node in a tree.
}
\description{
Given a rooted phylogenetic tree, calculate the phylogenetic depth of each node (mean distance to its descending tips).
}
\usage{
get_all_node_depths(tree, as_edge_count=FALSE)
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
A numeric vector of size Nnodes, with the i-th element being the mean distance of the i-th node to all of its tips.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_all_distances_to_root}}
}

\examples{
# generate a random tree
Ntips = 1000
tree = generate_random_tree(list(birth_rate_intercept=1, 
                            	 death_rate_intercept=0.5),
                            max_tips=Ntips)$tree

# calculate node phylogenetic depths
node_depths = get_all_node_depths(tree)

# plot histogram of node depths
hist(node_depths, xlab="phylogenetic depth", ylab="# nodes", prob=FALSE);
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{phylogenetic distance}
