\name{find_nearest_tips}
\alias{find_nearest_tips}
\title{
Find nearest tip to each tip & node of a tree.
}
\description{
Given a rooted phylogenetic tree and a subset of potential target tips, for each tip and node in the tree find the nearest target tip. The set of target tips can also be taken as the whole set of tips in the tree.
}
\usage{
find_nearest_tips(tree, 
                  only_descending_tips = FALSE, 
                  target_tips          = NULL, 
                  as_edge_counts       = FALSE, 
                  check_input          = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{only_descending_tips}{
A logical indicating whether the nearest tip to a node or tip should be chosen from its descending tips only. If FALSE, then the whole set of possible target tips is considered.
}
\item{target_tips}{
Optional integer vector or character vector listing the subset of target tips to restrict the search to. If an integer vector, this should list tip indices (values in 1,..,Ntips). If a character vector, it should list tip names (in this case \code{tree$tip.label} must exist). If \code{target_tips} is \code{NULL}, then all tips of the tree are considered as target tips.
}
\item{as_edge_counts}{
Logical, specifying whether to count phylogenetic distance in terms of edge counts instead of cumulative edge lengths. This is the same as setting all edge lengths to 1.
}
\item{check_input}{
Logical, whether to perform basic validations of the input data. If you know for certain that your input is valid, you can set this to \code{FALSE} to reduce computation time.
}
}


\details{
Langille et al. (2013) introduced the Nearest Sequenced Taxon Index (NSTI) as a measure for how well a set of microbial operational taxonomic units (OTUs) is represented by a set of sequenced genomes of related organisms. Specifically, the NSTI of a microbial community is the average phylogenetic distance of any OTU in the community, to the closest relative with an available sequenced genome ("target tips"). In analogy to the NSTI, the function \code{find_nearest_tips} provides a means to find the nearest tip (from a subset of target tips) to each tip and node in a phylogenetic tree, together with the corresponding phylogenetic ("patristic") distance. 

If \code{only_descending_tips} is \code{TRUE}, then only descending target tips are considered when searching for the nearest target tip of a node/tip. In that case, if a node/tip has no descending target tip, its nearest target tip is set to NA. If \code{tree$edge.length} is missing or NULL, then each edge is assumed to have length 1. The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child).

The asymptotic time complexity of this function is O(Nedges), where Nedges is the number of edges in the tree.
}


\value{
A list with the following elements:
\item{nearest_tip_per_tip}{
An integer vector of size Ntips, listing the nearest target tip for each tip in the tree. Hence, \code{nearest_tip_per_tip[i]} is the index of the nearest tip (from the set of target tips), with respect to tip i (where i=1,..,Ntips). Some values may appear multiple times in this vector, if multiple tips share the same nearest target tip.
}
\item{nearest_tip_per_node}{
An integer vector of size Nnodes, listing the index of the nearest target tip for each node in the tree. Hence, \code{nearest_tip_per_node[i]} is the index of the nearest tip (from the set of target tips), with respect to node i (where i=1,..,Nnodes). Some values may appear multiple times in this vector, if multiple nodes share the same nearest target tip.
}
\item{nearest_distance_per_tip}{
Integer vector of size Ntips. Phylogenetic ("patristic") distance of each tip in the tree to its nearest target tip. If \code{only_descending_tips} was set to \code{TRUE}, then \code{nearest_distance_per_tip[i]} will be set to infinity for any tip i that is not a target tip.
}
\item{nearest_distance_per_node}{
Integer vector of size Nnodes. Phylogenetic ("patristic") distance of each node in the tree to its nearest target tip. If \code{only_descending_tips} was set to \code{TRUE}, then \code{nearest_distance_per_node[i]} will be set to infinity for any node i that has no descending target tips.
}
}


\author{Stilianos Louca}

\references{
M. G. I. Langille, J. Zaneveld, J. G. Caporaso et al (2013). Predictive functional profiling of microbial communities using 16S rRNA marker gene sequences. Nature Biotechnology. 31:814-821.
}

\seealso{
\code{\link{find_farthest_tips}}
}

\examples{
# generate a random tree
Ntips = 1000
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# pick a random set of "target" tips
target_tips = sample.int(n=Ntips, size=as.integer(Ntips/10), replace=FALSE)

# find nearest target tip to each tip & node in the tree
results = find_nearest_tips(tree, target_tips=target_tips)

# plot histogram of distances to target tips (across all tips of the tree)
distances = results$nearest_distance_per_tip
hist(distances, breaks=10, xlab="nearest distance", ylab="number of tips", prob=FALSE);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }
