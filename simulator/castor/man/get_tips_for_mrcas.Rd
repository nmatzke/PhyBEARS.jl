\name{get_tips_for_mrcas}
\alias{get_tips_for_mrcas}
\title{
Find tips with specific most recent common ancestors.
}
\description{
Given a rooted phylogenetic tree and a list of nodes ("MRCA nodes"), for each MRCA node find a set of descending tips ("MRCA-defining tips") such that their most recent common ancestor (MRCA) is that node. This may be useful for cases where nodes need to be described as MRCAs of tip pairs for input to certain phylogenetics algorithms (e.g., for tree dating).
}
\usage{
get_tips_for_mrcas(tree, mrca_nodes, check_input=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{mrca_nodes}{
Either an integer vector or a character vector, listing the nodes for each of which an MRCA-defining set of tips is to be found. If an integer vector, it should list node indices (i.e. from 1 to Nnodes). If a character vector, it should list node names; in that case \code{tree$node.label} must exist.
}
\item{check_input}{
Logical, whether to perform basic validations of the input data. If you know for certain that your input is valid, you can set this to \code{FALSE} to reduce computation time.
}
}
\details{
At most 2 MRCA-defining tips are assigned to each MRCA node.
This function assumes that each of the \code{mrca_nodes} has at least two children or has a child that is a tip (otherwise the problem is not well-defined). The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child).

The asymptotic time complexity of this function is O(Ntips+Nnodes) + O(Nmrcas), where Ntips is the number of tips, Nnodes is the number of nodes in the tree and Nmrcas is equal to \code{length(mrca_nodes)}.
}
\value{
A list of the same size as \code{mrca_nodes}, whose n-th element is an integer vector of tip indices (i.e. with values in 1,..,Ntips) whose MRCA is the n-th node listed in \code{mrca_nodes}.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_pairwise_mrcas}}, \code{\link{get_mrca_of_set}}
}

\examples{
# generate a random tree
Ntips = 1000
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# pick random nodes
focal_nodes = sample.int(n=tree$Nnode, size=3, replace=FALSE)

# get tips for mrcas
tips_per_focal_node = get_tips_for_mrcas(tree, focal_nodes);

# check correctness (i.e. calculate actual MRCAs of tips)
for(n in 1:length(focal_nodes)){
  mrca = get_mrca_of_set(tree, tips_per_focal_node[[n]])
  cat(sprintf("Focal node = \%d, should match mrca of tips = \%d\n",focal_nodes[n],mrca-Ntips))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MRCA}
