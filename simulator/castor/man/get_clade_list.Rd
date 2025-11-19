\name{get_clade_list}
\alias{get_clade_list}
\title{
Get a representation of a tree as a table listing tips/nodes.
}
\description{
Given a tree in standard "phylo" format, calculate an alternative representation of the tree structure as a list of tips/nodes with basic information on parents, children and incoming edge lengths. This function is analogous to the function \code{read.tree.nodes} in the R package \code{phybase}.
}
\usage{
get_clade_list(tree, postorder=FALSE, missing_value=NA)
}
\arguments{
\item{tree}{
A tree of class "phylo". If \code{postorder==TRUE}, then the tree must be rooted.
}
\item{postorder}{
Logical, specifying whether nodes should be ordered and indexed in postorder traversal, i.e. with the root node listed last. Note that regardless of the value of \code{postorder}, tips will always be listed first and indexed in the order in which they are listed in the input tree.
}
\item{missing_value}{
Value to be used to denote missing information in the returned arrays, for example to denote the (non-existing) parent of the root node.
}
}
\details{
This function is analogous to the function \code{read.tree.nodes} in the R package \code{phybase} v1.4, but becomes multiple orders of magnitude faster than the latter for large trees (i.e. with 1000-1000,000 tips). Specifically, calling \code{get_clade_list} with \code{postorder=TRUE} and \code{missing_value=-9} on a bifurcating tree yields a similar behavior as calling \code{read.tree.nodes} with the argument ``\code{name}'' set to the tree's tip labels.
 
The input tree can include monofurcations, bifurcations and multifurcations.
The asymptotic average time complexity of this function is O(Nedges), where Nedges is the number of edges in the tree.
}

\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether model fitting succeeded. If \code{FALSE}, the returned list will include an additional ``error'' element (character) providing a description of the error; in that case all other return variables may be undefined.
}
\item{Nsplits}{
The maximum number of children of any node in the tree. For strictly bifurcating trees this will be 2.
}
\item{clades}{
2D integer matrix of size Nclades x (Nsplits+1), with every row representing a specific tip/node in the tree. If \code{postorder==FALSE}, then rows are in the same order as tips/nodes in the original tree, otherwise nodes (but not tips) will be re-ordered and re-indexed in postorder fashion, with the root being the last row. The first column lists the parent node index, the remaining columns list the child tip/node indices. For the root, the parent index will be set to \code{missing_value}; for the tips, the child indices will be set to \code{missing_value}. For nodes with fewer than Nsplits children, superfluous column entries will also be \code{missing_value}.
}
\item{lengths}{
Numeric vector of size Nclades, listing the lengths of the incoming edges at each tip/node in \code{clades}. For the root, the value will be \code{missing_value}. If the tree's \code{edge_length} was \code{NULL}, then \code{lengths} will be \code{NULL} as well.
}
\item{old2new_clade}{
Integer vector of size Nclades, mapping old tip/node indices to tip/node indices in the returned \code{clades} and \code{lengths} arrays. If \code{postorder==FALSE}, this will simply be \code{c(1:Nclades)}.
}
}


\author{Stilianos Louca}

%\references{
%}

%\seealso{
%}

\examples{
# generate a random bifurcating tree
tree = generate_random_tree(list(birth_rate_intercept=1),
                            max_tips=100)$tree

# get tree structure as clade list
# then convert into a similar format as would be
# returned by phybase::read.tree.nodes v1.4
results = get_clade_list(tree,postorder=TRUE,missing_value=-9)
nodematrix = cbind( results$clades, 
                    results$lengths, 
                    matrix(-9,nrow=nrow(results$clades),ncol=3))
phybaseformat = list(   nodes = nodematrix, 
                        names = tree$tip.label, 
                        root  = TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{phylogenetic distance}
