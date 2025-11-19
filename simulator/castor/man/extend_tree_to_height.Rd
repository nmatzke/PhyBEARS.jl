\name{extend_tree_to_height}
\alias{extend_tree_to_height}
\title{
Extend a rooted tree up to a specific height.
}
\description{
Given a rooted phylogenetic tree and a specific distance from the root (``new height''), elongate terminal edges (i.e. leading into tips) as needed so that all tips have a distance from the root equal to the new height. If a tip already extends beyond the specified new height, its incoming edge remains unchanged.
}
\usage{
extend_tree_to_height(tree, new_height=NULL)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{new_height}{
Numeric, specifying the phylogenetic distance from the root to which tips are to be extended. If \code{NULL} or negative, then it is set to the maximum distance of any tip from the root.
}
}


\details{
The input tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). All tip, edge and node indices remain unchanged. This function provides a quick-and-dirty way to make a tree ultrametric, or to correct small numerical inaccuracies in supposed-to-be ultrametric trees.
}


\value{
A list with the following elements:
\item{tree}{
A new rooted tree of class "phylo", representing the extended tree.
}
\item{max_extension}{
Numeric. The largest elongation added to a terminal edge.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{trim_tree_at_height}}
}

\examples{
# generate a random non-ultrametric tree
tree = generate_random_tree(list(birth_rate_intercept=1,death_rate_intercept=0.5),
                            max_time=1000,
                            coalescent=FALSE)$tree
                            
# print min & max distance from root
span = get_tree_span(tree)
cat(sprintf("Min & max tip height = \%g & \%g\n",span$min_distance,span$max_distance))

# make tree ultrametric by extending terminal edges
extended = extend_tree_to_height(tree)$tree

# print new min & max distance from root
span = get_tree_span(extended)
cat(sprintf("Min & max tip height = \%g & \%g\n",span$min_distance,span$max_distance))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ultrametric}
