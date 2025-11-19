\name{find_farthest_tip_pair}
\alias{find_farthest_tip_pair}
\title{
Find the two most distant tips in a tree.
}
\description{
Given a phylogenetic tree, find the two most phylogenetically distant tips (to each other) in the tree.
}
\usage{
find_farthest_tip_pair(tree, as_edge_counts = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo". While the tree must be rooted for technical reasons, the outcome does not actually depend on the rooting.
}
\item{as_edge_counts}{
Logical, specifying whether to count phylogenetic distance in terms of edge counts instead of cumulative edge lengths. This is the same as setting all edge lengths to 1.
}
}


\details{
If \code{tree$edge.length} is missing or NULL, then each edge is assumed to have length 1. The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child).

The asymptotic time complexity of this function is O(Nedges), where Nedges is the number of edges in the tree.
}


\value{
A named list with the following elements:
\item{tip1}{
An integer between 1 and Ntips, specifying the first of the two most distant tips.
}
\item{tip2}{
An integer between 1 and Ntips, specifying the second of the two most distant tips.
}
\item{distance}{
Numeric, specifying the phylogenetic (patristic) distance between the \code{farthest_tip1} and \code{farthest_tip2}.
}
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{find_nearest_tips}}, \code{\link{find_farthest_tips}}
}

\examples{
# generate a random tree
Ntips = 1000
parameters = list(birth_rate_intercept=1,death_rate_intercept=0.9)
tree = generate_random_tree(parameters,Ntips,coalescent=FALSE)$tree

# find farthest pair of tips
results = find_farthest_tip_pair(tree)

# print results
cat(sprintf("Tip \%d and \%d have distance \%g\n",
            results$tip1,results$tip2,results$distance))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }
