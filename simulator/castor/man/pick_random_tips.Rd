\name{pick_random_tips}
\alias{pick_random_tips}
\title{
Pick random subsets of tips on a tree.
}
\description{
Given a rooted phylogenetic tree, this function picks random subsets of tips by traversing the tree from root to tips, choosing a random child at each node until reaching a tip. Multiple random independent subsets can be generated if needed.
}
\usage{
pick_random_tips( tree, 
                  size              = 1, 
                  Nsubsets          = 1, 
                  with_replacement  = TRUE, 
                  drop_dims         = TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{size}{
Integer. The size of each random subset of tips.
}
\item{Nsubsets}{
Integer. Number of independent subsets to pick.
}
\item{with_replacement}{
Logical. If \code{TRUE}, each tip can be picked multiple times within a subset (i.e. are "replaced" in the urn). If \code{FALSE}, tips are picked without replacement in each subset. In that case, \code{size} must not be greater than the number of tips in the tree.
}
\item{drop_dims}{
Logical, specifying whether to return a vector (instead of a matrix) if \code{Nsubsets==1}.
}
}


\details{
If \code{with_replacement==TRUE}, then each child of a node is equally probable to be traversed and each tip can be included multiple times in a subset. If \code{with_replacement==FALSE}, then only children with at least one descending tip not included in the subset remain available for traversal; each available child of a node has equal probability to be traversed. In any case, it is always possible for separate subsets to include the same tips.

This random sampling algorithm differs from a uniform sampling of tips at equal probabilities; instead, this algorithm ensures that sister clades have equal probabilities to be picked (if \code{with_replacement==TRUE} or if \code{size}<<Ntips).

The time required by this function per random subset decreases with the number of subsets requested.
}


\value{
A 2D integer matrix of size Nsubsets x size, with each row containing indices of randomly picked tips (i.e. in 1,..,Ntips) within a specific subset. If \code{drop_dims==TRUE} and \code{Nsubsets==1}, then a vector is returned instead of a matrix.
}

\author{Stilianos Louca}

%\references{
%}

%\seealso{
%}

\examples{
# generate random tree
Ntips = 1000
tree  = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# pick random tip subsets
Nsubsets = 100
size     = 50
subsets = pick_random_tips(tree, size, Nsubsets, with_replacement=FALSE)

# count the number of times each tip was picked in a subset ("popularity")
popularities = table(subsets)

# plot histogram of tip popularities
hist(popularities,breaks=20,xlab="popularity",ylab="# tips",main="tip popularities")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{random}
