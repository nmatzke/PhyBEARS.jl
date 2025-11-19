\name{tree_from_sampling_branching_ages}
\alias{tree_from_sampling_branching_ages}
\title{
Generate a random timetree with specific tip/sampling and node/branching ages.
}
\description{
Generate a random bifurcating timetree based on specific sampling (tip) ages and branching (node) ages, by randomly connecting tips and nodes. Age refers to time before present, i.e., measured in reverse chronological direction. The tree's root will have the greatest age provided. The tree thus corresponds to a homogenous birth-death-sampling model, i.e. where at any given time point all lineages were equally likely to split, be sampled or go extinct.
}
\usage{
tree_from_sampling_branching_ages(sampling_ages,
                                  branching_ages,
                                  tip_basename  = "",
                                  node_basename = NULL,
                                  edge_basename = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{sampling_ages}{
Numeric vector of size Ntips, listing sampling ages (time before present) in ascending order.
}
\item{branching_ages}{
Numeric vector of size Nnodes, listing branching ages (time before present) in ascending order. The last entry will be the root age. Note that Nnodes must be equal to Ntips-1.
}
\item{tip_basename}{
Character. Prefix to be used for tip labels (e.g. "tip."). If empty (""), then tip labels will be integers "1", "2" and so on.
}
\item{node_basename}{
Character. Prefix to be used for node labels (e.g. "node."). If \code{NULL}, no node labels will be included in the tree.
}
\item{edge_basename}{
Character. Prefix to be used for edge labels (e.g. "edge."). Edge labels (if included) are stored in the character vector \code{edge.label}. If \code{NULL}, no edge labels will be included in the tree.
}
}


\details{
Tips and nodes will be indexed in chronological order (i.e. in order of decreasing age). In particular, node 0 will be the root. Note that not all choices of \code{sampling_ages} and \code{branching_ages} are permissible. Specifically, at any given age T, the number of sampling events with age equal or smaller than T must be greater than the number of branching events with age equal or smaller than T. If this requirement is not satisfied, the function will return with an error.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the tree was successfully generated. If \code{FALSE}, the only other value returned is \code{error}.
}
\item{tree}{
A rooted, ultrametric bifurcating tree of class "phylo", with the requested tip and node ages.
}
\item{error}{
Character, containing an explanation of the error that occurred. Only included if \code{success==FALSE}.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{tree_from_branching_ages}}
}

\examples{
sampling_ages   = c(0, 0.1, 0.15, 0.25, 0.9, 1.9, 3)
branching_ages  = c(0.3, 0.35, 0.4, 1.1, 2.5, 3.5)
tree = tree_from_sampling_branching_ages(sampling_ages, branching_ages)$tree
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{tree model}
\keyword{simulation}
