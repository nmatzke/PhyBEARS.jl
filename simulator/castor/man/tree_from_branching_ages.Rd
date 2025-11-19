\name{tree_from_branching_ages}
\alias{tree_from_branching_ages}
\title{
Generate a random timetree with specific branching ages.
}
\description{
Generate a random timetree based on specific branching ages (time before present), by randomly connecting tips and nodes. The tree's root will have the greatest age provided. The tree thus corresponds to a homogenous birth-death model, i.e. where at any given time point all lineages were equally likely to split or go extinct.
}
\usage{
tree_from_branching_ages(   branching_ages,
                            tip_basename       = "",
                            node_basename      = NULL,
                            edge_basename      = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{branching_ages}{
Numeric vector of size Nnodes, listing branching ages (time before present) in ascending order. The last entry will be the root age.
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
Tips in the generated tree are guaranteed to be connected in random order, i.e. this function can also be used to connect a random set of labeled tips into a tree. Nodes will be indexed in chronological order (i.e. in order of decreasing age). In particular, node 0 will be the root.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the tree was successfully generated. If \code{FALSE}, the only other value returned is \code{error}.
}
\item{tree}{
A rooted, ultrametric bifurcating tree of class "phylo", with the requested branching ages.
}
\item{error}{
Character, containing an explanation of the error that occurred. Only included if \code{success==FALSE}.
}
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{tree_from_sampling_branching_ages}}
}


\examples{
Nnodes              = 100
branching_intervals = rexp(n=Nnodes, rate=1)
branching_ages      = cumsum(branching_intervals)
tree                = castor::tree_from_branching_ages(branching_ages)$tree
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{tree model}
\keyword{simulation}
