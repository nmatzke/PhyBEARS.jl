\name{is_bifurcating}
\alias{is_bifurcating}
\title{
Determine if a tree is bifurcating.
}
\description{
This function determines if a tree is strictly bifurcating, i.e. each node has exactly 2 children. If a tree has monofurcations or multifurcations, this function returns \code{FALSE}.
}
\usage{
is_bifurcating(tree)
}
\arguments{
\item{tree}{
A tree of class "phylo".
}
}
\details{
This functions accepts rooted and unrooted trees, that may include monofurcations, bifurcations and multifurcations.
}
\value{
A logical, indicating whether the input tree is strictly bifurcating.
}


\author{Stilianos Louca}

%\references{
%}

%\seealso{
%}

\examples{
# generate random tree
Ntips = 10
tree  = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# check if the tree is bifurcating (as expected)
is_bifurcating(tree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bifurcation}
