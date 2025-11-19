\name{tree_imbalance}
\alias{tree_imbalance}
\title{
Calculate various imbalance statistics for a tree.
}
\description{
Given a rooted phylogenetic tree, calculate various "imbalance" statistics of the tree, such as Colless' Index or Sackin's Index.
}
\usage{
tree_imbalance(tree, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
\item{type}{
Character, specifying the statistic to be calculated. Must be one of "Colless" (Shao 1990), "Colless_normalized" (Colless normalized by the maximum possible value in the case of a bifurcating tree), "Sackin" (Sackin 1972) or "Blum" (Blum and Francois 2006, Eq. 5).
}
}


\details{
The tree may include multifurcations and monofurcations. Note that the Colless Index is traditionally only defined for bifurcating trees. For non-bifurcating trees this function calculates a generalization of the index, by summing over all children pairs at each node.
}

\value{
Numeric, the requested imbalance statistic of the tree.
}


\author{Stilianos Louca}

\references{
M. J. Sackin (1972). "Good" and "Bad" Phenograms. Systematic Biology. 21:225-226.

K.T. Shao, R. R. Sokal (1990). Tree Balance. Systematic Biology. 39:266-276.

M. G. B. Blum and O. Francois (2006). Which random processes describe the Tree of Life? A large-scale study of phylogenetic tree imbalance. Systematic Biology. 55:685-691.
}

%\seealso{
%}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# calculate Colless statistic
colless_index = tree_imbalance(tree, type="Colless")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }
