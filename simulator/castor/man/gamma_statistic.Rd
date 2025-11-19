\name{gamma_statistic}
\alias{gamma_statistic}
\title{
Calculate the gamma-statistic of a tree.
}
\description{
Given a rooted ultrametric phylogenetic tree, calculate the gamma-statistic (Pybus and Harevy, 2000).
}
\usage{
gamma_statistic(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted tree of class "phylo". The tree is assumed to be ultrametric; any deviations from ultrametricity are ignored.
}
}


\details{
The tree may include multifurcations and monofurcations. If edge lengths are missing (i.e. \code{edge.length=NULL}), then each edge is assumed to have length 1.

This function is similar to the function \code{gammaStat} in the R package \code{ape} v5.3.
}

\value{
Numeric, the gamma-statistic of the tree.
}


\author{Stilianos Louca}

\references{
O. G. Pybus and P. H. Harvey (2000). Testing macro-evolutionary models using incomplete molecular phylogenies. Proceedings of the Royal Society of London. Series B: Biological Sciences. 267:2267-2272.
}

%\seealso{
%}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# calculate & print gamma statistic
gammastat = gamma_statistic(tree)
cat(sprintf("Tree has gamma-statistic \%g\n",gammastat))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }
