\name{castor-package}
\alias{castor-package}
\alias{castor}
\docType{package}
\title{
Efficient computations on large phylogenetic trees.
}

\description{
This package provides efficient tree manipulation functions including pruning, rerooting, calculation of most-recent common ancestors, calculating distances from the tree root and calculating pairwise distance matrices. Calculation of phylogenetic signal and mean trait depth (trait conservatism). Efficient ancestral state reconstruction and hidden character prediction of discrete characters on phylogenetic trees, using Maximum Likelihood and Maximum Parsimony methods. Simulating models of trait evolution, and generating random trees.
}

\details{
The most important data unit is a phylogenetic tree of class "phylo", with the tree topology encoded in the member variable \code{tree.edge}. See the \code{ape} package manual for details on the "phylo" format. The castor package was designed to be efficient for large phylogenetic trees (>10,000 tips), and scales well to trees with millions of tips. Most functions have asymptotically linear time complexity O(N) in the number of edges N. This efficiency is achived via temporary auxiliary data structures, use of dynamic programing, heavy use of C++, and integer-based indexing instead of name-based indexing of arrays. All functions support trees that include monofurcations (nodes with a single child) as well as multifurcations (nodes with more than 2 children). See the associated paper by Louca et al. for a comparison with other packages.

Throughout this manual, "Ntips" refers to the number of tips, "Nnodes" to the number of nodes and "Nedges" to the number of edges in a tree. In the context of discrete trait evolution/reconstruction, "Nstates" refers to the number of possible states of the trait. In the context of multivariate trait evolution, "Ntraits" refers to the number of traits.
}
\author{
Stilianos Louca

Maintainer: Stilianos Louca <louca@zoology.ubc.ca>
}
\references{
S. Louca and M. Doebeli (2017). Efficient comparative phylogenetics on large trees. Bioinformatics. DOI:10.1093/bioinformatics/btx701
}
%\examples{
%  \dontrun{
%     ## Optional simple examples of the most important functions
%     ## These can be in \dontrun{} and \donttest{} blocks.   
%  }
%}
