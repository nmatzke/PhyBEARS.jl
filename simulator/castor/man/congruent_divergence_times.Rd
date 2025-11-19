\name{congruent_divergence_times}
\alias{congruent_divergence_times}
\title{
Extract dating anchors for a target tree, using a dated reference tree
}
\description{
Given a reference tree and a target tree, this function maps target nodes to concordant reference nodes when possible, and extracts divergence times of the mapped reference nodes from the reference tree. This function can be used to define secondary dating constraints for a larger target tree, based on a time-calibrated smaller reference tree (Eastman et al. 2013). This only makes sense if the reference tree is time-calibrated. A provided mapping specifies which and how tips in the target tree correspond to tips in the reference tree.
}
\usage{
congruent_divergence_times(reference_tree, target_tree, mapping)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{reference_tree}{
A rooted tree object of class "phylo". Usually this tree will be time-calibrated (i.e. edge lengths represent time intervals).
}
\item{target_tree}{
A rooted tree object of class "phylo".
}
\item{mapping}{
A table mapping a subset of target tips to a subset of reference tips, as described by Eastman et al (2013). Multiple target tips may map to the same reference tip, but not vice versa (i.e. every target tip can appear at most once in the mapping). In general, a tip mapped to in the reference tree is assumed to represent a monophyletic group of tips in the target tree, although this assumption may be violated in practice (Eastman et al. 2013). 

The \code{mapping} must be in one of the following formats:

Option 1: A 2D integer array of size NM x 2 (with NM being the number of mapped target tips), listing target tip indices mapped to reference tip indices (mapping[m,1] (target tip) --> mapping[m,2] (reference tip)). 

Option 2: A 2D character array of size NM x 2, listing target tip labels mapped to reference tip labels. 

Option 3: A data frame of size NM x 1, whose row names are target tip labels and whose entries are either integers (reference tip indices) or characters (reference tip labels). This is the format used by \code{geiger::congruify.phylo} (v.206). 

Option 4: A vector of size NM, whose names are target tip labels and whose entries are either integers (reference tip indices) or characters (reference tip labels).
}
}

\details{
Both the reference and target tree may include monofurcations and/or multifurcations. In principle, neither of the two trees needs to be ultrametric, although in most applications \code{reference_tree} will be ultrametric. 

In special cases each reference tip may be found in the target tree, i.e. the reference tree is a subtree of the target tree. This may occur e.g. if a smaller subtree of the target tree has been extracted and dated, and subsequently the larger target tree is to be dated using secondary constraints inferred from the dated subtree.

The function returns a table that maps a subset of target nodes to an equally sized subset of concordant reference nodes. Ages (divergence times) of the mapped reference nodes are extracted and associated with the concordant target nodes.

For bifurcating trees the average time complexity of this function is O(TNtips x log(RNtips) x NM), where TNtips and RNtips are the number of tips in the target and reference tree, respectively. This function is similar to \code{geiger::congruify.phylo} (v.206). For large trees, this function tends to be much faster than \code{geiger::congruify.phylo}.
}

\value{
A named list with the following elements:
\item{Rnodes}{
Integer vector of length NC (where NC is the number of concordant node pairs found) and with values in 1,..,RNnodes, listing indices of reference nodes that could be matched with (i.e. were concordant to) a target node. Entries in \code{Rnodes} will correspond to entries in \code{Tnodes} and \code{ages}.
}
\item{Tnodes}{
Integer vector of length NC and with values in 1,..,TNnodes, listing indices of target nodes that could be matched with (i.e. were concordant to) a reference node. Entries in \code{Tnodes} will correspond to entries in \code{Rnodes} and \code{ages}.
}
\item{ages}{
Numeric vector of length NC, listing divergence times (ages) of the reference nodes listed in \code{Rnodes}. These ages can be used as fixed anchors for time-calibrating the target tree using a separate program (such as \code{PATHd8}).
}
}

\author{Stilianos Louca}

\references{
J. M. Eastman, L. J. Harmon, D. C. Tank (2013). Congruification: support for time scaling large phylogenetic trees. Methods in Ecology and Evolution. 4:688-691.
}

\seealso{
\code{\link{extend_tree_to_height}},
\code{\link{date_tree_red}},
\code{\link{get_tips_for_mrcas}},
\code{\link{tree_distance}}
}

\examples{
# generate random tree (target tree)
Ntips = 10000
tree  = castor::generate_random_tree(parameters=list(birth_rate_intercept=1), max_tips=Ntips)$tree

# extract random subtree (reference tree)
Nsubtips    = 10
subtips     = sample.int(n=Ntips,size=Nsubtips,replace=FALSE)
subtreeing  = castor::get_subtree_with_tips(tree, only_tips=subtips)
subtree     = subtreeing$subtree

# map subset of target tips to reference tips
mapping = matrix(c(subtreeing$new2old_tip,(1:Nsubtips)),ncol=2,byrow=FALSE)

# extract divergence times by congruification
congruification = congruent_divergence_times(subtree, tree, mapping)

cat("Concordant target nodes:\n")
print(congruification$target_nodes)

cat("Ages of concordant nodes:\n")
print(congruification$ages)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dating}
\keyword{ultrametric}
\keyword{tree comparison}

