\name{hsp_nearest_neighbor}
\alias{hsp_nearest_neighbor}
\title{
Hidden state prediction based on nearest neighbor.
}
\description{
Predict unknown (hidden) character states of tips on a tree using nearest neighbor matching.
}
\usage{
hsp_nearest_neighbor(tree, tip_states, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}

\item{tip_states}{
A vector of length Ntips, specifying the state of each tip in the tree. Tip states can be any valid data type (e.g., characters, integers, continuous numbers, and so on). \code{NA} values denote unknown (hidden) tip states to be predicted.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
For each tip with unknown state, this function seeks the closest tip with known state, in terms of patristic distance. The state of the closest tip is then used as a prediction of the unknown state. In the case of multiple equal matches, the precise outcome is unpredictable (this is unlikely to occur if edge lengths are continuous numbers, but may happen frequently if e.g. edge lengths are all of unit length). This algorithm is arguably one of the crudest methods for predicting character states, so use at your own discretion.

Any \code{NA} entries in \code{tip_states} are interpreted as unknown states.
If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. The tree may include multifurcations (i.e. nodes with more than 2 children) as well as monofurcations (i.e. nodes with only one child). Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. \code{tip_states} need not include names; if names are included, however, they are checked for consistency with the tree's tip labels (if \code{check_input==TRUE}).
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether HSP was successful. If \code{FALSE}, some return values may be \code{NULL}.
}
\item{states}{
Vector of length Ntips, listing the known and predicted state for each tip.
}
\item{nearest_neighbors}{
Integer vector of length Ntips, listing for each tip the index of the nearest tip with known state. Hence, \code{nearest_neighbors[n]} specifies the tip from which the unknown state of tip n was inferred. If tip n had known state, \code{nearest_neighbors[n]} will be n.
}
\item{nearest_distances}{
Numeric vector of length Ntips, listing for each tip the patristic distance to the nearest tip with known state. For tips with known state, distances will be zero.
}
}


\author{Stilianos Louca}

\references{
J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
}

\seealso{
\code{\link{hsp_max_parsimony}}, 
\code{\link{hsp_mk_model}}, 
}

\examples{\dontrun{
# generate random tree
Ntips = 20
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# simulate a binary trait
Q = get_random_mk_transition_matrix(2, rate_model="ER")
tip_states = simulate_mk_model(tree, Q)$tip_states

# print tip states
print(tip_states)

# set half of the tips to unknown state
tip_states[sample.int(Ntips,size=as.integer(Ntips/2),replace=FALSE)] = NA

# reconstruct all tip states via nearest neighbor
predicted_states = hsp_nearest_neighbor(tree, tip_states)$states

# print predicted tip states
print(predicted_states)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hidden state prediction}
