\name{asr_subtree_averaging}
\alias{asr_subtree_averaging}
\title{
Ancestral state reconstruction via subtree averaging.
}
\description{
Reconstruct ancestral states in a phylogenetic tree for a continuous (numeric) trait by averaging trait values over descending subtrees. That is, for each node the reconstructed state is set to the arithmetic average state of all tips descending from that node.
}
\usage{
asr_subtree_averaging(tree, tip_states, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
A numeric vector of size Ntips, specifying the known state of each tip in the tree.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
The function returns the estimated ancestral states (=averages) as well as the corresponding standard deviations. Note that reconstructed states are local estimates, i.e. they only take into account the tips descending from the reconstructed node.

The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). Edge lengths and distances between tips and nodes are not taken into account. All tip states are assumed to be known, and \code{NA} or \code{NaN} are not allowed in \code{tip_states}.

Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. The vector \code{tip_states} need not include item names; if it does, however, they are checked for consistency (if \code{check_input==TRUE}).
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether ASR was sucessful. If all input data are valid then this will always be \code{TRUE}, but it is provided for consistency with other ASR functions.
}
\item{ancestral_states}{
A numeric vector of size Nnodes, listing the reconstructed state (=average over descending tips) for each node. The entries in this vector will be in the order in which nodes are indexed in the tree.
}
\item{ancestral_stds}{
A numeric vector of size Nnodes, listing the standard deviations corresponding to \code{ancestral_stds}.
}
\item{ancestral_counts}{
A numeric vector of size Nnodes, listing the number of (descending) tips used to reconstruct the state of each node.
}
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{asr_independent_contrasts}},
\code{\link{asr_squared_change_parsimony}}
}

\examples{
# generate random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# simulate a continuous trait
tip_states = simulate_ou_model(tree, stationary_mean=0, spread=1, decay_rate=0.001)$tip_states

# reconstruct node states by averaging simulated tip states
node_states = asr_subtree_averaging(tree, tip_states)$ancestral_states
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ancestral state reconstruction}
