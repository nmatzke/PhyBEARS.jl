\name{asr_empirical_probabilities}
\alias{asr_empirical_probabilities}
\title{
Empirical ancestral state probabilities.
}
\description{
Given a rooted phylogenetic tree and the states of a discrete trait for each tip, calculate the empirical state frequencies/probabilities for each node in the tree, i.e. the frequencies/probabilities of states across all tips descending from that node. This may be used as a very crude estimate of ancestral state probabilities.

}
\usage{
asr_empirical_probabilities(tree, tip_states, Nstates=NULL, 
                            probabilities=TRUE, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
An integer vector of size Ntips, specifying the state of each tip in the tree as an integer from 1 to Nstates, where Nstates is the possible number of states (see below).
}
\item{Nstates}{
Either \code{NULL}, or an integer specifying the number of possible states of the trait. If \code{NULL}, then it will be computed based on the maximum value encountered in \code{tip_states}
}
\item{probabilities}{
Logical, specifying whether empirical frequencies should be normalized to represent probabilities. If \code{FALSE}, then the raw occurrence counts are returned.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
For this function, the trait's states must be represented by integers within 1,..,Nstates, where Nstates is the total number of possible states. If the states are originally in some other format (e.g., characters or factors), you should map them to a set of integers 1,..,Nstates. You can easily map any set of discrete states to integers using the function \code{\link{map_to_state_space}}.

The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). The function has asymptotic time complexity O(Nedges x Nstates).

Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. The vector \code{tip_states} need not include names; if it does, however, they are checked for consistency (if \code{check_input==TRUE}).
}


\value{
A list with the following elements:
\item{ancestral_likelihoods}{
A 2D integer (if \code{probabilities==FALSE}) or numeric (if \code{probabilities==TRUE}) matrix, listing the frequency or probability of each state for each node. This matrix will have size Nnodes x Nstates, where Nstates was either explicitly provided as an argument or inferred from \code{tip_states}. The rows in this matrix will be in the order in which nodes are indexed in the tree, i.e. the [n,s]-th entry will be the frequency or probability of the s-th state for the n-th node.
Note that the name was chosen for compatibility with other ASR functions.
}
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{asr_max_parsimony}}, 
\code{\link{asr_squared_change_parsimony}}
\code{\link{asr_mk_model}}, 
\code{\link{map_to_state_space}}
}

\examples{
\dontrun{
# generate a random tree
Ntips = 100
tree  = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# create a random transition matrix
Nstates = 3
Q = get_random_mk_transition_matrix(Nstates, rate_model="ER", max_rate=0.01)
cat(sprintf("Simulated ER transition rate=\%g\n",Q[1,2]))

# simulate the trait's evolution
simulation = simulate_mk_model(tree, Q)
tip_states = simulation$tip_states

# calculate empirical probabilities of tip states
asr_empirical_probabilities(tree, tip_states=tip_states, Nstates=Nstates)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum parsimony}
\keyword{ancestral state reconstruction}
\keyword{Sankoff's algorithm}
