\name{asr_max_parsimony}
\alias{asr_max_parsimony}
\title{
Maximum-parsimony ancestral state reconstruction.
}
\description{
Reconstruct ancestral states for a discrete trait using maximum parsimony. Transition costs can vary between transitions, and can optionally be weighted by edge length.
}
\usage{
asr_max_parsimony(tree, tip_states, Nstates=NULL, 
                  transition_costs="all_equal",
                  edge_exponent=0, weight_by_scenarios=TRUE,
                  check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
An integer vector of size Ntips, specifying the state of each tip in the tree as an integer from 1 to Nstates, where Nstates is the possible number of states (see below).
}
\item{Nstates}{
Either \code{NULL}, or an integer specifying the number of possible states of the trait. If \code{NULL}, then Nstates will be computed based on the maximum value encountered in \code{tip_states}
}
\item{transition_costs}{
Either "all_equal","sequential", "proportional", "exponential", or a quadratic non-negatively valued matrix of size Nstates x Nstates, specifying the transition costs between all possible states (which can include 0 as well as \code{Inf}). The [r,c]-th entry of the matrix is the cost of transitioning from state r to state c. The option "all_equal" specifies that all transitions are permitted and are equally costly. "sequential" means that only transitions between adjacent states are permitted from a node to its children, and all permitted transitions are equally costly. "proportional" means that all transitions are permitted, but the cost increases proportional to the distance between states. "exponential" means that all transitions are permitted, but the cost increases exponentially with the distance between states. The options "sequential" and "proportional" only make sense if states exhibit an order relation (as reflected in their integer representation).
}
\item{edge_exponent}{
Non-negative real-valued number. Optional exponent for weighting transition costs by the inverse length of edge lengths. If 0, edge lengths do not influence the ancestral state reconstruction (this is the conventional max-parsimony). If >0, then at each edge the transition costs are multiplied by \eqn{1/L^e}, where \eqn{L} is the edge length and \eqn{e} is the edge exponent. This parameter is mostly experimental; modify at your own discretion.
}
\item{weight_by_scenarios}{
Logical, indicating whether to weight each optimal state of a node by the number of optimal maximum-parsimony scenarios in which the node is in that state. If FALSE, then all optimal states of a node are weighted equally (i.e. are assigned equal probabilities).
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
For this function, the trait's states must be represented by integers within 1,..,Nstates, where Nstates is the total number of possible states. If the states are originally in some other format (e.g. characters or factors), you should map them to a set of integers 1,..,Nstates. The order of states (if relevant) should be reflected in their integer representation. For example, if your original states are "small", "medium" and "large" and \code{transition_costs=="sequential"}, it is advised to represent these states as integers 1,2,3. You can easily map any set of discrete states to integers using the function \code{\link{map_to_state_space}}.

This function utilizes Sankoff's (1975) dynamic programming algorithm for determining the smallest number (or least costly if transition costs are uneven) of state changes along edges needed to reproduce the observed tip states. The function has asymptotic time complexity O(Ntips+Nnodes x Nstates).

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. If \code{edge_exponent} is 0, then edge lengths do not influence the result. If \code{edge_exponent!=0}, then all edges must have non-zero length. The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). 

Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. None of the input vectors or matrixes need include row or column names; if they do, however, they are checked for consistency (if \code{check_input==TRUE}).

This function is meant for reconstructing ancestral states in all nodes of a tree, when the state of each tip is known. If some of the tips have unknown state, consider either pruning the tree to keep only tips with known states, or using the function \code{\link{hsp_max_parsimony}}.

Not all datasets are consistent with all possible transition cost models, i.e., it could happen that for some peculiar datasets some rather constrained models (e.g. "sequential") cannot possibly produce the data. In this case, castor will most likely return non-sensical ancestral state estimates and \code{total_cost=Inf}, although this has not thoroughly been tested.
}


\value{
A list with the following elements:
\item{success}{Boolean, indicating whether ASR was successful. If \code{FALSE}, the remaining returned elements may be undefined.}
\item{ancestral_likelihoods}{
A 2D numeric matrix, listing the probability of each node being in each state. This matrix will have size Nnodes x Nstates, where Nstates was either explicitly provided as an argument or inferred from \code{tip_states}. The rows in this matrix will be in the order in which nodes are indexed in the tree, i.e. the [n,s]-th entry will be the probability of the s-th state for the n-th node. These probabilities are calculated based on \code{scenario_counts} (see below), assuming that every maximum parsimony scenario is equally likely. Note that the name was chosen for compatibility with other ASR functions.
}
\item{scenario_counts}{
A 2D numeric matrix of size Nnodes x Nstates, listing for each node and each state the number of maximum parsimony scenarios in which the node was in the specific state. If only a single maximum parsimony scenario exists for the whole tree, then the sum of entries in each row will be one.
}
\item{total_cost}{Real number, specifying the total transition cost across the tree for the most parsimonious scenario. In the classical case where \code{transition_costs="all_equal"}, the \code{total_cost} equals the total number of state changes in the tree under the most parsimonious scenario. Under some constrained transition models (e.g., "sequential"), \code{total_cost} may sometimes be Inf, which basically means that the data violates the model.}
}


\author{Stilianos Louca}

\references{
D. Sankoff (1975). Minimal mutation trees of sequences. SIAM Journal of Applied Mathematics. 28:35-42.

J. Felsenstein (2004). Inferring Phylogenies. Sinauer Associates, Sunderland, Massachusetts.
}

\seealso{
\code{\link{hsp_max_parsimony}}, 
\code{\link{asr_squared_change_parsimony}}
\code{\link{asr_mk_model}}, 
\code{\link{hsp_mk_model}}, 
\code{\link{map_to_state_space}}
}

\examples{\dontrun{
# generate random tree
Ntips = 10
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# simulate a discrete trait
Nstates = 5
Q = get_random_mk_transition_matrix(Nstates, rate_model="ER")
tip_states = simulate_mk_model(tree, Q)$tip_states

# reconstruct node states via MPR
results = asr_max_parsimony(tree, tip_states, Nstates)
node_states = max.col(results$ancestral_likelihoods)

# print reconstructed node states
print(node_states)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum parsimony}
\keyword{ancestral state reconstruction}
\keyword{Sankoff's algorithm}
