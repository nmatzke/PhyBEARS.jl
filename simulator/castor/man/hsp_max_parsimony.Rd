\name{hsp_max_parsimony}
\alias{hsp_max_parsimony}
\title{
Hidden state prediction via maximum parsimony.
}
\description{
Reconstruct ancestral discrete states of nodes and predict unknown (hidden) states of tips on a tree using maximum parsimony. Transition costs can vary between transitions, and can optionally be weighted by edge length.
}
\usage{
hsp_max_parsimony(tree, tip_states, Nstates=NULL, 
                  transition_costs="all_equal",
                  edge_exponent=0.0, weight_by_scenarios=TRUE, 
                  check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}

\item{tip_states}{
An integer vector of size Ntips, specifying the state of each tip in the tree as an integer from 1 to Nstates, where Nstates is the possible number of states (see below). \code{tip_states} can include \code{NA} to indicate an unknown tip state that is to be predicted.
}
\item{Nstates}{
Either \code{NULL}, or an integer specifying the number of possible states of the trait. If \code{NULL}, then it will be computed based on the maximum non-\code{NA} value encountered in \code{tip_states}
}
\item{transition_costs}{
Same as for the function \code{\link{asr_max_parsimony}}.
}
\item{edge_exponent}{
Same as for the function \code{\link{asr_max_parsimony}}.
}
\item{weight_by_scenarios}{
Logical, indicating whether to weight each optimal state of a node by the number of optimal maximum-parsimony scenarios in which the node is in that state. If FALSE, then all possible states of a node are weighted equally (i.e. are assigned equal probabilities).
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
For this function, the trait's states must be represented by integers within 1,..,Nstates, where Nstates is the total number of possible states. If the states are originally in some other format (e.g. characters or factors), you should map them to a set of integers 1,..,Nstates. The order of states (if relevant) should be reflected in their integer representation. For example, if your original states are "small", "medium" and "large" and \code{transition_costs=="sequential"}, it is advised to represent these states as integers 1,2,3. You can easily map any set of discrete states to integers using the function \code{\link{map_to_state_space}}.

Any \code{NA} entries in \code{tip_states} are interpreted as unknown states. Prior to ancestral state reconstruction, the tree is temporarily prunned, keeping only tips with known state. The function then applies Sankoff's (1975) dynamic programming algorithm for ancestral state reconstruction, which determines the smallest number (or least costly if transition costs are uneven) of state changes along edges needed to reproduce the known tip states. The state probabilities of tips with unknown state are set to those of the most recent ancestor with reconstructed states, as described by Zaneveld and Thurber (2014). This function has asymptotic time complexity O(Ntips+Nnodes x Nstates).

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. If \code{edge_exponent} is 0, then edge lengths do not influence the result. If \code{edge_exponent!=0}, then all edges must have non-zero length. The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). 

Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. None of the input vectors or matrixes need include row or column names; if they do, however, they are checked for consistency (if \code{check_input==TRUE}).

This function is meant for reconstructing ancestral states in all nodes of a tree as well as predicting the states of tips with an a priory unknown state. If the state of all tips is known and only ancestral state reconstruction is needed, consider using the function \code{\link{asr_max_parsimony}} for improved efficiency.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether HSP was successful. If \code{FALSE}, some return values may be \code{NULL}.
}
\item{likelihoods}{
A 2D numeric matrix, listing the probability of each tip and node being in each state. This matrix will have (Ntips+Nnodes) rows and Nstates columns, where Nstates was either explicitly provided as an argument or inferred based on the number of unique values in \code{tip_states} (if \code{Nstates} was passed as NULL). In the latter case, the column names of this matrix will be the unique values found in \code{tip_states}. The rows in this matrix will be in the order in which tips and nodes are indexed in the tree, i.e. the rows 1,..,Ntips store the probabilities for tips, while rows (Ntips+1),..,(Ntips+Nnodes) store the probabilities for nodes. Each row in this matrix will sum up to 1. Note that the return value is named this way for compatibility with other HSP functions.
}
}


\author{Stilianos Louca}

\references{
D. Sankoff (1975). Minimal mutation trees of sequences. SIAM Journal of Applied Mathematics. 28:35-42.

J. Felsenstein (2004). Inferring Phylogenies. Sinauer Associates, Sunderland, Massachusetts.

J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
}

\seealso{
\code{\link{asr_max_parsimony}}, 
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

# print tip states
print(tip_states)

# set half of the tips to unknown state
tip_states[sample.int(Ntips,size=as.integer(Ntips/2),replace=FALSE)] = NA

# reconstruct all tip states via MPR
likelihoods = hsp_max_parsimony(tree, tip_states, Nstates)$likelihoods
estimated_tip_states = max.col(likelihoods[1:Ntips,])

# print estimated tip states
print(estimated_tip_states)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum parsimony}
\keyword{ancestral state reconstruction}
\keyword{Sankoff's algorithm}
\keyword{hidden state prediction}
