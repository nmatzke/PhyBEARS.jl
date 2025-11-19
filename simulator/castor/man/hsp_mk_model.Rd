\name{hsp_mk_model}
\alias{hsp_mk_model}
\title{
Hidden state prediction with Mk models and rerooting
}
\description{
Reconstruct ancestral states of a discrete trait and predict unknown (hidden) states of tips using a fixed-rates continuous-time Markov model (a.k.a. "Mk model"). This function can fit the model (i.e. estimate the transition matrix) using maximum likelihood, or use a specified transition matrix. The function can optionally calculate marginal ancestral state likelihoods for each node in the tree, using the rerooting method by Yang et al. (1995). A subset of the tips may have completely unknown states; in this case the fitted Markov model is used to predict their state likelihoods based on their most recent reconstructed ancestor, as described by Zaneveld and Thurber (2014). The function can account for biases in which tips have known state (``reveal bias'').
}
\usage{
hsp_mk_model( tree, 
              tip_states, 
              Nstates = NULL, 
              reveal_fractions = NULL,
              tip_priors = NULL, 
              rate_model = "ER", 
              transition_matrix = NULL, 
              include_likelihoods = TRUE,
              root_prior = "empirical", 
              Ntrials = 1, 
              optim_algorithm = "nlminb",
              optim_max_iterations = 200,
              optim_rel_tol = 1e-8,
              store_exponentials = TRUE,
              check_input = TRUE, 
              Nthreads = 1)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
An integer vector of size Ntips, specifying the state of each tip in the tree in terms of an integer from 1 to Nstates, where Nstates is the possible number of states (see below). Can also be \code{NULL}, in which case \code{tip_priors} must not be \code{NULL} (see below). \code{tip_states} can include \code{NA} to indicate an unknown (hidden) tip state that is to be predicted.
}
\item{Nstates}{
Either NULL, or an integer specifying the number of possible states of the trait. If \code{Nstates==NULL}, then it will be computed based on the maximum non-\code{NA} value encountered in \code{tip_states} or based on the number of columns in \code{tip_priors} (whichever is non-\code{NULL}).
}
\item{reveal_fractions}{
Either NULL, or a numeric vector of size Nstates, specifying the fraction of tips with revealed (i.e., non-hidden) state, depending on the tip state. That is, \code{reveal_fractions[s]} is the probability that a given tip at state \code{s} will have known (i.e., non-hidden) state, conditional upon being included in the tree. If the tree only contains a random subset of species (sampled independently of each species' state), then \code{reveal_fractions[s]} is the probability of knowing the state of a species (regardless of whether it is included in the tree), if its state is \code{s}. This variable can be used to account for biases in which tips have known state, depending on their state. Only the relative ratios among reveal fractions matter, i.e. multiplying \code{reveal_fractions} with a constant factor has no effect.
}
\item{tip_priors}{
A 2D numeric matrix of size Ntips x Nstates, where Nstates is the possible number of states for the character modelled. Can also be \code{NULL}. Each row of this matrix must be a probability vector, i.e. it must only contain non-negative entries and must sum up to 1. The [i,s]-th entry should be the prior probability of tip i being in state s. If you know for certain that tip i is in some state s, you can set the corresponding entry to 1 and all other entries in that row to 0. A row can include \code{NA} to indicate that neither the state nor the probability distribution of a state are known for that tip. If for all tips you either know the exact state or have no information at all, you can also use \code{tip_states} instead. If \code{tip_priors==NULL}, then \code{tip_states} must not be \code{NULL} (see above).
}
\item{rate_model}{
Rate model to be used for fitting the transition rate matrix. Similar to the \code{rate_model} option in the function \code{\link{asr_mk_model}}. See the details of \code{\link{asr_mk_model}} on the assumptions of each \code{rate_model}.
}
\item{transition_matrix}{
Either a numeric quadratic matrix of size Nstates x Nstates containing fixed transition rates, or \code{NULL}. The [r,c]-th entry in this matrix should store the transition (probability) rate from the state r to state c. Each row in this matrix must have sum zero.
If \code{NULL}, then the transition rates will be estimated using maximum likelihood, based on the \code{rate_model} specified.
}
\item{include_likelihoods}{
Boolean, specifying whether to include the marginal state likelihoods for all tips and nodes, as returned variables. Setting this to \code{TRUE} can substantially increase computation time. If \code{FALSE}, the Mk model is merely fitted, but ancestral states and hidden tip states are not reconstructed.
}
\item{root_prior}{
Prior probability distribution of the root's states. Similar to the \code{root_prior} option in the function \code{\link{asr_mk_model}}.
}
\item{Ntrials}{
Number of trials (starting points) for fitting the transition matrix. Only relevant if \code{transition_matrix=NULL}. A higher number may reduce the risk of landing in a local non-global optimum of the likelihood function, but will increase computation time during fitting.
}
\item{optim_algorithm}{
Either "optim" or "nlminb", specifying which optimization algorithm to use for maximum-likelihood estimation of the transition matrix. Only relevant if \code{transition_matrix==NULL}.
}
\item{optim_max_iterations}{
Maximum number of iterations (per fitting trial) allowed for optimizing the likelihood function.
}
\item{optim_rel_tol}{
Relative tolerance (stop criterion) for optimizing the likelihood function.
}
\item{store_exponentials}{
Logical, specifying whether to pre-calculate and store exponentials of the transition matrix during calculation of ancestral likelihoods. This may reduce computation time because each exponential is only calculated once, but will use up more memory since all exponentials are stored. Only relevant if \code{include_ancestral_likelihoods} is \code{TRUE}, otherwise exponentials are never stored.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
\item{Nthreads}{
Number of parallel threads to use for running multiple fitting trials simultaneously. This only makes sense if your computer has multiple cores/CPUs and \code{Ntrials>1}, and is only relevant if \code{transition_matrix==NULL}.
}
}


\details{
For this function, the trait's states must be represented by integers within 1,..,Nstates, where Nstates is the total number of possible states. Note that Nstates can be chosen to be larger than the number of states observed in the tips of the present tree, to account for potential states not yet observed. If the trait's states are originally in some other format (e.g. characters or factors), you should map them to a set of integers 1,..,Nstates. The order of states (if applicable) should be reflected in their integer representation. For example, if your original states are "small", "medium" and "large" and \code{rate_model=="SUEDE"}, it is advised to represent these states as integers 1,2,3. You can easily map any set of discrete states to integers using the function \code{\link{map_to_state_space}}.

This function allows the specification of the precise tip states (if these are known) using the vector \code{tip_states}. Alternatively, if some tip states are only known in terms of a probability distribution, you can pass these probability distributions using the matrix \code{tip_priors}. Note that exactly one of the two arguments, \code{tip_states} or \code{tip_priors}, must be non-\code{NULL}. In either case, the presence of \code{NA} in \code{tip_states} or in a row of \code{tip_priors} is interpreted as an absence of information about the tip's state (i.e. the tip has "hidden state").

Tips must be represented in \code{tip_states} or \code{tip_priors} in the same order as in \code{tree$tip.label}. None of the input vectors or matrixes need include row or column names; if they do, however, they are checked for consistency (if \code{check_input==TRUE}).

This method assumes that the tree is either complete (i.e. includes all species), or that the tree's tips represent a random subset of species that have been sampled independent of their state. The function does not require that tip state knowledge is independent of tip state, provided that the associated biases are known (provided via \code{reveal_fractions}). The rerooting method by Yang et al (2015) is used to reconstruct the marginal ancestral state likelihoods for each node by treating the node as a root and calculating its conditional scaled likelihoods. The state likelihoods of tips with hidden states are calculated from those of the most recent ancestor with previously calculated state likelihoods, using the exponentiated transition matrix along the connecting edges (essentially using the rerooting method). Attention: The state likelihoods for tips with known states or with provided priors are not modified, i.e. they are as provided in the input. In other words, for those tips the returned state likelihoods should not be considered as posteriors in a Bayesian sense.

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child).
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether HSP was successful. If \code{FALSE}, some return values may be \code{NULL}.
}
\item{Nstates}{
Integer, specifying the number of modeled trait states.
}
\item{transition_matrix}{
A numeric quadratic matrix of size Nstates x Nstates, containing the transition rates of the Markov model. The [r,c]-th entry is the transition rate from state r to state c. Will be the same as the input \code{transition_matrix}, if the latter was not \code{NULL}.
}
\item{loglikelihood}{
Log-likelihood of the Markov model. If \code{transition_matrix} was \code{NULL} in the input, then this will be the log-likelihood maximized during fitting.
}
\item{likelihoods}{
A 2D numeric matrix, listing the probability of each tip and node being in each state. Only included if \code{include_likelihoods} was \code{TRUE}. This matrix will have (Ntips+Nnodes) rows and Nstates columns, where Nstates was either explicitly provided as an argument, or inferred from \code{tip_states} or \code{tip_priors} (whichever was non-\code{NULL}). The rows in this matrix will be in the order in which tips and nodes are indexed in the tree, i.e. rows 1,..,Ntips store the probabilities for tips, while rows (Ntips+1),..,(Ntips+Nnodes) store the probabilities for nodes. For example, \code{likelihoods[1,3]} will store the probability that tip 1 is in state 3. Each row in this matrix will sum up to 1. Note that for tips with known state or fully provided prior, the likelihoods will be unchanged, i.e. these are not the posteriors in a Bayesian sense.
}
}



\author{Stilianos Louca}

\references{
Z. Yang, S. Kumar and M. Nei (1995). A new method for inference of ancestral nucleotide and amino acid sequences. Genetics. 141:1641-1650.

J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
}

\seealso{
\code{\link{hsp_max_parsimony}}, 
\code{\link{hsp_squared_change_parsimony}}, 
\code{\link{asr_mk_model}}, 
\code{\link{map_to_state_space}}
}

\examples{\dontrun{
# generate random tree
Ntips = 1000
tree  = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# simulate a discrete trait
Nstates = 5
Q = get_random_mk_transition_matrix(Nstates, rate_model="ER", max_rate=0.01)
tip_states = simulate_mk_model(tree, Q)$tip_states
cat(sprintf("Simulated ER transition rate=\%g\n",Q[1,2]))

# print states for first 20 tips
print(tip_states[1:20])

# set half of the tips to unknown state
# chose tips randomly, regardless of their state (no biases)
tip_states[sample.int(Ntips,size=as.integer(Ntips/2),replace=FALSE)] = NA

# reconstruct all tip states via Mk model max-likelihood
results = hsp_mk_model(tree, tip_states, Nstates, rate_model="ER", Ntrials=2, Nthreads=2)
estimated_tip_states = max.col(results$likelihoods[1:Ntips,])

# print Mk model fitting summary
cat(sprintf("Mk model: log-likelihood=\%g\n",results$loglikelihood))
cat(sprintf("Universal (ER) transition rate=\%g\n",results$transition_matrix[1,2]))

# print estimated states for first 20 tips
print(estimated_tip_states[1:20])
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Mk model}
\keyword{ancestral state reconstruction}
\keyword{rerooting}
\keyword{maximum likelihood}
\keyword{hidden state prediction}
