\name{asr_mk_model}
\alias{asr_mk_model}
\title{
Ancestral state reconstruction with Mk models and rerooting
}
\description{
Ancestral state reconstruction of a discrete trait using a fixed-rates continuous-time Markov model (a.k.a. "Mk model"). This function can estimate the (instantaneous) transition matrix using maximum likelihood, or take a specified transition matrix. The function can optionally calculate marginal ancestral state likelihoods for each node in the tree, using the rerooting method by Yang et al. (1995).
}
\usage{
asr_mk_model( tree, 
              tip_states, 
              Nstates               = NULL, 
              tip_priors            = NULL, 
              rate_model            = "ER", 
              transition_matrix     = NULL, 
              include_ancestral_likelihoods = TRUE, 
              reroot                = TRUE,
              root_prior            = "auto", 
              Ntrials               = 1, 
              optim_algorithm       = "nlminb",
              optim_max_iterations  = 200,
              optim_rel_tol         = 1e-8,
              store_exponentials    = TRUE,
              check_input           = TRUE, 
              Nthreads              = 1)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{tip_states}{
An integer vector of size Ntips, specifying the state of each tip in the tree in terms of an integer from 1 to Nstates, where Ntips is the number of tips and Nstates is the number of possible states (see below). Can also be \code{NULL}. If \code{tip_states==NULL}, then \code{tip_priors} must not be \code{NULL} (see below).
}
\item{Nstates}{
Either NULL, or an integer specifying the number of possible states of the trait. If \code{Nstates==NULL}, then it will be computed based on the maximum value encountered in \code{tip_states} or based on the number of columns in \code{tip_priors} (whichever is non-\code{NULL}).
}
\item{tip_priors}{
A 2D numeric matrix of size Ntips x Nstates, where Nstates is the possible number of states for the character modelled. Hence, \code{tip_priors[i,s]} is the likelihood of the observed state of tip i, if the tip's true state was in state s. For example, if you know for certain that a tip is in state k, then set \code{tip_priors[i,s]=1} for s=k and \code{tip_priors[i,s]=0} for all other s.
}
\item{rate_model}{
Rate model to be used for fitting the transition rate matrix. Can be "ER" (all rates equal), "SYM" (transition rate i-->j is equal to transition rate j-->i), "ARD" (all rates can be different), "SUEDE" (only stepwise transitions i-->i+1 and i-->i-1 allowed, all 'up' transitions are equal, all 'down' transitions are equal) or "SRD" (only stepwise transitions i-->i+1 and i-->i-1 allowed, and each rate can be different). Can also be an index matrix that maps entries of the transition matrix to the corresponding independent rate parameter to be fitted. Diagonal entries should map to 0, since diagonal entries are not treated as independent rate parameters but are calculated from the remaining entries in the transition matrix. All other entries that map to 0 represent a transition rate of zero. The format of this index matrix is similar to the format used by the \code{ace} function in the \code{ape} package. \code{rate_model} is only relevant if \code{transition_matrix==NULL}.
}
\item{transition_matrix}{
Either a numeric quadratic matrix of size Nstates x Nstates containing fixed transition rates, or \code{NULL}. The [r,c]-th entry in this matrix should store the transition rate from state r to state c. Each row in this matrix must have sum zero.
If \code{NULL}, then the transition rates will be estimated using maximum likelihood, based on the \code{rate_model} specified.
}
\item{root_prior}{
Prior probability distribution of the root's states, used to calculate the model's overall likelihood from the root's marginal ancestral state likelihoods. Can be "\code{flat}" (all states equal), "\code{empirical}" (empirical probability distribution of states across the tree's tips), "\code{stationary}" (stationary probability distribution of the transition matrix), "\code{likelihoods}" (use the root's state likelihoods as prior) or "\code{max_likelihood}" (put all weight onto the state with maximum likelihood). If "\code{stationary}" and \code{transition_matrix==NULL}, then a transition matrix is first fitted using a flat root prior, and then used to calculate the stationary distribution. \code{root_prior} can also be a non-negative numeric vector of size Nstates and with total sum equal to 1.
}
\item{include_ancestral_likelihoods}{
Include the marginal ancestral likelihoods for each node (conditional scaled state likelihoods) in the return values. Note that this may increase the computation time and memory needed, so you may set this to \code{FALSE} if you don't need marginal ancestral states.
}
\item{reroot}{
Reroot tree at each node when computing marginal ancestral likelihoods, according to Yang et al. (1995). This is the default and recommended behavior, but leads to increased computation time. If \code{FALSE}, ancestral likelihoods at each node are computed solely based on the subtree descending from that node, without rerooting. Caution: Rerooting is strictly speaking only valid if the Mk model is time-reversible (for example, if the transition matrix is symmetric). Do not use the rerooting method if \code{rate_model="ARD"}.
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
Logical, specifying whether to pre-calculate and store exponentials of the transition matrix during calculation of ancestral likelihoods. This may reduce computation time because each exponential is only calculated once, but requires more memory since all exponentials are stored. 

Only relevant if \code{include_ancestral_likelihoods==TRUE}, otherwise exponentials are never stored.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
\item{Nthreads}{
Number of parallel threads to use for running multiple fitting trials simultaneously. This only makes sense if your computer has multiple cores/CPUs and if \code{Ntrials>1}, and is only relevant if \code{transition_matrix==NULL}. This option is ignored on Windows, because Windows does not support forking.
}
}


\details{
For this function, the trait's states must be represented by integers within 1,..,Nstates, where Nstates is the total number of possible states. If the states are originally in some other format (e.g. characters or factors), you should map them to a set of integers 1,..,Nstates. The order of states (if relevant) should be reflected in their integer representation. For example, if your original states are "small", "medium" and "large" and \code{rate_model=="SUEDE"}, it is advised to represent these states as integers 1,2,3. You can easily map any set of discrete states to integers using the function \code{\link{map_to_state_space}}.

This function allows the specification of the precise tip states (if these are known) using the vector \code{tip_states}. Alternatively, if some tip states are only known in terms of a probability distribution, you can pass these probability distributions using the matrix \code{tip_priors}. Note that exactly one of the two arguments, \code{tip_states} or \code{tip_priors}, must be non-\code{NULL}. 

Tips must be represented in \code{tip_states} or \code{tip_priors} in the same order as in \code{tree$tip.label}. None of the input vectors or matrixes need include row or column names; if they do, however, they are checked for consistency (if \code{check_input==TRUE}).

The tree is either assumed to be complete (i.e. include all possible species), or to represent a random subset of species chosen independently of their states. The rerooting method by Yang et al (1995) is used to calculate the marginal ancestral state likelihoods for each node by treating the node as a root and calculating its conditional scaled likelihoods. Note that the re-rooting algorithm is strictly speaking only valid for reversible Mk models, that is, satisfying the criterion
\deqn{
\pi_i Q_{ij}=\pi_j Q_{ji},\quad\forall i,j,
}
where \eqn{Q} is the transition rate matrix and \eqn{\pi} is the stationary distribution of the model. The rate models ``ER'', `SYM'', ``SUEDE'' and ``SRD'' are reversible. For example, for ``SUEDE'' or ``SRD'' choose \eqn{\pi_{i+1}=\pi_{i}Q_{i,i+1}/Q_{i+1,i}}. In contrast, ``ARD'' models are generally not reversible.

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). This function is similar to \code{rerootingMethod} in the \code{phytools} package (v0.5-64) and similar to \code{ape::ace} (v4.1) with options \code{method="ML", type="discrete"} and \code{marginal=FALSE}, but tends to be much faster than \code{rerootingMethod} and \code{ace} for large trees.

Internally, this function uses \code{\link{fit_mk}} to estimate the transition matrix if the latter is not provided. If you only care about estimating the transition matrix but not the ancestral states, consider using the more versatile function \code{\link{fit_mk}}.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether ASR was successful. If \code{FALSE}, all other return values may be \code{NULL}.
}
\item{Nstates}{
Integer, specifying the number of modeled trait states.
}
\item{transition_matrix}{
A numeric quadratic matrix of size Nstates x Nstates, containing the transition rates of the Markov model. The [r,c]-th entry is the transition rate from state r to state c. Will be the same as the input \code{transition_matrix}, if the latter was not \code{NULL}.
}
\item{loglikelihood}{
Log-likelihood of the observed tip states under the fitted (or provided) Mk model. If \code{transition_matrix} was \code{NULL} in the input, then this will be the log-likelihood maximized during fitting.
}
\item{AIC}{
Numeric, the Akaike Information Criterion for the fitted Mk model (if applicable), defined as \eqn{2k-2\log(L)}, where \eqn{k} is the number of independent fitted parameters and \eqn{L} is the maximized likelihood. If the transition matrix was provided as input, then no fitting was performed and hence AIC will be \code{NULL}.
}
\item{ancestral_likelihoods}{
Optional, only returned if \code{include_ancestral_likelihoods} was \code{TRUE}. A 2D numeric matrix, listing the likelihood of each state at each node (marginal ancestral likelihoods). This matrix will have size Nnodes x Nstates, where Nstates was either explicitly provided as an argument, or inferred from \code{tip_states} or \code{tip_priors} (whichever was non-\code{NULL}). The rows in this matrix will be in the order in which nodes are indexed in the tree, i.e. the [n,s]-th entry will be the likelihood of the s-th state at the n-th node. For example, \code{likelihoods[1,3]} will store the likelihood of observing the tree's tip states (if \code{reroot=TRUE}) or the descending subtree's tip states (if \code{reroot=FALSE}), if the first node was in state 3. Note that likelihoods are rescaled (normalized) to sum to 1 for convenience and numerical stability. The marginal likelihoods at a node should not, however, be interpreted as a probability distribution among states.
}
}



\author{Stilianos Louca}

\references{
Z. Yang, S. Kumar and M. Nei (1995). A new method for inference of ancestral nucleotide and amino acid sequences. Genetics. 141:1641-1650.
}

\seealso{
\code{\link{hsp_mk_model}}, 
\code{\link{asr_max_parsimony}},
\code{\link{asr_squared_change_parsimony}},
\code{\link{hsp_max_parsimony}},
\code{\link{map_to_state_space}},
\code{\link{fit_mk}}
}

\examples{\dontrun{
# generate random tree
Ntips = 1000
tree  = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# create random transition matrix
Nstates = 5
Q = get_random_mk_transition_matrix(Nstates, rate_model="ER", max_rate=0.01)
cat(sprintf("Simulated ER transition rate=\%g\n",Q[1,2]))

# simulate the trait's evolution
simulation = simulate_mk_model(tree, Q)
tip_states = simulation$tip_states
cat(sprintf("Simulated states for last 20 nodes:\n"))
print(tail(simulation$node_states,20))

# reconstruct node states from simulated tip states
# at each node, pick state with highest marginal likelihood
results = asr_mk_model(tree, tip_states, Nstates, rate_model="ER", Ntrials=2)
node_states = max.col(results$ancestral_likelihoods)

# print Mk model fitting summary
cat(sprintf("Mk model: log-likelihood=\%g\n",results$loglikelihood))
cat(sprintf("Fitted ER transition rate=\%g\n",results$transition_matrix[1,2]))

# print reconstructed node states for last 20 nodes
print(tail(node_states,20))
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Mk model}
\keyword{ancestral state reconstruction}
\keyword{rerooting}
\keyword{maximum likelihood}
