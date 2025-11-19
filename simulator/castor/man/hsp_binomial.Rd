\name{hsp_binomial}
\alias{hsp_binomial}
\title{
Hidden state prediction for a binary trait based on the binomial distribution.
}
\description{
Estimate the state probabilities for a binary trait at ancestral nodes and tips with unknown (hidden) state, by fitting the probability parameter of a binomial distribution to empirical state frequencies. For each node, the states of its descending tips are assumed to be drawn randomly and independently according to some a prior unknown probability distribution. The probability P1 (probability of any random descending tip being in state 1) is estimated separately for each node based on the observed states in the descending tips via maximum likelihood. 

This function can account for potential state-measurement errors, hidden states and reveal biases (i.e., tips in one particular state being more likely to be measured than in the other state). Only nodes with a number of non-hidden tips above a certain threshold are included in the ML-estimation phase. All other nodes and hidden tips are then assigned the probabilities estimated for the most closely related ancestral node with estimated probabilities. This function is a generalization of \code{\link{hsp_empirical_probabilities}} that can account for potential state-measurement errors and reveal biases.
}
\usage{
hsp_binomial( tree, 
              tip_states,
              reveal_probs  = NULL,
              state1_probs  = NULL,
              min_revealed  = 1,
              max_STE       = Inf,
              check_input   = TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}

\item{tip_states}{
Integer vector of length Ntips, specifying the state of each tip in the tree (either 1 or 2). \code{tip_states} can include \code{NA} to indicate a hidden (non-measured) tip state.
}
\item{reveal_probs}{
2D numeric matrix of size Ntips x 2, listing tip-specific reveal probabilities at each tip conditional on the tip's true state. Hence \code{reveal_probs[n,s]} is the probability that tip \code{n} would have a measured (non-hidden) state if its true state was \code{s}.
May also be a vector of length 2 (same \code{reveal_probs} for all tips) or \code{NULL} (unbiased reveal probs).
}
\item{state1_probs}{
2D numeric matrix of size Ntips x 2, listing the probability of measuring state 1 (potentially erroneously) at each tip conditional upon its true state and conditional upon its state having been measured (i.e., being non-hidden). For example, for an incompletely sequenced genome with completion level \code{C_n} and state 1 indicating presence and state 2 indicating absence of a gene, and assuming error-free detection of genes within the covered regions, one has \code{state1_probs[n,1] = C_n} and \code{state1_probs[n,2]=0}. \code{state1_probs} may also be a vector of length 2 (same probabilities for all tips) or \code{NULL}. If \code{NULL}, state measurements are assumed error-free, and hence this is the same as \code{c(1,0)}.
}
\item{min_revealed}{
Non-negative integer, specifying the minimum number of tips with non-hidden state that must descend from a node for estimating its P1 via maximum likelihood. For nodes with too few descending tips with non-hidden state, the probability P1 will not be estimated via maximum likelihood, and instead will be set to the P1 estimated for the nearest possible ancestral node. It is advised to set this threshold greater than zero (typical values are 2--10).
}
\item{max_STE}{
Non-negative numeric, specifying the maximum acceptable estimated standard error (STE) for the estimated probability P1 for a node. If the STE for a node exceeds this threshold, the P1 for that node is set to the P1 of the nearest ancestor with STE below that threshold. Setting this to \code{Inf} disables this functionality. The STE is estimated based on the Observed Fisher Information Criterion (which, strictly speaking, only provides a lower bound for the STE).
}
\item{check_input}{
Logical, specifying whether to perform some additional time-consuming checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
This function currently only supports binary traits, and states must be represented by integers 1 or 2. Any \code{NA} entries in \code{tip_states} are interpreted as hidden (non-revealed) states.

The algorithm proceeds in two phases ("ASR" phase and "HSP" phase). In the ASR phase the state probability P1 is estimated separately for every node and tip satisfying the thresholds \code{min_revealed} and \code{max_STE}, via maximum-likelihood. In the HSP phase, the P1 of nodes and tips not included in the ASR phase is set to the P1 of the nearest ancestral node with estimated P1, as described by Zaneveld and Thurber (2014).

This function yields estimates for the state probabilities P1 (note that P2=1-P1). In order to obtain point estimates for tip states one needs to interpret these probabilities in a meaningful way, for example by choosing as point estimate for each tip the state with highest probability P1 or P2; the closest that probability is to 1, the more reliable the point estimate will be.

The tree may include multi-furcations (i.e. nodes with more than 2 children) as well as mono-furcations (i.e. nodes with only one child). This function has asymptotic time complexity O(Nedges x Nstates).
Tips must be represented in \code{tip_states} in the same order as in \code{tree$tip.label}. The vector \code{tip_states} need not include names; if it does, however, they are checked for consistency (if \code{check_input==TRUE}).
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether HSP was successful. If \code{FALSE}, an additional element \code{error} (character) will be returned describing the error, while all other return values may be \code{NULL}.
}
\item{P1}{
Numeric vector of length Ntips+Nnodes, listing the estimated probability of being in state 1 for each tip and node. A value of P1[n]=0 or P1[n]=1 means that the n-th tip/node is in state 2 or state 1 with absolute certainty, respectively. Note that even tips with non-hidden state may have have a P1 that is neither 0 or 1, if state measurements are erroneous (i.e., if \code{state1_probs[n,]} differs from \code{(1,0)}). 
}
\item{STE}{
Numeric vector of length Ntips+Nnodes, listing the standard error of the estimated P1 at each tip and node, according to the Observed Fisher Information Criterion. Note that the latter strictly speaking only provides a lower bound on the standard error.
}
\item{reveal_counts}{
Integer vector of length Ntips+Nnodes, listing the number of tips with non-hidden state descending from each tip and node.
}
\item{inheritted}{
Logical vector of length Ntips+Nnodes, specifying for each tip or node whether its returned P1 was directly maximum-likelihood estimated duirng the ASR phase (\code{inheritted[n]==FALSE}) or set to the P1 estimated for an ancestral node during the HSP phase (\code{inheritted[n]==TRUE}).
}
}


\author{Stilianos Louca}

\references{
J. R. Zaneveld and R. L. V. Thurber (2014). Hidden state prediction: A modification of classic ancestral state reconstruction algorithms helps unravel complex symbioses. Frontiers in Microbiology. 5:431.
}

\seealso{
\code{\link{hsp_max_parsimony}}, 
\code{\link{hsp_mk_model}}, 
\code{\link{hsp_empirical_probabilities}}
}

\examples{\dontrun{
# generate random tree
Ntips =50
tree = generate_random_tree(list(birth_rate_factor=1),max_tips=Ntips)$tree

# simulate a binary trait on the tips
Q = get_random_mk_transition_matrix(Nstates=2, rate_model="ER", min_rate=0.1, max_rate=0.5)
tip_states = simulate_mk_model(tree, Q)$tip_states

# print tip states
cat(sprintf("True tip states:\n"))
print(tip_states)

# hide some of the tip states
# include a reveal bias
reveal_probs = c(0.8, 0.3)
revealed = sapply(1:Ntips, FUN=function(n) rbinom(n=1,size=1,prob=reveal_probs[tip_states[n]]))
input_tip_states = tip_states
input_tip_states[!revealed] = NA

# predict state probabilities P1 and P2
hsp = hsp_binomial(tree, input_tip_states, reveal_probs=reveal_probs, max_STE=0.2)
probs = cbind(hsp$P1,1-hsp$P1)

# pick most likely state as a point estimate
# only accept point estimate if probability is sufficiently high
estimated_tip_states = max.col(probs[1:Ntips,])
estimated_tip_states[probs[cbind(1:Ntips,estimated_tip_states)]<0.8] = NA
cat(sprintf("ML-predicted tip states:\n"))
print(estimated_tip_states)

# calculate fraction of correct predictions
predicted = which((!revealed) & (!is.na(estimated_tip_states)))
if(length(predicted)>0){
  Ncorrect  = sum(tip_states[predicted]==estimated_tip_states[predicted])
  cat(sprintf("\%.2g\%\% of predictions are correct\n",(100.0*Ncorrect)/length(predicted)))
}else{
  cat(sprintf("None of the tip states could be reliably predicted\n"))
}
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ancestral state reconstruction}
\keyword{hidden state prediction}
