\name{fit_mk}
\alias{fit_mk}
\title{
Fit a Markov (Mk) model for discrete trait evolution.
}
\description{
Estimate the transition rate matrix of a continuous-time Markov model for discrete trait evolution ("Mk model") via maximum-likelihood, based on one or more phylogenetic trees and its tips' states.
}
\usage{
fit_mk( trees,
        Nstates,
        tip_states              = NULL,
        tip_priors              = NULL,
        rate_model              = "ER",
        root_prior              = "auto",
        oldest_ages             = NULL,
        guess_transition_matrix = NULL,
        Ntrials                 = 1,
        max_model_runtime       = NULL,
        optim_algorithm         = "nlminb",
        optim_max_iterations    = 200,
        optim_rel_tol           = 1e-8,
        check_input             = TRUE,
        Nthreads                = 1,
        Nbootstraps             = 0,	
        Ntrials_per_bootstrap   = NULL,
        verbose                 = FALSE,
        verbose_prefix          = "")
        
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{trees}{
Either a single phylogenetic tree of class "phylo", or a list of phylogenetic trees. Edge lengths should correspond (or be analogous) to time. The trees don't need to be ultrametric.
}
\item{Nstates}{
Integer, specifying the number of possible discrete states that the trait can have.
}
\item{tip_states}{
Either an integer vector of size Ntips (only permitted if trees[] is a single tree) or a list containing Ntrees such integer vectors (if trees[] is a list of trees), listing the state of each tip in each tree. Note that tip_states cannot include NAs or NaNs; if the states of some tips are uncertain, you should use the option \code{tip_priors} instead. Can also be \code{NULL}, in which case \code{tip_priors} must be provided.
}
\item{tip_priors}{
Either a numeric matrix of size Ntips x Nstates (only permitted if trees[] is a single tree), or a list containing Ntrees such matrixes (if trees[] is a list of trees), listing the likelihood of each state at each tip in each tree. Can also be \code{NULL}, in which case \code{tip_states} must be provided. Hence, \code{tip_priors[t][i,s]} is the likelihood of the observed state of tip i in tree t, if the tip's true state was in state s. For example, if you know for certain that a tip is in state k, then set \code{tip_priors[t][i,s]=1} for s=k and \code{tip_priors[t][i,s]=0} for all other s.
}
\item{rate_model}{
Rate model to be used for the transition rate matrix. Can be "ER" (all rates equal), "SYM" (transition rate i-->j is equal to transition rate j-->i), "ARD" (all rates can be different), "SUEDE" (only stepwise transitions i-->i+1 and i-->i-1 allowed, all 'up' transitions are equal, all 'down' transitions are equal) or "SRD" (only stepwise transitions i-->i+1 and i-->i-1 allowed, and each rate can be different). Can also be an index matrix that maps entries of the transition matrix to the corresponding independent rate parameter to be fitted. Diagonal entries should map to 0, since diagonal entries are not treated as independent rate parameters but are calculated from the remaining entries in the transition rate matrix. All other entries that map to 0 represent a transition rate of zero. The format of this index matrix is similar to the format used by the \code{ace} function in the \code{ape} package. \code{rate_model} is only relevant if \code{transition_matrix==NULL}.
}
\item{root_prior}{
Prior probability distribution of the root's states, used to calculate the model's overall likelihood from the root's marginal ancestral state likelihoods. Can be "\code{flat}" (all states equal), "\code{empirical}" (empirical probability distribution of states across the tree's tips), "\code{stationary}" (stationary probability distribution of the transition matrix), "\code{likelihoods}" (use the root's state likelihoods as prior), "\code{max_likelihood}" (put all weight onto the state with maximum likelihood) or ``\code{auto}'' (will be chosen automatically based on some internal logic). If "\code{stationary}" and \code{transition_matrix==NULL}, then a transition matrix is first fitted using a flat root prior, and then used to calculate the stationary distribution. \code{root_prior} can also be a non-negative numeric vector of size Nstates and with total sum equal to 1.
}
\item{oldest_ages}{
Optional numeric or numeric vector of size Ntrees, specifying the oldest age (time before present) for each tree to consider when fitting the Mk model. If \code{NULL}, the entire trees are considered from the present all the way to their root. If non-\code{NULL}, then each tree is ``cut'' at the corresponding oldest age, yielding multiple subtrees, each of which is assumed to be an independent realization of the Mk process. If \code{oldest_ages} is a single numeric, then all trees are cut at the same oldest age. This option may be useful if temporal variation is suspected in the Mk rates, and only data near the present are to be used for fitting to avoid violating the assumptions of a constant-rates Mk model.
}
\item{guess_transition_matrix}{
Optional 2D numeric matrix, specifying a reasonable first guess for the transition rate matrix. May contain \code{NA}. May also be \code{NULL}, in which case a reasonable first guess is automatically generated.
}
\item{Ntrials}{
Number of trials (starting points) for fitting the transition rate matrix. A higher number may reduce the risk of landing in a local non-global optimum of the likelihood function, but will increase computation time during fitting.
}
\item{max_model_runtime}{
Optional positive numeric, specifying the maximum time (in seconds) allowed for a single evaluation of the likelihood function. If a specific Mk model takes longer than this threshold to evaluate, then its likelihood is set to -Inf. This option can be used to avoid badly parameterized models during fitting and can thus reduce fitting time. If NULL or <=0, this option is ignored.
}
\item{optim_algorithm}{
Either "optim" or "nlminb", specifying which optimization algorithm to use for maximum-likelihood estimation of the transition matrix.
}
\item{optim_max_iterations}{
Maximum number of iterations (per fitting trial) allowed for optimizing the likelihood function.
}
\item{optim_rel_tol}{
Relative tolerance (stop criterion) for optimizing the likelihood function.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
\item{Nthreads}{
Number of parallel threads to use for running multiple fitting trials simultaneously. This only makes sense if your computer has multiple cores/CPUs and if \code{Ntrials>1}. This option is ignored on Windows, because Windows does not support forking.
}
\item{Nbootstraps}{
Integer, specifying the number of parametric bootstraps to perform for estimating standard errors and confidence intervals of estimated rate parameters. Set to 0 for no bootstrapping.
}
\item{Ntrials_per_bootstrap}{
Integer, specifying the number of fitting trials to perform for each bootstrap sampling. If \code{NULL}, this is set equal to \code{max(1,Ntrials)}. Decreasing \code{Ntrials_per_bootstrap} will reduce computation time, at the expense of potentially inflating the estimated confidence intervals; in some cases (e.g., for very large trees) this may be useful if fitting takes a long time and confidence intervals are very narrow anyway. Only relevant if \code{Nbootstraps>0}.
}
\item{verbose}{
Logical, specifying whether to print progress reports and warnings to the screen.
}
\item{verbose_prefix}{
Character, specifying the line prefix for printing progress reports to the screen.
}
}


\details{
The trait's states must be represented by integers within 1,..,Nstates, where Nstates is the total number of possible states. If the states are originally in some other format (e.g. characters or factors), you should map them to a set of integers 1,..,Nstates. The order of states (if relevant) should be reflected in their integer representation. For example, if your original states are "small", "medium" and "large" and \code{rate_model=="SUEDE"}, it is advised to represent these states as integers 1,2,3. You can easily map any set of discrete states to integers using the function \code{\link{map_to_state_space}}.

This function allows the specification of the precise tip states (if these are known) using the vector \code{tip_states}. Alternatively, if some tip states are not fully known, you can pass the state likelihoods using the matrix \code{tip_priors}. Note that exactly one of the two arguments, \code{tip_states} or \code{tip_priors}, must be non-\code{NULL}. 

Tips must be represented in \code{tip_states} or \code{tip_priors} in the same order as in \code{tree$tip.label}. None of the input vectors or matrixes need include row or column names; if they do, however, they are checked for consistency (if \code{check_input==TRUE}).

The tree is either assumed to be complete (i.e. include all possible species), or to represent a random subset of species chosen independently of their states. If the tree is not complete and tips are not chosen independently of their states, then this method will not be valid.

\code{fit_Mk} uses maximum-likelihood to estimate each free parameter of the transition rate matrix. The number of free parameters depends on the \code{rate_model} considered; for example, \code{ER} implies a single free parameter, while \code{ARD} implies Nstates x (Nstates-1) free parameters. If multiple trees are provided as input, the likelihood is the product of likelihoods for each tree, i.e. as if each tree was an independent realization of the same Markov process.

This function is similar to \code{\link{asr_mk_model}}, but focused solely on fitting the transition rate matrix (i.e., without estimating ancestral states) and with the ability to utilize multiple trees at once.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the fitting was successful. If \code{FALSE}, an additional element \code{error} (of type character) is included containing an explanation of the error; in that case the value of any of the other elements is undetermined.
}
\item{Nstates}{
Integer, the number of states assumed for the model.
}
\item{transition_matrix}{
A matrix of size Nstates x Nstates, the fitted transition rate matrix of the model. The [r,c]-th entry is the transition rate from state r to state c.
}
\item{loglikelihood}{
Numeric, the log-likelihood of the observed tip states under the fitted model.
}
\item{Niterations}{
Integer, the number of iterations required to reach the maximum log-likelihood. Depending on the optimization algorithm used (see \code{optim_algorithm}), this may be \code{NA}.
}
\item{Nevaluations}{
Integer, the number of evaluations of the likelihood function required to reach the maximum log-likelihood. Depending on the optimization algorithm used (see \code{optim_algorithm}), this may be \code{NA}.
}
\item{converged}{
Logical, indicating whether the fitting algorithm converged. Note that \code{fit_Mk} may return successfully even if convergence was not achieved; if this happens, the fitted transition matrix may not be reasonable. In that case it is recommended to change the optimization options, for example increasing \code{optim_max_iterations}.
}
\item{guess_rate}{
Numeric, the initial guess used for the average transition rate, prior to fitting.
}
\item{AIC}{
Numeric, the Akaike Information Criterion for the fitted model, defined as \eqn{2k-2\log(L)}, where \eqn{k} is the number of independent fitted parameters and \eqn{L} is the maximized likelihood.
}
\item{standard_errors}{
Numeric matrix of size Nstates x Nstates, estimated standard error of the fitted transition rates, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50lower}{
Numeric matrix of size Nstates x Nstates, lower bounds of the 50\% confidence intervals (25-75\% percentile) for the fitted transition rates, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50upper}{
Numeric matrix of size Nstates x Nstates, upper bounds of the 50\% confidence intervals for the fitted transition rates, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95lower}{
Numeric matrix of size Nstates x Nstates, lower bounds of the 95\% confidence intervals (2.5-97.5\% percentile) for the fitted transition rates, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95upper}{
Numeric matrix of size Nstates x Nstates, upper bounds of the 95\% confidence intervals for the fitted transition rates, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
}

\author{Stilianos Louca}

\references{
Z. Yang, S. Kumar and M. Nei (1995). A new method for inference of ancestral nucleotide and amino acid sequences. Genetics. 141:1641-1650.

M. Pagel (1994). Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters. Proceedings of the Royal Society of London B: Biological Sciences. 255:37-45.
}

\seealso{
\code{\link{asr_mk_model}},
\code{\link{simulate_mk_model}},
\code{\link{fit_musse}}
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

# fit Mk transition matrix
results = fit_mk(tree, Nstates, tip_states, rate_model="ER", Ntrials=2)

# print Mk model fitting summary
cat(sprintf("Mk model: log-likelihood=\%g\n",results$loglikelihood))
cat(sprintf("Fitted ER transition rate=\%g\n",results$transition_matrix[1,2]))
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Mk model}
\keyword{ancestral state reconstruction}
\keyword{rerooting}
\keyword{maximum likelihood}
