\name{model_adequacy_hbds}
\alias{model_adequacy_hbds}
\title{
Check if a birth-death-sampling model adequately explains a timetree.
}
\description{
Given a rooted timetree and a homogenous birth-death-sampling model (e.g., as used in molecular epidemiology), check if the model adequately explains various aspects of the tree, such as the branch length and node age distributions and other test statistics. The function uses bootstrapping to simulate multiple hypothetical trees according to the model and then compares the distribution of those trees to the original tree. This function may be used to quantify the "goodness of fit" of a birth-death-sampling model to a timetree. For background on the HBDS model see the documentation for \code{\link{generate_tree_hbds}}.
}
\usage{
model_adequacy_hbds(tree,
                    models,
                    splines_degree      = 1,
                    extrapolate         = FALSE,
                    Nbootstraps         = 1000,
                    max_sim_attempts    = 1000,
                    Nthreads            = 1,
                    max_extant_tips     = NULL,
                    max_model_runtime   = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
A rooted timetree of class "phylo".
}
\item{models}{
\itemize{Either a single HBDS model or a list of HBDS models, specifying the pool of models from which to randomly draw bootstraps. Every model should itself be a named list with some or all of the following elements: 
\item \code{stem_age}: Numeric, the age (time before present) at which the HBDS process started. If \code{NULL}, this is automatically set to the input tree's root age.
\item \code{end_age} : Numeric, the age (time before present) at which the HBDS process halted. This will typically be 0 (i.e., at the tree's youngest tip), however it may also be negative if the process actually halted after the youngest tip was sampled.
\item \code{ages}: Numeric vector, specifying discrete ages (times before present) in ascending order, on which all model variables (e.g., \eqn{\lambda}, \eqn{\mu} and \eqn{\psi}) will be specified. Age increases from tips to root; the youngest tip in the input tree has age 0. The age grid must cover \code{stem_age} and \code{end_age}.
\item \code{lambda}: Numeric vector of the same length as \code{ages}, listing the speciation rate (\eqn{\lambda}) of the HBDS model at the corresponding \code{ages}. Between grid points, the speciation rate is assumed to either be constant (if \code{splines_degree=0}), or linearly (if \code{splines_degree=1}) or quadratically (if \code{splines_degree=2}) or cubically (if \code{splines_degree=3}).
\item \code{mu}: Numeric vector of the same length as \code{ages}, listing the extinction rate (\eqn{\mu}) of the HBDS model at the corresponding \code{ages}. Between grid points, the extinction rate is assumed to either be constant (if \code{splines_degree=0}), or linearly (if \code{splines_degree=1}) or quadratically (if \code{splines_degree=2}) or cubically (if \code{splines_degree=3}). Note that in epidemiological models \eqn{\mu} usually corresponds to the recovery rate plus the death rate of infected hosts. If \code{mu} is not included, it is assumed to be zero.
\item \code{psi}: Optional numeric vector of the same length as \code{ages}, listing the Poissonian sampling rate (\eqn{\mu}) of the HBDS model at the corresponding \code{ages}. Between grid points, the sampling rate is assumed to either be constant (if \code{splines_degree=0}), or linearly (if \code{splines_degree=1}) or quadratically (if \code{splines_degree=2}) or cubically (if \code{splines_degree=3}). If \code{psi} is not included, it is assumed to be zero.
\item \code{kappa}: Optional numeric vector of the same length as \code{ages}, listing the retention probability upon sampling (\eqn{\kappa}) of the HBDS model at the corresponding \code{ages}. Between grid points, the retention probability is assumed to either be constant (if \code{splines_degree=0}), or linearly (if \code{splines_degree=1}) or quadratically (if \code{splines_degree=2}) or cubically (if \code{splines_degree=3}). Note that since \code{kappa} are actual probabilities, they must all be between 0 and 1. If \code{kappa} is not included, it is assumed to be zero.
\item \code{CSA_ages}: Numeric vector listing the ages (time before present) of concentrated sampling attempts, in ascending order. If empty or \code{NULL}, no concentrated sampling attempts are included, i.e. all sampling is assumed to be done according to the Poissonian rate \eqn{\psi}.
\item \code{CSA_probs}: Optional numeric vector, of the same length as \code{CSA_ages}, specifying the sampling probabilities for each concentrated sampling attempt listed in \code{CSA_ages}. Hence, a lineage extant at age \code{CSA_ages[k]} has probability \code{CSA_probs[k]} of being sampled. Note that since \code{CSA_probs} are actual probabilities, they must all be between 0 and 1. \code{CSA_probs} must be provided if and only if \code{CSA_ages} is provided.
\item \code{CSA_kappas}: Optional numeric vector, of the same length as \code{CSA_ages}, specifying the retention probability upon sampling for each concentrated sampling attempt listed in \code{CSA_ages}. Note that since \code{CSA_kappas} are actual probabilities, they must all be between 0 and 1. \code{CSA_kappas} must be provided if and only if \code{CSA_ages} is provided.
}
If you are assessing the adequacy of a single model with specific parameters, then \code{models} can be a single model. If you want to assess the adequacy of a distribution of models, such as sampled from the posterior distribution during a Bayesian analysis, \code{models} should list those posterior models.
}
\item{splines_degree}{
Integer, one of 0, 1, 2 or 3, specifying the polynomial degree of the model parameters \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} between age-grid points. For example, \code{splines_degree=0} means piecewise constant, \code{splines_degree=1} means piecewise linear and so on.
}
\item{extrapolate}{
Logical, specifying whether to extrapolate the model variables \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} (as constants) beyond the provided age grid all the way to \code{stem_age} and \code{end_age} if needed.
}
\item{Nbootstraps}{
Integer, the number of bootstraps (simulations) to perform for calculating statistical significances. A larger number will increase the accuracy of estimated statistical significances.
}
\item{max_sim_attempts}{
Integer, maximum number of simulation attempts per bootstrap, before giving up. Multiple attempts may be needed if the HBDS model has a high probability of leading to extinction early on.
}
\item{Nthreads}{
Integer, number of parallel threads to use for bootstrapping. Note that on Windows machines this option is ignored.
}
\item{max_extant_tips}{
Integer, optional maximum number of extant tips per simulation. A simulation is aborted (and that bootstrap iteration skipped) if the number of extant tips exceeds this threshold. Use this to avoid occasional explosions of runtimes, for example due to very large generated trees.
}
\item{max_model_runtime}{
Numeric, optional maximum computation time (in seconds) to allow for each HBDS model simulation (per bootstrap). Use this to avoid occasional explosions of runtimes, for example due to very large generated trees. Aborted simulations will be omitted from the bootstrap statistics. If \code{NULL} or <=0, this option is ignored.
}
}


\details{
In addition to model selection, the adequacy of any chosen model should also be assessed in absolute terms, i.e. not just relative to other competing models (after all, all considered models might be bad). This function essentially determines how probable it is for hypothetical trees generated by a candidate model (or a distribution of candidate models) to resemble the tree at hand, in terms of various test statistics. In particular, the function uses a Kolmogorov-Smirnov test to check whether the probability distributions of edge lengths and node ages in the tree resemble those expected under the provided models. All statistical significances are calculated using bootstrapping, i.e. by simulating trees from the provided models. For every bootstrap, a model is randomly chosen from the provided \code{models} list.

Note that even if an HBDS model appears to adequately explain a given timetree, this does not mean that the model even approximately resembles the true diversification history (i.e., the true speciation, extinction and sampling rates) that generated the tree (Louca and Pennell 2020). Hence, it is generally more appropriate to say that a given model "congruence class" rather than a specific model explains the tree.

Note that here "age" refers to time before present, i.e. age increases from tips to roots and the youngest tip in the input tree has age 0. In some situations the process that generated the tree (or which is being compared to the tree) might have halted after the last tip was sampled, in which case \code{end_age} should be negative. Similarly, the process may have started prior to the tree's root (e.g., sampled tips coalesce at a later time than when the monitoring started), in which case \code{stem_age} should be greater than the root's age. 

For convenience, it is possible to specify a model without providing an explicit age grid (i.e., omitting \code{ages}); in such a model \eqn{\lambda}, \eqn{\mu}, \eqn{\psi} and \eqn{\kappa} are assumed to be time-independent, and hence \code{lambda}, \code{mu}, \code{psi} and \code{kappa} must be provided as single numerics (or not provided at all).
}

\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the model evaluation was successful. If \code{FALSE}, then an additional return variable, \code{error}, will contain a description of the error; in that case all other return variables may be undefined. Note that \code{success} does not say whether the model explains the tree, but rather whether the computation was performed without errors.
}
\item{Nbootstraps}{
Integer, the number of bootstraps used.
}
\item{tree_Ntips}{
Integer, the number of tips in the original tree.
}
\item{bootstrap_mean_Ntips}{
Numeric, mean number of tips in the bootstrap trees.
}
\item{PNtips}{
Numeric, two-sided statistical significance of the tree's number of tips under the provided null model, i.e. the probability that \code{abs(bootstrap_mean_Ntips-tree_Ntips)} would be as large as observed.
}
\item{tree_Colless}{
Numeric, Colless imbalance statistic (Shao and Sokal, 1990) of the original tree.
}
\item{bootstrap_mean_Colless}{
Numeric, mean Colless statistic across all bootstrap trees.
}
\item{PColless}{
Numeric, two-sided statistical significance of the tree's Colless statistic under the provided null model, i.e. the probability that \code{abs(bootstrap_mean_Colless-tree_Colless)} would be as large as observed.
}
\item{tree_Sackin}{
Numeric, Sackin statistic (Sackin, 1972) of the original tree.
}
\item{bootstrap_mean_Sackin}{
Numeric, median Sackin statistic across all bootstrap trees.
}
\item{PSackin}{
Numeric, two-sided statistical significance of the tree's Sackin statistic under the provided null model, i.e. the probability that \code{abs(bootstrap_mean_Sackin-tree_Sackin)} would be as large as observed.
}
\item{tree_edgeKS}{
Numeric, Kolmogorov-Smirnov (KS) statistic of the original tree's edge lengths, i.e. the estimated maximum difference between the tree's and the model's (estimated) cumulative distribution function of edge lengths.
}
\item{bootstrap_mean_edgeKS}{
Numeric, mean KS statistic of the bootstrap trees' edge lengths.
}
\item{PedgeKS}{
Numeric, the one-sided statistical significance of the tree's edge-length KS statistic, i.e. the probability that the KS statistic of any tree generated by the model would be larger than the original tree's KS statistic. A low value means that the probability distribution of edge lengths in the original tree differs strongly from that expected based on the model.
}
\item{tree_tipKS}{
Numeric, Kolmogorov-Smirnov (KS) statistic of the original tree's tip ages (sampling times before present), i.e. the estimated maximum difference between the tree's and the model's (estimated) cumulative distribution function of tip ages.
}
\item{bootstrap_mean_tipKS}{
Numeric, mean KS statistic of the bootstrap trees' tip ages.
}
\item{PtipKS}{
Numeric, the one-sided statistical significance of the tree's tip-age KS statistic, i.e. the probability that the KS statistic of any tree generated by the model would be larger than the original tree's KS statistic. A low value means that the probability distribution of tip ages in the original tree differs strongly from that expected based on the model.
}
\item{tree_nodeKS}{
Numeric, Kolmogorov-Smirnov (KS) statistic of the original tree's node ages (divergence times before present), i.e. the estimated maximum difference between the tree's and the model's (estimated) cumulative distribution function of node ages.
}
\item{bootstrap_mean_nodeKS}{
Numeric, mean KS statistic of the bootstrap trees' node ages.
}
\item{PnodeKS}{
Numeric, the one-sided statistical significance of the tree's node-age KS statistic, i.e. the probability that the KS statistic of any tree generated by the model would be larger than the original tree's KS statistic. A low value means that the probability distribution of node ages in the original tree differs strongly from that expected based on the model.
}
\item{statistical_tests}{
Data frame, listing the above statistical test results in a more compact format (one test statistic per row).
}
\item{LTT_ages}{
Numeric vector, listing ages (time before present) on which the tree's LTT will be defined.
}
\item{tree_LTT}{
Numeric vector of the same length as \code{LTT_ages}, listing the number of lineages in the tree at every age in \code{LTT_ages}.
}
\item{bootstrap_LTT_CI}{
Named list containing the elements \code{means}, \code{medians}, \code{CI50lower}, \code{CI50upper}, \code{CI95lower} and \code{CI95upper}. Each of these elements is a numeric vector of length equal to \code{LTT_ages}, listing the mean or a specific percentile of LTTs of bootstrap trees at every age in \code{LTT_ages}. For example, \code{bootstrap_LTT_CI$CI95lower[10]} and \code{bootstrap_LTT_CI$CI95upper[10]} define the lower and upper bound, respectively, of the 95\% confidence interval of LTTs generated by the model at age \code{LTT_ages[10]}.
}
\item{fraction_LTT_in_CI95}{
Numeric, fraction of the tree's LTT contained within the equal-tailed 95\%-confidence interval of the distribution of LTT values predicted by the model. For example, a value of 0.5 means that at half of the time points between the present-day and the root, the tree's LTT is contained with the 95\%-CI of predicted LTTs.
}
}


\author{Stilianos Louca}

\references{
S. Louca and M. W. Pennell (2020). Extant timetrees are consistent with a myriad of diversification histories. Nature. 580:502-505.

O. G. Pybus and P. H. Harvey (2000). Testing macro-evolutionary models using incomplete molecular phylogenies. Proceedings of the Royal Society of London. Series B: Biological Sciences. 267:2267-2272.

M. J. Sackin (1972). "Good" and "Bad" Phenograms. Systematic Biology. 21:225-226.

K.T. Shao, R. R. Sokal (1990). Tree Balance. Systematic Biology. 39:266-276.
}

\seealso{
\code{\link{simulate_deterministic_hbds}},
\code{\link{generate_tree_hbds}},
\code{\link{model_adequacy_hbd}}
}

\examples{
\dontrun{
# generate a tree based on a simple HBDS process
max_time = 10
gen = castor::generate_tree_hbds(max_time           = max_time,
                                 lambda             = 1,
                                 mu                 = 0.1,
                                 psi                = 0.1,
                                 no_full_extinction = TRUE)
if(!gen$success) stop(sprintf("Could not generate tree: \%s",gen$error))
tree     = gen$tree
root_age = castor::get_tree_span(tree)$max_distance

# determine age of the stem, i.e. when the HBDS process started
stem_age = gen$root_time + root_age

# determine age at which the HBDS simulation was halted.
# This might be slightly negative, e.g. if the process
# halted after the last sampled tip
end_age = root_age - (gen$final_time-gen$root_time)

# compare the tree to a slightly different model
model = list(stem_age    = stem_age,
             end_age     = end_age,
             lambda      = 1.2,
             mu          = 0.1,
             psi         = 0.2)
adequacy = model_adequacy_hbds( tree, 
                                models = model,
                                Nbootstraps = 100)
if(!adequacy$success){
    cat(sprintf("Adequacy test failed: \%s\n",adequacy$error))
}else{
    print(adequacy$statistical_tests)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{HBDS}
