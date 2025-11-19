\name{fit_musse}
\alias{fit_musse}
\title{
Fit a discrete-state-dependent diversification model via maximum-likelihood.
}
\description{
The Binary State Speciation and Extinction (BiSSE) model (Maddison et al. 2007) and its extension to Multiple State Speciation Extinction (MuSSE) models (FitzJohn et al. 2009, 2012), Hidden State Speciation Extinction (HiSSE) models (Beaulieu and O'meara, 2016) or Several Examined and Concealed States-dependent Speciation and Extinction (SecSSE) models (van Els et al. 2018), describe a Poissonian cladogenic process whose birth/death (speciation/extinction) rates depend on the states of an evolving discrete trait. Specifically, extant tips either go extinct or split continuously in time at Poissonian rates, and birth/death rates at each extant tip depend on the current state of the tip; lineages tansition stochastically between states acccording to a continuous-time Markov process with fixed transition rates. At the end of the simulation (i.e., at "present-day"), extant lineages are sampled according to some state-dependent probability ("sampling_fraction"), which may depend on proxy state. Optionally, tips may also be sampled continuously over time according to some Poissonian rate (which may depend on proxy state), in which case the resulting tree may not be ultrametric.

This function takes as main input a phylogenetic tree (ultrametric unless Poissonian sampling is included) and a list of tip proxy states, and fits the parameters of a BiSSE/MuSSE/HiSSE/SecSSE model to the data via maximum-likelihood. Tips can have missing (unknown) proxy states, and the function can account for biases in species sampling and biases in the identification of proxy states. The likelihood is calculated using a mathematically equivalent, but computationally more efficient variant, of the classical postorder-traversal BiSSE/MuSSE/HiSSE/SecSSE algorithm, as described by Louca (2019). This function has been optimized for large phylogenetic trees, with a relatively small number of states (i.e. Nstates<<Ntips); its time complexity scales roughly linearly with Ntips.
}
\usage{
fit_musse(tree, 
          Nstates,
          NPstates              = NULL,
          proxy_map             = NULL,
          state_names           = NULL,
          tip_pstates           = NULL,
          tip_priors            = NULL,
          sampling_fractions    = 1,
          reveal_fractions      = 1,
          sampling_rates        = 0,
          transition_rate_model	= "ARD",
          birth_rate_model      = "ARD",
          death_rate_model      = "ARD",
          transition_matrix     = NULL,
          birth_rates           = NULL,
          death_rates           = NULL,
          first_guess           = NULL,
          lower                 = NULL,
          upper                 = NULL,
          root_prior            = "auto",
          root_conditioning     = "auto",
          oldest_age            = NULL,
          Ntrials               = 1,
          max_start_attempts    = 10,
          optim_algorithm       = "nlminb",
          optim_max_iterations  = 10000,
          optim_max_evaluations = NULL,
          optim_rel_tol         = 1e-6,
          check_input           = TRUE,
          include_ancestral_likelihoods = FALSE,
          Nthreads              = 1,
          Nbootstraps           = 0,
          Ntrials_per_bootstrap = NULL,
          max_condition_number  = 1e4,
          relative_ODE_step     = 0.1,
          E_value_step          = 1e-4,
          D_temporal_resolution = 100,
          max_model_runtime     = NULL,
          verbose               = TRUE,
          diagnostics           = FALSE,
          verbose_prefix        = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{tree}{
Phylogenetic tree of class "phylo", representing the evolutionary relationships between sampled species/lineages. Unless Poissonian sampling is included in the model (option \code{sampling_rates}), this tree should be ultrametric.
}
\item{Nstates}{
Integer, specifying the number of possible discrete states a tip can have, influencing speciation/extinction rates. For example, if \code{Nstates==2} then this corresponds to the common Binary State Speciation and Extinction (BiSSE) model (Maddison et al., 2007). In the case of a HiSSE/SecSSE model, \code{Nstates} refers to the total number of diversification rate categories. For example, in the case of the HiSSE model described by Beaulieu and O'meara (2016), \code{Nstates=4}.
}
\item{NPstates}{
Integer, optionally specifying a number of "proxy-states" that are observed instead of the underlying speciation/extinction-modulating states. To fit a HiSSE/SecSSE model, \code{NPstates} should be smaller than \code{Nstates}. Each state corresponds to a different proxy-state, as defined using the variable \code{proxy_map} (see below). For BiSSE/MuSSE with no hidden states, \code{NPstates} can be set to either \code{NULL} or equal to \code{Nstates}; in either case, \code{NPstates} will be considered equal to \code{Nstates}. For example, in the case of the HiSSE model described by Beaulieu and O'meara (2016), \code{NPstates=2}.
}
\item{proxy_map}{
Integer vector of size \code{Nstates} and with values in 1,..\code{NPstates}, specifying the correspondence between states (i.e. diversification-rate categories) and proxy-states, in a HiSSE/SecSSE model. Specifically, \code{proxy_map[s]} indicates which proxy-state the state s is represented by. Each proxy-state can represent multiple states (i.e. proxies are ambiguous), but each state must be represented by exactly one proxy-state. For example, to setup the HiSSE model described by Beaulieu and O'meara (2016), use \code{proxy_map=c(1,2,1,2)}. For non-HiSSE models, set this to \code{NULL} or to \code{c(1:Nstates)}. See below for more details.
}
\item{state_names}{
Optional character vector of size \code{Nstates}, specifying a name/description for each state. This does not influence any of the calculations. It is merely used to add human-readable row/column names (rather than integers) to the returned vectors/matrices. If \code{NULL}, no row/column names are added.
}
\item{tip_pstates}{
Integer vector of size Ntips, listing the proxy state at each tip, in the same order as tips are indexed in the tree. The vector may (but need not) include names; if it does, these are checked for consistency with the tree (if \code{check_input==TRUE}). Values must range from 1 to \code{NPstates} (which is assumed equal to \code{Nstates} in the case of BiSSE/MuSSE). States may also be \code{NA}, corresponding to unknown tip proxy states (no information available).
}
\item{tip_priors}{
Numeric matrix of size Ntips x \code{Nstates} (or of size Ntips x \code{NPstates}), listing prior likelihoods of each state (or each proxy-state) at each tip. Can be provided as an alternative to \code{tip_pstates}. Thus, \code{tip_priors[i,s]} is the likelihood of observing the data (i.e., sampling tip i and observing the observed state) if the tip \code{i} was at state \code{s} (or proxy-state \code{s}). Hence, \code{tip_priors} should account for sampling fractions as well as reveal fractions. Either \code{tip_pstates} or \code{tip_priors} must be non-NULL, but not both.
}
\item{sampling_fractions}{
Numeric vector of size \code{NPstates}, with values between 0 and 1, listing the sampling fractions of extant species depending on proxy-state. That is, \code{sampling_fractions[p]} is the probability that an extant species, having proxy state \code{p}, is included in the phylogeny at present-day. If all extant species are included in the tree with the same probability (i.e., independent of state), this can also be a single number. If \code{NULL} (default), all extant species are assumed to be included in the tree. Irrelevant if \code{tip_priors} is provided and valid for all tips.
}
\item{reveal_fractions}{
Numeric vector of size \code{NPstates}, with values between 0 and 1, listing the probabilities of proxy-state identification depending on proxy-state. That is, \code{reveal_fractions[p]} is the probability that a species with proxy-state \code{p} will have a known ("revealed") state, conditional upon being included in the tree. This can be used to incorporate reveal biases for tips, depending on their proxy state. Can also be \code{NULL} or a single number (in which case reveal fractions are assumed to be independent of proxy-state). Note that only the relative values in \code{reveal_fractions} matter, for example \code{c(1,2,1)} has the same effect as \code{c(0.5,1,0.5)}, because \code{reveal_fractions} is normalized internally anyway. Irrelevant if \code{tip_priors} is provided and valid for all tips.
}
\item{sampling_rates}{
Numeric vector of size \code{NPstates}, listing Poissonian per-lineage sampling rates over time. Hence, \code{sampling_rates[p]} is the rate at which lineages are sampled over time when they are in proxy state \code{p}. Can also be a single numeric, in which case sampling rates are the same for all proxy states. If \code{NULL}, Poissonian sampling is assumed to not occur. Note that earlier MuSSE/HiSSE models (e.g., by Beaulieu and O'Meara, 2016) do not include Poissonian sampling (i.e., all tips are assumed to have been sampled at present-day). Poissonian sampling through time is common in epidemiological models but uncommon in macroevolution models.
}
\item{transition_rate_model}{
Either a character or a 2D integer matrix of size Nstates x Nstates, specifying the model for the transition rates between states. This option controls the parametric complexity of the state transition model, i.e. the number of independent rates and the correspondence between independent and dependent rates. If a character, then it must be one of "ER", "SYM", "ARD", "SUEDE" or "SRD", as used for Mk models (see the function \code{\link{asr_mk_model}} for details). For example, "ARD" (all rates different) specifies that all transition rates should be considered as independent parameters with potentially different values.

If an integer matrix, then it defines a custom parametric structure for the transition rates, by mapping entries of the transition matrix to a set of independent transition-rate parameters (numbered 1,2, and so on), similarly to the option \code{rate_model} in the function \code{\link{asr_mk_model}}, and as returned for example by the function \code{\link{get_transition_index_matrix}}. Entries must be between 1 and Nstates, however 0 may also be used to denote a fixed value of zero. For example, if \code{transition_rate_model[1,2]=transition_rate_model[2,1]}, then the transition rates 1->2 and 2->1 are assumed to be equal.
Entries on the diagonal are ignored, since the diagonal elements are always adjusted to ensure a valid Markov transition matrix. To construct a custom matrix with the proper structure, it may be convenient to first generate an "ARD" matrix using \code{\link{get_transition_index_matrix}}, and then modify individual entries to reduce the number of independent rates.
}
\item{birth_rate_model}{
Either a character or an integer vector of length Nstates, specifying the model for the various birth (speciation) rates. This option controls the parametric complexity of the possible birth rates, i.e. the number of independent birth rates and the correspondence between independent and dependent birth rates. If a character, then it must be either "ER" (equal rates) or "ARD" (all rates different). If an integer vector, it must map each state to an indepedent birth-rate parameter (indexed 1,2,..). For example, the vector \code{c(1,2,1)} specifies that the birth-rates \eqn{\lambda_1} and \eqn{\lambda_3} must be the same, but \eqn{\lambda_2} is independent.
}
\item{death_rate_model}{
Either a character or an integer vector of length Nstates, specifying the model for the various death (extinction) rates. Similar to \code{birth_rate_model}.
}
\item{transition_matrix}{
Either \code{NULL} or a 2D matrix of size Nstates x Nstates, specifying known (and thus fixed) transition rates between states. For example, setting some elements to 0 specifies that these transitions cannot occur directly. May also contain \code{NA}, indicating rates that are to be fitted. If \code{NULL} or empty, all rates are considered unknown and are therefore fitted. Note that, unless \code{transition_rate_model=="ARD"}, values in \code{transition_matrix} are assumed to be consistent with the rate model, that is, rates specified to be equal under the transition rate model are expected to also have equal values in \code{transition_matrix}.
}
\item{birth_rates}{
Either \code{NULL}, or a single number, or a numeric vector of length Nstates, specifying known (and thus fixed) birth rates for each state. May contain \code{NA}, indicating rates that are to be fitted. For example, the vector \code{c(5,0,NA)} specifies that \eqn{\lambda_1=5}, \eqn{\lambda_2=0} and that \eqn{\lambda_3} is to be fitted. If \code{NULL} or empty, all birth rates are considered unknown and are therefore fitted. If a single number, all birth rates are considered fixed at that given value.
}
\item{death_rates}{
Either \code{NULL}, or a single number, or a numeric vector of length Nstates, specifying known (and thus fixed) death rates for each state. Similar to \code{birth_rates}.
}
\item{first_guess}{
Either \code{NULL}, or a named list containing optional initial suggestions for various model parameters, i.e. start values for fitting. The list can contain any or all of the following elements:
\itemize{
\item\code{transition_matrix}: A single number or a 2D numeric matrix of size Nstates x Nstates, specifying suggested start values for the transition rates. May contain NA, indicating rates that should be guessed automatically by the function. If a single number, then that value is used as a start value for all transition rates.
\item\code{birth_rates}: A single number or a numeric vector of size Nstates, specifying suggested start values for the birth rates. May contain NA, indicating rates that should be guessed automatically by the function (by fitting a simple birth-death model, see \code{\link{fit_tree_model}}).
\item\code{death_rates}: A single number or a numeric vector of size Nstates, specifying suggested start values for the death rates. May contain NA, indicating rates that should be guessed automatically by the function (by fitting a simple birth-death model, see \code{\link{fit_tree_model}}).
}
Start values are only relevant for fitted (i.e., non-fixed) parameters.
}
\item{lower}{
Either \code{NULL} or a named list containing optional lower bounds for various model parameters. The list can contain any or all of the elements \code{transition_matrix}, \code{birth_rates} and \code{death_rates}, structured similarly to \code{first_guess}. For example, \code{list(transition_matrix=0.1, birth_rates=c(5,NA,NA))} specifies that all transition rates between states must be 0.1 or greater, that the birth rate \eqn{\lambda_1} must be 5 or greater, and that all other model parameters have unspecified lower bound. For parameters with unspecified lower bounds, zero is used as a lower bound. Lower bounds only apply to fitted (i.e., non-fixed) parameters.
}
\item{upper}{
Either \code{NULL} or a named list containing optional upper bounds for various model parameters. The list can contain any or all of the elements \code{transition_matrix}, \code{birth_rates} and \code{death_rates}, structured similarly to \code{upper}. For example, \code{list(transition_matrix=2, birth_rates=c(10,NA,NA))} specifies that all transition rates between states must be 2 or less, that the birth rate \eqn{\lambda_1} must be 10 or less, and that all other model parameters have unspecified upper bound. For parameters with unspecified upper bounds, infinity is used as an upper bound. Upper bounds only apply to fitted (i.e., non-fixed) parameters.
}
\item{root_prior}{
Either a character or a numeric vector of size Nstates, specifying the prior probabilities of states for the root, i.e. the weights for obtaining a single model likelihood by averaging the root's state likelihoods. If a character, then it must be one of "flat", "empirical", "likelihoods", "max_likelihood" or "auto". "empirical" means the root's prior is set to the proportions of (estimated) extant species in each state (correcting for sampling fractions and reveal fractions, if applicable). "likelihoods" means that the computed state-likelihoods of the root are used, after normalizing to obtain a probability distribution; this is the approach used in the package \code{hisse::hisse} v1.8.9 under the option \code{root.p=NULL}, and the approach in the package \code{diversitree::find.mle} v0.9-10 under the option \code{root=ROOT.OBS}. If "max_likelihood", then the root's prior is set to a Dirac distribution, with full weight given to the maximum-likelihood state at the root (after applying the conditioning). If a numeric vector, \code{root_prior} specifies custom probabilities (weights) for each state. Note that if \code{root_conditioning} is "\code{madfitz}" or "\code{herr_als}" (see below), then the prior is set before the conditioning and not updated afterwards for consistency with other R packages.
}
\item{root_conditioning}{
Character, specifying an optional modification to be applied to the root's state likelihoods prior to averaging. Can be "none" (no modification), "madfitz", "herr_als", "crown" or "stem". "madfitz" and "herr_als" (after van Els, Etiene and Herrera-Alsina 2018) are the options implemented in the package \code{hisse} v1.8.9, conditioning the root's state-likelihoods based on the birth-rates and the computed extinction probability (after or before averaging, respectively). See van Els (2018) for a comparison between "madfitz" and "herr_als". The option "stem" conditions the state likelihoods on the probability that the stem lineage would survive until the present. The option "crown" conditions the state likelihoods on the probability that a split occurred at \code{oldest_age} and that the two child lineages survived until the present; this option is only recommended if \code{oldest_age} is equal to the root age.
}
\item{oldest_age}{
Strictly positive numeric, specifying the oldest age (time before present) to consider for fitting. If this is smaller than the tree's root age, then the tree is split into multiple subtrees at \code{oldest_age}, and each subtree is considered as an independent realization of the same diversification/evolution process whose parameters are to be estimated. The \code{root_conditioning} and \code{root_prior} are applied separately to each subtree, prior to calculating the joint (product) likelihood of all subtrees. This option can be used to restrict the fitting to a small (recent) time interval, during which the MuSSE/BiSSE assumptions (e.g., time-independent speciation/extinction/transition rates) are more likely to hold. If \code{oldest_age} is \code{NULL}, it is automatically set to the root age. In principle \code{oldest_age} may also be older than the root age.
}
\item{Ntrials}{
Non-negative integer, specifying the number of trials for fitting the model, using alternative (randomized) starting parameters at each trial. A larger \code{Ntrials} reduces the risk of landing on a local non-global optimum of the likelihood function, and thus increases the chances of finding the truly best fit. If 0, then no fitting is performed, and only the first-guess (i.e., provided or guessed start params) is evaluated and returned. Hence, setting \code{Ntrials=0} can be used to obtain a reasonable set of start parameters for subsequent fitting or for Markov Chain Monte Carlo.
}
\item{max_start_attempts}{
Integer, specifying the number of times to attempt finding a valid start point (per trial) before giving up on that trial. Randomly chosen extreme start parameters may occasionally result in Inf/undefined likelihoods, so this option allows the algorithm to keep looking for valid starting points. For complex models (e.g., with >5 states), setting this to 10, 100 or even 1000 may be beneficial.
}
\item{optim_algorithm}{
Character, specifying the optimization algorithm for fitting. Must be one of either "optim", "nlminb" or "subplex" (requires the \code{nloptr} package).
}
\item{optim_max_iterations}{
Integer, maximum number of iterations allowed for fitting. Only relevant for "optim" and "nlminb".
}
\item{optim_max_evaluations}{
Integer, maximum number of function evaluations allowed for fitting. Only relevant for "nlminb" and "subplex" (requires the \code{nloptr} package).
}
\item{optim_rel_tol}{
Numeric, relative tolerance for the fitted log-likelihood.
}
\item{check_input}{
Logical, specifying whether to check the validity of input variables. If you are certain that all input variables are valid, you can set this to \code{FALSE} to reduce computation.
}
\item{include_ancestral_likelihoods}{
Logical, specifying whether to include the state likelihoods for each node, in the returned variables. These are the ``D'' variables calculated as part of the likelihood based on the subtree descending from each node, and may be used for "local" ancestral state reconstructions.
}
\item{Nthreads}{
Integer, specifying the number of threads for running multiple fitting trials in parallel. Only relevant if \code{Ntrials>1}. Should generally not exceed the number of CPU cores on a machine. Must be a least 1.
}
\item{Nbootstraps}{
Integer, specifying an optional number of bootstrap samplings to perform, for estimating standard errors and confidence intervals of maximum-likelihood fitted parameters. If 0, no bootstrapping is performed. Typical values are 10-100. At each bootstrap sampling, a simulation of the fitted MuSSE/HiSSE model is performed, the parameters are estimated anew based on the simulation, and subsequently compared to the original fitted parameters. Each bootstrap sampling will thus use roughly as many computational resources as the original maximum-likelihood fit (e.g., same number of trials, same optimization parameters etc).
}
\item{Ntrials_per_bootstrap}{
Integer, specifying the number of fitting trials to perform for each bootstrap sampling. If \code{NULL}, this is set equal to \code{max(1,Ntrials)}. Decreasing \code{Ntrials_per_bootstrap} will reduce computation time, at the expense of potentially inflating the estimated confidence intervals; in some cases (e.g., for very large trees) this may be useful if fitting takes a long time and confidence intervals are very narrow anyway. Only relevant if \code{Nbootstraps>0}.
}
\item{max_condition_number}{
Positive unitless number, specifying the maximum permissible condition number for the "G" matrix computed for the log-likelihood. A higher condition number leads to faster computation (roughly on a log-scale) especially for large trees, at the potential expense of lower accuracy. Typical values are 1e2-1e5. See Louca (2019) for further details on the condition number of the G matrix.
}
\item{relative_ODE_step}{
Positive unitless number, specifying the default relative time step for the ordinary differential equation solvers.
}
\item{E_value_step}{
Positive unitless number, specifying the relative difference between subsequent recorded and interpolated E-values, in the ODE solver for the extinction probabilities E (Louca 2019). Typical values are 1e-2 to 1e-5. A smaller \code{E_value_step} increases interpolation accuracy, but also increases memory requirements and adds runtime (scaling with the tree's age span, not Ntips).
}
\item{D_temporal_resolution}{
Positive unitless number, specifying the relative resolution for interpolating G-map over time (Louca 2019). This is relative to the typical time scales at which G-map varies. For example, a resolution of 10 means that within a typical time scale there will be 10 interpolation points. Typical values are 1-1000. A greater resolution increases interpolation accuracy, but also increases memory requirements and adds runtime (scaling with the tree's age span, not Ntips).
}
\item{max_model_runtime}{
Numeric, optional maximum number of seconds for evaluating the likelihood of a model, prior to cancelling the calculation and returning Inf. This may be useful if extreme model parameters (e.g., reached transiently during fitting) require excessive calculation time. Parameters for which the calculation of the likelihood exceed this threshold, will be considered invalid and thus avoided during fitting. For example, for trees with 1000 tips a time limit of 10 seconds may be reasonable.
If 0, no time limit is imposed.
}
\item{verbose}{
Logical, specifying whether to print progress reports and warnings to the screen. In any case, fatal errors are always reported.
}
\item{diagnostics}{
Logical, specifying whether to print detailed information (such as model likelihoods) at every iteration of the fitting routine. For debugging purposes mainly.
}
\item{verbose_prefix}{
Character, specifying the line prefix for printing progress reports, warnings and errors to the screen.
}
}


\details{
HiSSE/SecSSE models include two discrete traits, one trait that defines the rate categories of diversification rates (as in BiSSE/MuSSE), and one trait that does not itself influence diversification but whose states (here called "proxy states") each represent one or more of the diversity-modulating states. HiSSE models (Beaulieu and O'meara, 2016) and SecSSE models (van Els et al., 2018) are closely related to BiSSE/MuSSE models, the main difference being the fact that the actual diversification-modulating states are not directly observed. In essence, a HiSSE/SecSSE model is a BiSSE/MuSSE model, where the final tip states are replaced by their proxy states, thus "masking" the underlying diversity-modulating trait. This function is able to fit HiSSE/SecSSE models with appropriate choice of the input variables \code{Nstates}, \code{NPstates} and \code{proxy_map}.
Note that the terminology and setup of HiSSE/SecSSE models followed here differs from their description in the original papers by Beaulieu and O'meara (2016) and van Els et al. (2018), in order to achieve what we think is a more intuitive unification of BiSSE/MuSSE/HiSSE/SecSSE. For ease of terminology, when considering a BiSSE/MuSSE model, here we use the terms "states" and "proxy-states" interchangeably, since under BiSSE/MuSSE the proxy trait can be considered identical to the diversification-modulating trait. A distinction between "states" and "proxy-states" is only relevant for HiSSE/SecSSE models.

As an example of a HiSSE model, \code{Nstates=4}, \code{NPstates=2} and \code{proxy_map=c(1,2,1,2)} specifies that states 1 and 3 are represented by proxy-state 1, and states 2 and 4 are represented by proxy-state 2. This is the original case described by Beaulieu and O'Meara (2016); in their terminology, there would be 2 "hidden"" states ("0" and "1") and 2 "observed" states ("A" and "B"), and the 4 diversification rate categories (\code{Nstates=4}) would be called "0A", "1A", "0B" and "1B". The somewhat different terminology used here allows for easier generalization to an arbitrary number of diversification-modulating states and an arbitrary number of proxy states. For example, if there are 6 diversification modulating states, represented by 3 proxy-states as 1->A, 2->A, 3->B, 4->C, 5->C, 6->C, then one would set \code{Nstates=6}, \code{NPstates=3} and \code{proxy_map=c(1,1,2,3,3,3)}.

The run time of this function scales asymptotically linearly with tree size (Ntips), although run times can vary substantially depending on model parameters. As a rule of thumb, the higher the birth/death/transition rates are compared to the tree's overall time span, the slower the calculation becomes.

The following arguments control the tradeoff between accuracy and computational efficiency:
\itemize{
\item\code{max_condition_number}: A smaller value means greater accuracy, at longer runtime and more memory.
\item\code{relative_ODE_step}: A smaller value means greater accuracy, at longer runtime.
\item\code{E_value_step}: A smaller value means greater accuracy, at longer runtime and more memory.
\item\code{D_temporal_resolution}: A greater value means greater accuracy, at longer runtime and more memory.
}
Typically, the default values for these arguments should be fine. For smaller trees, where cladogenic and sampling stochasticity is the main source of uncertainty, these parameters can probably be made less stringent (i.e., leading to lower accuracy and faster computation), but then again for small trees computational efficiency may not be an issue anyway.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the fitting was successful. If \code{FALSE}, an additional element \code{error} (of type character) is included containing an explanation of the error; in that case the value of any of the other elements is undetermined.
}
\item{Nstates}{
Integer, the number of states assumed for the model.
}
\item{NPstates}{
Integer, the number of proxy states assumed for the model. Note that in the case of a BiSSE/MuSSE model, this will be the same as \code{Nstates}.
}
\item{root_prior}{
Character, or numeric vector of length Nstates, specifying the root prior used.
}
\item{parameters}{
Named list containing the final maximum-likelihood fitted model parameters. If \code{Ntrials>1}, then this contains the fitted parameters yielding the highest likelihood. Will contain the following elements:
\itemize{
\item\code{transition_matrix}: 2D numeric matrix of size Nstates x Nstates, listing the fitted transition rates between states.
\item\code{birth_rates}: Numeric vector of length Nstates, listing the fitted state-dependent birth rates.
\item\code{death_rates}: Numeric vector of length Nstates, listing the fitted state-dependent death rates.
}
}
\item{start_parameters}{
Named list containing the default start parameter values for the fitting. Structured similarly to \code{parameters}. Note that if \code{Ntrials>1}, only the first trial will have used these start values, all other trials will have used randomized start values. Will be defined even if \code{Ntrials==0}, and can thus be used to obtain a reasonable guess for the start parameters without actually fitting the model.
}
\item{loglikelihood}{
Numeric, the maximized log-likelihood of the model, if fitting succeeded.
}
\item{AIC}{
Numeric, the Akaike Information Criterion for the fitted model, defined as \eqn{2k-2\log(L)}, where \eqn{k} is the number of fitted parameters and \eqn{L} is the maximized likelihood.
}
\item{Niterations}{
The number of iterations needed for the best fit. Only relevant if the optimization method was "optim" or "nlminb".
}
\item{Nevaluations}{
Integer, the number of function evaluations needed for the best fit. Only relevant if the optimization method was "nlminb" or "subplex".
}
\item{converged}{
Logical, indicating whether convergence was successful during fitting. If convergence was not achieved, and the fitting was stopped due to one of the stopping criteria \code{optim_max_iterations} or \code{optim_max_evaluations}, the final likelihood will still be returned, but the fitted parameters may not be reasonable.
}
\item{warnings}{
Character vector, listing any warnings encountered during evaluation of the likelihood function at the fitted parameter values. For example, this vector may contain warnings regarding the differential equation solvers or regarding the rank of the G-matrix (Louca, 2019).
}
\item{subroots}{
Integer vector, listing indices of tips/nodes in the tree that were considered as starting points of independent MuSSE processes. If \code{oldest_age} was equal to or greater than the root age, then \code{subroots} will simply list the tree's root.
}
\item{ML_subroot_states}{
Integer vector, with values between 1 and Nstates, giving the maximum-likelihood estimate of each subroot's state.
}
\item{ML_substem_states}{
Integer vector, with values between 1 and Nstates, giving the maximum-likelihood estimate of the state at each subroot's stem (i.e., exactly at \code{oldest_age}).
}
\item{trial_start_loglikelihoods}{
Numeric vector of length \code{Ntrials}, listing the initial loglikelihoods (i.e., at the starting parameter values) for each fitting trial.
}
\item{trial_loglikelihoods}{
Numeric vector of length \code{Ntrials}, listing the maximized loglikelihoods for each fitting trial. These may be used for diagnosing the robustness of maximum-likelihood estimates and the assessing the needed for increasing \code{Ntrials}.
}
\item{trial_Nstart_attempts}{
Integer vector of length \code{Ntrials}, listing the number of random start attempts for each trial (see option \code{max_start_attempts}).
}
\item{trial_Niterations}{
Integer vector of length \code{Ntrials}, listing the number of iterations of each trial. Depending on the fitting algorithm used (option \code{optim_algorithm}), these may be \code{NA} (not available).
}
\item{trial_Nevaluations}{
Integer vector of length \code{Ntrials}, listing the number of likelihood evaluations of each trial. Depending on the fitting algorithm used (option \code{optim_algorithm}), these may be \code{NA} (not available).
}
\item{standard_errors}{
Named list containing the elements "transition_matrix" (numeric matrix of size Nstates x Nstates), "birth_rates" (numeric vector of size Nstates) and "death_rates" (numeric vector of size Nstates), listing standard errors of all model parameters estimated using parametric bootstrapping. Only included if \code{Nbootstraps>0}. Note that the standard errors of non-fitted (i.e., fixed) parameters will be zero.
}
\item{CI50lower}{
Named list containing the elements "transition_matrix" (numeric matrix of size Nstates x Nstates), "birth_rates" (numeric vector of size Nstates) and "death_rates" (numeric vector of size Nstates), listing the lower end of the 50\% confidence interval (i.e. the 25\% quantile) for each model parameter, estimated using parametric bootstrapping. Only included if \code{Nbootstraps>0}.
}
\item{CI50upper}{
Similar to \code{CI50lower}, but listing the upper end of the 50\% confidence interval (i.e. the 75\% quantile) for each model parameter. For example, the confidence interval for he birth-rate \eqn{\lambda_1} will be between \code{CI50lower$birth_rates[1]} and \code{CI50upper$birth_rates[1]}. Only included if \code{Nbootstraps>0}.
}
\item{CI95lower}{
Similar to \code{CI50lower}, but listing the lower end of the 95\% confidence interval (i.e. the 2.5\% quantile) for each model parameter. Only included if \code{Nbootstraps>0}.
}
\item{CI95upper}{
Similar to \code{CI50upper}, but listing the upper end of the 95\% confidence interval (i.e. the 97.5\% quantile) for each model parameter. Only included if \code{Nbootstraps>0}.
}
\item{CI}{
2D numeric matrix, listing maximum-likelihood estimates, standard errors and confidence intervals for all model parameters (one row per parameter, one column for ML-estimates, one column for standard errors, two columns per confidence interval). Standard errors and confidence intervals are as estimated using parametric bootstrapping. This matrix contains the same information as \code{parameters}, \code{standard_errors}, \code{CI50lower}, \code{CI50upper}, \code{CI95lower} and \code{CI95upper}, but in a more compact format. Only included if \code{Nbootstraps>0}.
}
\item{ancestral_likelihoods}{
2D matrix of size Nnodes x Nstates, listing the computed state-likelihoods for each node in the tree. These may be used for "local" ancestral state reconstructions, based on the information contained in the subtree descending from each node. Note that for each node the ancestral likelihoods have been normalized for numerical reasons, however they should not be interpreted as actual probabilities. For each node n and state s, \code{ancestral_likelihoods[n,s]} is proportional to the likelihood of observing the descending subtree and associated tip proxy states, if node n was at state s.
Only included if \code{include_ancestral_likelihoods==TRUE}.
}
}

\author{Stilianos Louca}

\references{
W. P. Maddison, P. E. Midford, S. P. Otto (2007). Estimating a binary character's effect on speciation and extinction. Systematic Biology. 56:701-710.

R. G. FitzJohn, W. P. Maddison, S. P. Otto (2009). Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology. 58:595-611

R. G. FitzJohn (2012). Diversitree: comparative phylogenetic analyses of diversification in R. Methods in Ecology and Evolution. 3:1084-1092

J. M. Beaulieu and B. C. O'Meara (2016). Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic Biology. 65:583-601.

D. Kuehnert, T. Stadler, T. G. Vaughan, A. J. Drummond (2016). Phylodynamics with migration: A computational framework to quantify population structure from genomic data. Molecular Biology and Evolution. 33:2102-2116.

P. van Els, R. S. Etiene, L. Herrera-Alsina (2018). Detecting the dependence of diversification on multiple traits from phylogenetic trees and trait data. Systematic Biology. syy057.

S. Louca and M. W. Pennell (2020). A general and efficient algorithm for the likelihood of diversification and discrete-trait evolutionary models. Systematic Biology. 69:545-556.
}

\seealso{
\code{\link{simulate_dsse}},
\code{\link{asr_mk_model}},
\code{\link{fit_tree_model}}
}

\examples{
# EXAMPLE 1: BiSSE model
# - - - - - - - - - - - - - -
# Choose random BiSSE model parameters
Nstates = 2
Q = get_random_mk_transition_matrix(Nstates, rate_model="ARD", max_rate=0.1)
parameters = list(birth_rates       = runif(Nstates,5,10),
                  death_rates       = runif(Nstates,0,5),
                  transition_matrix = Q)
rarefaction = 0.5 # randomly omit half of the tips

# Simulate a tree under the BiSSE model
simulation = simulate_musse(Nstates, 
                            parameters         = parameters, 
                            max_tips           = 1000,
                            sampling_fractions = rarefaction)
tree       = simulation$tree
tip_states = simulation$tip_states

\dontrun{
# fit BiSSE model to tree & tip data
fit = fit_musse(tree,
                Nstates            = Nstates,
                tip_pstates        = tip_states,
                sampling_fractions = rarefaction)
if(!fit$success){
  cat(sprintf("ERROR: Fitting failed"))
}else{
  # compare fitted birth rates to true values
  errors = (fit$parameters$birth_rates - parameters$birth_rates)
  relative_errors = errors/parameters$birth_rates
  cat(sprintf("BiSSE relative birth-rate errors:\n"))
  print(relative_errors)
}
}


# EXAMPLE 2: HiSSE model, with bootstrapping
# - - - - - - - - - - - - - -
# Choose random HiSSE model parameters
Nstates  = 4
NPstates = 2
Q = get_random_mk_transition_matrix(Nstates, rate_model="ARD", max_rate=0.1)
rarefaction = 0.5 # randomly omit half of the tips
parameters = list(birth_rates       = runif(Nstates,5,10),
                  death_rates       = runif(Nstates,0,5),
                  transition_matrix = Q)
                  
# reveal the state of 30\% & 60\% of tips (in state 1 & 2, respectively)
reveal_fractions = c(0.3,0.6)

# use proxy map corresponding to Beaulieu and O'Meara (2016)
proxy_map = c(1,2,1,2)

# Simulate a tree under the HiSSE model
simulation = simulate_musse(Nstates, 
                            NPstates            = NPstates,
                            proxy_map           = proxy_map,
                            parameters          = parameters, 
                            max_tips            = 1000,
                            sampling_fractions  = rarefaction,
                            reveal_fractions    = reveal_fractions)
tree       = simulation$tree
tip_states = simulation$tip_proxy_states

\dontrun{
# fit HiSSE model to tree & tip data
# run multiple trials to ensure global optimum
# also estimate confidence intervals via bootstrapping
fit = fit_musse(tree,
                Nstates            = Nstates,
                NPstates           = NPstates,
                proxy_map          = proxy_map,
                tip_pstates        = tip_states,
                sampling_fractions = rarefaction,
                reveal_fractions   = reveal_fractions,
                Ntrials            = 5,
                Nbootstraps        = 10,
                max_model_runtime  = 0.1)
if(!fit$success){
  cat(sprintf("ERROR: Fitting failed"))
}else{
  # compare fitted birth rates to true values
  errors = (fit$parameters$birth_rates - parameters$birth_rates)
  relative_errors = errors/parameters$birth_rates
  cat(sprintf("HiSSE relative birth-rate errors:\n"))
  print(relative_errors)
  
  # print 95\%-confidence interval for first birth rate
  cat(sprintf("CI95 for lambda1: \%g-\%g",
              fit$CI95lower$birth_rates[1],
              fit$CI95upper$birth_rates[1]))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{tree model}
\keyword{MuSSE}
\keyword{BiSSE}
\keyword{HiSSE}
