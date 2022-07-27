
#######################################################
# Compare BioGeoBEARS and diversitree-ClaSSE calculations
#
# This script will illustrate, for any input tree
# and geography file, and BioGeoBEARS parameters, the 
# equivalence of comparing BioGeoBEARS models (such as 
# DEC and DEC+J, but this applies to any 2 models) and 
# comparing 2 equivalent ClaSSE models, in the case where:
#
# * the claSSE lambdas = BGB_cladogenesis_probs * birthRate
# * the birthRate = the Maximum Likelihood estimate under Yule
# * the sampling is assumed to be 100%
#
# The key insight is that BioGeoBEARS calculates the 
# likelihood of the geography data (just a complex 
# character dataset). Adding a Yule process likelihood, 
# i.e. the probability density of a tree under a pure-birth
# process with an extinction rate of 0 and assuming 
# complete sampling, creates a special case of the 
# ClaSSE model.
#
# Therefore, comparing the likelihood difference between
# two BioGeoBEARS models is exactly equivalent to comparing
# two equivalent ClaSSE models. The likelihood difference
# will be identical, because the Yule-process likelihood
# is a constant across the different BioGeoBEARS of the 
# geography data.
#
# Note: This is only set up for non-time-stratified analyses.
# Also, the Yule-process assumption fails for a tree
# that includes fossils (non-contemporaneous tips).
# 
#######################################################

library(ape)
library(diversitree)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

source("/GitHub/BioGeoJulia.jl/Rsrc/ClaSSE_functions_v3.R")  # utility functions from diversitree
source("/GitHub/BioGeoJulia.jl/Rsrc/ClaSSE_mods_v2.R")       # helper functions in plain-R


wd = "/GitHub/BioGeoJulia.jl/Rsrc/"
setwd(wd)

# Load simple example tree (newick format, must be ultrametric, i.e. 
# all the tips come to the present)
trfn = "Psychotria_tree.newick"
tr = read.tree(trfn)

# Load geography data
geogfn = "Psychotria_geog.data"


#######################################################
# BioGeoBEARS Q and C matrices
#######################################################
library(BioGeoBEARS)
max_range_size = 4
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
include_null_range = BioGeoBEARS_run_object$include_null_range

# Run the Maximum Likelihood optimization
res = bears_optim_run(BioGeoBEARS_run_object)
res$total_loglikelihood
# -34.54313 for Psychotria, under DEC

# Extract the ML model's transition matrices
mats = get_Qmat_COOmat_from_res(res, numstates=ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node), include_null_range=res$inputs$include_null_range, max_range_size=res$inputs$max_range_size, timeperiod_i=1)
numstates = length(mats$states_list)

# Explanation:
# mats$Qmat = anagenetic transition matrix
# mats$COO_weights_columnar = Cmat, the cladogenetic transition matrix.
#      This is in a "table" format (COO-like), and includes only
#      cladogenetic range-inheritance scenarios allowed under the 
#      particular BioGeoBEARS model.
#      COO_weights_columnar codes 4 columns, giving 0-based state-indices
#        * for i, j, and k (ancestor, left descendant, right descendant)
#        * per-event weights for each cladogenetic range-inheritance scenario
#      The null_range is left out the Cmat, so to convert to 
#        standard 1-based state counts, add (1+sum(include_null_range)).
#      To convert the per-event weights to per-event probabilities:
#      per-event_prob_of_scenario = per-weight_prob_of_scenario / sum_of_all_weights_for_that_ancestor_i
# mats$Rsp_rowsums = the sum of the weights for each ancestor state i (again excluding null_range)
#
# To convert per-event probabilities to diversitree's lambdas, multiply by birthRate
#
# I.e. if
# A = state 2
# B = state 3
# AB = state 4
#
# And if 
# AB -> A, B  has probability 1/6 (vicariance)
# AB -> B, A  has probability 1/6 (vicariance)
# AB -> A, AB has probability 1/6 (subset sympatry)
# AB -> AB, A has probability 1/6 (subset sympatry)
# AB -> B, AB has probability 1/6 (subset sympatry)
# AB -> AB, B has probability 1/6 (subset sympatry)
# 
# And the ML birthRate under the Yule model is 0.9, then the lambda rates of
# each individual type of speciation event are:
#
# When birthRate = 0.9,
#
# AB -> A, B  has rate 1/6 * 0.9 = lambda423 = 0.15
# AB -> B, A  has rate 1/6 * 0.9 = lambda432 = 0.15
# AB -> A, AB has rate 1/6 * 0.9 = lambda424 = 0.15
# AB -> AB, A has rate 1/6 * 0.9 = lambda442 = 0.15
# AB -> B, AB has rate 1/6 * 0.9 = lambda434 = 0.15
# AB -> AB, B has rate 1/6 * 0.9 = lambda443 = 0.15
#
# Total speciation rate when in state AB = 0.15*6 = 0.9
# 
# In diversitree, the lambdas of equivalent left/right scenarios are combined, so
# diversitree will use:
#
# lambda423 = 0.3
# lambda424 = 0.3
# lambda434 = 0.3
# 


# Extract the parameters
birthRate = 0.3288164        # 0.3288164 for Psychotria tree
birthRate = yule(tr)$lambda  # The ML lambda from Yule. Equals (#speciations-1)/tree_length
birthRate
(tr$Nnode-1)/sum(tr$edge.length)

deathRate = 0.0     # Yule process means 0.0 extinction rate
d_val = 0.03505038	# ML estimate for Psychotria under DEC model
e_val = 0.02832370	# ML estimate for Psychotria under DEC model
j_val = 0.0         # Under DEC, j=0.0
d_val = res$output@params_table["d","est"] # Extract from ML result
e_val = res$output@params_table["e","est"] # Extract from ML result
j_val = 0.0



#######################################################
# Set up an equivalent ClaSSE model in diversitree
#######################################################

# This .R file contains a bunch of extract functions that diversitree
# uses behind the scenes. They are very hard to access normally, 
# because diversitree does tons of "functions writing other functions".
source('/GitHub/BioGeoJulia.jl/Rsrc/ClaSSE_mods_v2.R', chdir = TRUE)

# Convert the BioGeoBEARS cladogenesis matrix into a data.frame
# This handily shows the 
# * per-event weights (column "wt") and 
# * per-event probabilities (column "prob")
include_null_range = res$inputs$include_null_range
Carray_df = get_Cevent_probs_df_from_mats(mats, include_null_range=include_null_range)
head(Carray_df)
tail(Carray_df)

# Set up a ClaSSE model from diversitree

# Get the tip statenums
numtips = length(tr$tip.label)
numstates = ncol(res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
tip_statenums = rep(0, times=numtips)
for (i in 1:numtips)
	{ # Find the "1" (the observed state for each tip)
	TF = res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[i,] == 1
	tip_statenums[i] = (1:numstates)[TF]
	}
tip_statenums
names(tip_statenums) = tr$tip.label
states = tip_statenums

# Set the sampling rate to 1 for each state
sampling.f = rep(1, times=numstates)
k = numstates

# Create the ClaSSE likelihood function for k states
# (strict=FALSE means that some states in the state space can be absent from the tips)
classe_Kstates = make.classe(tree=tr, states=states, k=k, sampling.f=sampling.f, strict=FALSE)

# The names of the ClaSSE parameters:
# Note that diverstree creates ALL the possible parameters, which gets
# ridiculous quickly, e.g. 
# 4 areas = 16 geographic range states = 2432 parameters in param_names
param_names = argnames(classe_Kstates)
length(param_names)
length(param_names[grepl(pattern="lambda", x=param_names)]) # 2176 speciation rates
length(param_names[grepl(pattern="mu", x=param_names)])     #   16 extinction rates
length(param_names[grepl(pattern="q", x=param_names)])      #  240 Q transition rates

# Most parameters will be zero
classe_params = rep(0, times=length(param_names))
names(classe_params) = param_names
head(classe_params)
tail(classe_params)

# Make a data.frame to match up with the BioGeoBEARS Carray_df
lambda_ijk_df = classe_lambdas_to_df(classe_params, k=numstates)
head(lambda_ijk_df)

# Fill in the params from the BioGeoBEARS "res" results
classe_params = BGBres_into_classe_params(res, classe_params, birthRate=birthRate)
classe_params[153]
classe_params["lambda020202"]
classe_params["lambda060203"]
classe_params["q0206"]
classe_params["q1516"]
classe_params[classe_params != 0.0]

# Set up various assumptions about the root state probabilities
# All probabilities equal, except null range has prob=0
root_probs_equal = rep(1, times=numstates)
root_probs_equal[sum(include_null_range)] = 0
root_probs_equal = root_probs_equal / sum(root_probs_equal)

# Highly biased towards the last state
root_probs_biased = rep(0.01, times=numstates)
root_probs_biased[sum(include_null_range)] = 0
root_probs_biased[length(root_probs_biased)] = 0.01 * (numstates-include_null_range)
root_probs_biased = root_probs_biased / sum(root_probs_biased)

# All states, except null range, get "probability" 1
# (i.e., ignore root state frequencies, like DEC-type models)
root_probs_single = rep(1, times=numstates)
root_probs_single[sum(include_null_range)] = 0

# Do the ClaSSE calculation, under these different assumptions
res1 = classe_Kstates(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
res2 = classe_Kstates(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
root_probs = root_probs_equal
res3 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = root_probs_biased
res4 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = root_probs_single
res5 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
res6 = classe_Kstates(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)

res1t = classe_Kstates(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
res2t = classe_Kstates(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
root_probs = root_probs_equal
res3t = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = root_probs_biased
res4t = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = root_probs_single
res5t = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
res6t = classe_Kstates(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)

# get_classe_LnLs returns the total log-likelihood, and 
# the total of the branch likelihoods
LnLs1 = get_classe_LnLs(res1)
LnLs2 = get_classe_LnLs(res2)
LnLs3 = get_classe_LnLs(res3)
LnLs4 = get_classe_LnLs(res4)
LnLs5 = get_classe_LnLs(res5)
LnLs6 = get_classe_LnLs(res6)
LnLs1t = get_classe_LnLs(res1t)
LnLs2t = get_classe_LnLs(res2t)
LnLs3t = get_classe_LnLs(res3t)
LnLs4t = get_classe_LnLs(res4t)
LnLs5t = get_classe_LnLs(res5t)
LnLs6t = get_classe_LnLs(res6t)

LnLst = as.data.frame(rbind(LnLs1, LnLs2, LnLs3, LnLs4, LnLs5, LnLs6, LnLs1t, LnLs2t, LnLs3t, LnLs4t, LnLs5t, LnLs6t), stringsAsFactors=FALSE)
names(LnLst) = c("ttl_LnL", "branch_LnL")
ObsDiff = (LnLst$ttl_LnL - LnLst$branch_LnL)
exp_ObsDiff = exp((LnLst$ttl_LnL - LnLst$branch_LnL))
LnLdiff = round((LnLst$ttl_LnL - LnLst$branch_LnL - log(birthRate)), digits=4)
exp_LnLdiff = exp((LnLst$ttl_LnL - LnLst$branch_LnL - log(birthRate)))
LnLst2 = cbind(LnLst, ObsDiff, LnLdiff, exp_ObsDiff, exp_LnLdiff)

# Put the total and branch lnLs in a table; other columns were experimenting
# with various assumptions about constants (ignore except in ultra-simple cases)
all_lnLs = cft(LnLst2, numdigits_inbetween_have_fixed_digits=8)
all_lnLs$ttl_LnL = as.numeric(all_lnLs$ttl_LnL)
all_lnLs$branch_LnL = as.numeric(all_lnLs$branch_LnL)
all_lnLs$ObsDiff = as.numeric(all_lnLs$ObsDiff)
all_lnLs$LnLdiff = as.numeric(all_lnLs$LnLdiff)
all_lnLs$exp_ObsDiff = as.numeric(all_lnLs$exp_ObsDiff)
all_lnLs$exp_LnLdiff = as.numeric(all_lnLs$exp_LnLdiff)
all_lnLs

DEC_all_lnLs = all_lnLs
dput(DEC_all_lnLs)


# The diversitree outputs, with intermediates stored, are useful for 
# branch-by-branch comparison. However, they have to be transposed to 
# compare to BioGeoBEARS (so that rows = nodes).
# The columns are the Es for each state, then the Ds for each state
init = t(attr(res2, "intermediates")$init)
init

base = t(attr(res2, "intermediates")$base)
base

# lq is the Ds log-likelihood summed at each branch bottom, and extracted
Ds_cols = (numstates+1):(2*numstates)
lq = t(attr(res2, "intermediates")$lq)
lq
rowSums(base[,Ds_cols])
base[,Ds_cols] * exp(c(lq))
rowSums(base[,Ds_cols] * exp(c(lq)))
log(rowSums(base[,Ds_cols] * exp(c(lq))))

# Store the likelihoods at branch-bottoms for comparison
base_likes = apply(X=base[,Ds_cols], MARGIN=2, FUN="*", exp(lq))
base_normlikes = base_likes / rowSums(base_likes)

# Diversitree normalized likelihoods at branch bottoms match BioGeoBEARS
round(base_normlikes - res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS, 5)

# These match the lqs, but this is because base_likes = base_normlikes * exp(lq)
tmp_rowSums = (base_likes / res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS)[,numstates]
tmp_rowSums
log(tmp_rowSums)
sum(log(tmp_rowSums), na.rm=TRUE) # matches sum(lq)
sum(lq)
log(tmp_rowSums) - attr(res4,"intermediates")$lq


# Diversitree birthdeath calculation
lik.bd <- make.bd(tree=tr, sampling.f=NULL, unresolved=NULL, times=NULL, control=list(method="ode"))
diversitree_bd = lik.bd(pars=c(birthRate=birthRate, deathRate=deathRate), intermediates=TRUE)
c(diversitree_bd)    # log-likelihood = 0.4870967
yule(tr)$loglik      # log-likelihood = 0.4870968

# Likelihood equation in the birthdeath function
# (derived by pulling apart the birthdeath() function from ape)
# This version stores all of the pieces, for comparison
# bd_ape$lnL = 0.4870968
bd_ape = bd_liks(tr, birthRate=birthRate, deathRate=deathRate)
bd_ape

# Note how this equals -(tr$Nnode-1)
bd_ape$lnl_Births_above_root + bd_ape$lnl_branching_times
-(tr$Nnode-1)

# The diversitree birth-death function also stores 
# branch-bottom likelihoods in "lq"
bd_lq = attr(diversitree_bd,"intermediates")$lq
sum(bd_lq)
sum(bd_lq) - -(tr$Nnode-1)
bd_ape$lnl_numBirths

# Compare bd_lq and -birthRate * trtable$edge.length
trtable = prt(tr, printflag=FALSE) # prints the tree to node-order table
bd_lq
-birthRate * trtable$edge.length

# Differences
round(bd_lq - (-birthRate * trtable$edge.length), digits=4)

# What is that -1.1123?
log(birthRate) # i.e., a log(birthRate) for every internal node


#######################################################
# Matching diversitree and BioGeoBEARS
#######################################################
# BioGeoBEARS stores the likelihoods at each node, 
# including the root node.
#
# But diversitree stores the likelihoods at branch bottoms,
# and treats the root node differently, depending on 
# various user options.
#
# So, let's start by summing the BioGeoBEARS likelihoods
# but exclude the root node.
root_nodenum = length(tr$tip.label) + 1
sumBGBlike_not_root = sum(log(res$computed_likelihoods_at_each_node[-root_nodenum]))
sumBGBlike_not_root

# Let's take the sum of the branch-bottom likelihoods from the birth-death
# process
sum(bd_lq)
bd_ape$lnl_numBirths + bd_ape$lnl_Births_above_root + bd_ape$lnl_branching_times
bd_ape$lnl_numBirths + -(tr$Nnode-1)
sum(-birthRate * trtable$edge.length, na.rm=TRUE) + (tr$Nnode-1)*log(birthRate) 
sum_branchBot_likes = sum(-birthRate * trtable$edge.length, na.rm=TRUE) + (tr$Nnode-1)*log(birthRate) 

# Add the lnL of root speciation event, -1 for extra node
all_lnLs
sumBGBlike_not_root + bd_ape$lnl_numBirths + bd_ape$lnl_Births_above_root + bd_ape$lnl_branching_times + (1-log(1/birthRate))
sumBGBlike_not_root + bd_ape$lnL - bd_ape$lnl_topology + (1-log(1/birthRate))
sumBGBlike_not_root + sum_branchBot_likes + (1-log(1/birthRate))
# For DEC, with anagenetic change, matches the branch_lnLs!


# We can also add the root state likelihoods, if desired
BGBlnL_at_root = log(res$computed_likelihoods_at_each_node[root_nodenum]) - 1
d_root_orig_BGB = res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[root_nodenum,] * exp(BGBlnL_at_root)
d_root_orig_BGB
sum(d_root_orig_BGB)

vals = t(attr(res1, "intermediates")$vals)	# Es and Ds at the root
E_indices = 1:numstates
d_root_orig_diversitree = vals[-E_indices]
d_root_orig_diversitree
sum(d_root_orig_diversitree)

# Match BioGeoBEARS to diverstree res1
root.p = d_root_orig_BGB/sum(d_root_orig_BGB)
rootlikes = log(sum(root.p * d_root_orig_BGB))
rootlikes

sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
c(res1)



# Match BioGeoBEARS to diverstree res2 (ROOT.FLAT state frequencies)
root.p = rep(1/numstates, times=numstates)
rootlikes = log(sum(root.p * d_root_orig_BGB))
rootlikes

sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
c(res2)

# Match BioGeoBEARS to diverstree res3 (all equal, except null range -- ROOT.GIVEN)
# Set up various assumptions about the root state probabilities
# All probabilities equal, except null range has prob=0
root_probs_equal = rep(1, times=numstates)
root_probs_equal[sum(include_null_range)] = 0
root_probs_equal = root_probs_equal / sum(root_probs_equal)
root.p = root_probs_equal
rootlikes = log(sum(root.p * d_root_orig_BGB))
rootlikes

sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
c(res3)



# Match BioGeoBEARS to diverstree res5 (all 1s, except null range -- ROOT.GIVEN)
# All states, except null range, get "probability" 1
# (i.e., ignore root state frequencies, like DEC-type models)
root_probs_single = rep(1, times=numstates)
root_probs_single[sum(include_null_range)] = 0
root.p = root_probs_single
rootlikes = log(sum(root.p * d_root_orig_BGB))
rootlikes

sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
c(res5)



#######################################################
# Matching diversitree with BioGeoBEARS+Yule+SFs
# (simpler!)
#######################################################


bd_ape = bd_liks(tr, birthRate=birthRate, deathRate=deathRate)

bd_ape$lnl_topology
bd_ape$lnl_numBirths
bd_ape$lnl_Births_above_root
bd_ape$lnl_numtips_wOneMinusDeathRate
bd_ape$lnl_branching_times
bd_ape$lnL

bgb1 = sum(log(res$computed_likelihoods_at_each_node[-root_nodenum]))
bgb2 = sum(log(res$computed_likelihoods_at_each_node))
bgb_root_lnL = sum(log(res$computed_likelihoods_at_each_node[root_nodenum]))
equal_root_prob = log(1/numstates)
equal_root_prob2 = log(1/(numstates-include_null_range)) 


all_lnLs
bgb1 + bd_ape$lnL - bd_ape$lnl_topology

# Matches classe branch_lnL
bgb1 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate)
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL

# res1 match
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL +  log(sum(d_root_orig_BGB*d_root_orig_BGB/sum(d_root_orig_BGB)))

# res2 match
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob
bgb2 + bd_ape$lnL - bd_ape$lnl_topology - log(1/birthRate) + equal_root_prob

# res3 match
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2
bgb2 + bd_ape$lnL - bd_ape$lnl_topology - log(1/birthRate) + equal_root_prob2


# res5 match
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) 
bgb2 + bd_ape$lnL - bd_ape$lnl_topology - log(1/birthRate) 


# res1t match
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL +  log(sum(d_root_orig_BGB*d_root_orig_BGB/sum(d_root_orig_BGB))) + log(1/birthRate)

# res2t match
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/birthRate)
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2



# res3t match
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/birthRate)
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2


# res5t match
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/(birthRate))
bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 
bgb1 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 + bgb_root_lnL
(bgb1 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate)) + equal_root_prob2 + bgb_root_lnL - (1-log(1/birthRate))



#######################################################
# NOTES
#######################################################

# To get the "cache" object from diversitree:

# Get the constant part of bd_lq
cache <- diversitree:::make.cache.bd(tree=tr, sampling.f=NULL, unresolved=NULL, times=NULL, control=list(method="ode"))
cache$const             # 36.39545
bd_ape$lnl_topology     # 36.39545
lfactorial(tr$Nnode)    # 36.39545







#######################################################
# Compare to Julia
#######################################################
base = t(attr(res2, "intermediates")$base)
base
lq = t(attr(res2, "intermediates")$lq)
lq
all_lnLs
R_result_branch_lnL = sum(lq)
R_result_total_LnLs1 = c(res1)
R_result_total_LnLs1t = c(res1t)

# Does the total of branch likelihoods (lq) + node likelihoods match R?
computed_likelihoods_at_each_node_x_lambda = rep(0.0, times=tr$Nnode + length(tr$tip.label))

computed_likelihoods_at_each_node_just_before_speciation = get_sum_log_computed_likes_at_each_node(tr, base, lq, classe_params)
computed_likelihoods_at_each_node_just_before_speciation
rowSums(computed_likelihoods_at_each_node_just_before_speciation)
log(rowSums(computed_likelihoods_at_each_node_just_before_speciation))
TF = is.finite(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)))
sum(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)[TF]))
R_result_sum_log_computed_likelihoods_at_each_node = sum(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)[TF]))
# > rowSums(computed_likelihoods_at_each_node_just_before_speciation)
#  [1] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#  [8] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
# [15] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.019560684 0.031880406
# [22] 0.001104266 0.026160645 0.047882797 0.052434844 0.052130315 0.287861948 0.277916166
# [29] 0.256410560 0.281575666 0.018628353 0.017883012 0.024094055 0.049780977 0.236860054
# [36] 0.049857979 0.240287100
# > log(rowSums(computed_likelihoods_at_each_node_just_before_speciation))
#  [1]      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf
# [10]      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf      -Inf
# [19]      -Inf -3.934234 -3.445764 -6.808575 -3.643499 -3.038999 -2.948184 -2.954009 -1.245274
# [28] -1.280436 -1.360975 -1.267354 -3.983070 -4.023904 -3.725790 -3.000122 -1.440286 -2.998577
# [37] -1.425921
# > TF = is.finite(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)))
# > sum(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)[TF]))
# [1] -52.52497
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = R_result_sum_log_computed_likelihoods_at_each_node + sum(lq)
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
# -120.1545

# res5 uses root.p of e.g.:
# 0,1,1,1,1,1,1,1,1,1....
c(res5)

# DEC
DEC_lnL = -34.54313
DEC_R_result_branch_lnL = -67.6295
DEC_R_result_total_LnLs1 = -72.60212
DEC_R_result_total_LnLs1t = -71.48986
DEC_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -120.1545
DEC_R_result_total_LnLs5 = -71.56373

# DEC+J
DECj_lnL = -20.94759
DECj_R_result_branch_lnL = -55.37332
DECj_R_result_total_LnLs1 = -58.83758
DECj_R_result_total_LnLs1t = -57.72533
DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -96.34151
DECj_R_result_total_LnLs5 = -57.96819

# Differences between DEC and DEC+J
DEC_lnL - DECj_lnL
DEC_R_result_branch_lnL - DECj_R_result_branch_lnL
DEC_R_result_total_LnLs1 - DECj_R_result_total_LnLs1
DEC_R_result_total_LnLs1t - DECj_R_result_total_LnLs1t
DEC_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda - DECj_R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
DEC_R_result_total_LnLs5 - DECj_R_result_total_LnLs5 # EXACT MATCH

# EXACT MATCH
DEC_lnL - DECj_lnL
-13.59554
DEC_R_result_total_LnLs5 - DECj_R_result_total_LnLs5 # EXACT MATCH
-13.59554



DEC_all_lnLs = structure(list(ttl_LnL = c(-72.602116, -74.33632, -74.271782, 
-74.605042, -71.563731, -74.263777, -71.48986, -73.159526, -73.159526, 
-73.492786, -73.159526, -73.079568), branch_LnL = c(-67.629498, 
-67.629498, -67.629498, -67.629498, -67.629498, -67.629498, -67.629498, 
-67.629498, -67.629498, -67.629498, -67.629498, -67.629498), 
    ObsDiff = c(-4.972618, -6.706822, -6.642284, -6.975544, -3.934234, 
    -6.634279, -3.860362, -5.530028, -5.530028, -5.863288, -5.530028, 
    -5.45007), LnLdiff = c(-3.8604, -5.5946, -5.53, -5.8633, 
    -2.822, -5.522, -2.7481, -4.4178, -4.4178, -4.751, -4.4178, 
    -4.3378), exp_ObsDiff = c(0.00692499, 0.00122254, 0.00130405, 
    0.00093446, 0.0195607, 0.00131453, 0.0210604, 0.00396588, 
    0.00396588, 0.00284188, 0.00396588, 0.004296), exp_LnLdiff = c(0.0210604, 
    0.00371801, 0.00396588, 0.00284188, 0.0594882, 0.00399775, 
    0.064049, 0.0120611, 0.0120611, 0.00864277, 0.0120611, 0.0130651
    )), .Names = c("ttl_LnL", "branch_LnL", "ObsDiff", "LnLdiff", 
"exp_ObsDiff", "exp_LnLdiff"), row.names = c("LnLs1", "LnLs2", 
"LnLs3", "LnLs4", "LnLs5", "LnLs6", "LnLs1t", "LnLs2t", "LnLs3t", 
"LnLs4t", "LnLs5t", "LnLs6t"), class = "data.frame")



DECj_all_lnLs = 
structure(list(ttl_LnL = c(-58.837585, -60.740782, -60.676244, 
-61.330746, -57.968194, -59.96035, -57.725329, -59.563988, -59.563988, 
-60.218491, -59.563988, -58.848094), branch_LnL = c(-55.37332, 
-55.37332, -55.37332, -55.37332, -55.37332, -55.37332, -55.37332, 
-55.37332, -55.37332, -55.37332, -55.37332, -55.37332), ObsDiff = c(-3.464265, 
-5.367462, -5.302924, -5.957426, -2.594874, -4.58703, -2.352009, 
-4.190668, -4.190668, -4.845171, -4.190668, -3.474774), LnLdiff = c(-2.352, 
-4.2552, -4.1907, -4.8452, -1.4826, -3.4748, -1.2398, -3.0784, 
-3.0784, -3.7329, -3.0784, -2.3625), exp_ObsDiff = c(0.031296, 
0.00466596, 0.00497702, 0.00258656, 0.0746553, 0.0101831, 0.0951777, 
0.0151362, 0.0151362, 0.00786628, 0.0151362, 0.0309688), exp_LnLdiff = c(0.0951777, 
0.0141902, 0.0151362, 0.00786628, 0.227043, 0.0309688, 0.289456, 
0.0460323, 0.0460323, 0.023923, 0.0460323, 0.0941827)), .Names = c("ttl_LnL", 
"branch_LnL", "ObsDiff", "LnLdiff", "exp_ObsDiff", "exp_LnLdiff"
), row.names = c("LnLs1", "LnLs2", "LnLs3", "LnLs4", "LnLs5", 
"LnLs6", "LnLs1t", "LnLs2t", "LnLs3t", "LnLs4t", "LnLs5t", "LnLs6t"
), class = "data.frame")


# BioGeoBEARS differences
DEC_lnL - DECj_lnL

# Diversitree differences
DEC_R_result_total_LnLs5 - DECj_R_result_total_LnLs5
DEC_all_lnLs$branch_LnL - DECj_all_lnLs$branch_LnL
DEC_all_lnLs$ttl_LnL - DECj_all_lnLs$ttl_LnL



# BioGeoBEARS
# DEC_lnL - DECj_lnL
# -13.59554

DEC_R_result_total_LnLs5 - DECj_R_result_total_LnLs5 # EXACT MATCH
# -13.59554

# Diversitree branch likelihoods
DEC_all_lnLs$branch_LnL - DECj_all_lnLs$branch_LnL
#  [1] -12.25618 -12.25618 -12.25618 -12.25618 -12.25618 -12.25618
#  [7] -12.25618 -12.25618 -12.25618 -12.25618 -12.25618 -12.25618

# Diversitree total likelihood differences
DEC_all_lnLs$ttl_LnL - DECj_all_lnLs$ttl_LnL
#  [1] -13.76453 -13.59554 -13.59554 -13.27430 -13.59554 -14.30343
#  [7] -13.76453 -13.59554 -13.59554 -13.27429 -13.59554 -14.23147

# Exact matches in lnL difference occur with:
# res2: If root=ROOT.FLAT, root.p=NULL, condition.surv=FALSE
root.p = rep(1/numstates, times=numstates)

# res3: If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root_probs = root_probs_equal

# res5: If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root_probs = root_probs_single   # 0,1,1,1,1,1...







# Key parts of the root calculation
lq = t(attr(res2, "intermediates")$lq)			# Branch likelihoods
vals = t(attr(res2, "intermediates")$vals)	# Es and Ds at the root
numstates = length(vals) / 2
E_indices = 1:numstates
d_root_orig = vals[-E_indices]							# Original D likelihoods at root

# res1: If root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
root.p = d_root_orig/sum(d_root_orig)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# res2: If root=ROOT.FLAT, root.p=NULL, condition.surv=FALSE
root.p = rep(1/numstates, times=numstates)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# res3: If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root_probs = root_probs_equal
root.p = root_probs
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# res4: If root=ROOT.GIVEN, root.p=c(0.75,0.25), condition.surv=FALSE
root_probs = root_probs_biased
root.p = root_probs
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# res5: If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root_probs = root_probs_single
root.p = root_probs
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# res6: Equilibrium root frequencies
# If root=ROOT.EQUI, condition.surv=FALSE

# Project the ClaSSE model onto an instantaneous rate matrix, A
A = projection.matrix.classe(pars=classe_params, k) 

# Calculate equilibrium frequencies by eigenvectors
evA <- eigen(A)
i <- which(evA$values == max(evA$values))
equilibrium_root_freqs = evA$vectors[, i]/sum(evA$vectors[, i])
equilibrium_root_freqs
# 0.2652666 0.2285983 0.2285983 0.2775368

loglik = log(sum(equilibrium_root_freqs * d_root_orig)) + sum(lq)
loglik
# -12.269765 matches!



# If root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
root.p = d_root_orig/sum(d_root_orig)
#lambda <- classe_params[E_indices]
e.root <- vals[E_indices]

# BiSSE
#d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)

# MuSSE/ClaSSE
pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)

# res1t: If root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# res2t: If root=ROOT.FLAT, root.p=NULL, condition.surv=TRUE
root.p = rep(1/numstates, times=numstates)
pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# res3t: If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root_probs = root_probs_equal
root.p = root_probs
pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# res4t: If root=ROOT.GIVEN, root.p=c(0.75,0.25), condition.surv=TRUE
root_probs = root_probs_biased
root.p = root_probs
pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# res5t: If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root_probs = root_probs_biased
root.p = root_probs
pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik


# res6t: Equilibrium root frequencies
# If root=ROOT.EQUI, condition.surv=TRUE
# Project the ClaSSE model onto an instantaneous rate matrix, A
A = projection.matrix.classe(pars=classe_params, k) 

# Calculate equilibrium frequencies by eigenvectors
evA <- eigen(A)
i <- which(evA$values == max(evA$values))
equilibrium_root_freqs = evA$vectors[, i]/sum(evA$vectors[, i])
equilibrium_root_freqs
# 0.2652666 0.2285983 0.2285983 0.2775368

pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(equilibrium_root_freqs * lambda * (1 - e.root)^2)
loglik = log(sum(equilibrium_root_freqs * d.root)) + sum(lq)
loglik
# -12.94599 matches!


cft(LnLst2, numdigits_inbetween_have_fixed_digits=8)

init = t(attr(res2, "intermediates")$init)
init

base = t(attr(res2, "intermediates")$base)
base

apply(X=base[,4:6], MARGIN=2, FUN="*", exp(lq))

# Get Es,Ds matrix
Dindexes = (numstates+1):(numstates*2)
EsDs_branch_bottoms = base
EsDs_branch_bottoms[,Dindexes] = EsDs_branch_bottoms[,Dindexes] * exp(attr(res2, "intermediates")$lq)
EsDs_branch_bottoms[1,]
