
#######################################################
# Compare ClaSSE and BiSSE calculations
#
# E.g.:
# diversitree versus plain-R
# diversitree versus BioGeoBEARS+Yule+BFs
# 
#
# Run this script with:
# 
# source("/GitHub/PhyBEARS.jl/Rsrc/_compare_ClaSSE_calcs_n4_compare2julia.R")
# 
#######################################################

library(ape)
library(diversitree)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

source("/GitHub/PhyBEARS.jl/Rsrc/ClaSSE_functions_v3.R")  # utility functions from diversitree
source("/GitHub/PhyBEARS.jl/Rsrc/ClaSSE_mods_v2.R") # assisting functions for ClaSSE models
source("/GitHub/PhyBEARS.jl/Rsrc/ClaSSE_pureR_v1.R") # simple implementations in plain-R


# Load simple example tree
wd = "/GitHub/PhyBEARS.jl/data/"
setwd(wd)
trfn = "treeorang.newick"
tr = read.tree(trfn)

trstr = "(((chimp:1,human:1):1,gorilla:2):1,orang:3);"
tr = read.tree(file="", text=trstr)

# Run a BiSSE model from diversitree

# Setup
states = c(2,1,2,2)		# Tip states
names(states) = tr$tip.label
states

sampling.f = c(1,1,1)		# Proportion of species in each state; for 3 states
											# (Let's assume we have all species)
k = length(sampling.f)

# Create the BiSSE likelihood function. 
# (strict=FALSE means that some states in the state space can be absent from the tips)
classe_3states = make.classe(tree=tr, states=states, k=k, sampling.f=sampling.f, strict=FALSE)

# Input some parameters
birthRate = 0.2
deathRate = 0.1
d_val = 0.0
e_val = 0.0
j_val = 0.0

# The names of the parameters:
param_names = argnames(classe_3states)
param_names

# Most parameters will be zero
classe_params = rep(0, times=length(param_names))
names(classe_params) = param_names


# This is basically a DEC model for 3 states

# All extinction rates are the same (state-independent)
# Here, deathRate is 0 for all states
#classe_params[grepl(pattern="lambda", x=param_names)] = birthRate
classe_params[grepl(pattern="mu", x=param_names)] = deathRate
classe_params[grepl(pattern="q", x=param_names)] = 0

# For DEC
classe_params[param_names == "q31"] = e_val
classe_params[param_names == "q32"] = e_val
classe_params[param_names == "q13"] = d_val
classe_params[param_names == "q23"] = d_val
classe_params



# The birthRate (lambda) is state-independent.  However, 
# only certain cladogenesis scenarios are allowed under DEC.
#
# Disallowed cladogenesis scenarios have a rate of 0.
#
# If there is more than one cladogenesis scenario conditional 
# on a certain ancestor, DEC assigns each a weight of 1, and 
# then divides by the sum of the weights. I.e., if there are
# six possible cladogenetic range-inheritance events, they 
# each get a conditional probability of 1/6.
# 
# To translate to ClaSSE, if the speciation rate for a lineage 
# in a certain state is lambda, then the rate of each individual 
# allowed scenario would be lambda * 1/6
# 
y_val = (3-j_val)/3
total_of_weights = y_val + j_val + j_val
yprob = y_val / total_of_weights
jprob = j_val / total_of_weights


# Specifying the nonzero lambdas
# Narrow sympatry (ancestor A or B; rangesize of 1 area)
classe_params[param_names=="lambda111"] = yprob * birthRate
classe_params[param_names=="lambda222"] = yprob * birthRate
classe_params[param_names=="lambda333"] = 0

# Jump dispersal speciation
classe_params[param_names=="lambda112"] = jprob * birthRate
classe_params[param_names=="lambda121"] = jprob * birthRate
classe_params[param_names=="lambda212"] = jprob * birthRate
classe_params[param_names=="lambda221"] = jprob * birthRate

# Subset sympatry for state AB
# classe_params[param_names=="lambda312"] = 1/6 * birthRate
# classe_params[param_names=="lambda321"] = 1/6 * birthRate
# classe_params[param_names=="lambda313"] = 1/6 * birthRate
# classe_params[param_names=="lambda331"] = 1/6 * birthRate
# classe_params[param_names=="lambda323"] = 1/6 * birthRate
# classe_params[param_names=="lambda332"] = 1/6 * birthRate

# For diversitree ClaSSE, you have to lump lambda312 and lambda321
classe_params[param_names=="lambda312"] = 1/3 * birthRate
classe_params[param_names=="lambda313"] = 1/3 * birthRate
classe_params[param_names=="lambda323"] = 1/3 * birthRate

classe_params_DEC = classe_params




# To see the function:
dput(classe_3states)

classe_3states_default <- function(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL, intermediates=FALSE) 
	{
	## Note that this uses MuSSE's cache...
	pars2 <- f.pars(pars)
	ans <- all.branches(pars2, intermediates)
	ans$branchLnL = sum(ans$lq)
	rootfunc.classe(ans, pars, condition.surv, root, root.p, intermediates)
	}


# Do the ClaSSE calculation, under many different assumptions
res1 = classe_3states(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
res2 = classe_3states(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(1/3, 1/3, 1/3)
res3 = classe_3states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.1, 0.1, 0.8)
res4 = classe_3states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0, 0, 1)
res5 = classe_3states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
res6 = classe_3states(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)

res1t = classe_3states(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
res2t = classe_3states(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(1/3, 1/3, 1/3)
res3t = classe_3states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.1, 0.1, 0.8)
res4t = classe_3states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0, 0, 1)
res5t = classe_3states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
res6t = classe_3states(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)

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
cft(LnLst2, numdigits_inbetween_have_fixed_digits=8)



init = t(attr(res2, "intermediates")$init)
init

lq = t(attr(res2, "intermediates")$lq)
lq

base = t(attr(res2, "intermediates")$base)
base


base_likes = apply(X=base[,4:6], MARGIN=2, FUN="*", exp(lq))
base_normlikes = base_likes / rowSums(base_likes)


# Matches diversitree BiSSE "lq"
lnls = log(rowSums(base_likes))
lnls[!is.finite(lnls)] = NA
lnls
# -0.1903792 -0.1903792 -0.3632521 -0.5213903         NA -5.0488890 -5.3927962
sum(lnls, na.rm=TRUE)
# -11.70709



EsDs = t(attr(res1t,"intermediates")$init)
Ds = EsDs[,((ncol(EsDs)/2)+1):ncol(EsDs)]
sum_log_Ds = sum(log(rowSums(Ds))); sum_log_Ds
branch_lqs = attr(res1t,"intermediates")$lq; branch_lqs
sum(attr(res1t,"intermediates")$lq)
sum(branch_lqs)

# This corresponds to:
# Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
# R_sum_lq_nodes = R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
# OLD, for 2 states, 4 EsDs columns:
# sum(log(rowSums(EsDs[,3:4]))) + sum(attr(res1t,"intermediates")$lq)
sum_log_Ds + sum(branch_lqs)

# ...but is double-counting lnLs



# Sum of the branch likelihoods
lq
sum(lq)
# -0.2757956 -0.2757956 -0.511628 -0.717602    0 -3.698663 -3.685024
# -9.164509

# Add the root probabilities
# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
LnLs1
# -12.609029  -9.164509


# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
LnLs1t
# -10.538667  -9.164509

# Does the total of branch likelihoods (lq) + node likelihoods match R?
computed_likelihoods_at_each_node_x_lambda = rep(0.0, times=tr$Nnode + length(tr$tip.label))
likes_at_node5 = rep(0.0, times=ncol(base_likes))
likes_at_node4 = rep(0.0, times=ncol(base_likes))

# Internal node
state_i = 1
likes_at_node5[state_i] = sum(base_normlikes[1,state_i] * base_normlikes[2,state_i] * birthRate) # sympatry
state_i = 2
likes_at_node5[state_i] = sum(base_normlikes[1,state_i] * base_normlikes[2,state_i] * birthRate) # sympatry
state_i = 3
likes_at_node5[state_i] = sum(base_normlikes[1,2] * base_normlikes[2,1] * 1/6*birthRate) # vicariance
# Add small probs for subset sympatry
likes_at_node5[state_i] = likes_at_node5[state_i] + sum(base_normlikes[1,2] * base_normlikes[2,3] * 1/6*birthRate)
likes_at_node5[state_i] = likes_at_node5[state_i] + sum(base_normlikes[2,1] * base_normlikes[1,3] * 1/6*birthRate)
likes_at_node5

# Root node
state_i = 1
likes_at_node4[state_i] = sum(base_normlikes[3,state_i] * base_normlikes[5,state_i] * birthRate) # sympatry
state_i = 2
likes_at_node4[state_i] = sum(base_normlikes[3,state_i] * base_normlikes[5,state_i] * birthRate) # sympatry
state_i = 3
likes_at_node4[state_i] = sum(base_normlikes[3,2] * base_normlikes[5,3] * 1/6*birthRate) # subset sympatry
likes_at_node4



computed_likelihoods_at_each_node_just_before_speciation = get_sum_log_computed_likes_at_each_node(tr, base, lq, classe_params)
computed_likelihoods_at_each_node_just_before_speciation
rowSums(computed_likelihoods_at_each_node_just_before_speciation)
log(rowSums(computed_likelihoods_at_each_node_just_before_speciation))
TF = is.finite(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)))
sum(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)[TF]))
# 0.000000 0.000000 0.000000 0.000000 2.127192 2.063961 1.538840
# -Inf      -Inf      -Inf      -Inf 0.7548027 0.7246270 0.4310292
# 1.910459

R_result_sum_log_computed_likelihoods_at_each_node = sum(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)[TF]))
R_result_sum_log_computed_likelihoods_at_each_node
log(R_result_sum_log_computed_likelihoods_at_each_node)
sum(log(R_result_sum_log_computed_likelihoods_at_each_node))
# [1] NA 0.004184687
# log(R_result_sum_log_computed_likelihoods_at_each_node)
# NA -5.476323
# sum(log(R_result_sum_log_computed_likelihoods_at_each_node))
# NA

R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = R_result_sum_log_computed_likelihoods_at_each_node + sum(lq)
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
# -7.25405



R_result_branch_lnL = -9.164509
R_result_total_LnLs1 = -12.609029
R_result_total_LnLs1t = -10.538667
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -7.25405


















#######################################################
# Connect to BioGeoBEARS
#######################################################
numstates = k

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
# This version stores all of the piece, for comparison
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
sumBGBlike_not_root + bd_ape$lnl_numBirths + bd_ape$lnl_Births_above_root + bd_ape$lnl_branching_times
sumBGBlike_not_root + bd_ape$lnL - bd_ape$lnl_topology
sumBGBlike_not_root + sum_branchBot_likes

# Matches the branch_lnLs!
# For DEC+J, without anagenetic change, matches the branch_lnLs!


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



# Match BioGeoBEARS to diverstree res2
root.p = rep(1/numstates, times=numstates)
rootlikes = log(sum(root.p * d_root_orig_BGB))
rootlikes

sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
c(res2)

# Match BioGeoBEARS to diverstree res3 (all equal, except null range)
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



# Match BioGeoBEARS to diverstree res5 (all 1s, except null range)
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
# OLDER
#######################################################

# Key parts of the calculation
lq = t(attr(res2, "intermediates")$lq)			# Branch likelihoods
vals = t(attr(res2, "intermediates")$vals)	# Es and Ds at the root
nstates = length(vals) / 2
E_indices = 1:nstates
d_root_orig = vals[-E_indices]							# Original D likelihoods at root

# If root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
root.p = d_root_orig/sum(d_root_orig)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# If root=ROOT.FLAT, root.p=NULL, condition.surv=FALSE
root.p = rep(1/nstates, times=nstates)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root.p = c(0.3333333, 0.3333333, 0.3333333)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.75,0.25), condition.surv=FALSE
root.p = c(0.1, 0.1, 0.8)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root.p = c(0, 0, 1)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# If root=ROOT.EQUI, condition.surv=FALSE

# Project the ClaSSE model onto an instantaneous rate matrix, A
A = projection.matrix.classe(pars=classe_params, k) 

# Calculate equilibrium frequencies by eigenvectors
evA <- eigen(A)
 		# https://www-sciencedirect-com.ezproxy.auckland.ac.nz/topics/mathematics/dominant-eigenvalue
		# Predicting Population Growth: Modeling with Projection Matrices
		# Janet Steven, James Kirkwood, in Mathematical Concepts and Methods in Modern Biology, 2013"
		# 7.8.4 Finding the Stable Distribution
		#
		# Suppose that A is a projection matrix that meets the assumptions of the Perron-Frobenius 
		# theorem and that va is any vector.
		# 
		# ...so the equilibrium state is the normalized eigenvector for the dominant eigenvalue.
# OLD: i <- which(evA$values == max(evA$values))
i <- which(abs(evA$values) == max(abs(evA$values)))
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

loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.FLAT, root.p=NULL, condition.surv=TRUE
root.p = rep(1/nstates, times=nstates)
pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root.p = c(0.3333333, 0.3333333, 0.3333333)
pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.75,0.25), condition.surv=TRUE
root.p = c(0.1, 0.1, 0.8)
pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root.p = c(0, 0, 1)
pars = classe_params
nsum <- k * (k + 1)/2
lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
i <- seq_len(k)
e.root <- vals[i]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik


yule(tr)$loglik - log(birthRate)


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
Dindexes = (nstates+1):(nstates*2)
EsDs_branch_bottoms = base
EsDs_branch_bottoms[,Dindexes] = EsDs_branch_bottoms[,Dindexes] * exp(attr(res2, "intermediates")$lq)
EsDs_branch_bottoms[1,]


# Yay! Figured it out!

#projection.matrix.classe <- function(pars, k) 









LnLs = rbind(LnLs1, LnLs2, LnLs3, LnLs4, LnLs5)
print(LnLs)

init = t(attr(res2, "intermediates")$init)
init

base = t(attr(res2, "intermediates")$base)
base
# Columns 4-6, the Ds, sum to 1
rowSums(base[,4:6])



lq = attr(res2, "intermediates")$lq
lq
sum(lq)

rowSums(base[,4:6]) * exp(lq)
base[,4] * exp(lq)
base[,5] * exp(lq)
base[,6] * exp(lq)
