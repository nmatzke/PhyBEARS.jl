
#######################################################
# Compare ClaSSE and BiSSE calculations
#
# E.g.:
# diversitree versus plain-R
# diversitree versus BioGeoBEARS+Yule+BFs
# 
#######################################################

library(ape)
library(diversitree)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

source("/GitHub/BioGeoJulia.jl/Rsrc/ClaSSE_functions_v3.R")  # utility functions from diversitree
source("/GitHub/BioGeoJulia.jl/Rsrc/ClaSSE_pureR_v1.R") # simple implementations in plain-R
source("/GitHub/BioGeoJulia.jl/Rsrc/ClaSSE_mods_v2.R") # simple implementations in plain-R


# Load simple example tree
wd = "/GitHub/BioGeoJulia.jl/Rsrc/"
setwd(wd)
trfn = "Psychotria_tree.newick"
tr = read.tree(trfn)

trstr = "((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);"
tr = read.tree(file="", text=trstr)

# Run a ClaSSE model from diversitree

# Setup
# Island numbers (KOMH = 1234) in Rnodenums order:
island_nums = c(3, 3, 2, 2, 3, 3, 2, 1, 1, 3, 4, 2, 1, 1, 1, 1, 1, 1, 2)
states = c(island_nums + 1)		# Tip states
names(states) = tr$tip.label
states

sampling.f = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)		# Proportion of species in each state; for 2 states
											# (Let's assume we have all species)
k = length(sampling.f)

# Create the BiSSE likelihood function. 
# (strict=FALSE means that some states in the state space can be absent from the tips)
classe_16states = make.classe(tree=tr, states=states, k=k, sampling.f=sampling.f, strict=FALSE)

# The names of the parameters:
param_names = argnames(classe_16states)
param_names

# Most parameters will be zero
classe_params = rep(0, times=length(param_names))
names(classe_params) = param_names



#######################################################
# BioGeoBEARS Q and C matrices
#######################################################
library(BioGeoBEARS)
trfn = "Psychotria_tree.newick"
geogfn = "Psychotria_geog.data"
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

res = bears_optim_run(BioGeoBEARS_run_object)
mats = get_Qmat_COOmat_from_res(res, numstates=ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node), include_null_range=TRUE, max_range_size=res$inputs$max_range_size, timeperiod_i=1)
numstates = length(mats$states_list)

source('/GitHub/BioGeoBEARS/R/BioGeoBEARS_extract_Qmat_COOmat_v1.R', chdir = TRUE)
source('/GitHub/BioGeoJulia.jl/Rsrc/ClaSSE_mods_v2.R', chdir = TRUE)
Carray_df = get_Cevent_probs_df_from_mats(mats, include_null_range=include_null_range)
head(Carray_df)
tail(Carray_df)

lambda_ijk_df = classe_lambdas_to_df(classe_params, k=numstates)
head(lambda_ijk_df)

# Put the Qmat and res parameter values from BioGeoBEARS "mats" into classe_params
# Input some parameters
yule(tr)

birthRate = 0.3288164
deathRate = 0.0
d_val = 0.03505038
e_val = 0.02832370
j_val = 0.0

# Blank out the params
classe_params = rep(0, times=length(param_names))
names(classe_params) = param_names
head(classe_params)

# Fill in the params from the BioGeoBEARS "res" results
classe_params = BGBres_into_classe_params(res, classe_params, birthRate=birthRate)
classe_params[153]
classe_params["lambda020202"]
classe_params["lambda060203"]
classe_params["q0206"]
classe_params["q1516"]
classe_params[classe_params != 0.0]

root_probs_equal = rep(1, times=numstates)
root_probs_equal[sum(include_null_range)] = 0
root_probs_equal = root_probs_equal / sum(root_probs_equal)

root_probs_biased = rep(0.01, times=numstates)
root_probs_biased[sum(include_null_range)] = 0
root_probs_biased[length(root_probs_biased)] = 0.01 * (numstates-include_null_range)
root_probs_biased = root_probs_biased / sum(root_probs_biased)

root_probs_single = rep(1, times=numstates)
root_probs_single[sum(include_null_range)] = 0
#root_probs_single[1+include_null_range] = 1
#root_probs_single = root_probs_single / sum(root_probs_single)

# Do the ClaSSE calculation, under many different assumptions
res1 = classe_16states(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
res2 = classe_16states(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
root_probs = root_probs_equal
res3 = classe_16states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = root_probs_biased
res4 = classe_16states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = root_probs_single
res5 = classe_16states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
res6 = classe_16states(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)

res1t = classe_16states(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
res2t = classe_16states(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
root_probs = root_probs_equal
res3t = classe_16states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = root_probs_biased
res4t = classe_16states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = root_probs_single
res5t = classe_16states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
res6t = classe_16states(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)

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
all_lnLs = cft(LnLst2, numdigits_inbetween_have_fixed_digits=8)
all_lnLs$ttl_LnL = as.numeric(all_lnLs$ttl_LnL)
cft(LnLst2, numdigits_inbetween_have_fixed_digits=8)


init = t(attr(res2, "intermediates")$init)
init

lq = t(attr(res2, "intermediates")$lq)
lq

base = t(attr(res2, "intermediates")$base)
base


Ds_cols = (numstates+1):(2*numstates)
base_likes = apply(X=base[,Ds_cols], MARGIN=2, FUN="*", exp(lq))
base_normlikes = base_likes / rowSums(base_likes)

base_likes / res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS

# These match!!!
round(base_normlikes - res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS, 5)

# Not match
round(base_normlikes - res$condlikes_of_each_state / rowSums(res$condlikes_of_each_state), 5)

tmp_rowSums = (base_likes / res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS)[,numstates]
tmp_rowSums
log(tmp_rowSums)
sum(log(tmp_rowSums), na.rm=TRUE) # matches sum(lq)
# Log-likelihoods

# Equal
log(tmp_rowSums) - attr(res4,"intermediates")$lq


# Diversitree birthdeath calculation
lik.bd <- make.bd(tree=tr, sampling.f=NULL, unresolved=NULL, times=NULL, control=list(method="ode"))
diversitree_bd = lik.bd(pars=c(birthRate=0.3288164, deathRate=0.0), intermediates=TRUE)
bd_lq = attr(diversitree_bd,"intermediates")$lq

# FIGURED IT OUT
bd_lq
log(exp(-birthRate * trtable$edge.length))
-birthRate * trtable$edge.length

sum(-birthRate * trtable$edge.length, na.rm=TRUE)
# -17

bd_lq[21:37] - (-birthRate * trtable$edge.length)[21:37]
# -1.112256 -1.112256 -1.112256....

exp(-1.112256)
#[1] 0.3288163

# 19 tips, 18 internal nodes
# birthRate = (17 speciation events above root)/total_brlen
sum(-birthRate * trtable$edge.length, na.rm=TRUE) + 17*log(birthRate)
sum(bd_lq)

sum(log(res$computed_likelihoods_at_each_node[-20]))
branches_above_node = (-birthRate * trtable$edge.length)[c(22,31)]
branches_above_node

sum(log(res$computed_likelihoods_at_each_node[-20])) + sum(branches_above_node)

res$total_loglikelihood + sum(bd_lq)
res$total_loglikelihood + sum(bd_lq)
sum(lq) - sum(bd_lq)

all_lnLs
sum(log(res$computed_likelihoods_at_each_node[-20])) +sum(bd_lq)
-67.62949 - -67.51729
-0.1122
exp(-0.1122)
# 0.8938655


# GOT IT F YEAH
all_lnLs
# -67.629498
sum(log(res$computed_likelihoods_at_each_node[-20])) +sum(bd_lq) + log(birthRate) + 1
# -67.62954



sum(rowSums(base_likes)[c(22,31)])


x = -birthRate * trtable$edge.length
exp(x[c(22,31)])
x[c(22,31)]
sum(x[c(22,31)])

# rows Sum to 1, so just summing "x", equals -17
y = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS * x
rowSums(y)
sum(rowSums(y), na.rm=TRUE)
# -17



rowSums(base_likes) /exp(x)
(rowSums(base_likes) /exp(x))[c(22,31)]

# This should be the character change process, times the tree process
x = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS * exp(attr(diversitree_bd,"intermediates")$lq)

rowSums(x)
log(rowSums(x))
sum(log(rowSums(x)), na.rm=TRUE)

log(tmp_rowSums) - attr(diversitree_bd,"intermediates")$lq


res$total_loglikelihood




attr(diversitree_bd,"intermediates")$lq


all_lnLs$ttl_LnL - res$total_loglikelihood - log(1/16)
LnLs4[1] - res$total_loglikelihood - log(1/16)

all_lnLs$ttl_LnL - res$total_loglikelihood
LnLs2[2] - res$total_loglikelihood



# Get the constant part of bd_lq
cache <- diversitree:::make.cache.bd(tree=tr, sampling.f=NULL, unresolved=NULL, times=NULL, control=list(method="ode"))
cache$const
# In diversitree:::rootfunc.bd.ode
# loglik <- log(d.root) + sum(lq) + const
diversitree_bd_lnL_minus_const = c(diversitree_bd) - cache$const - log(birthRate)
diversitree_bd_lnL_minus_const
# -35.90835

LnLs4[1] + cache$const + log(birthRate)
-74.33632

# diversitree res5:
# If root=ROOT.GIVEN, root.p=c(0.75,0.25), condition.surv=FALSE
root_probs_single

# Key parts of the calculation
lq = t(attr(res5, "intermediates")$lq)			# Branch likelihoods
vals = t(attr(res5, "intermediates")$vals)	# Es and Ds at the root
nstates = length(vals) / 2
E_indices = 1:nstates
d_root_orig = vals[-E_indices]							# Original D likelihoods at root

root.p = root_probs_single
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
log(sum(root.p * d_root_orig)) 
log(sum(d_root_orig))
loglik
c(res5)

d_root_orig / res$condlikes[20,]

log(sum(d_root_orig))
log(sum(res$condlikes[20,])) - 1
c(res5) - sum(lq)

res$total_loglikelihood - (log(sum(res$condlikes[20,])) - 1)

res$total_loglikelihood - log(sum(res$condlikes[20,]))

LnLs4
d_root_orig =


sum(lq) - cache$const

diversitree_bd_lnL_minus_const - res$total_loglikelihood

all_lnLs$ttl_LnL - sum(bd_lq)
LnLs2[2] - sum(bd_lq)

bd = bd_liks(tr, birthRate=birthRate, deathRate=deathRate)
bd

bd$lnl_Births_above_root + bd$lnl_branching_times
-17

bd$lnL + lnl_numBirths + bd$lnl_Births_above_root + bd$lnl_branching_times
bd$lnl_branching_times + bd$lnl_numBirths + bd$lnl_Births_above_root
bd$lnl_branching_times + bd$lnl_Births_above_root
bd$lnl_topology + bd$lnL

# Matches diversitree BiSSE "lq"
lnls = log(rowSums(base_likes))
lnls[!is.finite(lnls)] = NA
lnls
# -0.2757956 -0.2757956 -0.5116280         NA -3.6850244
sum(lnls, na.rm=TRUE)
# -5.085475



res = res1
claSSE_res_to_prt(res1, tr, classe_params)




# Sum of the branch likelihoods
lq
sum(lq)
# [1,] -1.092256 -1.092256 -1.998899    0 -3.154265
# [1] -7.337676

# Add the root probabilities
# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
LnLs1
# [1] -11.105443  -7.337676


# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
LnLs1t
# [1] -5.932152 -7.337676

# Does the total of branch likelihoods (lq) + node likelihoods match R?
computed_likelihoods_at_each_node_x_lambda = rep(0.0, times=tr$Nnode + length(tr$tip.label))

computed_likelihoods_at_each_node_just_before_speciation = get_sum_log_computed_likes_at_each_node(tr, base, lq, classe_params)
computed_likelihoods_at_each_node_just_before_speciation
rowSums(computed_likelihoods_at_each_node_just_before_speciation)
log(rowSums(computed_likelihoods_at_each_node_just_before_speciation))
TF = is.finite(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)))
sum(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)[TF]))
R_result_sum_log_computed_likelihoods_at_each_node = sum(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)[TF]))
# [1] 0.00000000 0.00000000 0.00000000 0.06912288 0.10564300
# [1]      -Inf      -Inf      -Inf -2.671869 -2.247690
# [1] -4.919559

R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = R_result_sum_log_computed_likelihoods_at_each_node + sum(lq)
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
# -9.861079



R_result_branch_lnL = -7.337676
R_result_total_LnLs1 = -11.105443
R_result_total_LnLs1t = -5.932152
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -12.25723











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
