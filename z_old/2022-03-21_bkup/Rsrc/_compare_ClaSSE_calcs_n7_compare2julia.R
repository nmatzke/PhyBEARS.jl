
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
trfn = "tree_small.newick"
tr = read.tree(trfn)

trstr = "((chimp:1,human:1):1,gorilla:2);"
tr = read.tree(file="", text=trstr)

# Run a BiSSE model from diversitree

# Setup
states = c(2,1,2)		# Tip states
names(states) = tr$tip.label
states

sampling.f = c(1,1,1)		# Proportion of species in each state; for 2 states
											# (Let's assume we have all species)
k = length(sampling.f)

# Create the BiSSE likelihood function. 
# (strict=FALSE means that some states in the state space can be absent from the tips)
classe_3states = make.classe(tree=tr, states=states, k=k, sampling.f=sampling.f, strict=FALSE)

# Input some parameters
birthRate = 0.2
deathRate = 1.0
d_val = 0.5
e_val = 0.4
j_val = 1.5

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
# We do x2, because lambda112 covers lambda121
classe_params[param_names=="lambda112"] = jprob * birthRate * 2
classe_params[param_names=="lambda121"] = jprob * birthRate * 2
classe_params[param_names=="lambda212"] = jprob * birthRate * 2
classe_params[param_names=="lambda221"] = jprob * birthRate * 2

# Vicariance etc. -- still just 1/6 as there are 
# no j events from AB in a 2-area system
# We do x2, because lambda312 covers lambda321
classe_params[param_names=="lambda312"] = 1/6 * birthRate * 2
classe_params[param_names=="lambda321"] = 1/6 * birthRate * 2
classe_params[param_names=="lambda313"] = 1/6 * birthRate * 2
classe_params[param_names=="lambda331"] = 1/6 * birthRate * 2
classe_params[param_names=="lambda323"] = 1/6 * birthRate * 2
classe_params[param_names=="lambda332"] = 1/6 * birthRate * 2

# For diversitree ClaSSE, you have to lump lambda312 and lambda321
# classe_params[param_names=="lambda312"] = 1/3 * birthRate
# classe_params[param_names=="lambda313"] = 1/3 * birthRate
# classe_params[param_names=="lambda323"] = 1/3 * birthRate
classe_lambdas_to_df(classe_params, k=3)

classe_params_DEC = classe_params


classe_lambdas_to_df(classe_params, k=3)

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
root_probs = c(0.3333333, 0.3333333, 0.3333333)
res3 = classe_3states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.1, 0.1, 0.8)
res4 = classe_3states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0, 0, 1)
res5 = classe_3states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
res6 = classe_3states(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)

res1t = classe_3states(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
res2t = classe_3states(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.3333333, 0.3333333, 0.3333333)
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
