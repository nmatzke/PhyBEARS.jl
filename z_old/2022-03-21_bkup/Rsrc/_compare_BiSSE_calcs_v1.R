
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

source("/GitHub/BioGeoJulia.jl/Rsrc/ClaSSE_mods_v2.R")  # utility functions from diversitree
source("/GitHub/BioGeoJulia.jl/Rsrc/ClaSSE_pureR_v1.R") # simple implementations in plain-R


# Load simple example tree
wd = "/GitHub/BioGeoJulia.jl/Rsrc/"
setwd(wd)
trfn = "tree_small.newick"
tr = read.tree(trfn)

trstr = "((chimp:1,human:1):1,gorilla:2);"

# Run a BiSSE model from diversitree

# Setup
states = c(2,2,2,3)		# Tip states
names(states) = tr$tip.label
states

sampling.f = c(1,1,1,1)		# Proportion of species in each state; for 2 states
											# (Let's assume we have all species)
k = length(sampling.f)

# Create the BiSSE likelihood function. 
# (strict=FALSE means that some states in the state space can be absent from the tips)
classe_4states = make.classe(tree=tr, states=states, k=k, sampling.f=sampling.f, strict=FALSE)

# Input some parameters
birthRate = 0.222222
deathRate = 0.0
d_val = 0.1
e_val = 0.02
j_val = 0.0

# The names of the parameters:
param_names = argnames(classe_4states)
param_names

# Most parameters will be zero
classe_params = rep(0, times=length(param_names))
names(classe_params) = param_names


# This is basically a DEC model for 4 states

# All extinction rates are the same (state-independent)
# Here, deathRate is 0 for all states
classe_params[grepl(pattern="lambda", x=param_names)] = birthRate
classe_params[grepl(pattern="mu", x=param_names)] = deathRate
classe_params[grepl(pattern="q", x=param_names)] = 0
classe_params[param_names == "q21"] = e
classe_params[param_names == "q31"] = e
classe_params[param_names == "q24"] = d
classe_params[param_names == "q34"] = d
classe_params[param_names == "q42"] = e
classe_params[param_names == "q43"] = e
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
# Null range cannot speciate (doesn't seem to matter anyway,
# as "null" cannot be an ancestor anyway)
classe_params[param_names=="lambda111"] = birthRate

# Narrow sympatry (ancestor A or B; rangesize of 1 area)
classe_params[param_names=="lambda222"] = yprob * birthRate
classe_params[param_names=="lambda333"] = yprob * birthRate

# Jump dispersal speciation
classe_params[param_names=="lambda223"] = jprob * birthRate
classe_params[param_names=="lambda323"] = jprob * birthRate

# Subset sympatry for state AB
classe_params[param_names=="lambda424"] = 1/6 * birthRate
classe_params[param_names=="lambda434"] = 1/6 * birthRate

# Vicariance for state AB
classe_params[param_names=="lambda423"] = 1/6 * birthRate

classe_params_DEC = classe_params




# To see the function:
dput(classe_4states)

classe_4states_default <- function(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL, intermediates=FALSE) 
	{
	## Note that this uses MuSSE's cache...
	pars2 <- f.pars(pars)
	ans <- all.branches(pars2, intermediates)
	ans$branchLnL = sum(ans$lq)
	rootfunc.classe(ans, pars, condition.surv, root, root.p, intermediates)
	}


# Do the BiSSE calculation, under many different assumptions
res1 = classe_4states(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
res2 = classe_4states(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.25,0.25,0.25,0.25)
res3 = classe_4states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.1, 0.1, 0.1, 0.7)
res4 = classe_4states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0, 0, 0, 1)
res5 = classe_4states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)

res1t = classe_4states(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
res2t = classe_4states(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.25,0.25,0.25,0.25)
res3t = classe_4states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.1, 0.1, 0.1, 0.7)
res4t = classe_4states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0, 0, 0, 1)
res5t = classe_4states(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)

# get_classe_LnLs returns the total log-likelihood, and 
# the total of the branch likelihoods
LnLs1 = get_classe_LnLs(res1)
LnLs2 = get_classe_LnLs(res2)
LnLs3 = get_classe_LnLs(res3)
LnLs4 = get_classe_LnLs(res4)
LnLs5 = get_classe_LnLs(res5)
LnLs1t = get_classe_LnLs(res1t)
LnLs2t = get_classe_LnLs(res2t)
LnLs3t = get_classe_LnLs(res3t)
LnLs4t = get_classe_LnLs(res4t)
LnLs5t = get_classe_LnLs(res5t)

LnLst = as.data.frame(rbind(LnLs1, LnLs2, LnLs3, LnLs4, LnLs5, LnLs1t, LnLs2t, LnLs3t, LnLs4t, LnLs5t), stringsAsFactors=FALSE)
names(LnLst) = c("ttl_LnL", "branch_LnL")
Ldiff = exp((LnLst$ttl_LnL - log(0.5)) - LnLst$branch_LnL)
LnLst2 = cbind(LnLst, Ldiff)
cft(LnLst2, numdigits_inbetween_have_fixed_digits=8)

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
root.p = c(0.25,0.25,0.25,0.25)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.75,0.25), condition.surv=FALSE
root.p = c(0.1, 0.1, 0.1, 0.7)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root.p = c(0, 0, 0, 1)
loglik = log(sum(root.p * d_root_orig)) + sum(lq)
loglik






# If root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
root.p = d_root_orig/sum(d_root_orig)
lambda <- classe_params[E_indices]
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
lambda <- classe_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root.p = c(0.5,0.5)
lambda <- classe_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.75,0.25), condition.surv=TRUE
root.p = c(0.75,0.25)
lambda <- classe_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root.p = c(0.5,0.5)
lambda <- classe_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root.p = c(1.0,0.0)
lambda <- classe_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik











LnLs = rbind(LnLs1, LnLs2, LnLs3, LnLs4, LnLs5)
print(LnLs)

init = t(attr(res2, "intermediates")$init)
init

base = t(attr(res2, "intermediates")$base)
base
# Columns 3-4 sum to 1
rowSums(base[,3:4])

lq = attr(res2, "intermediates")$lq
lq
sum(lq)
