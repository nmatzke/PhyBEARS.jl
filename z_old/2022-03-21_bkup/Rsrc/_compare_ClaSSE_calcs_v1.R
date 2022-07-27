
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
states = c(0, 0, 0)		# Tip states
names(states) = tr$tip.label
states

sampling.f = c(1,1)		# Proportion of species in each state; for 2 states
											# (Let's assume we have all species)

# Create the BiSSE likelihood function. 
# (strict=FALSE means that some states in the state space can be absent from the tips)
bisse_2areas = make.bisse(tree=tr, states=states, sampling.f=sampling.f, strict=FALSE)

# Input some parameters
birthRate = 0.222222
deathRate = 0.0
transitionRate = 0.01

# The names of the parameters:
param_names = argnames(bisse_2areas)
param_names

# Most parameters will be zero
bisse_params = rep(0, times=length(param_names))
names(bisse_params) = param_names

# All extinction rates are the same (state-independent)
# Here, deathRate is 0 for all states
bisse_params[grepl(pattern="lambda", x=param_names)] = birthRate
bisse_params[grepl(pattern="mu", x=param_names)] = deathRate
bisse_params[param_names == "q01"] = transitionRate
bisse_params[param_names == "q10"] = transitionRate

# This is basically an Mk, 2-state equal rates model, with a pure-birth (Yule) process.

# To see the function:
dput(bisse_2areas)

bisse_2areas_default <- function(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL, intermediates=FALSE) 
	{
	check.pars.bisse(pars)
	preset <- branches.unresolved.bisse(pars, unresolved)
	ans <- all.branches(pars, intermediates, preset)
	rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
	}


# Do the BiSSE calculation, under many different assumptions
res1 = bisse_2areas(pars=bisse_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
res2 = bisse_2areas(pars=bisse_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.5,0.5)
res3 = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.75, 0.25)
res4 = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(1, 0)
res5 = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)

res1t = bisse_2areas(pars=bisse_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
res2t = bisse_2areas(pars=bisse_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.5,0.5)
res3t = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.75, 0.25)
res4t = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(1, 0)
res5t = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)

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
root.p = d.root/sum(d.root)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.FLAT, root.p=NULL, condition.surv=FALSE
root.p = rep(1/nstates, times=nstates)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root.p = c(0.5,0.5)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.75,0.25), condition.surv=FALSE
root.p = c(0.75,0.25)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root.p = c(0.5,0.5)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE
root.p = c(1.0,0.0)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik





# If root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
root.p = d_root_orig/sum(d_root_orig)
lambda <- bisse_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.FLAT, root.p=NULL, condition.surv=TRUE
root.p = rep(1/nstates, times=nstates)
lambda <- bisse_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root.p = c(0.5,0.5)
lambda <- bisse_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.75,0.25), condition.surv=TRUE
root.p = c(0.75,0.25)
lambda <- bisse_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root.p = c(0.5,0.5)
lambda <- bisse_params[E_indices]
e.root <- vals[E_indices]
d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
loglik = log(sum(root.p * d.root)) + sum(lq)
loglik

# If root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE
root.p = c(1.0,0.0)
lambda <- bisse_params[E_indices]
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
