
library(BioGeoBEARS)
library(diversitree)
source("/GitHub/PhyBEARS.jl/Rsrc/ClaSSE_mods_v2.R")
source("/GitHub/PhyBEARS.jl/Rsrc/ClaSSE_functions_v3.R")

## Start with a simple tree evolved under a BiSSE with all rates
## asymmetric:
# lambda0, lambda1, mu0, mu1, q01, q10

# Gives weird result -- green nodes, blue truth
set.seed(47) # Rare to have a transition, but observed here
pars <- c(0.222222222, 0.222222222, 0.111111111, 0.05, 0.1, 0.15)
orig_pars = pars
bisse_params = pars
birthRate = pars[1]
deathRate = pars[3]

set.seed(48) # Rare to have a transition, but observed here
phy <- trees(pars, "bisse", max.taxa=4, max.t=Inf, x0=0)[[1]]
write.tree(phy, file="")

cols = c("blue", "green3")
lower = rep(0.0, times=length(pars))
upper = rep(10.0, times=length(pars))

## Here is the true history
h <- history.from.sim.discrete(phy, 0:1)
diversitree:::plot.history(h, phy, main="True history", cols=cols)

lik <- make.bisse(phy, phy$tip.state)
phy$tip.state

bisse_2areas = lik
#fit <- find.mle(lik, pars, method="subplex")
# MLE doesn't make much sense with tiny data
fit <- find.mle(func=lik, x.init=pars, method="subplex", fail.value=-1e10, lower=lower, upper=upper)
coef(fit)
st <- asr.marginal(lik, coef(fit))
nodelabels(thermo=t(st), piecol=cols, cex=0.5)
t(st)



res1 = bisse_2areas(pars=bisse_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
res2 = bisse_2areas(pars=bisse_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0, 1)
res3 = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.25, 0.75)
res4 = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.5,0.5)
res5 = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
res6 = bisse_2areas(pars=bisse_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)

# Res1t is BiSSE default
# "lik" matches res1t
res1t = bisse_2areas(pars=bisse_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
res2t = bisse_2areas(pars=bisse_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0, 1)
res3t = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.25, 0.75)
res4t = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.5,0.5)
res5t = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
res6t = bisse_2areas(pars=bisse_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)


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

# Compare to Julia:
# > LnLs1
# [1] -9.574440 -6.670978
# > LnLs1t
# [1] -7.464283 -6.670978


EsDs = t(attr(res1t,"intermediates")$init)
sum(log(rowSums(EsDs[,3:4])))
attr(res1t,"intermediates")$lq
sum(attr(res1t,"intermediates")$lq)

# This corresponds to:
# Julia_sum_lq_nodes = sum(log.(sum.(res.likes_at_each_nodeIndex_branchTop))) + Julia_sum_lq
# R_sum_lq_nodes = R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
sum(log(rowSums(EsDs[,3:4]))) + sum(attr(res1t,"intermediates")$lq)
# ...but is double-counting lnLs




#######################################################
# Look at ASRs
#######################################################
tmpmat = rbind(attr(res1,"intermediates")$root.p, 
attr(res2,"intermediates")$root.p, 
attr(res3,"intermediates")$root.p, 
attr(res4,"intermediates")$root.p, 
attr(res5,"intermediates")$root.p, 
attr(res6,"intermediates")$root.p, 
attr(res1t,"intermediates")$root.p, 
attr(res2t,"intermediates")$root.p, 
attr(res3t,"intermediates")$root.p, 
attr(res4t,"intermediates")$root.p, 
attr(res5t,"intermediates")$root.p, 
attr(res6t,"intermediates")$root.p
)

tmpnames = paste0("ancstate", seq(1,ncol(tmpmat),1))

tmprownames = c("res1", "res2", "res3", "res4", "res5", "res6", "res1t", "res2t", "res3t", "res4t", "res5t", "res6t")

rootstates_df = adf2(tmpmat)
names(rootstates_df) = tmpnames
row.names(rootstates_df) = tmprownames
rootstates_df

# rootstates_df
#       ancstate1 ancstate2
# res1  0.4303571 0.5696429 # <- this is what asr.marginal gives
# res2  0.5000000 0.5000000
# res3  0.0000000 1.0000000
# res4  0.2500000 0.7500000
# res5  0.5000000 0.5000000
# res6  0.5392658 0.4607342
# res1t 0.4303571 0.5696429 # <- this is what asr.marginal gives
# res2t 0.5000000 0.5000000
# res3t 0.0000000 1.0000000
# res4t 0.2500000 0.7500000
# res5t 0.5000000 0.5000000
# res6t 0.5392658 0.4607342

# Same for everything:
rootnode = length(phy$tip.label)+1
EsDs = t(attr(res1,"intermediates")$init)
EsDs[rootnode, 3:4] / sum(EsDs[rootnode, 3:4])


#######################################################
# TRYING TO FIT TO SIMULATED DATA
# This **FAILS** with such tiny data
#######################################################



#######################################################
# Look at find.mle
#######################################################
methods("find.mle")
diversitree:::find.mle.default
diversitree:::find.mle.dtlik
diversitree:::find.mle.mixed
diversitree:::do.mle.search


# change to create some ambiguity

phy2 = phy

plot(phy2, label.offset=0.1)
tiplabels(text=phy2$tip.state, tip=1:4, col="black", bg=cols[phy2$tip.state+1])
#diversitree:::plot.history(h, phy, main="True history", cols=cols)

lik2 <- make.bisse(phy2, phy2$tip.state)
fit2 <- find.mle(func=lik2, x.init=coef(fit), method="subplex", fail.value=-1e10, lower=lower, upper=upper)
coef(fit2 )
#fit2 <- find.mle(lik2, pars, method="subplex", fail.value=-1e10)
st2 <- asr.marginal(lik2, coef(fit2))
nodelabels(thermo=t(st2), piecol=cols, cex=0.5)
t(st2)

orig_pars = pars
names(pars) = names(coef(fit2))
st2 <- asr.marginal(lik2, pars)
t(st2)
st2 <- asr.marginal(lik2, pars, root=ROOT.OBS, root.p=NULL, condition.surv=FALSE)
t(st2)

# t(st2)
#             [,1]         [,2]
# [1,] 0.430357148 0.5696428522  # <- matches res1 or res1t
# [2,] 0.005145312 0.9948546878
# [3,] 0.999681941 0.0003180585
st2 <- asr.marginal(lik2, pars, root=ROOT.OBS, condition.surv=TRUE)
t(st2)
# > t(st2)
#           [,1]        [,2]
# [1,] 0.49926227 0.500737731
# [2,] 0.00570688 0.994293120
# [3,] 0.99963639 0.000363611







# Turn off extinction; the speciation/extinction thing dominates on large branches I guess
pars["mu0"] = 0
pars["mu1"] = 0
st2 <- asr.marginal(lik2, pars, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
t(st2)
st2 <- asr.marginal(lik2, pars)
t(st2)

#######################################################
# Edit to clarify node structure
#######################################################
st2[,1] = c(0.5, 0.5)
st2[,2] = c(0.75, 0.25)
st2[,3] = c(1.0, 0.0)

nodelabels(thermo=t(st2), piecol=cols, cex=0.5)


# So, this shows nodes 5,6,7:
t(st2)
