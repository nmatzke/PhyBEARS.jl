
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
fit <- find.mle(func=lik, x.init=pars, method="subplex", fail.value=-1e10, lower=lower, upper=upper)
coef(fit)
st <- asr.marginal(lik, coef(fit))
nodelabels(thermo=t(st), piecol=cols, cex=0.5)


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



EsDs = t(attr(res1t,"intermediates")$init)
sum(log(rowSums(EsDs[,3:4])))
attr(res1t,"intermediates")$lq

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
phy2$edge.length[c(1,4)] = c(1-phy2$edge.length[2], 1-phy2$edge.length[6])
phy2$tip.state[1:4] = c(0, 1, 0, 1)

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

#######################################################
# Edit to clarify node structure
#######################################################
st2[,1] = c(0.5, 0.5)
st2[,2] = c(0.75, 0.25)
st2[,3] = c(1.0, 0.0)

nodelabels(thermo=t(st2), piecol=cols, cex=0.5)


# So, this shows nodes 5,6,7:
t(st2)
