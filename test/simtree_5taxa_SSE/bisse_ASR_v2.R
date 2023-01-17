

library(diversitree)

## Start with a simple tree evolved under a BiSSE with all rates
## asymmetric:
# lambda0, lambda1, mu0, mu1, q01, q10

# Gives weird result -- green nodes, blue truth
set.seed(47) # Rare to have a transition, but observed here
pars <- c(0.222222222, 0.222222222, 0.111111111, 0.222222222, 0.1, 0.05)

# Seems accurate
pars <- c(0.222222222, 0.222222222, 0.111111111, 0.05, 0.1, 0.15)

set.seed(59) # Rare to have a transition, but observed here
phy <- trees(pars, "bisse", max.taxa=5, max.t=Inf, x0=0)[[1]]
cols = c("blue", "green3")
lower = rep(0.0, times=length(pars))
upper = rep(10.0, times=length(pars))

## Here is the true history
h <- history.from.sim.discrete(phy, 0:1)
diversitree:::plot.history(h, phy, main="True history", cols=cols)

lik <- make.bisse(phy, phy$tip.state)
#fit <- find.mle(lik, pars, method="subplex")
fit <- find.mle(func=lik, x.init=pars, method="subplex", fail.value=-1e10, lower=lower, upper=upper)
parsdf = cft(adf2(matrix(coef(fit),nrow=1))); names(parsdf)=names(coef(fit))
parsdf

st <- asr.marginal(lik, coef(fit))
states_df = cft(adf2(t(st)))
states_df

nodelabels(thermo=t(st), piecol=cols, cex=0.5)




#######################################################
# Look at find.mle
#######################################################
methods("find.mle")
diversitree:::find.mle.default
diversitree:::find.mle.dtlik
diversitree:::find.mle.mixed
diversitree:::do.mle.search


# Gives weird result -- green nodes, blue truth
set.seed(47) # Rare to have a transition, but observed here
pars <- c(0.222222222, 0.222222222, 0.111111111, 0.222222222, 0.1, 0.05)

# Seems accurate
pars <- c(0.222222222, 0.222222222, 0.111111111, 0.05, 0.1, 0.15)

set.seed(59) # Rare to have a transition, but observed here
phy <- trees(pars, "bisse", max.taxa=10, max.t=Inf, x0=0)[[1]]
cols = c("blue", "green3")
lower = rep(0.0, times=length(pars))
upper = rep(10.0, times=length(pars))

## Here is the true history
h <- history.from.sim.discrete(phy, 0:1)
diversitree:::plot.history(h, phy, main="True history", cols=cols)

plot(phy, label.offset=0.5, main="Inferred history")
lik <- make.bisse(phy, phy$tip.state)
#fit <- find.mle(lik, pars, method="subplex")
fit <- find.mle(func=lik, x.init=pars, method="subplex", fail.value=-1e10, lower=lower, upper=upper)
parsdf = cft(adf2(matrix(coef(fit),nrow=1))); names(parsdf)=names(coef(fit))
parsdf

st <- asr.marginal(lik, coef(fit))
states_df = cft(adf2(t(st)))
states_df

nodelabels(thermo=t(st), piecol=cols, cex=0.5)




#######################################################
# Small tree
#######################################################

set.seed(6534) # Rare to have a transition, but observed here
pars <- c(0.222222222, 0.222222222, 0.111111111, 0.0111111111, 0.06, 0.05)

phy <- trees(pars, "bisse", max.taxa=5, max.t=Inf, x0=0)[[1]]
cols = c("blue", "green3")
lower = rep(0.0, times=length(pars))
upper = rep(10.0, times=length(pars))

## Here is the true history
h <- history.from.sim.discrete(phy, 0:1)
diversitree:::plot.history(h, phy, main="True history", cols=cols)

lik <- make.bisse(phy, phy$tip.state)
lik(pars)
# -13.08989

reslnls = lik(pars, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# -15.63256
sum(attr(reslnls,"intermediates")$lq)
# -12.87511


reslnls = lik(pars, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
# -12.87511
sum(attr(reslnls,"intermediates")$lq)
# -12.87511


reslnls = lik(pars, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# -15.63256
sum(attr(reslnls,"intermediates")$lq)
# -12.87511


reslnls = lik(pars, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
# -12.87511
sum(attr(reslnls,"intermediates")$lq)
# -12.87511




#fit <- find.mle(lik, pars, method="subplex")
fit <- find.mle(func=lik, x.init=pars, method="subplex", fail.value=-1e10, lower=lower, upper=upper)

# ML on this tiny dataset is meaningless
fit$par[1:length(fit$par)] = pars

plot(phy, label.offset=0.2)
axisPhylo()
tiplabels(text=phy$tip.state, tip=1:5, col="black", bg=cols[phy$tip.state+1])
#diversitree:::plot.history(h, phy, main="True history", cols=cols)

#fit2 <- find.mle(lik2, pars, method="subplex", fail.value=-1e10)
st <- asr.marginal(lik, pars)
nodelabels(thermo=t(st), piecol=cols, cex=0.5)
t(st)

states_df = cft(adf2(t(st)))
states_df

nodelabels(thermo=t(st), piecol=cols, cex=0.5)


# Write out tree and data
write.tree(phy, file="")

# (sp5:9.201374725,(sp6:3.944281182,(sp7:1.688329525,(sp8:0.07452322119,sp9:0.07452322119)nd8:1.613806304)nd7:2.255951656)nd5:5.257093543)nd2;"

#  > phy$tip.state
#  sp5 sp6 sp7 sp8 sp9 
#   0   0   0   1   1 

# parameters: lambda, mu, q
0.22222222 0.22222222 0.11111111 0.01111111 0.06000000 0.05000000


