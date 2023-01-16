

library(diversitree)

## Start with a simple tree evolved under a BiSSE with all rates
## asymmetric:
# lambda0, lambda1, mu0, mu1, q01, q10

# Gives weird result -- green nodes, blue truth
set.seed(47) # Rare to have a transition, but observed here
pars <- c(0.222222222, 0.222222222, 0.111111111, 0.222222222, 0.1, 0.05)

# Seems accurate
pars <- c(0.222222222, 0.222222222, 0.111111111, 0.05, 0.1, 0.15)

set.seed(48) # Rare to have a transition, but observed here
phy <- trees(pars, "bisse", max.taxa=4, max.t=Inf, x0=0)[[1]]
cols = c("blue", "green3")


## Here is the true history
h <- history.from.sim.discrete(phy, 0:1)
diversitree:::plot.history(h, phy, main="True history", cols=cols)

lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars, method="subplex")
st <- asr.marginal(lik, coef(fit))
nodelabels(thermo=t(st), piecol=cols, cex=0.5)




# change to create some ambiguity

phy2 = phy
phy2$edge.length[c(1,4)] = c(1-phy2$edge.length[2], 1-phy2$edge.length[6])
phy2$tip.state[1:4] = c(1, 1, 1, 0)

plot(phy2, label.offset=0.1)
tiplabels(text=phy2$tip.state, tip=1:4, col="black", bg=cols[phy2$tip.state+1])
#diversitree:::plot.history(h, phy, main="True history", cols=cols)

lik2 <- make.bisse(phy2, phy2$tip.state)
fit2 <- find.mle(lik2, pars, method="subplex")
st2 <- asr.marginal(lik2, coef(fit2))
nodelabels(thermo=t(st2), piecol=cols, cex=0.5)
