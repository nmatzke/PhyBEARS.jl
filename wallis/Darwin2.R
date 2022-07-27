
#######################################################
# 
# This script demonstrates how BiSSE is calculated
# 
#######################################################



#######################################################
# OUTLINE
# 
# 1. Calculation of tip state likelihoods under DEC, in 
#    BioGeoBEARS
#
# 2. Calculation of tree likelihood under a Yule process
# 
# 3. Set up of equivalent BiSSE model
#
# 4. BiSSE likelihoods and comparison to DEC
# 
#######################################################

library(ape)

nexfn = "https://www.r-phylo.org/w/images/0/02/Geospiza.nex"
tr = read.nexus(file=nexfn)
# plot(tr)
# axisPhylo()
# title("Tree with 3 taxa")



#######################################################
#######################################################
# 1. Calculation of tip state likelihoods under DEC, in 
#    BioGeoBEARS
#######################################################
#######################################################
library(diversitree)
library(deSolve)		# for lsoda

# After installation, load the package, dependencies, updates
library(optimx)
library(FD)
library(snow)
library(parallel)
library(BioGeoBEARS)



#######################################################
# BiSSE directly, using lsoda
#######################################################

#######################################################
# Example lsoda run
#######################################################

# Define the function

# INPUTS
# t = current time point in integration
# y = current estimate of the variables in the system
#     (if names(y) exists, those variables will be available)
# parms = parameters (named)
#
# OUTPUTS:
# The return value of func should be a list, whose first element 
# is a vector containing the derivatives of y with respect to time, 
# and whose next elements are global values that are required at each 
# point in times. The derivatives must be specified in the same order 
# as the state variables y.
# 
SPCmod <- function(t, x, parms)
{
  with(data=
         as.list(c(parms, x)), 
       { # expr
         import <- sigimp(t)
         dS <- import - b*S*P + g*C     # substrate
         dP <- c*S*P  - d*C*P           # producer
         dC <- e*P*C  - f*C             # consumer
         res <- c(dS, dP, dC)
         list(res)
       }
  )
}

## Parameters 
parms  <- c(b = 0.0, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0.0)

## vector of timesteps
times  <- seq(0, 100, length = 101)

## external signal with rectangle impulse
signal <- as.data.frame(list(times = times,
                             import = rep(0,length(times))))

signal$import[signal$times >= 10 & signal$times <= 11] <- 0.2

sigimp <- approxfun(signal$times, signal$import, rule = 2)


## Start values for steady state
y <- xstart <- c(S = 1, P = 1, C = 1)

## Solving
out <-  lsoda(xstart, times, SPCmod, parms) 






#######################################################
#######################################################
# 2. Calculation of tree likelihood under a BirthDeath process
#######################################################
#######################################################
library(ape)
nexfn = "https://www.r-phylo.org/w/images/0/02/Geospiza.nex"
tr = read.nexus(file=nexfn)
plot(tr)
axisPhylo()

# Get the ML estimates of birthRate (speciation rate) 
# and deathRate (extinction rate)
# ...using ape::birthdeath
BD =  birthdeath(tr)
BD
names(BD)

# Calculate the birthRate and deathRate from the outputs
x1 = unname(BD$para["d/b"])
x2 = unname(BD$para["b-d"])
deathRate = (x2*x1) / (1-x1)
birthRate = deathRate+x2
c(birthRate, deathRate)

# You should get:
# c(birthRate, deathRate)
# [1] 3.682184 2.263549

# Get the log-likelihood of the tree under the ML parameters
# Convert the deviance to likelihood
BD_LnL = -1 * BD$dev / 2
BD_LnL
# -3.216395


# Set birthRate
#birthRate = 0.3.682184


#######################################################
# Likelihood equation in the birthdeath function
#######################################################
N <- length(tr$tip.label)
nb_node = tr$Nnode - 1
sptimes <- c(NA, branching.times(tr)) # NA so the number of times equals number of tips?
# a = "d/b"
# r = "b-d"

dev <- function(a=0.1, r=0.2, N, x, return_deviance=FALSE)
{
  if (r < 0 || a > 1) 
    return(1e+100)
  
  lnl_topology = lfactorial(tr$Nnode)
  lnl_numBirths = nb_node * log(r)
  lnl_Births_above_root = r * sum(sptimes[3:N])
  
  lnl_numtips_wOneMinusDeathRate = N * log(1 - a)
  # Interpretation: more tips are less likely, if relativeDeathRate is >0
  # If relativeDeathRate = 1, a=0, and lnl=-Inf... 
  #    CLEARLY WRONG EXCEPT IN A MAXIMUM LIKELIHOOD CONTEXT!!!
  # If relativeDeathRate = 0, a=0, and lnl=0, i.e. any number of tips is equiprobable
  
  lnl_branching_times = -2 * sum(log(exp(r * sptimes[2:N]) - a))
  # For each observed branching,
  # netDiversificationRate * timeOfEvent <- take exponent of that ; this means recorded events are less likely in the past
  # <- subtract "a", a constant (relativeDeathRate)
  #
  # This would be a straight likelihood as:
  # 1/
  # (exp(r*branching_time)-a)^2
  #
  # Sum the logs of these
  #
  # If deathRate = 0
  # lnl_branching_times = -2 * sum(log(exp(birthRate*sptimes[2:N]) - 0))
  # lnl_branching_times = -2 * sum(log( exp(birthRate*sptimes[2:N]) )
  # lnl_branching_times = -2 * sum( birthRate*sptimes[2:N] )
  #
  # Note: sum(X) = 9 = total branchlength of tr
  # In BD:
  # -2*sum(sptimes[2:N]) = -12
  # sum(sptimes[3:N]) = 3
  # So, lnl_branching_times + lnl_Births_above_root = yule's -lambda * X
  
  
  
  if (return_deviance == TRUE)
  {
    result = -2 * (lnl_topology + lnl_numBirths + lnl_Births_above_root + 
                     lnl_numtips_wOneMinusDeathRate + lnl_branching_times)
  } else {
    result = lnl_topology + lnl_numBirths + lnl_Births_above_root + 
      lnl_numtips_wOneMinusDeathRate + lnl_branching_times
  }
  return(result)
}
dev(a=0.1, r=0.2, N=N, x=x)
# 3.430055

dev(a=x1, r=x2, N=N, x=x)
# 22.36105

dev(a=deathRate/birthRate, r=birthRate-deathRate, N=N, x=x)
# 22.36105

dev(a=0/birthRate, r=birthRate-0, N=N, x=x)

birthRate = 3.682184
dev(a=0/birthRate, r=birthRate-0, N=N, x=x)


S <- rep(1, times=length(tr$tip.label))
bd.ext(phy=tr, S=S, conditional=TRUE)
#       d / b = 2.883955e-07   StdErr = NaN 
#       b - d = 0.2656457   StdErr = NaN 

bd.ext(phy=tr, S=S, conditional=FALSE) # same than older versions
#      d / b = 4.949791e-06   StdErr = NaN 
#      b - d = 0.242636   StdErr = NaN 




#######################################################
#######################################################
# 3. Calculation of tree likelihood under a Yule process
#######################################################
#######################################################
yule(phy=tr)
# $lambda
# [1] 2.61642
# 
# $se
# [1] 0.7552953
# 
# $loglik
# [1] 22.09385
# 
# attr(,"class")
# [1] "yule"


#######################################################
# Yule likelihood
#######################################################
yule_lik <- function(tr)
{
  # Total length of tree
  X <- sum(tr$edge.length)
  
  # Number of internal nodes
  tr$Nnode
  
  # Number of internal nodes, minus root
  nb_node <- tr$Nnode - 1
  
  # ML estimate of lambda (birthrate)
  lambda <- nb_node/X
  
  # Standard error of ML estimate
  se <- lambda/sqrt(nb_node)
  
  # Log-likelihood
  lnl_topology = lfactorial(tr$Nnode)
  lnl_numBirths = nb_node * log(lambda)
  loglik = -lambda * X + lnl_topology + lnl_numBirths
  loglik
  # -3.216395
  
  res = NULL
  res$lambda = lambda
  res$se = se
  res$loglik = loglik
  
  return(res)
}

yule_lik(tr=tr)
# 22.09385






#######################################################
#######################################################
# 3. Set up of equivalent BiSSE model for DEC, starting values
#######################################################
#######################################################

# ClaSSE helper functions
library(diversitree)
source("/GitHub/PhyBEARS.jl/Rsrc/ClaSSE_mods_v2.R")
source("/GitHub/PhyBEARS.jl/Rsrc/ClaSSE_functions_v3.R")


##################################################
# Set up states for BiSSE to match DEC states
##################################################

# States in the BiSSE model, and how they correspond
# to the geographic ranges in DEC
# 
# ====================
# (statenum = range)
# ====================
# 1 = null range
# 2 = A = Africa
# 3 = B = Asia
# 4 = AB = both
# ====================

# States at the tips of the tree
# (ranges A, A, A, B)
#states = c(1, 1, 1, 1)
states = rep(1, times=length(tr$tip.label))  # 14 binary states for 14 tips
names(states) = tr$tip.label
states


# Proportion of species in each state; for 2 states
# (Let's assume we have all species)
sampling.f = c(1,1)  # length(sampling.f) = 2 states

# Number of states
k = 2

# Make the BiSSE likelihood function. 
# (strict=FALSE means that some states in the state space can be absent from the tips)
bisse_2areas = diversitree::make.bisse(tree=tr, states=states, sampling.f=sampling.f, strict=FALSE)

# Look at all the parameters of this model!
# lambdas = speciation rates
# mus = extinction rates
# qs = anagenetic transition rates
birthRate = 3.682184
deathRate = 2.263549

lambda0 = birthRate
lambda1 = birthRate
mu0 = deathRate
mu1 = deathRate
q01 = 0.0
q10 = 0.0
parms = c(lambda0, lambda1, mu0, mu1, q01, q10)
names(parms) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
parms
bisse_params = parms

# Constraint parameters so you are only fitting 1 birthRate
#constraints = list(lambda0~lambda1, mu0~0.0, mu1~0.0, q01~0.0, q10~0.0)
#bisse_2areas_constrained = diversitree::constrain(f=bisse_2areas, formulae=constraints)

# Wait 1 seconds
#fit <- find.mle(func=bisse_2areas_constrained, x.init=bisse_params, method="subplex")
#fit$par.full
#fit$lnLik

# Compare to Yule
yule(tr)

bisse_params_orig = bisse_params
bisse_params = fit$par.full

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


init = t(attr(res2, "intermediates")$init)
init

base = t(attr(res2, "intermediates")$base)
base
# Columns 3-4 sum to 1
rowSums(base[,3:4])

lq = attr(res2, "intermediates")$lq
lq
sum(lq)


#######################################################
#######################################################
# 5. BiSSE likelihoods and comparison to yule
#######################################################
#######################################################
birthRate


# Speciation
lambda0 = bisse_params["lambda0"]
lambda1 = bisse_params["lambda1"]
# Extinction
mu0 = deathRate
mu1 = deathRate
# Character transition
q01 = 0.0 # ML
q10 = 0.0

parms = c(lambda0, lambda1, mu0, mu1, q01, q10)
names(parms) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
parms

bisse_params = parms

# Starting values of state variables
E0t = 0
E1t = 0
D0t = 0
D1t = 1
y = c(E0t, E1t, D0t, D1t)
names(y) = c("E0t", "E1t", "D0t", "D1t")
y

# t = current time point of integration
# y = state variable we are tracking (named)  MUST HAVE NAMES!!!
# parms = model parameters (named)            MUST HAVE NAMES!!!

define_BiSSE_eqns_in_R <- function(t, y, parms)
{
  with(data=
         as.list(c(y, parms)), 
       { # expr
         # When the limit is taken as deltaT goes to 0, the
         # change in E0 in dt is:
         # probs of: 
         # - extinction
         # - no change but later extinction
         # - state change, then extinction
         # - no change, speciation, extinction of both
         dE0t <- mu0 - (mu0 + q01 + lambda0)*E0t + q01*E1t + lambda0*(E0t)^2
         dE1t <- mu1 - (mu1 + q10 + lambda1)*E1t + q10*E0t + lambda1*(E1t)^2
         
         # probs of:
         # - no change
         # - character change, no speciation
         # - speciation followed by extinction of 1
         # - extinction (prob of observed clade = 0, since the clade is extant)
         dD0t <- -1*(lambda0 + mu0 + q01)*D0t + q01*D1t + 2*lambda0*E0t*D0t + 0
         dD1t <- -1*(lambda1 + mu1 + q10)*D1t + q10*D0t + 2*lambda1*E1t*D1t + 0
         
         # Return the list of coupled differential equations
         res <- c(dE0t, dE1t, dD0t, dD1t)
         return(list(res))
         #return(list_of_diff_eqns)
       }
  )
}

# One step test:
one_step_result = define_BiSSE_eqns_in_R(t=t, y=y, parms=parms)
one_step_result
# [[1]]
# [1]   0.000000  0.000000  0.000000 -2.616461

# LSODA inputs:
# y = initial state values
# times = times at which you want estimates
# func
times = seq(0,1,1/50)
out <- lsoda(y=y, times=times, func=define_BiSSE_eqns_in_R, parms=parms) 
out

parms2 = parms
parms2["mu0"] = 0.1
parms2["mu1"] = 0.1
times = seq(0,1,1/50)
out <- lsoda(y=y, times=times, func=define_BiSSE_eqns_in_R, parms=parms2) 
out





# Downpass
numsteps = 100000
numstates = 2
num_internal_nodes = tr$Nnode
numtips = length(tr$tip.label)
num_internal_nodes = tr$Nnode
numnodes = numtips + num_internal_nodes
tipnums <- 1:numtips


# Reorder the edge matrix into pruningwise order
# This is CRUCIAL!!
tr2 <- reorder(tr, "pruningwise")
edgelengths = tr2$edge.length


# Define matrices to store data
# We have Es for each state, and Ds for each state
# But we only need to record the Ds
condlikes_of_each_treeState_BRANCHTOP_AT_NODE = matrix(data=0, nrow=numnodes, ncol=numstates)

# Fill in the likelihoods of tip nodes manually
#tip_states_Ds = y[c(((length(y)/2)+1), length(y))]
fakedata = rep(1, times=numnodes)

#fakedata = rep(2, times=length(phy$tip.label))
#fakedata[1:5] = 1

for (i in 1:length(fakedata))
{
  condlikes_of_each_treeState_BRANCHTOP_AT_NODE[i,1] =1
}

zeros = matrix(data=0, nrow=numnodes, ncol=numstates)
condlikes_of_each_treeState_BRANCHTOP_AT_NODE = cbind(zeros, condlikes_of_each_treeState_BRANCHTOP_AT_NODE)
condlikes_of_each_treeState_x_lambda_BRANCHTOP_AT_NODE = condlikes_of_each_treeState_BRANCHTOP_AT_NODE

relative_probs_of_each_state_BRANCHTOP_AT_NODE = condlikes_of_each_treeState_BRANCHTOP_AT_NODE[,3:4]
relative_probs_of_each_state_BRANCHTOP_AT_NODE[1:numtips,] = relative_probs_of_each_state_BRANCHTOP_AT_NODE[1:numtips,] / rowSums(relative_probs_of_each_state_BRANCHTOP_AT_NODE[1:numtips,])

relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = matrix(data=0, nrow=numnodes, ncol=numstates)

# Store both the Es and the Ds likelihoods
condlikes_treeStates_BRANCHTOP_AT_NODE_DOWNPASS = condlikes_of_each_treeState_BRANCHTOP_AT_NODE
condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS = matrix(data=0, nrow=numnodes, ncol=numstates*2)
condlikesProbs_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS = matrix(data=0, nrow=numnodes, ncol=numstates*2)

computed_likelihoods_at_each_node = rep(0, numnodes)
computed_likelihoods_at_each_node[1:numtips] = 1
computed_likelihoods_at_each_node_x_lambda = computed_likelihoods_at_each_node



# DEFINE DOWNPASS THROUGH THE BRANCHES	
i = 1
edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)

# Speciation
lambda0 = bisse_params["lambda0"]
lambda1 = bisse_params["lambda1"]
# Extinction
mu0 = deathRate
mu1 = deathRate
# Character transition
q01 = 0.0 # ML
q10 = 0.0


parms = c(lambda0, lambda1, mu0, mu1, q01, q10)
names(parms) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
parms
numsteps=100

names_y = c("E0t", "E1t", "D0t", "D1t")

i = 1
for (i in edges_to_visit)
{
  # First edge visited is i
  #print(i)
  
  # Its sister is j 
  j <- i + 1
  
  # Get the node numbers at the tips of these two edges		
  left_desc_nodenum <- tr2$edge[i, 2]
  right_desc_nodenum <- tr2$edge[j, 2]
  
  left_edgenum_TF = tr2$edge[, 2] == left_desc_nodenum
  right_edgenum_TF = tr2$edge[, 2] == right_desc_nodenum
  left_edgenum = (1:length(edgelengths))[left_edgenum_TF]
  right_edgenum = (1:length(edgelengths))[right_edgenum_TF]
  
  
  
  # And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
  anc <- tr2$edge[i, 1]
  
  
  # Calculate the downpass on two branches
  # The Es are 0 (no extinctions above, already accounted for)
  
  # Left branch
  # The Ds are relative state probabilities
  edgelength_Left = edgelengths[left_edgenum]
  times = seq(from=0, to=edgelength_Left, by=edgelength_Left/numsteps)
  y = condlikes_treeStates_BRANCHTOP_AT_NODE_DOWNPASS[left_desc_nodenum,]
  names(y) = names_y
  
  out_matrixL = lsoda(y=y, times=times, func=define_BiSSE_eqns_in_R, parms=parms) 
  condlikes_tempLeft = out_matrixL[nrow(out_matrixL),][2:5]
  condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[left_desc_nodenum,] = condlikes_tempLeft
  condExtinct_Left = condlikes_tempLeft[1:2]
  condlikes_Left = condlikes_tempLeft[3:4]
  
  # 
  # R, 100 steps, mu=0
  # 100 1.00   0   0 0.7687171 0.0338017454
  
  # HiSSE, 2 steps, tiny mu
  # 2    1 1.992624e-07 1.992624e-07 0.7666849 0.03405248
  
  # HiSSE, 100 steps
  # 100 1.00 1.974810e-07 1.974810e-07 0.7687171 0.0338017546
  
  # Right branch
  edgelength_Right = edgelengths[right_edgenum]
  times = seq(from=0, to=edgelength_Right, by=edgelength_Right/numsteps)
  # The Ds are relative state probabilities
  y = condlikes_treeStates_BRANCHTOP_AT_NODE_DOWNPASS[right_desc_nodenum,]
  names(y) = names_y
  
  out_matrixR = lsoda(y=y, times=times, func=define_BiSSE_eqns_in_R, parms=parms) 
  condlikes_tempRight = out_matrixR[nrow(out_matrixR),][2:5]
  condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[right_desc_nodenum,] = condlikes_tempRight
  condExtinct_Right = condlikes_tempRight[1:2]
  condlikes_Right = condlikes_tempRight[3:4]
  
  
  # Conditional likelihoods of states at the bottom of right branch
  #condlikes_Right = independent_likelihoods_on_each_branch[[j]] %*% relative_probs_of_each_state_BRANCHTOP_AT_NODE[right_desc_nodenum,]
  
  
  
  # Every node (except maybe the root) has a branch below it, and there is also a 
  # relative_probs_of_each_state_BRANCHTOP_AT_NODE at the bottom of this branch
  condprobs_Left = condlikes_Left / sum(condlikes_Left)
  condprobs_Right = condlikes_Right / sum(condlikes_Right)
  relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_desc_nodenum,] = condprobs_Left
  relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_desc_nodenum,] = condprobs_Right
  
  
  
  # If there is no speciational model, you are assuming 100% sympatry (range duplication)
  # at each speciation event
  #
  # In this case, you can just multiply the two conditional likelihood matrices together
  #
  # Also, if a branch is extremely short (a "hook"), this is essentially a zero-length
  # branch, we are assuming that this represents the range of a lineage at that 
  # point.  There is no speciation event here -- both "lineages" inherit
  # the same range.  This allows fossils to closely influence ancestral states.
  #
  # This was developed with Kaitlin Maguire over several years of screwing around.
  
  # Check for a short "hook" branch; if found, use just allopatric speciational model
  
  # get the correct edge
  left_edge_TF = tr2$edge[,2] == left_desc_nodenum
  right_edge_TF = tr2$edge[,2] == right_desc_nodenum
  
  node_likelihood = condlikes_Left * condlikes_Right
  lambda0 = parms["lambda0"]
  lambda1 = parms["lambda1"]
  node_likelihood_x_lambda = node_likelihood * c(lambda0, lambda1)
  
  E_average = colMeans(rbind(condExtinct_Left, condExtinct_Right))
  
  # birthRate*base[1,3:4]^2
  D_relprobs_multiplied = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_desc_nodenum,] * relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_desc_nodenum,]
  D_combined = c(lambda0, lambda1) * D_relprobs_multiplied
  
  # Store the various options
  condlikes_of_each_treeState_BRANCHTOP_AT_NODE[anc,] = c(E_average, node_likelihood)
  condlikes_of_each_treeState_x_lambda_BRANCHTOP_AT_NODE[anc,] = c(E_average, node_likelihood_x_lambda)
  condlikes_treeStates_BRANCHTOP_AT_NODE_DOWNPASS[anc, ] = c(E_average, D_combined)
  condlikesProbs_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[left_desc_nodenum,] = c(E_average, condprobs_Left)
  condlikesProbs_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[right_desc_nodenum,] = c(E_average, condprobs_Right)
  
  # node_likelihood_x_lambda
  #          D0t          D1t 
  # 0.1313168545 0.0002539018 
  
  
  # 	> v*c(0.333333, 0.333333)
  # 			  D0           D1 
  # 	0.1306233807 0.0002576823 
  # 	> v
  # 			 D0          D1 
  # 	0.587805801 0.001159572 
  # 	>  sequence(length(desRows))
  # 	[1] 1 2
  # 	> compD[focal,]
  # 	[1] 0.1306233895 0.0002576823
  
  
  total_likelihood_for_node = sum(node_likelihood)
  total_likelihood_for_node_x_lambda = sum(node_likelihood_x_lambda)
  
  computed_likelihoods_at_each_node[anc] = total_likelihood_for_node
  computed_likelihoods_at_each_node_x_lambda[anc] = sum(total_likelihood_for_node_x_lambda)
  
  #print(total_likelihood_for_node)
  relative_probs_of_each_state_BRANCHTOP_AT_NODE[anc, ] = node_likelihood_x_lambda / total_likelihood_for_node_x_lambda
} # END for (i in edges_to_visit)


#######################################################
#######################################################
# START PROOF OF MATCHING THIS CODE TO DIVERSITREE
#######################################################
#######################################################
# Best, matches diversitree BiSSE "init"
init
condlikes_treeStates_BRANCHTOP_AT_NODE_DOWNPASS
# [,1] [,2]     [,3] [,4]
# [1,]    0    0 1.000000    0
# [2,]    0    0 1.000000    0
# [3,]    0    0 1.000000    0
# [4,]    0    0 1.000000    0
# [5,]    0    0 1.000000    0
# [6,]    0    0 1.000000    0
# [7,]    0    0 1.000000    0
# [8,]    0    0 1.000000    0
# [9,]    0    0 1.000000    0
# [10,]    0    0 1.000000    0
# [11,]    0    0 1.000000    0
# [12,]    0    0 1.000000    0
# [13,]    0    0 1.000000    0
# [14,]    0    0 1.000000    0
# [15,]    0    0 2.616461    0
# [16,]    0    0 2.616461    0
# [17,]    0    0 2.616461    0
# [18,]    0    0 2.616461    0
# [19,]    0    0 2.616461    0
# [20,]    0    0 2.616461    0
# [21,]    0    0 2.616461    0
# [22,]    0    0 2.616461    0
# [23,]    0    0 2.616461    0
# [24,]    0    0 2.616461    0
# [25,]    0    0 2.616461    0
# [26,]    0    0 2.616461    0
# [27,]    0    0 2.616461    0


# Matches matches diversitree BiSSE "base" -- but a weird combination of 
# - left columns: average E's passed down
# - right columns: normalized conditional likelihoods
base
condlikesProbs_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
# [,1] [,2] [,3] [,4]
# [1,]    0    0    1    0
# [2,]    0    0    1    0
# [3,]    0    0    1    0
# [4,]    0    0    1    0
# [5,]    0    0    1    0
# [6,]    0    0    1    0
# [7,]    0    0    1    0
# [8,]    0    0    1    0
# [9,]    0    0    1    0
# [10,]    0    0    1    0
# [11,]    0    0    1    0
# [12,]    0    0    1    0
# [13,]    0    0    1    0
# [14,]    0    0    1    0
# [15,]    0    0    0    0
# [16,]    0    0    1    0
# [17,]    0    0    1    0
# [18,]    0    0    1    0
# [19,]    0    0    1    0
# [20,]    0    0    1    0
# [21,]    0    0    1    0
# [22,]    0    0    1    0
# [23,]    0    0    1    0
# [24,]    0    0    1    0
# [25,]    0    0    1    0
# [26,]    0    0    1    0
# [27,]    0    0    1    0

# lq and likelihoods
base[,3:4] * exp(lq)
condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
# [,1] [,2]       [,3] [,4]
# [1,]    0    0 0.86596973    0
# [2,]    0    0 0.86596973    0
# [3,]    0    0 0.74990357    0
# [4,]    0    0 0.61898408    0
# [5,]    0    0 0.60430957    0
# [6,]    0    0 0.55070663    0
# [7,]    0    0 0.79710518    0
# [8,]    0    0 0.94901639    0
# [9,]    0    0 0.94901639    0
# [10,]    0    0 0.91249186    0
# [11,]    0    0 0.29583224    0
# [12,]    0    0 0.24723253    0
# [13,]    0    0 0.21734689    0
# [14,]    0    0 0.09980853    0
# [15,]    0    0 0.00000000    0
# [16,]    0    0 1.20151390    0
# [17,]    0    0 2.30018128    0
# [18,]    0    0 2.18662536    0
# [19,]    0    0 1.84247924    0
# [20,]    0    0 1.99595662    0
# [21,]    0    0 2.38437755    0
# [22,]    0    0 2.55443103    0
# [23,]    0    0 2.15967425    0
# [24,]    0    0 2.26577559    0
# [25,]    0    0 1.37897302    0
# [26,]    0    0 2.28560312    0
# [27,]    0    0 2.51576151    0

# Matches diversitree BiSSE "lq"
lnls = log(rowSums(condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[,3:4]))
lnls[!is.finite(lnls)] = NA
lnls
sum(lnls, na.rm=TRUE)
# [1] -0.14390533 -0.14390533 -0.28781066 -0.47967573 -0.50366867 -0.59655304 -0.22676863
# [8] -0.05232921 -0.05232921 -0.09157612 -1.21796274 -1.39742595 -1.52626063 -2.30450158
# [15]          NA  0.18358234  0.83298794  0.78235942  0.61111208  0.69112344  0.86893810
# [22]  0.93782951  0.76995740  0.81791712  0.32133903  0.82662994  0.92257554
sum(lnls, na.rm=TRUE)
# -0.458321


# Sum of the branch likelihoods
lq
sum(lq)
# [1] -0.14390533 -0.14390533 -0.28781066 -0.47967570 -0.50366865 -0.59655300 -0.22676863
# [8] -0.05232921 -0.05232921 -0.09157612 -1.21796237 -1.39742540 -1.52625991 -2.30449993
# [15]  0.00000000  0.18358244  0.83298794  0.78235943  0.61111209  0.69112345  0.86893810
# [22]  0.93782951  0.76995740  0.81791712  0.32133908  0.82662994  0.92257554
# -0.4583174

# Add the root probabilities
# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)

LnLs1
LnLs1t
#[1] -0.4075469 -1.7110530
LnLs1t
#[1] -0.191117 -1.711053

# Does the total of branch likelihoods (lq) + node likelihoods match R?
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = sum(log(computed_likelihoods_at_each_node_x_lambda))
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
# 12.04537


R_result_branch_lnL = -1.711053
R_result_total_LnLs1 = -0.4075469
R_result_total_LnLs1t = -0.191117
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = 15.23452


#######################################################
#######################################################
# END PROOF OF MATCHING
#######################################################
#######################################################



