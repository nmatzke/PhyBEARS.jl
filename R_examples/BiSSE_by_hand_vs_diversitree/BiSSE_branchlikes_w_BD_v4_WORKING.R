
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

# Set your working directory to wherever you unzipped this directory
wd = "/GitHub/phyBEARS/ex/BiSSE_branchlikes_w_BD/"
setwd(wd)



library(ape)

par(mfrow=c(2,1))
tree_high_deathrate_fn = "tree_high_deathrate.newick"
trstr = "(((chimp:1,human:1):1,gorilla:2):3,orang:5);"
tree_high_deathrate = read.tree(file="", text=trstr)
plot(tree_high_deathrate)
axisPhylo()
title("Tree with high deathRate")


trfn = "tree.newick"
trstr = "(((chimp:1,human:1):1,gorilla:2):1,orang:3);"
tr = read.tree(file="", text=trstr)
plot(tr)
axisPhylo()
title("Tree with 0 deathRate")



trfn = "tree_small.newick"
tr = read.tree(file=trfn)
plot(tr)
axisPhylo()
title("Tree with 3 taxa")



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
source("/drives/Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp.R")
sourceall("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)




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
trfn = "tree_small.newick"
tr = read.tree(trfn)
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
# [1] 0.2222222 0.0000000

# Get the log-likelihood of the tree under the ML parameters
# Convert the deviance to likelihood
BD_LnL = -1 * BD$dev / 2
BD_LnL
# -3.216395


# Set birthRate
#birthRate = 0.22222222222222222222


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
# -4.063987

dev(a=x1, r=x2, N=N, x=x)
# -3.216395

dev(a=deathRate/birthRate, r=birthRate-deathRate, N=N, x=x)
# -3.216395

dev(a=0/birthRate, r=birthRate-0, N=N, x=x)

birthRate = 0.222222222222
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
# [1] 0.2222222
# 
# $se
# [1] 0.1571348
# 
# $loglik
# [1] -3.216395
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
# -3.216395 






#######################################################
#######################################################
# 3. Set up of equivalent BiSSE model for DEC, starting values
#######################################################
#######################################################

# ClaSSE helper functions
library(diversitree)
source("/GitHub/phyBEARS/R/ClaSSE_mods_v1.R")

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
states = c(1, 1, 1)
names(states) = tr$tip.label
states


# Proportion of species in each state; for 2 states
# (Let's assume we have all species)
sampling.f = c(1,1)

# Number of states
k = 2

# Make the BiSSE likelihood function. 
# (strict=FALSE means that some states in the state space can be absent from the tips)
bisse_2areas = make.bisse(tree=tr, states=states, sampling.f=sampling.f, strict=FALSE)
# Look at all the parameters of this model!
# lambdas = speciation rates
# mus = extinction rates
# qs = anagenetic transition rates
birthRate = 0.222222222
deathRate = 0.1

lambda0 = birthRate
lambda1 = birthRate
mu0 = deathRate
mu1 = deathRate
q01 = 0.1
q10 = 0.1
parms = c(lambda0, lambda1, mu0, mu1, q01, q10)
names(parms) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
parms
bisse_params = parms

res1t = bisse_2areas(pars=bisse_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
res2t = bisse_2areas(pars=bisse_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0, 1)
res3t = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.25, 0.75)
res4t = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
root_probs = c(0.5,0.5)
res5t = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)

LnLs1t = get_classe_LnLs(res1t)
LnLs2t = get_classe_LnLs(res2t)
LnLs3t = get_classe_LnLs(res3t)
LnLs4t = get_classe_LnLs(res4t)
LnLs5t = get_classe_LnLs(res5t)

LnLst = rbind(LnLs1t, LnLs2t, LnLs3t, LnLs4t, LnLs5t)
print(LnLst)




res1 = bisse_2areas(pars=bisse_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
res2 = bisse_2areas(pars=bisse_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0, 1)
res3 = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.25, 0.75)
res4 = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
root_probs = c(0.5,0.5)
res5 = bisse_2areas(pars=bisse_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)

LnLs1 = get_classe_LnLs(res1)
LnLs2 = get_classe_LnLs(res2)
LnLs3 = get_classe_LnLs(res3)
LnLs4 = get_classe_LnLs(res4)
LnLs5 = get_classe_LnLs(res5)

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


#######################################################
#######################################################
# 5. BiSSE likelihoods and comparison to yule
#######################################################
#######################################################
birthRate


# Speciation
lambda0 = birthRate
lambda1 = birthRate
# Extinction
mu0 = 0
mu1 = 0
# Character transition
q01 = 0 # ML
q10 = 0

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
# [1]  0.0000000  0.0000000  0.0000000 -0.2222222

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
condlikes_of_each_treeState_BRANCHTOP_AT_NODE[1,numstates] = 1	# chimp
condlikes_of_each_treeState_BRANCHTOP_AT_NODE[2,numstates] = 1	# human
condlikes_of_each_treeState_BRANCHTOP_AT_NODE[3,numstates] = 1	# gorilla
# condlikes_of_each_treeState_BRANCHTOP_AT_NODE[4,numstates] = 1	# orang if a different state
# condlikes_of_each_treeState_BRANCHTOP_AT_NODE[4,numstates] = 1	# orang same states for all tips

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

birthRate = 0.222222222
deathRate = 0.1

lambda0 = birthRate
lambda1 = birthRate
mu0 = deathRate
mu1 = deathRate
q01 = 0.1
q10 = 0.1

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
#            [,1]       [,2]        [,3]      [,4]
# [1,] 0.00000000 0.00000000 0.000000000 1.0000000
# [2,] 0.00000000 0.00000000 0.000000000 1.0000000
# [3,] 0.00000000 0.00000000 0.000000000 1.0000000
# [4,] 0.15069356 0.15069356 0.003615038 0.1672756
# [5,] 0.08603219 0.08603219 0.001825474 0.1837656
#
#            [,1]       [,2]        [,3]      [,4]
# [1,] 0.00000000 0.00000000 0.000000000 1.0000000
# [2,] 0.00000000 0.00000000 0.000000000 1.0000000
# [3,] 0.00000000 0.00000000 0.000000000 1.0000000
# [4,] 0.15026899 0.15026899 0.003582702 0.1674394
# [5,] 0.08566205 0.08566205 0.001808968 0.1839313


# Matches matches diversitree BiSSE "base" -- but a weird combination of 
# - left columns: average E's passed down
# - right columns: normalized conditional likelihoods
base
condlikesProbs_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
#            [,1]       [,2]       [,3]      [,4]
# [1,] 0.08603219 0.08603219 0.09063462 0.9093654
# [2,] 0.08603219 0.08603219 0.09063462 0.9093654
# [3,] 0.15069356 0.15069356 0.16483998 0.8351600
# [4,]         NA         NA         NA        NA
# [5,] 0.15069356 0.15069356 0.09868766 0.9013123
# 
#            [,1]       [,2]       [,3]      [,4]
# [1,] 0.08566205 0.08566205 0.08981508 0.9101849
# [2,] 0.08566205 0.08566205 0.09063463 0.9093654
# [3,] 0.15026899 0.15026899 0.16484005 0.8351599
# [4,] 0.00000000 0.00000000 0.00000000 0.0000000
# [5,] 0.15026899 0.15026899 0.09780486 0.9021951

# lq and likelihoods
condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
#            [,1]       [,2]       [,3]      [,4]
# [1,] 0.08603219 0.08603219 0.06700011 0.6722329
# [2,] 0.08603219 0.08603219 0.06700011 0.6722329
# [3,] 0.15069358 0.15069358 0.09311736 0.4717779
# [4,] 0.00000000 0.00000000 0.00000000 0.0000000
# [5,] 0.15069357 0.15069357 0.01399609 0.1278260

# Matches diversitree BiSSE "lq"
lnls = log(rowSums(condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[,3:4]))
lnls[!is.finite(lnls)] = NA
lnls
sum(lnls, na.rm=TRUE)
# [1] -0.2993006 -0.3021421 -0.5711150         NA -1.9499680
# [1] -3.122526

lq
sum(lq)
# [1] -0.3021421 -0.3021421 -0.5711149  0.0000000 -1.9531821
# [1] -3.128581

#######################################################
#######################################################
# END PROOF OF MATCHING
#######################################################
#######################################################








#######################################################
#######################################################
# Other notes/experiments - disregard mostly...
#######################################################
#######################################################

condlikes_of_each_treeState_BRANCHTOP_AT_NODE
condlikes_of_each_treeState_x_lambda_BRANCHTOP_AT_NODE

condlikes_treeStates_BRANCHTOP_AT_NODE_DOWNPASS
condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
relative_probs_of_each_state_BRANCHTOP_AT_NODE

relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS

init
base


log(rowSums(init[,3:4]))
sum(log(rowSums(init[,3:4])))
# [1] -3.450941		# NO

lq
sum(lq)
# -0.3021421 -0.3021421 -0.5711149  0.0000000 -1.9531821
# -3.128581
# YES

log(rowSums(condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[,3:4]))



log(computed_likelihoods_at_each_node)
# 0.000  0.000  0.000  0.000 -0.888 -0.666 -0.444
sum(log(computed_likelihoods_at_each_node))  # corresponds to -lambda * sum(branchlengths)
# -1.998

# Log-likelihood
sum(log(rowSums(condlikes_of_each_treeState_BRANCHTOP_AT_NODE)))  # corresponds to -lambda * sum(branchlengths)

sum(log(rowSums(condlikes_of_each_treeState_x_lambda_BRANCHTOP_AT_NODE)))  # corresponds to -lambda * sum(branchlengths)


nb_node = tr$Nnode - 1
lambda = birthRate
sum(log(computed_likelihoods_at_each_node)) + lfactorial(tr$Nnode) + (nb_node * log(lambda))
# -3.214395


# Compare to Yule:
# -3.216395 



times = seq(0,1,1/50)
top_y = init[5,]
names(top_y) = names_y
out <- lsoda(y=top_y, times=times, func=define_BiSSE_eqns_in_R, parms=parms) 
tail(out)
brbot_likes = out[nrow(out),]
log(sum(brbot_likes[4:5]))
brbot_likes[4:5] / sum(brbot_likes[4:5])







times = seq(0,1,1/50)
#y = c(0, 0, 0, 1)
out <- lsoda(y=y, times=times, func=define_BiSSE_eqns_in_R, parms=parms) 
tail(out)

out <- ode(y=y, times=times, func=define_BiSSE_eqns_in_R, parms=parms) 
tail(out)

out <- ode(y=y, times=times, func=define_BiSSE_eqns_in_R, parms=parms, method="bdf") 
tail(out)


# EVERYTHING MATCHES


likes_of_treeStates_above = out[nrow(out), 4:5]
likes_of_treeStates_above
normalized_likes_of_treeStates_above = likes_of_treeStates_above / sum(likes_of_treeStates_above)
normalized_likes_of_treeStates_above

init
base

# Log-likelihood of data at the base of the tip-branches of length 1
log(sum(likes_of_treeStates_above))

# branchlengths + states log-likelihood
attr(res2, "intermediates")$lq

# Total (branchlengths + states) log-likelihood
attr(res2, "intermediates")$lq
log(rowSums(condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[,3:4]))


condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS


condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
log(rowSums(condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[,3:4]))
sum(.Last.value)

relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS

likes_of_treeStates_above * likes_of_treeStates_above
normalized_likes_of_treeStates_above * normalized_likes_of_treeStates_above

likes_of_treeStates_above * birthRate
normalized_likes_of_treeStates_above * birthRate







base[1,3:4]
# [1] 0.09063462 0.90936538
base[1,3:4]^2
# [1] 0.008214635 0.826945386
base[1,3:4] / init[5,3:4]
# [1] 49.649899  4.948506
(base[1,3:4] / init[5,3:4])^0.5
# [1] 7.046268 2.224524
(base[1,3:4] / init[5,3:4])^0.222
# [1] 2.379545 1.426176
(base[1,3:4] / init[5,3:4])^(1/0.222222)
# [1] 42819373.869     1333.941
(base[1,3:4]^2 / init[5,3:4])
# 4.5 4.5
log(1/4.5)
# [1] -1.504077
1/0.2222222
# [1] 4.5
node_likelihood_x_lambda
#          D0t          D1t 
# 0.0006872412 0.0321185677 



init
base
condlikes_of_each_treeState_BRANCHTOP_AT_NODE
condlikes_of_each_treeState_x_lambda_BRANCHTOP_AT_NODE







colSums(log(mat[3:4,]))


condlikes_treeStates_BRANCHTOP_AT_NODE_DOWNPASS
condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS

condlikes_of_each_treeState_BRANCHTOP_AT_NODE
relative_probs_of_each_state_BRANCHTOP_AT_NODE
attr(res2, "intermediates")$lq





mat = attr(res1, "intermediates")$base
class(mat)

t(condlikes_treeStates_BRANCHTOP_AT_NODE_DOWNPASS)
t(condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS)





#sum(log(computed_likelihoods_at_each_node))
log(computed_likelihoods_at_each_node_x_lambda)
# 0.000000  0.000000  0.000000  0.000000 -2.392077 -2.170077 -1.948077
sum(log(computed_likelihoods_at_each_node_x_lambda))
# -6.510232

sum(log(computed_likelihoods_at_each_node_x_lambda)) + lfactorial(tr$Nnode)
# -4.718473

root.p = c(0.5, 0.5)
rootnode = length(tr$tip.label)+1
root_likes = relative_probs_of_each_state_BRANCHTOP_AT_NODE[rootnode,] * root.p
root_likes
log(sum(root_likes))
# -0.6931472

sum(log(computed_likelihoods_at_each_node_x_lambda))
sum(log(computed_likelihoods_at_each_node_x_lambda)) + log(sum(root_likes))
# -7.203379


