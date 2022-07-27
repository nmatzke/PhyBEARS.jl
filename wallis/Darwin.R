
#Example likelihood calculation on a single branch
#(binary discrete trait)
library(ape)
library(expm)

# NM2: load the "Darwin's Finches" tree (genus Geospiza)
nexfn = "https://www.r-phylo.org/w/images/0/02/Geospiza.nex"
phy = read.nexus(file=nexfn)
phy

# NM2: create some fake data
fakedata = rep(2, times=length(phy$tip.label))
fakedata[1:5] = 1


# Specifying the parameters
p12 = 0.1
p21 = 0.01
param_vector = c(p12, p21)

discrete_lnL_2x2 <- function(param_vector, fakedata, phy)
{

  # Specifying the Q transition matrix
  Qmat = matrix(data=0, nrow=2, ncol=2)
  Qmat

  Qmat[1,2] = param_vector[1]
  Qmat[2,1] = param_vector[2]

  # (the diagonals should be -1 * the sum of the rows)
  diag(Qmat) = -1 * rowSums(Qmat)
  Qmat

  # NM2: Store the likelihoods during the downpass
  number_of_discrete_states = 2
  number_of_tip_nodes = length(phy$tip.label)
  number_of_internal_nodes = phy$Nnode
  total_number_of_nodes = number_of_tip_nodes + number_of_internal_nodes

  TIPS = 1:number_of_tip_nodes

  # Table storing the likelihoods
  liks = matrix(data=0, nrow=total_number_of_nodes, ncol=number_of_discrete_states)

  # comp stores the likelihood computation results
  comp = rep(0, times=total_number_of_nodes)


  # Load the tip data into the liks table
  for (i in 1:length(fakedata))
  {
    liks[i,fakedata[i]] =1
  }
  liks

  # reorder edge matrix
  phy = reorder(phy, "postorder")

  # Load the edge matrix
  e1 = phy$edge[,1] # ancestor nodes
  e2 = phy$edge[,2] # descendant nodes
  EL = phy$edge.length # edge lengths

  for (i in seq(from=1, by=2, length.out=number_of_internal_nodes))
  {
    j <- i + 1L
    anc <- e1[i]
    des1 <- e2[i]
    des2 <- e2[j]
    v.l <- expm(Qmat * EL[i]) %*% liks[des1, ]
    v.r <- expm(Qmat * EL[j]) %*% liks[des2, ]
    v <- v.l * v.r
    txt = paste0("i=", i, ", j=", j)
    #	cat(txt, "\n")
    #	print(v)

    comp[anc] <- sum(v)
    liks[anc, ] <- v/comp[anc]
  }


  lnL = sum(log(comp[-TIPS]))

  # Error check
  if (is.finite(lnL) == FALSE)
  {
    lnL = -1e50
  }

  # print to screen during each optimization step
  txt = paste0("p1=", param_vector[1], ", p2=", param_vector[2], ", lnL=", lnL)
  cat(txt, "\n")

  return(lnL)
}

discrete_lnL_2x2(param_vector, fakedata, phy)


# Compare to ace's ML result for a 2-rate modle
"ARD" = "all rates different"
ace(x=fakedata, phy=phy, type="discrete", method="ML", model="ARD")
Log-likelihood: -5.057634
p2 = 0.4688
p1 = 2.1555


discrete_lnL_2x2(fakedata, phy, param_vector=c(2.1555,0.4688) )
# -5.057634
# matches, yay!

# Optimization with optim
# par = starting parameter guesses
# fn = the function we are optimizing
# method = L-BFGS-B (this searches within
#         lower/upper limits for params)
control=list(fnscale=-1) # change from minimizer
                         # to maximizer
# fakedata, phy : inputs to discrete_lnL_2x2
starting_parameters = c(0.1, 0.1)
optim(par=starting_parameters, fn=discrete_lnL_2x2, method="L-BFGS-B", lower=c(0,0), upper=c(100,100), control=list(fnscale=-1), fakedata=fakedata, phy=phy)



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
# 22.36105


# Set birthRate
#birthRate = 3.682184


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
# 2.869079

dev(a=x1, r=x2, N=N, x=x)
# 18.38194

dev(a=deathRate/birthRate, r=birthRate-deathRate, N=N, x=x)
# 18.38195

dev(a=0/birthRate, r=birthRate-0, N=N, x=x)

birthRate = 3.682184
dev(a=0/birthRate, r=birthRate-0, N=N, x=x)


S <- rep(1, times=length(tr$tip.label))
bd.ext(phy=tr, S=S, conditional=TRUE)
# d / b = 1.298031e-06   StdErr = 0.5292301 
# b - d = 1.96011   StdErr = 0.8721752 

bd.ext(phy=tr, S=S, conditional=FALSE) # same than older versions
# d / b = 4.355594e-06   StdErr = 0.5004528 
# b - d = 1.541568   StdErr = 0.7703756 




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
  
  
  res = NULL
  res$lambda = lambda
  res$se = se
  res$loglik = loglik
  
  return(res)
  }

yule_lik(tr=tr)
# $lambda
# [1] 2.61642
# 
# $se
# [1] 0.7552953
# 
# $loglik
# [1] 22.09385






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
# (ranges A, A, A, A...)
# states =  rep(2, times=length(tr$tip.label))
# states[1:5] = 1
states = rep(1, times=length(tr$tip.label)) # 14 states for 14 tips
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

lambda0 = birthRate*2
lambda1 = birthRate
mu0 = deathRate
mu1 = deathRate*2
q01 = 0.4
q10 = 0.2
parms = c(lambda0, lambda1, mu0, mu1, q01, q10)
names(parms) = c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
parms
bisse_params = parms

# Constraint parameters so you are only fitting 1 birthRate
# constraints = list(lambda0~lambda1, mu0~mu1, q01~q10)
# bisse_2areas_constrained = constrain(f=bisse_2areas, formulae=constraints)
# 
# Wait 1 seconds
# fit <- find.mle(func=bisse_2areas_constrained, x.init=bisse_params, method="subplex")
# fit$par.full
# fit$lnLik
# 
# Compare to Yule
# yule(tr)
# 
# bisse_params_orig = bisse_params
# bisse_params = fit$par.full

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

# ttl_LnL branch_LnL    ObsDiff   LnLdiff exp_ObsDiff exp_LnLdiff
# LnLs1   -8.071427   -9.12651   1.055082 -0.248400    2.872212    0.780029
# LnLs2   -8.712652   -9.12651   0.413858 -0.889600    1.512642    0.410800
# LnLs3   -8.045821   -9.12651   1.080688 -0.222800    2.946707    0.800261
# LnLs4   -8.324654   -9.12651   0.801856 -0.501700    2.229675    0.605531
# LnLs5   -8.712652   -9.12651   0.413858 -0.889600    1.512642    0.410800
# LnLs6  -10.455174   -9.12651  -1.328664   -2.6322    0.264831   0.0719222
# LnLs1t  -6.422516   -9.12651   2.703994    1.4005   14.939281    4.057179
# LnLs2t  -9.218484   -9.12651 -0.0919740   -1.3955    0.912129    0.247714
# LnLs3t  -5.855906   -9.12651   3.270603    1.9671   26.327221    7.149893
# LnLs4t  -8.202651   -9.12651   0.923859 -0.379600    2.518992    0.684103
# LnLs5t  -9.218484   -9.12651 -0.0919740   -1.3955    0.912129    0.247714
# LnLs6t -11.555104   -9.12651  -2.428594   -3.7321   0.0881607   0.0239425

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
mu0 = bisse_params["mu0"]
mu1 = bisse_params["mu1"]
# Character transition
q01 = bisse_params["q01"]
q10 = bisse_params["q10"]

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
# [1]  2.263549  4.527098  0.400000 -8.409282

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
  #   liks

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

# Speciation
lambda0 = bisse_params["lambda0"]
lambda1 = bisse_params["lambda1"]
# Extinction
mu0 = bisse_params["mu0"]
mu1 = bisse_params["mu1"]
# Character transition
q01 = bisse_params["q01"]
q10 = bisse_params["q10"]

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
# [,1]       [,2]         [,3]         [,4]
# [1,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [2,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [3,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [4,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [5,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [6,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [7,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [8,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [9,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [10,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [11,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [12,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [13,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [14,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [15,] 0.34032057 0.82565349 0.0785771910 2.9467074078
# [16,] 0.32296750 0.75484827 0.0686757722 2.9934403512
# [17,] 0.31781527 0.73677237 0.1941292826 2.5023823440
# [18,] 0.30871508 0.70679161 0.4556108127 1.3788904994
# [19,] 0.28090491 0.62494223 0.4946605287 0.2637140538
# [20,] 0.24328043 0.52701302 0.5082949255 0.0251604970
# [21,] 0.22481301 0.48206122 7.0267572837 0.0003531790
# [22,] 0.21940564 0.46918031 6.9495651423 0.0023754136
# [23,] 0.16335245 0.34117351 7.1082332562 0.0009958782
# [24,] 0.09896247 0.20257089 7.1977326413 0.0004767213
# [25,] 0.13910329 0.28814378 0.0045052886 3.4972487132
# [26,] 0.06810815 0.13830153 0.0005943981 3.6109085383
# [27,] 0.04142494 0.08358918 0.0004492111 3.6248920362
condlikes_treeStates_BRANCHTOP_AT_NODE_DOWNPASS
# [,1]       [,2]         [,3]         [,4]
# [1,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [2,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [3,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [4,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [5,] 0.00000000 0.00000000 1.0000000000 0.0000000000
# [6,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [7,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [8,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [9,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [10,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [11,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [12,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [13,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [14,] 0.00000000 0.00000000 0.0000000000 1.0000000000
# [15,] 0.34032057 0.82565349 0.0785772099 2.9467073685
# [16,] 0.32296750 0.75484826 0.0686757474 2.9934404908
# [17,] 0.31781526 0.73677236 0.1941292675 2.5023822075
# [18,] 0.30871509 0.70679161 0.4556107675 1.3788896493
# [19,] 0.28090493 0.62494230 0.4946609443 0.2637139148
# [20,] 0.24328046 0.52701310 0.5082953414 0.0251604964
# [21,] 0.22481305 0.48206130 7.0267568719 0.0003531796
# [22,] 0.21940567 0.46918038 6.9495647030 0.0023754174
# [23,] 0.16335247 0.34117354 7.1082330936 0.0009958792
# [24,] 0.09896248 0.20257090 7.1977326096 0.0004767214
# [25,] 0.13910330 0.28814380 0.0045052911 3.4972486560
# [26,] 0.06810816 0.13830153 0.0005943981 3.6109085351
# [27,] 0.04142494 0.08358918 0.0004492111 3.6248920353


# Matches matches diversitree BiSSE "base" -- but a weird combination of 
# - left columns: average E's passed down
# - right columns: normalized conditional likelihoods
base
condlikesProbs_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
# [,1]       [,2]        [,3]        [,4]
# [1,] 0.09896248 0.20257090 0.988621642 0.011378358
# [2,] 0.09896248 0.20257090 0.988621642 0.011378358
# [3,] 0.16335247 0.34117354 0.976418218 0.023581782
# [4,] 0.21940567 0.46918038 0.958618550 0.041381450
# [5,] 0.22481305 0.48206130 0.956252711 0.043747289
# [6,] 0.24328046 0.52701310 0.069531520 0.930468480
# [7,] 0.13910330 0.28814380 0.031275471 0.968724529
# [8,] 0.04142494 0.08358918 0.007810116 0.992189884
# [9,] 0.04142494 0.08358918 0.007810116 0.992189884
# [10,] 0.06810816 0.13830153 0.013425793 0.986574207
# [11,] 0.30871509 0.70679161 0.106504301 0.893495699
# [12,] 0.31781526 0.73677236 0.112544180 0.887455820
# [13,] 0.32296750 0.75484826 0.115939668 0.884060332
# [14,] 0.34032057 0.82565349 0.125156634 0.874843366
# [15,] 0.00000000 0.00000000 0.000000000 0.000000000
# [16,] 0.34032057 0.82565349 0.085252524 0.914747476
# [17,] 0.32296750 0.75484826 0.080433295 0.919566705
# [18,] 0.31781526 0.73677236 0.234224593 0.765775407
# [19,] 0.30871509 0.70679161 0.580886588 0.419113412
# [20,] 0.28090493 0.62494230 0.922758549 0.077241451
# [21,] 0.24328046 0.52701310 0.992656349 0.007343651
# [22,] 0.22481305 0.48206130 0.997807503 0.002192497
# [23,] 0.21940567 0.46918038 0.984410624 0.015589376
# [24,] 0.16335247 0.34117354 0.988531026 0.011468974
# [25,] 0.28090493 0.62494230 0.072792069 0.927207931
# [26,] 0.13910330 0.28814380 0.019560661 0.980439339
# [27,] 0.06810816 0.13830153 0.006011765 0.993988235
condlikesProbs_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
# [,1]       [,2]        [,3]        [,4]
# [1,] 0.09896248 0.20257090 0.988621642 0.011378358
# [2,] 0.09896248 0.20257090 0.988621642 0.011378358
# [3,] 0.16335247 0.34117354 0.976418218 0.023581782
# [4,] 0.21940567 0.46918038 0.958618550 0.041381450
# [5,] 0.22481305 0.48206130 0.956252711 0.043747289
# [6,] 0.24328046 0.52701310 0.069531520 0.930468480
# [7,] 0.13910330 0.28814380 0.031275471 0.968724529
# [8,] 0.04142494 0.08358918 0.007810116 0.992189884
# [9,] 0.04142494 0.08358918 0.007810116 0.992189884
# [10,] 0.06810816 0.13830153 0.013425793 0.986574207
# [11,] 0.30871509 0.70679161 0.106504301 0.893495699
# [12,] 0.31781526 0.73677236 0.112544180 0.887455820
# [13,] 0.32296750 0.75484826 0.115939668 0.884060332
# [14,] 0.34032057 0.82565349 0.125156634 0.874843366
# [15,] 0.00000000 0.00000000 0.000000000 0.000000000
# [16,] 0.34032057 0.82565349 0.085252524 0.914747476
# [17,] 0.32296750 0.75484826 0.080433295 0.919566705
# [18,] 0.31781526 0.73677236 0.234224593 0.765775407
# [19,] 0.30871509 0.70679161 0.580886588 0.419113412
# [20,] 0.28090493 0.62494230 0.922758549 0.077241451
# [21,] 0.24328046 0.52701310 0.992656349 0.007343651
# [22,] 0.22481305 0.48206130 0.997807503 0.002192497
# [23,] 0.21940567 0.46918038 0.984410624 0.015589376
# [24,] 0.16335247 0.34117354 0.988531026 0.011468974
# [25,] 0.28090493 0.62494230 0.072792069 0.927207931
# [26,] 0.13910330 0.28814380 0.019560661 0.980439339
# [27,] 0.06810816 0.13830153 0.006011765 0.993988235


# lq and likelihoods
condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS
# [1,] 0.09896248 0.20257090 0.601554008 0.006923475
# [2,] 0.09896248 0.20257090 0.601554008 0.006923475
# [3,] 0.16335248 0.34117356 0.386193121 0.009327071
# [4,] 0.21940568 0.46918041 0.228373385 0.009858376
# [5,] 0.22481306 0.48206131 0.214689102 0.009821741
# [6,] 0.24328047 0.52701312 0.018900198 0.252921816
# [7,] 0.13910330 0.28814381 0.017244206 0.534120979
# [8,] 0.04142494 0.08358918 0.006695332 0.850568780
# [9,] 0.04142494 0.08358918 0.006695332 0.850568780
# [10,] 0.06810816 0.13830153 0.010329973 0.759082537
# [11,] 0.30871507 0.70679157 0.012372479 0.103796344
# [12,] 0.31781526 0.73677234 0.010667106 0.084114390
# [13,] 0.32296750 0.75484825 0.009568187 0.072959105
# [14,] 0.34032057 0.82565349 0.004918531 0.034380470
# [15,] 0.00000000 0.00000000 0.000000000 0.000000000
# [16,] 0.34032057 0.82565349 0.130872973 1.404248417
# [17,] 0.32296750 0.75484826 0.189765397 2.169523726
# [18,] 0.31781527 0.73677237 0.342696258 1.120413378
# [19,] 0.30871510 0.70679166 0.239334953 0.172681709
# [20,] 0.28090493 0.62494229 0.269940202 0.022595914
# [21,] 0.24328046 0.52701309 5.564029159 0.041162574
# [22,] 0.22481304 0.48206128 6.532107017 0.014353091
# [23,] 0.21940565 0.46918034 4.200263734 0.066516442
# [24,] 0.16335245 0.34117352 4.619673712 0.053597627
# [25,] 0.28090494 0.62494230 0.084743377 1.079440830
# [26,] 0.13910329 0.28814379 0.050684149 2.540442463
# [27,] 0.06810816 0.13830153 0.019565251 3.234928386

# Matches diversitree BiSSE "lq"
lnls = log(rowSums(condlikes_treeStates_BRANCHBOTTOM_BELOW_NODE_DOWNPASS[,3:4]))
lnls[!is.finite(lnls)] = NA
lnls
sum(lnls, na.rm=TRUE)
# [1] -0.4967954 -0.4967954 -0.9275534 -1.4345113 -1.4938313 -1.3026078 -0.5953579 -0.1540092
# [9] -0.1540092 -0.2621280 -2.1527108 -2.3561811 -2.4946262 -3.2365562         NA  0.4286095
# [17]  0.8583604  0.3805641 -0.8866915 -1.2291671  1.7236933  1.8789245  1.4508595  1.5418593
# [25]  0.1520206  0.9520928  1.1800367
sum(lnls, na.rm=TRUE)
# -9.126511


# Sum of the branch likelihoods
lq
sum(lq)
# [1] -0.4967953 -0.4967953 -0.9275529 -1.4345101 -1.4938301 -1.3026068 -0.5953577 -0.1540092
# [9] -0.1540092 -0.2621280 -2.1527119 -2.3561821 -2.4946271 -3.2365566  0.0000000  0.4286096
# [17]  0.8583604  0.3805646 -0.8866917 -1.2291679  1.7236933  1.8789245  1.4508596  1.5418594
# [25]  0.1520213  0.9520928  1.1800367

# -9.12651

# Add the root probabilities
# Assuming diversitree options:
# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
# i.e., the root state probs are just the root_Ds/sum(root_Ds)
LnLs1
# -8.071427 -9.126510

# root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
LnLs1t
# -6.422516 -9.126510

# Does the total of branch likelihoods (lq) + node likelihoods match R?
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = sum(log(computed_likelihoods_at_each_node_x_lambda))
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
# 5.442733


R_result_branch_lnL = sum(lnls, na.rm=TRUE)
R_result_total_LnLs1 = LnLs1
R_result_total_LnLs1t = LnLs1t


R_result_branch_lnL
R_result_total_LnLs1
R_result_total_LnLs1t
R_result_sum_log_computed_likelihoods_at_each_node_x_lambda

#R_result_branch_lnL = -6.496757
#R_result_total_LnLs1 = -4.991969
#R_result_total_LnLs1t = -3.310171
#R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = 18.0192

# When some are type 0 and some are type 1:
#R_result_sum_log_computer_likelihoods_at_each_node_x_lambda = 5.442733 

#######################################################
#######################################################
# END PROOF OF MATCHING
#######################################################
#######################################################



