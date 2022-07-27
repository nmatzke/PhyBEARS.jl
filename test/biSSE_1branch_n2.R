#######################################################
#
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, pure birth, no Q transitions, branchlength=1)
# 
#
# Run with:
#
# source("/GitHub/PhyBEARS.jl/test/biSSE_1branch_n2.R")
# 
#######################################################


# ClaSSE helper functions
library(deSolve)  # for lsoda
library(ape)
library(diversitree)
# library(rexpokit)
# library(cladoRcpp)
# library(BioGeoBEARS)

source("/GitHub/PhyBEARS.jl/Rsrc/ClaSSE_functions_v3.R")  # utility functions from diversitree
source("/GitHub/PhyBEARS.jl/Rsrc/ClaSSE_pureR_v1.R") # simple implementations in plain-R



#######################################################
#######################################################
# 5. BiSSE likelihoods and comparison to yule
#######################################################
#######################################################
# Look at all the parameters of this model!
# lambdas = speciation rates
# mus = extinction rates
# qs = anagenetic transition rates
birthRate = 0.222222222
deathRate = 0.4


# Speciation
lambda0 = birthRate
lambda1 = birthRate
# Extinction
mu0 = deathRate
mu1 = deathRate
# Character transition
q01 = 0.1 # ML
q10 = 0.1

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
times = seq(from=0, to=1, by=1/50)
out <- lsoda(y=y, times=times, func=define_BiSSE_eqns_in_R, parms=parms) 
tail(out)

# The last values of "out" are:
result_EsDs = out[nrow(out), ]
result_EsDs

txt = paste0(result_EsDs, collapse=" ")
juliatxt = paste0("R_result_EsDs = [", txt, "]")
cat("\n")
cat("Output of this R script:")
cat(juliatxt)
cat("\n")



