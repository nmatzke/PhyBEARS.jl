

# t = current time point of integration
# y = state variable we are tracking (named)  MUST HAVE NAMES!!!
# parms = model parameters (named)            MUST HAVE NAMES!!!

define_BiSSE_eqns_in_R <- function(t, y, parms)
	{
	defaults='
	# Starting values of state variables
	E0t = 0
	E1t = 0
	D0t = 0
	D1t = 1
	y = c(E0t, E1t, D0t, D1t)
	names(y) = c("E0t", "E1t", "D0t", "D1t")
	y
	# One step test:
	one_step_result = define_BiSSE_eqns_in_R(t=t, y=y, parms=parms)
	one_step_result
# [[1]]
# [1]  0.0000000  0.0000000  0.0000000 -0.2222222

	'
	
	
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
	} # END define_BiSSE_eqns_in_R <- function(t, y, parms)



#######################################################
# Add ClaSSE/tree likelihoods to BioGeoBEARS result
# (assumes lineage extinction rate = 0)
#######################################################
add_ClaSSE_to_BGBres <- function(res)
	{
	# Get the tree
	tr = read.tree(res$inputs$trfn)
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=res$inputs$geogfn)
	
	birthRate = yule(tr)$lambda  # The ML lambda from Yule. Equals (#speciations-1)/tree_length
birthRate
(tr$Nnode-1)/sum(tr$edge.length)

deathRate = 0.0     # Yule process means 0.0 extinction rate

	}




