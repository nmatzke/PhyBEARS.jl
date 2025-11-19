# Predict various deterministic features of a homogenous birth-death-sampling cladogenic model, backward in time.
# The speciation rate lambda, extinction rate mu, sampling rate psi and retention probability kappa are specified on a discrete age-grid, and assumed to vary linearly (or polynomially, as splines) between grid points (see "splines_degree" argument).
# This function calculates, among others, the following features over time:
#	Deterministic LTT curve
#	Deterministic total population size
#   Pulled diversification rate (PDR)
#	Basic reproduction ratio (R0)
simulate_deterministic_hbds = function(	age_grid						= NULL,		# numeric vector listing grid ages in ascending order, on which the various model parameters (lambda, mu, rho, kappa) are specified. The age grid must generally cover the maximum possible simulation period (i.e. from 0 to max(requested_ages)), unless splines_degree=0 (in which case everything is extrapolated as constant if needed)
										lambda							= NULL,		# numeric vector of the same length as age_grid[], listing per-capita birth rates (speciation rates) at each age_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
										mu								= NULL,		# numeric vector of the same length as age_grid[], listing per-capita death rates (extinction rates) at each age_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
										psi								= NULL,		# numeric vector of the same length as age_grid[], listing per-capita sampling rates (Poissonian detection rates) at each age_grid point. Can also be a single number. Can also be NULL, which is the same as being zero.
										kappa							= NULL,		# numeric vector of the same length as age_grid[], listing the retention probability (upon sampling) at each age_grid point, i.e. the probability that a sampled lineage remains in the pool. If 0, then every sampled lineage becomes a tip.
										splines_degree					= 1,		# polynomial degree of time-dependent model parameters (lambda, mu, psi, kappa) between time-grid points
										CSA_ages						= NULL,		# optional numeric vector listing concentrated sampling ages, in ascending order
										CSA_probs						= NULL,		# optional numeric vector listing concentrated sampling probabilities, corresponding to CSA_ages[]
										CSA_kappas						= NULL,		# optional numeric vector listing retention probabilities during concentrated sampling attempts, corresponding to CSA_ages[]
										requested_ages					= NULL,		# optional numeric vector listing ages (in ascending order), for which the various model variables are requested. If NULL, it will be set to age_grid.
										age0							= 0,		# non-negative numeric, age at which N0 or LTT0 are specified
										N0								= NULL, 	# number of extant species, i.e. total diversity (sampled or not) at age0. Used as a "boundary condition" that determines the scale of the LTT.
										LTT0							= NULL, 	# number of lineages represented in the tree at age0. Used as a "boundary condition" that determines the scale of the LTT.
										ODE_relative_dt					= 0.001,	# (positive unitless number) relative integration time step for the ODE solvers. Relative to the typical time scales of the dynamics, as estimated from the theoretically maximum possible rate of change. Typical values are 0.001 - 0.1.
										ODE_relative_dy					= 1e-4){	# (positive unitless mumber) unitless number, relative step for interpolating E over time. So a E_value_step of 0.001 means that E is recorded and interpolated between points between which E differs by roughy 0.001. Typical values are 0.01-0.0001. A smaller E_value_step increases interpolation accuracy, but also increases memory requirements and adds runtime (scales with the tree's age span, not Ntips).
	# check validity of input variables
	if(is.null(requested_ages) && is.null(age_grid)) return(list(success=FALSE, error="At least one of age_grid and requested_ages must be non-null"))
	if(is.null(requested_ages)){
		requested_ages = age_grid
	}else{
		if(any(diff(requested_ages)<0)) return(list(success=FALSE, error="requested_ages must be in ascending order"))
		if(any(requested_ages<0)) return(list(success=FALSE, error="requested_ages must be non-negative"))
	}
	oldest_age = tail(requested_ages,1)
	if(is.null(age_grid) || (length(age_grid)<=1)){
		# create dummy age grid
		NG 		 = 2;
		age_grid = seq(from=0,to=1.01*oldest_age,length.out=NG)
		if((!is.null(lambda)) && (length(lambda)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of lambda values (%d); since no age grid was provided, you must either provide a single (constant) birth rate or NULL",length(lambda))))
		if((!is.null(mu)) && (length(mu)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of mu values (%d); since no age grid was provided, you must either provide a single (constant) death rate or none",length(mu))))
		if((!is.null(psi)) && (length(psi)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of psi values (%d); since no age grid was provided, you must either provide a single (constant) sampling rate or none",length(psi))))
		if((!is.null(kappa)) && (length(kappa)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of kappa values (%d); since no age grid was provided, you must either provide a single (constant) retention probability or none",length(kappa))))
		if(!is.null(lambda)){
			lambda = rep(lambda,times=NG)
		}else{
			lambda = rep(0,times=NG)
		}
		if(!is.null(mu)){
			mu = rep(mu,times=NG)
		}else{
			mu = rep(0,times=NG)
		}
		if(!is.null(psi)){
			psi = rep(psi,times=NG)
		}else{
			psi = rep(0,times=NG)
		}
		if(!is.null(kappa)){
			kappa = rep(kappa,times=NG)
		}else{
			kappa = rep(0,times=NG)
		}
	}else{
		NG = length(age_grid);
		if(any(diff(age_grid)<0))return(list(success = FALSE, error = sprintf("Values in age_grid must be strictly increasing")))
		if((splines_degree>0) && ((age_grid[1]>oldest_age) || (age_grid[NG]<oldest_age))) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including oldest_age (%g), when splines_degree>0",oldest_age)))
		if((splines_degree>0) && ((age_grid[1]>0) || (age_grid[NG]<0))) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including present-day (age 0), when splines_degree>0")))
		if(is.null(lambda)){
			lambda = rep(0,times=NG)
		}else if(length(lambda)==1){
			lambda = rep(lambda,times=NG)
		}else if(length(lambda)!=NG){
			return(list(success=FALSE, error=sprintf("Expected either a single birth-rate lambda or exactly %d birth-rates (=age_grid length), but instead got %d",NG,length(lambda))))
		}
		if(is.null(mu)){
			mu = rep(0,times=NG)
		}else if(length(mu)==1){
			mu = rep(mu,times=NG)
		}else if(length(mu)!=NG){
			return(list(success=FALSE, error=sprintf("Expected either a single death-rate mu or exactly %d death-rates (=age_grid length), but instead got %d",NG,length(mu))))
		}
		if(is.null(psi)){
			psi = rep(0,times=NG)
		}else if(length(psi)==1){
			psi = rep(psi,times=NG)
		}else if(length(psi)!=NG){
			return(list(success=FALSE, error=sprintf("Expected either a single sampling-rate psi or exactly %d sampling-rates (=age_grid length), but instead got %d",NG,length(psi))))
		}
		if(is.null(kappa)){
			kappa = rep(0,times=NG)
		}else if(length(kappa)==1){
			kappa = rep(kappa,times=NG)
		}else if(length(kappa)!=NG){
			return(list(success=FALSE, error=sprintf("Expected either a single retention probability kappa or exactly %d probabilities (=age_grid length), but instead got %d",NG,length(kappa))))
		}
	}
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 0,1,2,3.",splines_degree)))
	NCE = (if(is.null(CSA_ages)) 0  else length(CSA_ages))
	if((NCE==0) && (!is.null(CSA_probs)) && (length(CSA_probs)>0)) return(list(success=FALSE, error="CSA_ages is missing while CSA_probs was provided; either provide both or none"))
	if((NCE==0) && (!is.null(CSA_kappas)) && (length(CSA_kappas)>0)) return(list(success=FALSE, error="CSA_ages is missing while CSA_kappas was provided; either provide both or none"))
	if((NCE>0) && is.null(CSA_probs)) return(list(success=FALSE, error="CSA_probs is missing while CSA_ages was provided; either provide both or none"))
	if((NCE>0) && is.null(CSA_kappas)) return(list(success=FALSE, error="CSA_kappas is missing while CSA_kappas was provided; either provide both or none"))
	if((NCE>0) && (!is.null(CSA_probs)) && (!is.null(CSA_kappas))){
		if(length(CSA_ages)!=length(CSA_probs)) return(list(success=FALSE, error="Number of CSA_ages (%d) differs from number of CSA_probs (%d)",length(CSA_ages),length(CSA_probs)))
		if(length(CSA_ages)!=length(CSA_kappas)) return(list(success=FALSE, error="Number of CSA_ages (%d) differs from number of CSA_kappas (%d)",length(CSA_ages),length(CSA_kappas)))
		if(any(diff(CSA_ages)<=0)) return(list(success=FALSE, error="CSA_ages must be in strictly increasing order"))
		if(any(CSA_probs<0) || any(CSA_probs>1)) return(list(success=FALSE, error="CSA_probs must be true probabilities, and thus between 0 and 1"))
		if(any(CSA_kappas<0) || any(CSA_kappas>1)) return(list(success=FALSE, error="CSA_kappas must be true probabilities, and thus between 0 and 1"))
	}
	if(is.null(LTT0) && is.null(N0)) return(list(success=FALSE, error="Either N0 or LTT0 must be specified"))
	if((!is.null(LTT0)) && (!is.null(N0))) return(list(success=FALSE, error="Either N0 or LTT0 must be specified, but not both"))
	if(age0<0) return(list(success=FALSE, error="The reference age (age0) must be non-negative"))
	if((!is.null(LTT0)) && (age0==0) && ((NCE==0) || CSA_ages[1]>0)) return(list(success=FALSE, error="LTT0 cannot be specified at age=0 since there is no concentrated sampling attempt at age 0. Either set age0>0, or use N0 instead of N0 for scaling"))
	
	if(splines_degree==0){
		# represent the time-curves on a time-grid as piecewise-linear functions, because simulate_deterministic_HBDS_CPP() cannot work directly with splines_degree 0
		refinement_factor	= 100
		refined_age_grid 	= c(unlist(lapply(seq_len(NG-1), FUN=function(g) seq(from=age_grid[g],to=age_grid[g+1]*(1-1/refinement_factor),length.out=refinement_factor))),age_grid[NG])
		if(refined_age_grid[1]>0) refined_age_grid = c(0,refined_age_grid)
		if(tail(refined_age_grid,1)<oldest_age) refined_age_grid = c(refined_age_grid,oldest_age)
		lambda 				= evaluate_spline(Xgrid = age_grid, Ygrid = lambda, splines_degree = 0, Xtarget = refined_age_grid, extrapolate="const", derivative = 0)
		mu 					= evaluate_spline(Xgrid = age_grid, Ygrid = mu, splines_degree = 0, Xtarget = refined_age_grid, extrapolate="const", derivative = 0)
		psi 				= evaluate_spline(Xgrid = age_grid, Ygrid = psi, splines_degree = 0, Xtarget = refined_age_grid, extrapolate="const", derivative = 0)
		kappa 				= evaluate_spline(Xgrid = age_grid, Ygrid = kappa, splines_degree = 0, Xtarget = refined_age_grid, extrapolate="const", derivative = 0)
		splines_degree 		= 1
		age_grid 			= refined_age_grid
	}

	simulation = simulate_deterministic_HBDS_CPP(	CSA_ages				= (if(NCE==0) numeric(0) else CSA_ages),
													CSA_probs				= (if(NCE==0) numeric(0) else CSA_probs),
													CSA_kappas				= (if(NCE==0) numeric(0) else CSA_kappas),
													age_grid				= age_grid,
													lambdas					= lambda,
													mus						= mu,
													psis					= psi,
													kappas					= kappa,
													splines_degree			= splines_degree,
													age0					= age0,
													N0						= (if(is.null(N0)) -1 else N0),
													LTT0					= (if(is.null(LTT0)) -1 else LTT0),
													requested_ages			= requested_ages,
													ODE_relative_dt			= ODE_relative_dt,
													ODE_relative_dy			= ODE_relative_dy,
													runtime_out_seconds		= -1)
	if(!simulation$success) return(list(success = FALSE, error = sprintf("Could not simulate model: %s",simulation$error)))
	NRA = length(requested_ages)
	
	# calculate area-under-the-curve for the LTT, in order to get a normalized version of the LTT
	# Note that the LTT's AUC is only calculated over the domain spanned by requested_ages!
	LTT_AUC = sum(0.5 * (simulation$LTT[2:NRA]+simulation$LTT[1:(NRA-1)]) * diff(requested_ages))
	simulation$nLTT = simulation$nonscaledLTT/sum(0.5 * (simulation$nonscaledLTT[2:NRA]+simulation$nonscaledLTT[1:(NRA-1)]) * diff(requested_ages)) # use the nonscaled-LTT to get the nLTT, since the scaled LTT may sometimes be NaN.

	return(list(success				= TRUE,
				ages				= requested_ages,
				total_diversity		= simulation$total_diversity,
				LTT					= simulation$LTT,
				nLTT				= simulation$nLTT,
				Pmissing			= simulation$Pmissing,
				lambda				= simulation$lambda,
				mu					= simulation$mu,
				psi					= simulation$psi,
				kappa				= simulation$kappa,
				PDR					= simulation$PDR,	# pulled diversification rate
				IPDR				= simulation$IPDR,	# age-integrated PDR, \int_0^t PDR(s) ds
				PSR					= simulation$PSR,	# pulled speciation rate
				PRP					= simulation$PRP,	# pulled retention probability
				diversification_rate= simulation$diversification, # net diversification rate, lambda-mu-psi
				branching_density	= simulation$PSR * simulation$nLTT, # non-normalized probability density of branching points over time (in units 1/time)
				sampling_density	= simulation$psi * simulation$total_diversity/LTT_AUC, # non-normalized probability density of sampling points over time (in units 1/time)
				lambda_psi			= simulation$lambda * simulation$psi,
				kappa_psi			= simulation$psi * simulation$kappa,
				Reff				= simulation$Reff, # effective reproduction ratio
				removal_rate		= simulation$mu+simulation$psi, # removal rate = mu+psi
				sampling_proportion = simulation$psi/(simulation$mu + simulation$psi),
				CSA_ages			= CSA_ages,
				CSA_probs			= CSA_probs,
				CSA_kappas			= CSA_kappas,
				CSA_pulled_probs 	= simulation$CSA_pulled_probs, # pulled_rho_k := rho_k/(1-Pmissing)
				CSA_psis			= simulation$CSA_psis,
				CSA_PSRs			= simulation$CSA_PSRs)) # PSRs exactly at the concentrated sampling attempts
}
