# Given some HBDS congruence class (in terms of PDR, PSR lambda_psi and CSA_pulled_probs), as well as some desired profile for a specific variable (e.g., the sampling rate psi), calculate the corresponding model in the congruence class
# The model in the congruence class must be specified through exactly one of the following variable choices:
# 	lambda (speciation rate)
#	mu (extinction rate) and lambda0 (present-day speciation rate)
#	psi (sampling rate)
#	Reff (effective reproduction ratio) and lambda0
#	removal_rate (aka. become-uninfectious rate, mu+psi) and lambda0
# ATTENTION: This function is currently only implemented for HBDS models without retention, i.e. where kappa=0
# ATTENTION: Not all psi are allowed. In particular, psi(CSA_ages[k]) must be the same for all models in the congruence class whenever CSA_ages[k]>0.
# ATTENTION: The input time profiles (e.g., for lambda, mu etc) should be at least continuous, i.e., piecewise constant with jumps may lead to incorrect outcomes because the routine does not properly capture Dirac-type derivatives.
congruent_hbds_model = function(age_grid,						# numeric vector of size NG, listing discrete ages in ascending order
								PSR,							# numeric vector of size NG, listing the PSR on the age_grid
								PDR,							# numeric vector of size NG, listing the PDR on the age_grid
								lambda_psi,						# numeric vector of size NG, listing lambda*psi on the age_grid
								lambda				= NULL,		# numeric vector of size NG. Exactly one of lambda, psi, mu, Reff or removal_rate must be provided.
								mu					= NULL,		# numeric vector of size NG. Exactly one of psi, mu, Reff or removal_rate must be provided.
								psi					= NULL,		# numeric vector of size NG. Exactly one of psi, mu, Reff or removal_rate must be provided.
								Reff				= NULL,		# numeric vector of size NG. Exactly one of psi, mu, Reff or removal_rate must be provided.
								removal_rate		= NULL,		# numeric vector of size NG. Exactly one of psi, mu, Reff or removal_rate must be provided.
								lambda0				= NULL,		# numeric, specifying lambda at age 0 (present-day). Only relevant if mu or Reff or removal_rate is provided.
								CSA_ages			= NULL,		# numeric vector of size NCSA, listing the ages of concentrated sampling attempts, in ascending order
								CSA_pulled_probs 	= NULL,		# numeric vector of size NCSA, listing the pulled probabilities of concentrated sampling attempts, in ascending order
								CSA_PSRs			= NULL, 	# numeric vector of size NCSA, listing the PSR at the concentrated sampling attempts, in ascending order
								splines_degree		= 1, 		# integer between 1 and 3. 0-splines are not allowed, because internally the 1st derivatives are needed.
								ODE_relative_dt		= 0.001,
								ODE_relative_dy		= 1e-4){
	# basic error checking
	NCSA = (if(is.null(CSA_ages)) 0 else length(CSA_ages))
	if((NCSA==0) && (!is.null(CSA_pulled_probs)) && (length(CSA_pulled_probs)>0)) return(list(success=FALSE, error="No CE ages were provided, but CSA_pulled_probs were"))
	if((NCSA>0) && is.null(CSA_pulled_probs)) return(list(success=FALSE, error="Missing CSA_pulled_probs"))
	if((NCSA>0) && (length(CSA_pulled_probs)!=NCSA)) return(list(success=FALSE, error=sprintf("Expected %d CSA_pulled_probs, but instead got %d",NCSA,length(CSA_pulled_probs))))
	if((NCSA>0) && is.null(CSA_PSRs)) return(list(success=FALSE, error="Missing CSA_PSRs"))
	if((NCSA>0) && (length(CSA_PSRs)!=NCSA)) return(list(success=FALSE, error=sprintf("Expected %d CSA_PSRs, but instead got %d",NCSA,length(CSA_PSRs))))
	if(is.null(CSA_ages)){
		CSA_ages 			= numeric(0)
		CSA_pulled_probs 	= numeric(0)
		CSA_PSRs 			= numeric(0)
	}
	NG = length(age_grid)
	if(is.null(PSR)){
		PSR = rep(0,times=NG)
	}else if(length(PSR)==1){
		PSR = rep(PSR,times=NG)
	}else if(length(PSR)!=NG){
		return(list(success=FALSE, error=sprintf("Expected %d PSR values, but instead got %d",NG,length(PSR))))
	}
	if(is.null(PDR)){
		PDR = rep(0,times=NG)
	}else if(length(PDR)==1){
		PDR = rep(PDR,times=NG)
	}else if(length(PDR)!=NG){
		return(list(success=FALSE, error=sprintf("Expected %d PDR values, but instead got %d",NG,length(PDR))))
	}
	if(is.null(lambda_psi)){
		lambda_psi = rep(0,times=NG)
	}else if(length(psi)==1){
		lambda_psi = rep(lambda_psi,times=NG)
	}else if(length(lambda_psi)!=NG){
		return(list(success=FALSE, error=sprintf("Expected %d lambda_psi values, but instead got %d",NG,length(lambda_psi))))
	}
	if(is.null(lambda) && is.null(mu) && is.null(psi) && is.null(Reff) && is.null(removal_rate)) return(list(success=FALSE, error=sprintf("Expecting either lambda, mu, psi, Reff or removal_rate")))
	if(sum(!c(is.null(lambda),is.null(mu),is.null(psi),is.null(Reff),is.null(removal_rate)))>1) return(list(success=FALSE, error=sprintf("Only one of lambda, mu, psi, Reff or removal_rate must be provided")))
	if((!is.null(lambda)) && (!is.null(lambda0))) return(list(success=FALSE, error=sprintf("lambda0 must not be provided if lambda is provided")))
	if((!is.null(mu)) && is.null(lambda0)) return(list(success=FALSE, error=sprintf("lambda0 must be provided when mu is provided")))
	if((!is.null(psi)) && (!is.null(lambda0))) return(list(success=FALSE, error=sprintf("lambda0 must not be provided if psi is provided")))
	if((!is.null(Reff)) && is.null(lambda0)) return(list(success=FALSE, error=sprintf("lambda0 must be provided when Reff is provided")))
	if((!is.null(removal_rate)) && is.null(lambda0)) return(list(success=FALSE, error=sprintf("lambda0 must be provided when removal_rate is provided")))
	if(!is.null(lambda)){
		if(length(lambda)==1){
			lambda = rep(lambda,times=NG)
		}else if(length(lambda)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d lambda values, but instead got %d",NG,length(lambda))))
		}
		if(NCSA>0) return(list(success=FALSE, error=sprintf("Providing lambda to define a model is only available in the absence of CSAs")))
	}
	if(!is.null(mu)){
		if(length(mu)==1){
			mu = rep(mu,times=NG)
		}else if(length(mu)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d mu values, but instead got %d",NG,length(mu))))
		}
		if(NCSA>0) return(list(success=FALSE, error=sprintf("Providing mu to define a model is only available in the absence of CSAs")))
	}
	if(!is.null(psi)){
		if(length(psi)==1){
			psi = rep(psi,times=NG)
		}else if(length(psi)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d psi values, but instead got %d",NG,length(psi))))
		}
	}
	if(!is.null(Reff)){
		if(length(Reff)==1){
			Reff = rep(Reff,times=NG)
		}else if(length(Reff)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d Reff values, but instead got %d",NG,length(Reff))))
		}
		if(NCSA>0) return(list(success=FALSE, error=sprintf("Providing Reff to define a model is only available in the absence of CSAs")))
	}
	if(!is.null(removal_rate)){
		if(length(removal_rate)==1){
			removal_rate = rep(removal_rate,times=NG)
		}else if(length(removal_rate)!=NG){
			return(list(success=FALSE, error=sprintf("Expected %d removal_rate values, but instead got %d",NG,length(removal_rate))))
		}
		if(NCSA>0) return(list(success=FALSE, error=sprintf("Providing removal_rate to define a model is only available in the absence of CSAs")))
	}
	if(!(splines_degree %in% c(1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree (%d): Expected one of 1,2,3.",splines_degree)))
	
	results = get_congruent_HBDS_CPP(	CSA_ages			= CSA_ages,
										CSA_pulled_probs	= CSA_pulled_probs,
										CSA_PSRs			= CSA_PSRs,
										age_grid			= age_grid,
										PSRs				= PSR,
										PDRs				= PDR,
										lambda_psis			= lambda_psi,
										lambdas				= (if(is.null(lambda)) numeric(0) else lambda),
										mus					= (if(is.null(mu)) numeric(0) else mu),
										psis				= (if(is.null(psi)) numeric(0) else psi),
										Reffs				= (if(is.null(Reff)) numeric(0) else Reff),
										removal_rates		= (if(is.null(removal_rate)) numeric(0) else removal_rate),
										lambda0				= (if(is.null(lambda0)) 0 else lambda0),
										splines_degree		= splines_degree,
										ODE_relative_dt		= ODE_relative_dt,
										ODE_relative_dy		= ODE_relative_dy,
										runtime_out_seconds	= -1)

	if(!results$success) return(list(success=FALSE, error=results$error))
	return(list(success			= TRUE,
				valid			= results$valid,
				ages			= age_grid,
				lambda			= results$lambdas,
				mu				= results$mus,
				psi				= results$psis,
				lambda_psi		= results$lambda_psis,
				Reff			= results$Reffs,
				removal_rate	= results$removal_rates,
				Pmissing		= results$Pmissings,
				CSA_probs		= results$CSA_probs,
				CSA_Pmissings	= results$CSA_Pmissings))
}