# WARNING: THIS FUNCTION IS EXPERIMENTAL, AND NOT YET WORKING FULLY
# THE PROBLEM IS THAT THE FORMULA USED TO CALCULATE THE RATE OF DECAY OF THE GENE-LTT WRONGLY ASSUMES ZERO VARIANCE IN THE NUMBER OF ALLELES IN EACH SPECIES-BRANCH
# THIS CAUSES MILD INACCURACIES WHEN COALESCENCE TIMES ARE LONG
# 
# Predict various deterministic features of a nested homogenous birth-death (HBD) multispecies coalescent (MSC) gene tree model, backward in time.
# The pulled speciation rate (PSR) PSR_p and coalescence_time (Ne*generation_time) are specified on a discrete age-grid, and assumed to vary linearly (or polynomially, as splines) between grid points (see "degree" argument).
# This function calculates, among others, the following features over time:
#	Deterministic species-LTT curve
#	Deterministic gene-LTT curve
#
simulate_deterministic_hbd_msc = function(	sLTT0, 						# positive numeric, number of extant species represented in the tree at present-day (age 0), i.e. after rarefaction..
											oldest_age,					# positive numeric, specifying how far back (time before present) to simulate the model
											gLTT0			= NULL, 	# positive numeric, number of extant alleles represented in the gene tree at present-day (age 0), i.e. after rarefaction. This is equal to the number of species sampled multiplied by the number of alleles sampled per species. If NULL, this is assumed to be equal to sLTT0 (i.e., one allele per species). gLTT0 must be at least as large as sLTT0.
											rho0			= 1,		# numeric within (0,1], specifying the fraction of extant species represented in the tree at present-day. Can also be NULL, which is equivalent to setting rho0=1.
											age_grid		= NULL,		# either NULL, or empty, or a numeric vector of size NG, listing ages in ascending order, on which birth/mu are specified. If NULL or empty, then PSR[] and coalescent_time[] must be a single scalar. The returned time series will be defined on an age-grid that may be finer than this grid. If of size >=2, the age_grid must cover oldest_age and 0.
											PSR				= 0,		# either a single numeric (constant PSR over time), or a numeric vector of size NG (listing PSRs at each age in grid_ages[]).
											CT				= 0,		# either a single non-negative numeric (constant coalescence-time over time), or a numeric vector of size NG (listing coalescence time at each age in grid_ages[]).
											splines_degree	= 1,		# integer, either 1 or 2 or 3, specifying the degree for the splines defined by PSR, mu and PDR on the age grid.
											relative_dt		= 1e-3){	# maximum relative time step allowed for integration. Smaller values increase integration accuracy but increase computation time. Typical values are 0.0001-0.001. The default is usually sufficient.	
	# check validity of input variables
	age0 = 0;
	if(is.null(gLTT0)) gLTT0 = sLTT0;
	if(gLTT0<sLTT0) return(list(success = FALSE, error = sprintf("gLTT0 (number of alleles at present-day) must not be smaller than sLTT0 (number of species at present-day).")))
	if(is.null(rho0)) rho0 = 1;
	if(is.null(PSR)) return(list(success = FALSE, error = sprintf("Missing PSR (pulled speciation rates)")))
	if(is.null(CT)) return(list(success = FALSE, error = sprintf("Missing CT (coalescence times)")))
	if(is.null(age_grid) || (length(age_grid)<=1)){
		if((!is.null(PSR)) && (length(PSR)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of PSR values; since no age grid was provided, you must either provide a single (constant) PSR or none")))
		if((!is.null(CT)) && (length(CT)!=1)) return(list(success = FALSE, error = sprintf("Invalid number of CT values; since no age grid was provided, you must provide a single (constant) CT")))
		# create dummy age grid
		NG 			= 2;
		age_grid	= seq(from=0,to=oldest_age,length.out=NG)
		if(!is.null(PSR)) PSR = rep(PSR,times=NG);
		if(!is.null(CT)) CT = rep(CT,times=NG);
	}else{
		NG = length(age_grid);
		if((age_grid[1]>oldest_age) || (age_grid[NG]<oldest_age)) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including oldest_age (%g)",oldest_age)))
		if((age_grid[1]>age0) || (age_grid[NG]<age0)) return(list(success = FALSE, error = sprintf("Age grid must cover the entire requested age interval, including age %g",age0)))
		if((!is.null(PSR)) && (length(PSR)!=1) && (length(PSR)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of PSR values; since an age grid of size %d was provided, you must either provide zero, one or %d PSRs",NG,NG)))
		if((!is.null(CT)) && (length(CT)!=1) && (length(CT)!=NG)) return(list(success = FALSE, error = sprintf("Invalid number of CT values; since an age grid of size %d was provided, you must either provide one or %d CTs",NG,NG)))
		if((!is.null(PSR)) && (length(PSR)==1)) PSR = rep(PSR,times=NG);
		if((!is.null(CT)) && (length(CT)==1)) CT = rep(CT,times=NG);
	}
	if(!(splines_degree %in% c(0,1,2,3))) return(list(success = FALSE, error = sprintf("Invalid splines_degree: Expected one of 0,1,2,3.")))

		
	# simulate model backward in time
	simulation = simulate_deterministic_HBD_MSC_CPP(	oldest_age			= oldest_age,
														age_grid 			= age_grid,
														PSRs 				= (if(is.null(PSR)) numeric() else PSR),
														CTs				 	= (if(is.null(CT)) numeric() else CT),
														rho0				= rho0,
														sLTT0				= sLTT0,
														gLTT0				= gLTT0,
														splines_degree	= splines_degree,
														relative_dt		= relative_dt);
	if(!simulation$success) return(list(success = FALSE, error = sprintf("Could not simulate model: %s",simulation$error)))

	return(list(success			= TRUE,
				ages			= simulation$refined_age_grid, # potentially refined ages grid, on which all returned variables are defined
				PSR				= simulation$refined_PSRs,	# PSRs defined (refined) on ages[]
				CT				= simulation$refined_CTs, 	# coalescence times (Ne*generation_time) defined (refined) on ages[]
				sLTT			= simulation$sLTT,
				gLTT			= simulation$gLTT));
}

