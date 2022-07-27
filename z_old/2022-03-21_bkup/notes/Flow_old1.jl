module Flow
print("\n\nStarting module 'Flow'...loading dependencies...\n")
using LinearAlgebra  # for mul! (matrix multiplication)
using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA
using DifferentialEquations
using Distributed
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloNetworks
#using Plots						# for plot
using DataFrames          # for DataFrame()
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs

export parameterized_ClaSSE_As_v5, check_linearDynamics_of_As, calc_Gs_SSE_condnums!, calc_Gs_SSE, calc_Gs_SSE!, run_Gs


# Construct interpolation function for calculating linear dynamics A, 
# at any timepoint t

#function get_LinearDynamics_A(inputs)

# 	modelD.getLinearDynamics(ages[a], dynamics); // get linear dynamics at this age
# 		// calculate largest singular value (sigma1) of the dynamics at this age
# 		// then kappa_rate <= 2*sigma1, since sigma1(exp(t*A))/sigma2(exp(t*A)) <= exp(t*2*sigma1(A)) [So and Thompson (2000) Singular values of Matrix Exponentials. Theorem 2.1]
#
# NOTE: In Julia, 
# When p=2, the operator norm is the spectral norm, equal to the largest singular value of A.
# 
# 
#
# This version excludes Xi (and Xj), just like castor's get_LinearDynamics_A
# calculation of A
parameterized_ClaSSE_As_v5 = (t, A, p) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	two = 1.0
	# Iterate through the ancestral states
  @inbounds for i in 1:n
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		A[i,i] = A[i,i]  + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i]) # *u[i]  
		
		# case 2: anagenetic change
		@inbounds for m in 1:length(Qi_sub_i)
			A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_sub_i[m]] #* u[Qj_sub_i[m]])
		end
		
		# case 34: change + eventual extinction
		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			# excluding the u[], i.e. the Ds, i.e. the Xs, just as is done in 
			# 2*speciation_rates[r]*current_E[r]
			#rate_sp_then_ex = Cijk_vals[Ci_sub_i] * ((u[Ck_sub_i] * uE[Cj_sub_i]) + (u[Cj_sub_i] * uE[Ck_sub_i]))
			rate_sp_then_ex = Cijk_vals[Ci_sub_i[m]] * (uE[Cj_sub_i[m]] + uE[Ck_sub_i[m]])
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex
		end
		
# 		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
# 			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
# 			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
# 				 (u[Ck_sub_i].*uE[Cj_sub_i] 
# 			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))

  end # End @inbounds for i in 1:n
 	return(A)
end


# Check the linear dynamics (matrix A) for timespans where the kappa rate
# exceeds log(max_condition_number). max_condition_number is usually between
# 1e4 (slower) and 1e8 (faster)
# 
# "exponential growth rate of the condition number ("kappa_rate") based on the 
#  largest singular value of the dynamics A" 

function check_linearDynamics_of_As(tvals, p_Ds_v5; max_condition_number=1e8)

	# build an A matrix to hold the linear dynamics
	n = p_Ds_v5.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))
	
	# build arrays to hold the output for each t
	upper_bound_kappa_rates_A = collect(repeat([0.0], length(tvals)))
	cond_num_too_big_TF = collect(repeat([false], length(tvals)))
	
	for i in 1:length(tvals)
		t = tvals[i]
		A_at_t = Flow.parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
		sigma1_of_A = opnorm(A,1)  # the 1-norm should be adequate here (fastest calc.)
		upper_bound_kappa_rates_A[i] = 2*sigma1_of_A
		if (upper_bound_kappa_rates_A[i] > log(max_condition_number))
			cond_num_too_big_TF[i] = true
		end
	end
	
	
	# Creat dataframe
	kappa_Arates_df = DataFrames.DataFrame(tvals=tvals, ub_kappa_ratesA=upper_bound_kappa_rates_A, cond_num_too_big=cond_num_too_big_TF)

	return(kappa_Arates_df)
end

# Calculate G flow matrix down time series, outputting various
# condition numbers and kappa rates
calc_Gs_SSE_condnums! = (dG, G, pG, t) -> begin
	#A = pG.A             # Initial A
	p_Ds_v5 = pG.p_Ds_v5 # Calculate Ds down a timespan
	
	# A as a function of time t
	A = parameterized_ClaSSE_As_v5(t, pG.A[:,:], p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)
#	mul!(dG, A, G)
	mul!(dG, parameterized_ClaSSE_As_v5(t, pG.A[:,:], p_Ds_v5), G)

	condG1 = cond(G,1)
	condG2 = cond(G,2)
	condGInf = cond(G,Inf)
	tmpstr = paste0(["\nAt time t=", string(round(t, digits=6)), ", Condition number of G=", string(round(condG1, digits=6)), ", ", string(round(condG2, digits=6)), ", ", string(round(condGInf, digits=6))])
	print(tmpstr)
	#display(cond(G))
	
	# From Louca & Pennell code:
	# phylogenetics_cpp_routines.cpp
	# max_condition_number,			// (INPUT) unitless number, 
	# the maximum acceptable condition number for the Gmap 
	# (as estimated from the linearized dynamics), when choosing 
	# the integration interval size. A larger max_condition number 
	# leads to fewer age-splits, thus faster computation but also 
	# lower accuracy. Hence, this number controls the trade-off 
	# between speed and accuracy. Typical values are 
	# 1e4 (slower, more accurate) up to 
	# 1e8 (faster, less accurate).
	
	# The condition number is kappa:
	# https://julia.quantecon.org/tools_and_techniques/iterative_methods_sparsity.html
	# The cond() operation can be expensive, so the inequality in Louca, Supp. Mat. Eq. 3 is useful
	# It looks *extremely* conservative, it blows up at e.g.
	# time = 0.00015
	#
	# upper_bound_condition_number: 18425.777249466213
	# 
	# At time t=0.00015, Condition number=2.175764, 1.658813, 1.70285
	# opnorms p=1 & p=2:
	# 32759.807937823098
	# 30274.76762619003
	# 30274.767626190034
	# 17479.145238634184
	#
	# vs.
	# 
	# upper_bound_condition_number: 1.0971658583069032e6
	# 
	# At time t=0.0002, Condition number=3.261165, 2.255482, 2.334215
	# opnorms p=1 & p=2:
	# 34805.07811454515
	# 32164.890247616975
	# 32164.89024761698
	# 18570.408042916435
	# 
	# Probably we could just use opnorm(A,1), or even opnorm(A,1)^2
	
	
	# Matrix norms
	# See: 
	# Lambers, Jim (2009). Vector Norms & Matrix Norms. pp. 1-16.
	# https://www.math.usm.edu/lambers/mat610/sum10/lecture2.pdf
	# These notes have the inequalities between norm forms
	
	# https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/
	# Note: operator norm = matrix norm
	# When p=1, ||A||1, much faster, seems to always be bigger
	# When p=2, ||A||2, the operator norm is the spectral norm, equal to the largest singular value of A
	# When p=Inf, ||A||Inf, this is just opnorm(t(A),1), and 
	#
	# ||A||2 <= sqrt(ncol(A))*||A||Inf
	# ||A||2 <= sqrt(nrow(A))*||A||Inf
	# 
	# (doesn't seem to be true, actually. But cond of opnorm(1) does seem conservative
	
	# Actually, instead of tracking the condition number kappa, they are tracking the 
	# *growth rate* of kappa:
	#
	# // determine linear dynamics (matrix form) of D at various representative 
	# ages, and keep track of the largest estimated exponential growth rate of 
	# the condition number of Gmap ("kappa_rate")
	# 
	# // calculate largest singular value (sigma1) of the dynamics at this age
	#	// then kappa_rate <= 2*sigma1, since 
	# // sigma1(exp(t*A))/sigma2(exp(t*A)) <= exp(t*2*sigma1(A))
	# [So and Thompson (2000) Singular values of Matrix Exponentials. Theorem 2.1]
	# 
	
	# From Louca & Pennell:
	# max_condition_number,	
	# // (INPUT) unitless number, the maximum acceptable condition number for the Gmap 
	# (as estimated from the linearized dynamics), when choosing the integration 
	# interval size. A larger max_condition number leads to fewer age-splits, thus 
	# faster computation but also lower accuracy. Hence, this number controls the 
	# trade-off between speed and accuracy. Typical values are 1e4 (slower, more accurate) 
	# up to 1e8 (faster, less accurate).
	#
	# The comparison is done in:
	# const double DeltaT 		= max(root_age/max_Nintervals, 
	# min(1.0000001*root_age,log(max_condition_number)/max_kappa_rate));
	# 
	# If the maximum observed kappa was 20, and the max condition number
	# was 1e8, log(1e8) would be 18.4
	# 

# Rescaling after Delta_int time interval:
#
# Just take mean of the vector:
#
# inline double vector_mean(const std::vector<double> &values){
# 	double S = 0;
# 	for(long i=0; i<values.size(); ++i) S += values[i];
# 	return (S/values.size());
# }
# 
# 		const double initial_mean = vector_mean(initial);
# 		if(initial_mean<=0) return false;
# 		scale = log(initial_mean);
# 		shape = initial/initial_mean;
# 		return true; 
# 	}
# 
# 	// record a new time series point, provided by the numerical solver
# 	void registerState(double age, const MuSSEstateD &state){
# 		trajectory.push_back(state); 
# 		ages.push_back(age); 
# 
# 		// make sure entries are in [0,1]
# 		const long i = trajectory.size()-1;
# 		for(long s=0; s<trajectory[i].size(); ++s) trajectory[i][s] = max(0.0, min(1.0, trajectory[i][s]));
# 	}
# 	
# 	
# 	// record a new trajectory point, provided by the numerical solver
# 	// The state X to be recorded is provided in rescaled format, i.e. state = exp(scale) * shape
# 	// You can either record shape and scale separately, or combine them to obtain the actual state
# 	void registerScaledState(double age, const MuSSEstateD &shape, const double scale){
# 		trajectory_shape.push_back(shape);
# 		trajectory_scale.push_back(scale);
# 		ages.push_back(age); 






	
	print("\nopnorms(A, p=1, p=2, p=Inf, sqrt(nrow(A))*||A||Inf:\n")
	display(opnorm(A,1))
	display(opnorm(A,2))
	display(opnorm(transpose(A),2))
	display(1/(sqrt(size(A)[1]))*opnorm(transpose(A),2))
	print("\n")
	#display(opnorm(A,Inf))

	sigma1_of_A = opnorm(A,1)
	upper_bound_condition_number = exp(2*t*sigma1_of_A)
	upper_bound_kappa_growth_rate = 2*sigma1_of_A
	tmpstr = paste0(["\nupper_bound_condition_number of A=", string(upper_bound_condition_number), "\nupper_bound_kappa_growth_rate=", string(upper_bound_kappa_growth_rate)])
	print(tmpstr)
	print("\n")
	


	# No return needed, what is returned is G (as .u)
	
	#display(dG)
	#return(dG)
end # End calc_Gs_SSE_condnums



# Calculate G flow matrix down time series
# (no printed outputs)
calc_Gs_SSE = (dG, G, pG, t; max_condition_number=1e8) -> begin
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)

	# The new dG is A %*% G
	mul!(dG, A, G)

	# No return needed, what is returned is G (as .u)
	return(dG)
end # End calc_Gs_SSE




# Doesn't match, risky
calc_Gs_SSE! = (dG, G, pG, t) -> begin
# 	n = pG.n
# 	tmpzero = repeat([0.0], n^2)
# 	A = reshape(tmpzero, (n,n))
	#tmpA = pG.A[:,:]

	p_Ds_v5 = pG.p_Ds_v5
	#A = parameterized_ClaSSE_As_v5(t, tmpA, p_Ds_v5)
	#display(A)
	#dG = A * G
	mul!(dG, parameterized_ClaSSE_As_v5(t, pG.A[:,:], p_Ds_v5), G)
	#display(dG)
	#return(dG)
end # End calc_Gs_SSE



run_Gs = (p_Ds_v5) -> begin
	n = p_Ds_v5.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	G0 = reshape(tmpzero, (n,n))
	for i in 1:n
		G0[i,i] = 1.0
	end
	G = G0
	
	for i in 1:100
		t = 0.01
		A = parameterized_ClaSSE_As_v5(t, A, p_Ds_v5)
		G = A * G
		#Base.showarray(G)
		display(G)
	end

end


# This version includes Xi, Xj in the A equation (but this can't be more efficient I don't think,
# since these will be different on every branch)
parameterized_ClaSSE_A_v5xx = (du,u,A,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	two = 1.0
	# Iterate through the ancestral states
  @inbounds for i in 1:n
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		A[i,i] = A[i,i] + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])#*u[i]  
		
		# case 2: anagenetic change
		@inbounds for m in 1:length(Qi_sub_i)
			A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_sub_i[m]]# * u[Qj_sub_i[m]]
		end
		
		# case 34: change + eventual extinction
		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			rate_sp_then_ex = (u[Ck_sub_i] * uE[Cj_sub_i]) + ( uE[Ck_sub_i])
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + rate_sp_then_ex
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + rate_sp_then_ex
		end
		
# 		du[i] = -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i])*u[i] +  # case 1: no event
# 			(sum(Qij_vals[Qi_sub_i] .* u[Qj_sub_i])) + 	# case 2	
# 			(sum(Cijk_vals[Ci_sub_i] .*                                               # case 34: change + eventual extinction
# 				 (u[Ck_sub_i].*uE[Cj_sub_i] 
# 			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))

  end # End @inbounds for i in 1:n
end


end # end module Flow

