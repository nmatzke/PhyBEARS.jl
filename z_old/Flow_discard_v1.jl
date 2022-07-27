#######################################################
# These functions use
# parameterized_ClaSSE_As_v7
#
# But, this seems to need Cijk_rates_sub_i, which has
# to be updated after every update of p_Ds_v5. It is 
# simpler to just update 1 params vector and then use the
# Ci_eq_i saved True/False variable!

# If we
# just use calc_Gs_SSE and parameterized_ClaSSE_As_v6,
# we get matching, even on half-matrix Cmodels. So just do that.
# 
#######################################################

export calc_Gs_SSE_v7, calc_Gs_SSE_v7!, parameterized_ClaSSE_As_v7

# These are wholly redundant with calc_Gs_SSE (i.e. v5)

# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
calc_Gs_SSE_v7! = (dG, G, pG, t) -> begin
	# Have to use pG.A[:,:] to avoid overwriting A, which screws everything up!
	#A = parameterized_ClaSSE_As_v5(t, pG.A[:,:], pG.p_Ds_v5)
	#display(A)
	#dG = A * G
	# Have to use pG.A[:,:] to avoid overwriting A, which screws everything up!
	#mul!(dG, A(t), G(t))  # Equation 20
	mul!(dG, parameterized_ClaSSE_As_v7(t, pG.A[:,:], pG.p_Ds_v5), G)
	#display(dG)
	#return(dG)
end # End calc_Gs_SSE


# Calculate G flow matrix down time series
# (no printed outputs)
# For this one, the Ds and Es are passed
calc_Gs_SSE_v7 = (dG, G, pG, t; max_condition_number=1e8) -> begin
	n = pG.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v7(t, A, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)

	# The new dG is A %*% G
	mul!(dG, A, G)

	# No return needed, what is returned is G (as .u)
	return(dG)
end # End calc_Gs_SSE




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
parameterized_ClaSSE_As_v7_NEEDS_UPDATER_FOR_Cijk_rates_sub_i_BUT_As_v6_WORKS_FINE(t, A, p; max_condition_number=1e8, print_warnings=true) -> begin

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
	# (u = c(uE, uD))
	
	two = 1.0
	# Iterate through the ancestral states
  @inbounds for i in 1:n
		# These are the (e.g.) j state-indices (left descendant) when the ancestor==state i
		#	Ci_sub_i = Any[]
		#	Cj_sub_i = Any[]
		#	Ck_sub_i = Any[]

		# The push! operation may get slow at huge n
		# This will have to change for non-Mk models
		#	for i in 1:n
		#		push!(Qi_eq_i, Qarray_ivals .== i)
		#		push!(Qi_sub_i, Qarray_ivals[Qarray_ivals .== i])
		#		push!(Qj_sub_i, Qarray_jvals[Qarray_ivals .== i])

		#		push!(Ci_eq_i, Carray_ivals .== i)
		#		push!(Ci_sub_i, Carray_ivals[Carray_ivals .== i])
		#		push!(Cj_sub_i, Carray_jvals[Carray_ivals .== i])
		#		push!(Ck_sub_i, Carray_kvals[Carray_ivals .== i])
		#	end

		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also)
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i[i]  # The Qij rates for anc==i
		
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		Cijk_pair_sub_i = p.p_TFs.Cijk_pair_sub_i[i] # The Cijk pairs (1 or 2 eventsX for anc==i
		Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i[i] # The Cijk rates for anc==i

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		#A[i,i] = A[i,i]  + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i]) # *u[i]  
		A[i,i] = A[i,i] + (-sum(Qij_vals_sub_i) - sum(Cijk_rates_sub_i) - mu[i]) #*uD[i] #+ 2*sum(Cijk_vals[Ci_sub_i])*uE[i]
		
		# case 2: anagenetic change (non-diagonal)
		@inbounds for m in 1:length(Qi_sub_i)
			#A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_sub_i[m]] #* uD[Qj_sub_i[m]]
			A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals_sub_i[m] #* uD[Qj_sub_i[m]]
		end
		
		# case 34: change + eventual extinction (non-diagonal)
		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			# excluding the u[], i.e. the Ds, i.e. the Xs, just as is done in 
			# 2*speciation_rates[r]*current_E[r]

			# case 34: change + eventual extinction
			#rate_sp_then_ex1 = Cijk_rates_sub_i[m] * uE[Ck_sub_i[m]] # rate of going to j, times k going extinct
			#rate_sp_then_ex2 = Cijk_rates_sub_i[m] * uE[Cj_sub_i[m]] # rate of going to k, times j going extinct
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + Cijk_rates_sub_i[m]/Cijk_pair_sub_i[m] * uE[Ck_sub_i[m]]
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + (Cijk_rates_sub_i[m]/Cijk_pair_sub_i[m] * (Cijk_pair_sub_i[m]-1.0)) * uE[Cj_sub_i[m]]
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + Cijk_rates_sub_i[m]/Cijk_pair_sub_i[m] * uE[Cj_sub_i[m]]
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + (Cijk_rates_sub_i[m]/Cijk_pair_sub_i[m] * (Cijk_pair_sub_i[m]-1.0)) * uE[Ck_sub_i[m]]
			# <== WORKS 2022-03-24

		end
  end # End @inbounds for i in 1:n
	
	# Error check; large growth rate in the condition number
	# (indicated by the norm of the "A" matrix giving the linear dynamics)
	# suggests that the G matrix is approaching singularity and numerical
	# errors will accumulate.
	# Check the maximum condition number of A; if it is exceeded, either raise the max cond number,
	# or divide the tree into smaller chunks
	# (after Louca & Pennell)
	if (print_warnings == true)
		sigma1_of_A = opnorm(A,1)  # the 1-norm should be adequate here (fastest calc.)
		if (2*sigma1_of_A > log(max_condition_number))
			warning_txt = join(["WARNING in parameterized_ClaSSE_As_v5 at t=", string(t), ": 2*opnorm(A,1)>log(max_condition_number) (", string(round(sigma1_of_A; digits=2)), " > ", string(round(log(max_condition_number); digits=2)), ")\n"], "")
			display(warning_txt)
		end # END if (2*sigma1_of_A > log(max_condition_number))
	end # END if (print_warnings == true)
 	return(A)
end # END parameterized_ClaSSE_As_v7


#######################################################
# parameterized_ClaSSE_As_v7
#######################################################
# Construct interpolation function for calculating linear dynamics A, 
# at any timepoint t
#
#	This version is more complex, but actually reduces down to 
# parameterized_ClaSSE_As_v6
#
#######################################################
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
function parameterized_ClaSSE_As_v7(t, A, p; max_condition_number=1e8, print_warnings=true)

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
	# (u = c(uE, uD))
	
	two = 1.0
	# Iterate through the ancestral states
  @inbounds for i in 1:n
		# These are the (e.g.) j state-indices (left descendant) when the ancestor==state i
		#	Ci_sub_i = Any[]
		#	Cj_sub_i = Any[]
		#	Ck_sub_i = Any[]

		# The push! operation may get slow at huge n
		# This will have to change for non-Mk models
		#	for i in 1:n
		#		push!(Qi_eq_i, Qarray_ivals .== i)
		#		push!(Qi_sub_i, Qarray_ivals[Qarray_ivals .== i])
		#		push!(Qj_sub_i, Qarray_jvals[Qarray_ivals .== i])

		#		push!(Ci_eq_i, Carray_ivals .== i)
		#		push!(Ci_sub_i, Carray_ivals[Carray_ivals .== i])
		#		push!(Cj_sub_i, Carray_jvals[Carray_ivals .== i])
		#		push!(Ck_sub_i, Carray_kvals[Carray_ivals .== i])
		#	end

		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also)
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		Cijk_pair_sub_i = p.p_TFs.Cijk_pair_sub_i[i] # The Cijk pairs (1 or 2 eventsX for anc==i

		# Calculation of "A" (the from-to probabilities between every pair of states)
		# Pull out the Q transitions - diagonal
		# case 1: no event
		#A[i,i] = A[i,i]  + -(sum(Cijk_vals[Ci_sub_i]) + sum(Qij_vals[Qi_sub_i]) + mu[i]) # *u[i]  
		A[i,i] = A[i,i] + (-sum(Qij_vals[Qi_eq_i]) - sum(Cijk_vals[Ci_eq_i]) - mu[i]) #*uD[i] #+ 2*sum(Cijk_vals[Ci_sub_i])*uE[i]
		
		# case 2: anagenetic change (non-diagonal)
		@inbounds for m in 1:length(Qi_sub_i)
			#A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_sub_i[m]] #* uD[Qj_sub_i[m]]
			A[Qi_sub_i[m],Qj_sub_i[m]] = A[Qi_sub_i[m],Qj_sub_i[m]] + Qij_vals[Qi_eq_i][m] #* uD[Qj_sub_i[m]]
		end
		
		# case 34: change + eventual extinction (non-diagonal)
		@inbounds for m in 1:length(Ci_sub_i)
			# each cladogenesis event puts probability in 2 places
			# excluding the u[], i.e. the Ds, i.e. the Xs, just as is done in 
			# 2*speciation_rates[r]*current_E[r]

			# case 34: change + eventual extinction
			#rate_sp_then_ex1 = Cijk_vals[Ci_eq_i][m] * uE[Ck_sub_i[m]] # rate of going to j, times k going extinct
			#rate_sp_then_ex2 = Cijk_vals[Ci_eq_i][m] * uE[Cj_sub_i[m]] # rate of going to k, times j going extinct
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + Cijk_vals[Ci_eq_i][m]/Cijk_pair_sub_i[m] * uE[Ck_sub_i[m]]
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + (Cijk_vals[Ci_eq_i][m]/Cijk_pair_sub_i[m] * (Cijk_pair_sub_i[m]-1.0)) * uE[Cj_sub_i[m]]
			A[Ci_sub_i[m],Ck_sub_i[m]] = A[Ci_sub_i[m],Ck_sub_i[m]] + Cijk_vals[Ci_eq_i][m]/Cijk_pair_sub_i[m] * uE[Cj_sub_i[m]]
			A[Ci_sub_i[m],Cj_sub_i[m]] = A[Ci_sub_i[m],Cj_sub_i[m]] + (Cijk_vals[Ci_eq_i][m]/Cijk_pair_sub_i[m] * (Cijk_pair_sub_i[m]-1.0)) * uE[Ck_sub_i[m]]
			# <== WORKS 2022-03-24

		end
  end # End @inbounds for i in 1:n
	
	# Error check; large growth rate in the condition number
	# (indicated by the norm of the "A" matrix giving the linear dynamics)
	# suggests that the G matrix is approaching singularity and numerical
	# errors will accumulate.
	# Check the maximum condition number of A; if it is exceeded, either raise the max cond number,
	# or divide the tree into smaller chunks
	# (after Louca & Pennell)
	if (print_warnings == true)
		sigma1_of_A = opnorm(A,1)  # the 1-norm should be adequate here (fastest calc.)
		if (2*sigma1_of_A > log(max_condition_number))
			warning_txt = join(["WARNING in parameterized_ClaSSE_As_v5 at t=", string(t), ": 2*opnorm(A,1)>log(max_condition_number) (", string(round(sigma1_of_A; digits=2)), " > ", string(round(log(max_condition_number); digits=2)), ")\n"], "")
			display(warning_txt)
		end # END if (2*sigma1_of_A > log(max_condition_number))
	end # END if (print_warnings == true)
 	return(A)
end # END parameterized_ClaSSE_As_v7




# Calculate G flow matrix down time series
# (no printed outputs)
# For this one, the Ds and Es are passed
calc_Gs_SSE_v7 = (dG, G, pG, t; max_condition_number=1e8) -> begin
	n = pG.n
	tmpzero = repeat([0.0], n^2)
	A = reshape(tmpzero, (n,n))

	p_Ds_v5 = pG.p_Ds_v5
	A = parameterized_ClaSSE_As_v7(t, A, p_Ds_v5)
	#display(A)
	#dG = A * G
	#display(G)

	# The new dG is A %*% G
	mul!(dG, A, G)

	# No return needed, what is returned is G (as .u)
	return(dG)
end # End calc_Gs_SSE




