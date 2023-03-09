"""
include("/Users/nmat471/HD/GitHub/PhyBEARS.jl/notes/parameterized_ClaSSE_As_v7.jl")
"""



parameterized_ClaSSE_As_v7 = (t, A, p) -> begin
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  psi = p.params.psi_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals

	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	# Zero-out A
	A .= 0.0;
	
	# loop through i=states
	@inbounds for i in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		
		Qi_eq_i_index = p.p_TFs.Qi_eq_i_index[i]
		Ci_eq_i_index = p.p_TFs.Ci_eq_i_index[i]
		# For a particular i/state, loop through all of the transitions from that state,
		# for the Q matrix and the C matrix
		# Q matrix
		@inbounds @simd for mi in 1:length(Qi_eq_i_index)
			# looping through mi, with +=, does a sum
			A[i,i] -= Qij_vals[Qi_eq_i_index[mi]] # term2 / part of case1
			# 
			A[Qi_sub_i[mi],Qj_sub_i[mi]] += Qij_vals[Qi_eq_i_index[mi]] # term3 / case2
		end
		
		# C matrix
		@inbounds @simd for mi in 1:length(Ci_eq_i_index)
			# term1, part of case1
			A[i,i] -= Cijk_vals[Ci_eq_i_index[mi]]
			
			# term4, Case 3/4
			A[Ci_sub_i[mi],Cj_sub_i[mi]] += Cijk_vals[Ci_eq_i_index[mi]] * uE[Ck_sub_i[mi]]
			A[Ci_sub_i[mi],Ck_sub_i[mi]] += Cijk_vals[Ci_eq_i_index[mi]] * uE[Cj_sub_i[mi]]
		end		
		
		# Additional parts of term1
		A[i,i] -= p.params.mu_vals[i]
		A[i,i] -= p.params.psi_vals[i]
	end

end # END parameterized_ClaSSE_As_v7 = (t, A, p)


parameterized_ClaSSE_As_v7_simd = (t, A, p) -> begin
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  psi = p.params.psi_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals

	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	# Zero-out A
	A .= 0.0;
	
	# loop through i=states
	@inbounds for i in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		
		Qi_eq_i_index = p.p_TFs.Qi_eq_i_index[i]
		Ci_eq_i_index = p.p_TFs.Ci_eq_i_index[i]
		# For a particular i/state, loop through all of the transitions from that state,
		# for the Q matrix and the C matrix
		# Q matrix
		sum_Qij_vals_inbounds_simd_A!(A, Qij_vals, Qi_sub_i, Qj_sub_i, Qi_eq_i_index)
		
		# C matrix
		sum_Cijk_vals_inbounds_simd_A!(A, Cijk_vals, Ci_eq_i_index, Ci_sub_i, Cj_sub_i, Ck_sub_i, uE)		
		
		# Additional parts of term1
		A[i,i] -= p.params.mu_vals[i]
		A[i,i] -= p.params.psi_vals[i]
	end

end # END parameterized_ClaSSE_As_v7 = (t, A, p)





function sum_Qij_vals_inbounds_simd_A!(A, Qij_vals, Qi_sub_i, Qj_sub_i, Qi_eq_i_index)
	@inbounds @simd for mi in 1:length(Qi_eq_i_index)
		# looping through mi, with +=, does a sum
		A[i,i] -= Qij_vals[Qi_eq_i_index[mi]] # term2 / part of case1
		# 
		A[Qi_sub_i[mi],Qj_sub_i[mi]] += Qij_vals[Qi_eq_i_index[mi]] # term3 / case2
	end
	return A
end;



function sum_Cijk_vals_inbounds_simd_A!(A, Cijk_vals, Ci_eq_i_index, Ci_sub_i, Cj_sub_i, Ck_sub_i, uE)
	# C matrix
	@inbounds @simd for mi in 1:length(Ci_eq_i_index)
		# term1, part of case1
		A[i,i] -= Cijk_vals[Ci_eq_i_index[mi]]
		
		# term4, Case 3/4
		A[Ci_sub_i[mi],Cj_sub_i[mi]] += Cijk_vals[Ci_eq_i_index[mi]] * uE[Ck_sub_i[mi]]
		A[Ci_sub_i[mi],Ck_sub_i[mi]] += Cijk_vals[Ci_eq_i_index[mi]] * uE[Cj_sub_i[mi]]
	end		
	return A
end;
