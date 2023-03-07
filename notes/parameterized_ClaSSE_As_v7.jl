parameterized_ClaSSE_As_v7 = (t, A, p)
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  psi = p.params.psi_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals

	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v7
	uE = p.uE
	uE = sol_Es(t)
	
	p.terms
	
	@inbounds for st in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[st]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[st]	# To get the is for Qi_eq_i (somehow this works for Qijs also
		Qj_sub_i = p.p_TFs.Qj_sub_i[st]
		Ci_eq_i = p.p_TFs.Ci_eq_i[st]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[st]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[st]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[st]	# To get the is for Ci_eq_i
		
		Qi_eq_i_index = p.p_TFs.Qi_eq_i_index[st]
		
		# Q loop
		@inbounds for mi in 1:length(Qi_eq_i_index)
			# looping through mi sums 
			A[i,i] += Qij_vals[Qi_eq_i_index[mi]] 
			# 
			A[Qi_sub_i[mi],Qj_sub_i[mi]] += Qij_vals[Qi_eq_i_index[mi]][mi]
		end
		A[i,i] -= p.params.mu_vals[i]
		A[i,i] -= p.params.psi_vals[i]
		
		# C loop
		
		
		
	end

end # END parameterized_ClaSSE_As_v6 = (t, A, p)





function sum_Qij_vals_inbounds_simd(Qij_vals_sub_i, tmp_u, Qj_sub_i; term2=Float64(0.0), term3=Float64(0.0))
    @inbounds @simd for it=1:length(Qij_vals_sub_i)
    	term2 += Qij_vals_sub_i[it]
    	term3 += Qij_vals_sub_i[it] * tmp_u[Qj_sub_i[it]]
    end
    return term2, term3
end;


function sum_Cijk_rates_Ds_inbounds_simd(Cijk_rates_sub_i, tmp_u, tmp_uE, Cj_sub_i, Ck_sub_i; term1=Float64(0.0), term4=Float64(0.0))
    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
    	term4 += Cijk_rates_sub_i[it] * (tmp_u[Ck_sub_i[it]] * tmp_uE[Cj_sub_i[it]] + tmp_u[Cj_sub_i[it]] * tmp_uE[Ck_sub_i[it]])
    end
    return term1, term4
end;
