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


end # END parameterized_ClaSSE_As_v6 = (t, A, p)





function sum_Cijk_rates_Ds_inbounds_simd(Cijk_rates_sub_i, tmp_u, tmp_uE, Cj_sub_i, Ck_sub_i; term1=Float64(0.0), term4=Float64(0.0))
    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
    	term4 += Cijk_rates_sub_i[it] * (tmp_u[Ck_sub_i[it]] * tmp_uE[Cj_sub_i[it]] + tmp_u[Cj_sub_i[it]] * tmp_uE[Ck_sub_i[it]])
    end
    return term1, term4
end;


function sum_Qij_vals_inbounds_simd(Qij_vals_sub_i, tmp_u, Qj_sub_i; term2=Float64(0.0), term3=Float64(0.0))
    @inbounds @simd for it=1:length(Qij_vals_sub_i)
    	term2 += Qij_vals_sub_i[it]
    	term3 += Qij_vals_sub_i[it] * tmp_u[Qj_sub_i[it]]
    end
    return term2, term3
end;
