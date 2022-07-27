

# Plain version
Ds_v5_tmp = (du,u,p,t) -> begin
#	let p=p;

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
	
  @inbounds for i in 1:n
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		du[i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 	# case 2	
			(sum(Cijk_vals[Ci_eq_i] .*                                               
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
 # end

	end # END let p=p;
end


# @simd-enhanced version; 10 times faster on a single core
function Ds_v7_simd_sums(du,u,p,t)
#	let p=p;
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
	uE = p.sol_Es_v5(t)
	terms = p.terms
  for i in 1:n
		terms[1] = 0.0
		terms[2] = 0.0
		terms[3] = 0.0
		terms[4] = 0.0

		terms[1], terms[4] = sum_Cijk_rates_Ds_inbounds_simd_tmp(p.p_TFs.Cijk_rates_sub_i[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i])
	
		terms[2], terms[3] = sum_Qij_vals_inbounds_simd_tmp(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i])
		
		du[i] = -(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4]
  end
#	end # END let p=p;
end




# Inner-loop functions
function sum_Cijk_rates_Ds_inbounds_simd_tmp(Cijk_rates_sub_i, tmp_u, tmp_uE, Cj_sub_i, Ck_sub_i)
    term1 = 0.0
    term4 = 0.0
    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
    	term4 += Cijk_rates_sub_i[it] * (tmp_u[Ck_sub_i[it]] * tmp_uE[Cj_sub_i[it]] + tmp_u[Cj_sub_i[it]] * tmp_uE[Ck_sub_i[it]])
    end
    return term1, term4
end;

function sum_Cijk_rates_Es_inbounds_simd_tmp(Cijk_rates_sub_i, tmp_u, Cj_sub_i, Ck_sub_i)
    term1 = 0.0
    term4 = 0.0
    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
    	term4 += Cijk_rates_sub_i[it] * tmp_u[Cj_sub_i[it]] * tmp_u[Ck_sub_i[it]]
    end
    return term1, term4
end;

function sum_Qij_vals_inbounds_simd_tmp(Qij_vals_sub_i, tmp_u, Qj_sub_i)
    term2 = 0.0
    term3 = 0.0
    @inbounds @simd for it=1:length(Qij_vals_sub_i)
    	term2 += Qij_vals_sub_i[it]
    	term3 += Qij_vals_sub_i[it] * tmp_u[Qj_sub_i[it]]
    end
    return term2, term3
end;


# Only needs to be calculated once, then interpolator re-used by Ds_v7_simd_sums via p.sol_Es_v5(t)
Es_v7_simd_sums = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
	terms = p.terms
  for i in 1:n
		terms[1] = 0.0
		terms[2] = 0.0
		terms[3] = 0.0
		terms[4] = 0.0

		terms[1], terms[4] = sum_Cijk_rates_Es_inbounds_simd_tmp(p.p_TFs.Cijk_rates_sub_i[i], u, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i])
	
		terms[2], terms[3] = sum_Qij_vals_inbounds_simd_tmp(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i])
		
		du[i] = mu[i] -(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4]
  end
end





