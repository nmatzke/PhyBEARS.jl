	Qij_vals_sub_i = Vector{Vector{Float64}}(undef, n)
	Qij_vals_sub_i_t = Vector{Vector{Float64}}(undef, n)


		Qij_vals_sub_i[i] = Qmat.Qij_vals[Qarray_ivals .== i]	# list of Qij rates lists for anc==i
		Qij_vals_sub_i_t[i] = Qmat.Qij_vals[Qarray_ivals .== i]
		Qji_vals_sub_j[i] = Qmat.Qij_vals[Qarray_jvals .== i]	# list of Qij rates lists for anc==i
		Qji_vals_sub_j_t[i] = Qmat.Qij_vals[Qarray_jvals .== i]
