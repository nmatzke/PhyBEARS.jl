module SSEs
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

print("PhyBEARS: loading SSEs.jl dependencies...")
using DataFrames  # for e.g. DataFrame()
using PhyloBits.TrUtils # for get_area_of_range etc.
using PhyBEARS.TimeDep # for get_area_of_range etc.

print("...done.\n")


export parameterized_ClaSSE, parameterized_ClaSSE_Es, parameterized_ClaSSE_Ds, parameterized_ClaSSE_v5, parameterized_ClaSSE_Es_v5, parameterized_ClaSSE_Ds_v5, parameterized_ClaSSE_Es_v5orig, parameterized_ClaSSE_Ds_v5orig, parameterized_ClaSSE_Es_v6, parameterized_ClaSSE_Ds_v6, parameterized_ClaSSE_Es_v6orig, parameterized_ClaSSE_Ds_v6orig, parameterized_ClaSSE_Es_v7orig, parameterized_ClaSSE_Ds_v7orig, parameterized_ClaSSE_Ds_v7_forloops_sucks, parameterized_ClaSSE_Es_v7_simd_sums, parameterized_ClaSSE_Ds_v7_simd_sums, parameterized_ClaSSE_Es_v10_simd_sums, parameterized_ClaSSE_Ds_v10_simd_sums, parameterized_ClaSSE_Es_v11_simd_sums, parameterized_ClaSSE_Ds_v11_simd_sums, parameterized_ClaSSE_Es_v12_simd_sums, parameterized_ClaSSE_Ds_v12_simd_sums, sum_range_inbounds_simd, sum_Cijk_rates_Ds_inbounds_simd, sum_Cijk_rates_Es_inbounds_simd, sum_Qij_vals_inbounds_simd



parameterized_ClaSSE_Ds_v5orig_return_du = (du,u,p,t) -> begin

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
  @inbounds for i in 1:n
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 	# case 2	
			(sum(Cijk_vals[Ci_eq_i] .*                                               # case 3/4: change + eventual extinction
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
  end
	return du
end








#######################################################
# Various versions of ClaSSE equations
# trying for more and more efficient
#######################################################


# Load g3 -- 100 states
#include("/drives/GDrive/__GDrive_projects/2018-01-22_Marsden/software/Julia_EPIRK/ClaSSE_25states_v1.jl")
#include("/drives/GDrive/__GDrive_projects/2018-01-22_Marsden/software/Julia_EPIRK/ClaSSE_25states_v1.jl")

# Seems to be the same as
# function parameterized_ClaSSE!(du,u,p,t)
# (in-place function with "!")
parameterized_ClaSSE = (du,u,p,t) -> begin

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
	
	two = 1.0
  @inbounds for i in 1:n
		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qarray_ivals .== i] .* u[Qarray_jvals[Qarray_ivals .== i]])) + 			# case 3: change + eventual extinction
			(two * sum(Cijk_vals[Carray_ivals .== i] .* u[Carray_jvals[Carray_ivals .== i]] .* u[Carray_kvals[Carray_ivals .== i]])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required

		# Calculation of "D" (likelihood of tip data)
		du[n+i] = -(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[n+i] +  # case 1: no event
			(sum(Qij_vals[Qarray_ivals .== i] .* u[(n.+Qarray_jvals[Qarray_ivals .== i])])) + 	# case 2	
			(sum(Cijk_vals[Carray_ivals .== i] .*                                               # case 34: change + eventual extinction
				 (u[(n.+Carray_kvals[Carray_ivals .== i])].*u[Carray_jvals[Carray_ivals .== i]] 
			 .+ u[(n.+Carray_jvals[Carray_ivals .== i])].*u[Carray_kvals[Carray_ivals .== i]]) ))
  end
end




"""
This function describes a system of ODEs for calculating the "E" part of the 
likelihood calculation (probability of a lineage going 
between time t and the present) of a ClaSSE model

du = time increment
u = vector of E values (of length numstates)
p = structure containing parameters
t = time

This function is input into `ODEProblem()`, along with:

u0    = starting vector of Es (typically, E=0 for all states, as any lineage alive at 
        t=0 mya has a 0 probability of becoming extinct at 0 mya.
tspan = vector of start- and end-times, e.g. tspan=(0.0, 10.0)
p     = structure containing parameters

# Example

```julia-repl

# Initialize evolving state vector
u = repeat([0.0], 2*n)

# Starting values
du = repeat([0.0], 2*n)
u0 = repeat([0.0], 2*n)
u0[n+1] = 1.0
tspan = (0.0, 10.0)

prob = ODEProblem(parameterized_ClaSSE, u0, tspan, p)
sol = solve(prob, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
```
"""
parameterized_ClaSSE_Es = (du,u,p,t) -> begin

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
	
	two = 1.0
  @inbounds for i in 1:n
		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qarray_ivals .== i] .* u[Qarray_jvals[Qarray_ivals .== i]])) + 			# case 3: change + eventual extinction
			(two * sum(Cijk_vals[Carray_ivals .== i] .* u[Carray_jvals[Carray_ivals .== i]] .* u[Carray_kvals[Carray_ivals .== i]])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required
  end
end

parameterized_ClaSSE_Ds = (du,u,p,t) -> begin

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
	sol_Es = p.sol_Es
	uE = p.uE
	uE = sol_Es(t)
	
	two = 1.0
  @inbounds for i in 1:n
		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qarray_ivals .== i] .* u[(Qarray_jvals[Qarray_ivals .== i])])) + 	# case 2	
			(sum(Cijk_vals[Carray_ivals .== i] .*                                               # case 34: change + eventual extinction
				 (u[(Carray_kvals[Carray_ivals .== i])].*uE[Carray_jvals[Carray_ivals .== i]] 
			 .+ u[(Carray_jvals[Carray_ivals .== i])].*uE[Carray_kvals[Carray_ivals .== i]]) ))
  end
end




 
 
# Modeled on:
# https://gitter.im/JuliaDiffEq/Lobby?at=5cbe49dbb4700e023db73d9e
# https://gist.github.com/nmatzke/b6332845747d7452b6e3b45564f460e8
# Pre-allocating the Carray_ivals .== i, Qarray_jvals[Qarray_ivals .== i
# Reduces GC (Garbage Collection) from 40% to ~5%
parameterized_ClaSSE_v5 = (du,u,p,t) -> begin

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
	
	two = 1.0
  @inbounds for i in 1:n
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
			(two * sum(Cijk_vals[Ci_eq_i] .* u[Cj_sub_i] .* u[Ck_sub_i])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required

		# Calculation of "D" (likelihood of tip data)
		du[n+i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[n+i] +  # case 1: no event
			(sum(Qij_vals[Qi_eq_i] .* u[(n.+Qj_sub_i)])) + 	# case 2	
			(sum(Cijk_vals[Ci_eq_i] .*                                               # case 34: change + eventual extinction
				 (u[(n.+Ck_sub_i)].*u[Cj_sub_i] 
			 .+ u[(n.+Cj_sub_i)].*u[Ck_sub_i]) ))
  end
end
 
 
parameterized_ClaSSE_Es_v5 = (du,u,p,t) -> begin

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
	
	two = 1.0
  @inbounds for i in 1:n
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
			(two * sum(Cijk_vals[Ci_eq_i] .* u[Cj_sub_i] .* u[Ck_sub_i])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required
  end
end

parameterized_ClaSSE_Ds_v5 = (du,u,p,t) -> begin

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
  @inbounds for i in 1:n
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 	# case 2	
			(sum(Cijk_vals[Ci_eq_i] .*                                               # case 3/4: change + eventual extinction
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
  end
end



# Assumes a C matrix where i,j,k=1,1,2 but not 1,2,1; such values will have doubled rates
# (These calcs match v6, which is more complex/explicit)
parameterized_ClaSSE_Es_v5orig = (du,u,p,t) -> begin

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
	
	two = 1.0
  @inbounds for i in 1:n
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
			(two * sum(Cijk_vals[Ci_eq_i] .* u[Cj_sub_i] .* u[Ck_sub_i])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required, if rates are doubled
  end
end




parameterized_ClaSSE_Es_v6 = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
  #Cijk_rates = p.params.Cijk_rates
  #Cijk_pair = p.params.Cijk_pair
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
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
	
  @inbounds for i in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also)
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i[i]  # The Qij rates for anc==i

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		Cijk_not_y_sub_i = p.p_TFs.Cijk_not_y_sub_i[i]
		Cijk_pair_sub_i = p.p_TFs.Cijk_pair_sub_i[i] # The Cijk pairs (1 or 2 events for anc==i
		Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i[i] # The Cijk rates for anc==i


		# Calculation of "E" (prob. of extinction)
		# Old:
		#	-(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
		#	(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
		#	(two * sum(Cijk_vals[Ci_eq_i] .* u[Cj_sub_i] .* u[Ck_sub_i])) 

		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
			sum((Cijk_vals[Ci_eq_i]./ Cijk_pair_sub_i) .* u[Cj_sub_i] .* u[Ck_sub_i]) + # case 4: both go extinct
			sum(((Cijk_vals[Ci_eq_i]./ Cijk_pair_sub_i) .* u[Cj_sub_i] .* u[Ck_sub_i]) .* (Cijk_pair_sub_i.-1.0)) # case 4, reversed
			
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required
  end
end

parameterized_ClaSSE_Ds_v6 = (du,u,p,t) -> begin

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
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also)
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i[i]  # The Qij rates for anc==i

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		Cijk_pair_sub_i = p.p_TFs.Cijk_pair_sub_i[i] # The Cijk pairs (1 or 2 events for anc==i
		Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i[i] # The Cijk rates for anc==i

		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 	# case 2	
			(sum((Cijk_vals[Ci_eq_i] ./ Cijk_pair_sub_i) .*                                               # case 3/4: change + eventual extinction
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]  
			 .+ u[Ck_sub_i].*uE[Cj_sub_i] .* (Cijk_pair_sub_i .- 1.0) # reversed for pair=2 events
			 .+ u[Cj_sub_i].*uE[Ck_sub_i] .* (Cijk_pair_sub_i .- 1.0) )
			 ))
  end
end

# v6 works for both:
# * a C matrix where i,j,k=1,1,2 but not 1,2,1; such values will have doubled rates
# * a C matrix where i,j,k=1,1,2 and 1,2,1; such values will have non-doubled rates
parameterized_ClaSSE_Es_v6orig = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
  #Cijk_rates = p.params.Cijk_rates
  #Cijk_pair = p.params.Cijk_pair
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
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
	
  @inbounds for i in 1:n
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also)
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i[i]  # The Qij rates for anc==i

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		Cijk_not_y_sub_i = p.p_TFs.Cijk_not_y_sub_i[i]
		Cijk_pair_sub_i = p.p_TFs.Cijk_pair_sub_i[i] # The Cijk pairs (1 or 2 events for anc==i
		Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i[i] # The Cijk rates for anc==i


		# Calculation of "E" (prob. of extinction)
		# Old:
		#	-(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
		#	(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
		#	(two * sum(Cijk_vals[Ci_eq_i] .* u[Cj_sub_i] .* u[Ck_sub_i])) 

		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
			sum((Cijk_vals[Ci_eq_i]./ Cijk_pair_sub_i) .* u[Cj_sub_i] .* u[Ck_sub_i]) + # case 4: both go extinct
			sum(((Cijk_vals[Ci_eq_i]./ Cijk_pair_sub_i) .* u[Cj_sub_i] .* u[Ck_sub_i]) .* (Cijk_pair_sub_i.-1.0)) # case 4, reversed
			
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required
  end
end

# v6 works for both:
# * a C matrix where i,j,k=1,1,2 but not 1,2,1; such values will have doubled rates
# * a C matrix where i,j,k=1,1,2 and 1,2,1; such values will have non-doubled rates
parameterized_ClaSSE_Ds_v6orig = (du,u,p,t) -> begin

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
		Qi_eq_i = p.p_TFs.Qi_eq_i[i]		# Use to get the Qij_vals for Qi_eq_i
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]	# To get the is for Qi_eq_i (somehow this works for Qijs also)
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i[i]  # The Qij rates for anc==i

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]		# Use to get the lambdas/Cijk_vals for Ci_eq-I: USE:  Cijk_vals[Ci_eq_i] not  Cijk_vals[Ci_sub_i]
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]  # To get the is for Ci_eq_i
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]	# To get the js for Ci_eq_i
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]	# To get the is for Ci_eq_i
		Cijk_pair_sub_i = p.p_TFs.Cijk_pair_sub_i[i] # The Cijk pairs (1 or 2 events for anc==i
		Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i[i] # The Cijk rates for anc==i

		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 	# case 2	
			(sum((Cijk_vals[Ci_eq_i] ./ Cijk_pair_sub_i) .*                                               # case 3/4: change + eventual extinction
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]  
			 .+ u[Ck_sub_i].*uE[Cj_sub_i] .* (Cijk_pair_sub_i .- 1.0) # reversed for pair=2 events
			 .+ u[Cj_sub_i].*uE[Ck_sub_i] .* (Cijk_pair_sub_i .- 1.0) )
			 ))
  end
end




# Assumes a C matrix where i,j,k=1,1,2 but not 1,2,1; such values will have doubled rates
# (These calcs match v6, which is more complex/explicit)
parameterized_ClaSSE_Es_v7orig = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
#  Qij_vals = p.params.Qij_vals
#  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
#	Qarray_ivals = p.p_indices.Qarray_ivals
#	Qarray_jvals = p.p_indices.Qarray_jvals
#	Carray_ivals = p.p_indices.Carray_ivals
#	Carray_jvals = p.p_indices.Carray_jvals
#	Carray_kvals = p.p_indices.Carray_kvals

	Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i
	Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i
	
  @inbounds for i in 1:n
#		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
#		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
#		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
#		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "E" (prob. of extinction)
		du[i] = mu[i] +                                         # case 1: lineage extinction
			-(sum(Cijk_rates_sub_i[i]) + sum(Qij_vals_sub_i[i]) + mu[i])*u[i] +  # case 2: no event + eventual extinction
			(sum(Qij_vals_sub_i[i] .* u[Qj_sub_i])) + 			# case 3: change + eventual extinction
			(sum(Cijk_rates_sub_i[i] .* u[Cj_sub_i] .* u[Ck_sub_i])) 
			# case 4 & 5: speciation from i producing j,k or k,j, eventually both daughters go extinct
			# Because Carray contains all nonzero i->j,k rates, iterating once through on a particular
			# "i" does the double-summation required, if rates are doubled
  end
end

# Assumes a C matrix where i,j,k=1,1,2 but not 1,2,1; such values will have doubled rates
# (These calcs match v6, which is more complex/explicit)
parameterized_ClaSSE_Ds_v7orig = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
#  Qij_vals = p.params.Qij_vals
#  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
#	Qarray_ivals = p.p_indices.Qarray_ivals
#	Qarray_jvals = p.p_indices.Qarray_jvals
#	Carray_ivals = p.p_indices.Carray_ivals
#	Carray_jvals = p.p_indices.Carray_jvals
#	Carray_kvals = p.p_indices.Carray_kvals

	Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i
	Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i
	
	# Pre-calculated solution of the Es
#	sol_Es = p.sol_Es_v5
#	uE = p.uE
	uE = p.sol_Es_v5(t)
	
  @inbounds for i in 1:n
		#Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		#Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		#Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		#Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_rates_sub_i[i]) + sum(Qij_vals_sub_i[i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals_sub_i[i] .* u[Qj_sub_i])) + 	# case 2	
			(sum(Cijk_rates_sub_i[i] .*                                               # case 3/4: change + eventual extinction
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
  end
end


# Assumes a C matrix where i,j,k=1,1,2 but not 1,2,1; such values will have doubled rates
# (These calcs match v6, which is more complex/explicit)

parameterized_ClaSSE_Ds_v5orig = (du,u,p,t) -> begin

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
  @inbounds for i in 1:n
		Qi_sub_i = p.p_TFs.Qi_sub_i[i]
		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Ci_eq_i]) + sum(Qij_vals[Qi_eq_i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qi_eq_i] .* u[Qj_sub_i])) + 	# case 2	
			(sum(Cijk_vals[Ci_eq_i] .*                                               # case 3/4: change + eventual extinction
				 (u[Ck_sub_i].*uE[Cj_sub_i] 
			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
  end
end





# Takes 2.5 seconds, still sucks
parameterized_ClaSSE_Ds_v7_forloops_sucks = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
#  Qij_vals = p.params.Qij_vals
#  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
#	Qarray_ivals = p.p_indices.Qarray_ivals
#	Qarray_jvals = p.p_indices.Qarray_jvals
#	Carray_ivals = p.p_indices.Carray_ivals
#	Carray_jvals = p.p_indices.Carray_jvals
#	Carray_kvals = p.p_indices.Carray_kvals

#	Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i
#	Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i
	
	# Pre-calculated solution of the Es
#	sol_Es = p.sol_Es_v5
#	uE = p.uE
	uE = p.sol_Es_v5(t)
	#terms = p.terms
	terms = Vector{Float64}(undef, 4)
  for i in 1:n
		#Qi_sub_i = p.p_TFs.Qi_sub_i[i]
#		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		#Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		#Ci_sub_i = p.p_TFs.Ci_sub_i[i]
#		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
#		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		#Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "D" (likelihood of tip data)
#		du[i] = -(sum(Cijk_rates_sub_i[i]) + sum(Qij_vals_sub_i[i]) + mu[i])*u[i] +  # case 1: no event
#			(sum(Qij_vals_sub_i[i] .* u[Qj_sub_i])) + 	# case 2	
#			(sum(Cijk_rates_sub_i[i] .*                                               # case 3/4: change + eventual extinction
#				 (u[Ck_sub_i].*uE[Cj_sub_i] 
#			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
		
#		Cijk_rates_temp = @view p.p_TFs.Cijk_rates_sub_i[i]
#		Qij_rates_temp = @view p.p_TFs.Qij_vals_sub_i[i]
		
		terms[1] = 0.0
		terms[2] = 0.0
		terms[3] = 0.0
		terms[4] = 0.0
		@inbounds for (it,tmp_Cijk_rate) in enumerate(p.p_TFs.Cijk_rates_sub_i[i])
			terms[1] += tmp_Cijk_rate
			terms[4] += p.p_TFs.Cijk_rates_sub_i[i][it] * (u[p.p_TFs.Ck_sub_i[i]][it] * uE[p.p_TFs.Cj_sub_i[i]][it] + u[p.p_TFs.Cj_sub_i[i]][it] * uE[p.p_TFs.Ck_sub_i[i]][it])
		end

		@inbounds for (it,tmp_Qij_val) in enumerate(p.p_TFs.Qij_vals_sub_i[i])
			terms[2] += tmp_Qij_val
			terms[3] += tmp_Qij_val * u[p.p_TFs.Qj_sub_i[i]][it]
		end
		
		du[i] = -(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4]
  end
end



# Damn! 2022-06-15
#julia> @time sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, 
#       abstol=solver_options.abstol, reltol=solver_options.reltol);
#  0.513253 seconds (3.17 M allocations: 195.377 MiB, 7.21% gc time)
#
#julia> @time sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, 
#       abstol=solver_options.abstol, reltol=solver_options.reltol);
#  0.039865 seconds (697.30 k allocations: 15.083 MiB)

# julia> @benchmark sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, 
#        abstol=solver_options.abstol, reltol=solver_options.reltol)
# BenchmarkTools.Trial: 11 samples with 1 evaluation.
#  Range (min … max):  478.275 ms … 521.345 ms  ┊ GC (min … max): 5.15% … 9.54%
#  Time  (median):     488.913 ms               ┊ GC (median):    5.16%
#  Time  (mean ± σ):   494.444 ms ±  13.522 ms  ┊ GC (mean ± σ):  6.39% ± 2.09%
# 
#   ▁    ▁    █  ▁ ▁      ▁  ▁                ▁   ▁             ▁  
#   █▁▁▁▁█▁▁▁▁█▁▁█▁█▁▁▁▁▁▁█▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
#   478 ms           Histogram: frequency by time          521 ms <
# 
#  Memory estimate: 195.38 MiB, allocs estimate: 3174625.
# 
# julia> @benchmark sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, 
#        abstol=solver_options.abstol, reltol=solver_options.reltol)
# BenchmarkTools.Trial: 114 samples with 1 evaluation.
#  Range (min … max):  39.619 ms … 67.727 ms  ┊ GC (min … max): 0.00% … 33.74%
#  Time  (median):     41.977 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   44.080 ms ±  7.037 ms  ┊ GC (mean ± σ):  4.91% ± 10.21%
# 
#    ▃▁▃▆█                                                       
#   ▅██████▆▆▃▃▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▃▄▃▁▃▁▃ ▃
#   39.6 ms         Histogram: frequency by time        67.4 ms <
# 
#  Memory estimate: 15.08 MiB, allocs estimate: 697302.
parameterized_ClaSSE_Es_v7_simd_sums = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
#  Qij_vals = p.params.Qij_vals
#  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
#	Qarray_ivals = p.p_indices.Qarray_ivals
#	Qarray_jvals = p.p_indices.Qarray_jvals
#	Carray_ivals = p.p_indices.Carray_ivals
#	Carray_jvals = p.p_indices.Carray_jvals
#	Carray_kvals = p.p_indices.Carray_kvals

#	Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i
#	Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i
	
	#terms = p.terms
	terms = Vector{Float64}(undef, 4)

  @inbounds @simd for i in 1:n
		#Qi_sub_i = p.p_TFs.Qi_sub_i[i]
#		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		#Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		#Ci_sub_i = p.p_TFs.Ci_sub_i[i]
#		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
#		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		#Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "D" (likelihood of tip data)
#		du[i] = -(sum(Cijk_rates_sub_i[i]) + sum(Qij_vals_sub_i[i]) + mu[i])*u[i] +  # case 1: no event
#			(sum(Qij_vals_sub_i[i] .* u[Qj_sub_i])) + 	# case 2	
#			(sum(Cijk_rates_sub_i[i] .*                                               # case 3/4: change + eventual extinction
#				 (u[Ck_sub_i].*uE[Cj_sub_i] 
#			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
		
#		Cijk_rates_temp = @view p.p_TFs.Cijk_rates_sub_i[i]
#		Qij_rates_temp = @view p.p_TFs.Qij_vals_sub_i[i]
		
		terms .= 0.0

		terms[1], terms[4] = sum_Cijk_rates_Es_inbounds_simd(p.p_TFs.Cijk_rates_sub_i[i], u, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=terms[1], term4=terms[4])
	
		terms[2], terms[3] = sum_Qij_vals_inbounds_simd(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i]; term2=terms[2], term3=terms[3])
		
		du[i] = mu[i] -(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4]
  end
end



# Takes 0.041236 seconds, holy moly
# 
# julia> @time sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, 
#        abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);
#   0.618784 seconds (3.35 M allocations: 231.215 MiB, 6.98% gc time)
# 
# julia> @time sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, 
#        abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);
#   0.042717 seconds (649.83 k allocations: 14.951 MiB)
#
#
# julia> @benchmark sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=false, 
#        abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
# BenchmarkTools.Trial: 9 samples with 1 evaluation.
#  Range (min … max):  580.423 ms … 671.283 ms  ┊ GC (min … max): 4.24% … 8.30%
#  Time  (median):     615.601 ms               ┊ GC (median):    8.21%
#  Time  (mean ± σ):   616.233 ms ±  24.939 ms  ┊ GC (mean ± σ):  6.53% ± 2.14%
# 
#   █         █      █  █  █ █ █   █                            █  
#   █▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁█▁▁█▁▁█▁█▁█▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
#   580 ms           Histogram: frequency by time          671 ms <
# 
#  Memory estimate: 231.22 MiB, allocs estimate: 3345985.
# 
# julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=false, 
#        abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
# BenchmarkTools.Trial: 113 samples with 1 evaluation.
#  Range (min … max):  40.291 ms … 69.749 ms  ┊ GC (min … max): 0.00% … 35.60%
#  Time  (median):     42.087 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   44.531 ms ±  7.518 ms  ┊ GC (mean ± σ):  4.95% ± 10.06%
# 
#    ▅█▄▁ ▂                                                      
#   ▇████▅█▆▄▃▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▃▃▃▄▁▁▄ ▃
#   40.3 ms         Histogram: frequency by time        68.5 ms <
# 
#  Memory estimate: 14.95 MiB, allocs estimate: 649827.
# 
parameterized_ClaSSE_Ds_v7_simd_sums = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
#  Qij_vals = p.params.Qij_vals
#  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
#	Qarray_ivals = p.p_indices.Qarray_ivals
#	Qarray_jvals = p.p_indices.Qarray_jvals
#	Carray_ivals = p.p_indices.Carray_ivals
#	Carray_jvals = p.p_indices.Carray_jvals
#	Carray_kvals = p.p_indices.Carray_kvals

#	Qij_vals_sub_i = p.p_TFs.Qij_vals_sub_i
#	Cijk_rates_sub_i = p.p_TFs.Cijk_rates_sub_i
	
	# Pre-calculated solution of the Es
#	sol_Es = p.sol_Es_v5
#	uE = p.uE
	uE = p.sol_Es_v5(t)
	terms = Vector{Float64}(undef, 4)
  @inbounds @simd for i in 1:n
		#Qi_sub_i = p.p_TFs.Qi_sub_i[i]
#		Qj_sub_i = p.p_TFs.Qj_sub_i[i]
		#Qi_eq_i  = p.p_TFs.Qi_eq_i[i]

  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		#Ci_sub_i = p.p_TFs.Ci_sub_i[i]
#		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
#		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		#Ci_eq_i  = p.p_TFs.Ci_eq_i[i]

		# Calculation of "D" (likelihood of tip data)
#		du[i] = -(sum(Cijk_rates_sub_i[i]) + sum(Qij_vals_sub_i[i]) + mu[i])*u[i] +  # case 1: no event
#			(sum(Qij_vals_sub_i[i] .* u[Qj_sub_i])) + 	# case 2	
#			(sum(Cijk_rates_sub_i[i] .*                                               # case 3/4: change + eventual extinction
#				 (u[Ck_sub_i].*uE[Cj_sub_i] 
#			 .+ u[Cj_sub_i].*uE[Ck_sub_i]) ))
		
#		Cijk_rates_temp = @view p.p_TFs.Cijk_rates_sub_i[i]
#		Qij_rates_temp = @view p.p_TFs.Qij_vals_sub_i[i]
		
		terms .= 0.0

		terms[1], terms[4] = sum_Cijk_rates_Ds_inbounds_simd(p.p_TFs.Cijk_rates_sub_i[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=terms[1], term4=terms[4])
	
		terms[2], terms[3] = sum_Qij_vals_inbounds_simd(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i]; term2=terms[2], term3=terms[3])
		
		du[i] = -(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4]
  end
end


# Time-varying areas & extinction rates
parameterized_ClaSSE_Es_v10_simd_sums = (du,u,p,t) -> begin
 	# The row of bmo that refers to "u", the effect of area on extinction rate
 	# u_row = (1:Rnrow(bmo))[bmo.rownames .== "u"][]
 	max_extinction_rate = p.setup.max_extinction_rate

 	# Get the area of areas at time t
	p.setup.area_of_areas .= p.interpolators.area_of_areas_interpolator(t)
 	
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals # base extinction rate of each range
  mu_t = p.params.mu_t_vals # mu_t = mu at time t
  
  # Populate changing mus with time
  @inbounds for i in 1:n
  	# total_area = get_area_of_range(tval, state_as_areas_list, area_of_areas_interpolator)
  	mu_t[i] = mu[i] * get_area_of_range(t, p.states_as_areas_lists[i], p.setup.area_of_areas)^p.bmo.est[p.setup.bmo_rows.u_e]
  end
  # Correct "Inf" max_extinction_rates
  mu_t[mu_t .> max_extinction_rate] .= max_extinction_rate
  
  
  # Get the e_vals for the Qij matrix, at time t
  update_Qij_e_vals!(p)
  # (updates p.params.Qij_vals)
  
  
 # Get the d_vals for the Qij matrix, at time t
  # 1. Update the distance matrices etc.
  p.setup.distmat .= p.interpolators.distances_interpolator(t)
  
  
  # Update the vicariance minimum distance 
  p.setup.vicdist_t .= p.interpolators.vicariance_mindists_interpolator(t)
    
  # Using the current t's distmat, etc. update the dmat_t, then 
  # propagate through the Q matrix
  update_Qij_d_vals!(p)
  
  # Using the current t's distmat, etc. update the jmat_t, then
  # propagate through the C matrix
  #update_Cijk_j_rates!(p)
  # update vicariance also
  update_Cijk_rates!(p)
  
  # Populate changing "e" with time
	#terms = Vector{Float64}(undef, 4)

  @inbounds for i in 1:n
		p.terms .= 0.0

		p.terms[1], p.terms[4] = sum_Cijk_rates_Es_inbounds_simd(p.p_TFs.Cijk_rates_sub_i_t[i], u, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=p.terms[1], term4=p.terms[4])
	
		#terms[2], terms[3] = sum_Qij_vals_inbounds_simd(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i]; term2=terms[2], term3=terms[3])
		p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals_t[p.p_TFs.Qi_eq_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		
		du[i] = mu_t[i] -(p.terms[1] + p.terms[2] + mu_t[i])*u[i] + p.terms[3] + p.terms[4]
  end
end # END parameterized_ClaSSE_Es_v10_simd_sums = (du,u,p,t) -> begin


# Time-varying areas & extinction rates
parameterized_ClaSSE_Ds_v10_simd_sums = (du,u,p,t) -> begin
 	# The row of bmo that refers to "u", the effect of area on extinction rate
 	# u_row = (1:Rnrow(bmo))[bmo.rownames .== "u"][]
 	max_extinction_rate = p.setup.max_extinction_rate
 	
 	# Get the area of areas at time t
	p.setup.area_of_areas .= p.interpolators.area_of_areas_interpolator(t)

  # NOT the slow step:
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals # base extinction rate of each range
  mu_t = p.params.mu_t_vals # mu_t = mu at time t
  
  # Populate changing mus with time
  @inbounds for i in 1:n
  	# total_area = get_area_of_range(tval, state_as_areas_list, area_of_areas_interpolator)
  	mu_t[i] = mu[i] * get_area_of_range(t, p.states_as_areas_lists[i], p.setup.area_of_areas)^p.bmo.est[p.setup.bmo_rows.u_e]
  end
  # Correct "Inf" max_extinction_rates
  mu_t[mu_t .> max_extinction_rate] .= max_extinction_rate


  # Get the e_vals for the Qij matrix, at time t
  # elist_actual = elist_base * area_of_area_lost^u_e
  # THIS IS THE SLOW STEP; pre-allocate e_rows
  ##prtQp(p)
  update_Qij_e_vals!(p);
  ##prtQp(p)
  # (updates p.params.Qij_vals)


  # Get the d_vals for the Qij matrix, at time t
  # 1. Update the distance matrices etc.
  p.setup.distmat .= p.interpolators.distances_interpolator(t)
  
  
  # Update the vicariance minimum distance 
  p.setup.vicdist_t .= p.interpolators.vicariance_mindists_interpolator(t)
  
  # ...others?
  
  # Using the current t's distmat, etc. update the dmat, then 
  # propagate through the 
  ##prtQp(p)
  update_Qij_d_vals!(p);
  ##prtQp(p)
  
  # Using the current t's distmat, etc. update the jmat_t, then
  # propagate through the C matrix
  ##x1 = prtCp(p).rates_t
  #update_Cijk_j_rates!(p)
  # update vicariance also
  update_Cijk_rates!(p)
  ##x2 = prtCp(p).rates_t
	##x1 .- x2
	
	
	# Pre-calculated solution of the Es
#	sol_Es = p.sol_Es_v5
#	uE = p.uE
	uE = p.sol_Es_v10(t)
	#terms = Vector{Float64}(undef, 4)
  @inbounds for i in 1:n
		p.terms .= 0.0

		p.terms[1], p.terms[4] = sum_Cijk_rates_Ds_inbounds_simd(p.p_TFs.Cijk_rates_sub_i_t[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=p.terms[1], term4=p.terms[4])
		
		# Is this the slow step?? -- NO!
		#terms[2], terms[3] = sum_Qij_vals_inbounds_simd(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i]; term2=terms[2], term3=terms[3])
		
		################################################
		# Caused a slight difference!! :
		# WRONG: sum_Qij_vals_inbounds_simd(p.params.Qij_vals_t[p.p_TFs.Qi_sub_i[i]]
		# RIGHT: sum_Qij_vals_inbounds_simd(p.p_TFs.Qij_vals_sub_i[i]
		#
		# Original definition in update:
		# p_Ds_v5.p_TFs.Qij_vals_sub_i[i] .= Qij_vals[p_Ds_v5.p_TFs.Qi_eq_i[i]]
		# eq rather than sub!!!
		#
		# Differences, dej 0.34, 0.28, 0.11, bd = 0.1, 0.1
		#  -39.96536521864187 - -40.092317498361034
		# 0.12695227971916268
		#
		# Differences, dej 0.34, 0.28, 0.11, bd = 0.34, 0.0
		# -26.775073112165167 - -26.663186300866577
		# -0.11188681129858935
		# 
		# Differences, dej 0.34, 0.28, 0.0, bd = 0.34, 0.0
		# -34.7239266954144 - -34.51134005431238
		# -0.2125866411020212
		#
		# # Differences, dej 0.0, 0.0, 0.0, bd = 0.34, 0.0
		# -20.921822175682088 - -20.921822175682088
		# 0.0
		
		# Broke
		#p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals_t[p.p_TFs.Qi_sub_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		# Works
		#p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		# Works
		p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals_t[p.p_TFs.Qi_eq_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		
		
		
		
		du[i] = -(p.terms[1] + p.terms[2] + mu_t[i])*u[i] + p.terms[3] + p.terms[4]
  end
end # END parameterized_ClaSSE_Ds_v10_simd_sums = (du,u,p,t) -> begin




# EXPERIMENTING/TESTING
# Time-varying areas & extinction rates
parameterized_ClaSSE_Es_v11_simd_sums = (du,u,p,t) -> begin
 	# The row of bmo that refers to "u", the effect of area on extinction rate
 	# u_row = (1:Rnrow(bmo))[bmo.rownames .== "u"][]
 	##max_extinction_rate = p.setup.max_extinction_rate

 	# Get the area of areas at time t
	##p.setup.area_of_areas .= p.area_of_areas_interpolator(t)
 	
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals # base extinction rate of each range
  ##mu_t = p.params.mu_t_vals # mu_t = mu at time t
  
  # Populate changing mus with time
  ##@inbounds for i in 1:n
  	# total_area = get_area_of_range(tval, state_as_areas_list, area_of_areas_interpolator)
  ##	mu_t[i] = mu[i] * get_area_of_range(t, p.states_as_areas_lists[i], p.setup.area_of_areas)^p.bmo.est[p.setup.bmo_rows.u_e]
  ##end
  # Correct "Inf" max_extinction_rates
  ##mu_t[mu_t .> max_extinction_rate] .= max_extinction_rate
  
  
  # Get the e_vals for the Qij matrix, at time t
  ##update_Qij_e_vals!(p)
  # (updates p.params.Qij_vals)
  
  
  # Get the d_vals for the Qij matrix, at time t
  # 1. Update the distance matrices etc.
  ##p.setup.distmat .= p.distances_interpolator(t)
  # ...others?
  
  # Using the current t's distmat, etc. update the dmat_t, then 
  # propagate through the Q matrix
  ##update_Qij_d_vals!(p)
  
  # Using the current t's distmat, etc. update the jmat_t, then
  # propagate through the C matrix
  ##update_Cijk_j_rates!(p)
  
  # Populate changing "e" with time
	#terms = Vector{Float64}(undef, 4)

  @inbounds for i in 1:n
		p.terms .= 0.0

		p.terms[1], p.terms[4] = sum_Cijk_rates_Es_inbounds_simd(p.p_TFs.Cijk_rates_sub_i[i], u, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=p.terms[1], term4=p.terms[4])
		
		# Broke
		#terms[2], terms[3] = sum_Qij_vals_inbounds_simd(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i]; term2=terms[2], term3=terms[3])
		# Works
		#p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals[p.p_TFs.Qi_sub_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		# Works
		p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals_t[p.p_TFs.Qi_eq_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		
		du[i] = mu[i] -(p.terms[1] + p.terms[2] + mu[i])*u[i] + p.terms[3] + p.terms[4]
  end
end # END parameterized_ClaSSE_Es_v11_simd_sums = (du,u,p,t) -> begin

# EXPERIMENTING/TESTING
# Time-varying areas & extinction rates
parameterized_ClaSSE_Ds_v11_simd_sums = (du,u,p,t) -> begin
 	# The row of bmo that refers to "u", the effect of area on extinction rate
 	# u_row = (1:Rnrow(bmo))[bmo.rownames .== "u"][]
 	##max_extinction_rate = p.setup.max_extinction_rate
 	
 	# Get the area of areas at time t
	p.setup.area_of_areas .= p.area_of_areas_interpolator(t)

  # NOT the slow step:
  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals # base extinction rate of each range
  ##mu_t = p.params.mu_t_vals # mu_t = mu at time t
  
  # Populate changing mus with time
  ##@inbounds for i in 1:n
  	# total_area = get_area_of_range(tval, state_as_areas_list, area_of_areas_interpolator)
  ##	mu_t[i] = mu[i] * get_area_of_range(t, p.states_as_areas_lists[i], p.setup.area_of_areas)^p.bmo.est[p.setup.bmo_rows.u_e]
  ##end
  # Correct "Inf" max_extinction_rates
  ##mu_t[mu_t .> max_extinction_rate] .= max_extinction_rate


  # Get the e_vals for the Qij matrix, at time t
  # elist_actual = elist_base * area_of_area_lost^u_e
  # THIS IS THE SLOW STEP; pre-allocate e_rows
  ##update_Qij_e_vals!(p)
  # (updates p.params.Qij_vals)


  # Get the d_vals for the Qij matrix, at time t
  # 1. Update the distance matrices etc.
  ##p.setup.distmat .= p.distances_interpolator(t)
  # ...others?
  
  # Using the current t's distmat, etc. update the dmat, then 
  # propagate through the 
  ##update_Qij_d_vals!(p)
  
  # Using the current t's distmat, etc. update the jmat_t, then
  # propagate through the C matrix
  ##update_Cijk_j_rates!(p)
	
	# Pre-calculated solution of the Es
#	sol_Es = p.sol_Es_v5
#	uE = p.uE
	uE = p.sol_Es_v10(t)
	#terms = Vector{Float64}(undef, 4)
  @inbounds for i in 1:n
		p.terms .= 0.0

		p.terms[1], p.terms[4] = sum_Cijk_rates_Ds_inbounds_simd(p.p_TFs.Cijk_rates_sub_i[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=p.terms[1], term4=p.terms[4])
		
		# Is this the slow step?? -- NO!
		# Works
		#terms[2], terms[3] = sum_Qij_vals_inbounds_simd(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i]; term2=terms[2], term3=terms[3])
		# Broke
		#p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals[p.p_TFs.Qi_sub_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals_t[p.p_TFs.Qi_eq_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		
		du[i] = -(p.terms[1] + p.terms[2] + mu[i])*u[i] + p.terms[3] + p.terms[4]
  end
end # END parameterized_ClaSSE_Ds_v11_simd_sums = (du,u,p,t) -> begin





# Time-varying areas & extinction rates
parameterized_ClaSSE_Es_v12_simd_sums = (du,u,p,t) -> begin
  # Update all rates in p.params
  
  # We don't even need this, since we have the interpolator
  # update_QC_mats_time_t!(p, t)
  
  # Interpolate the current Q_vals_t and C_rates_t
  p.params.Qij_vals_t .= p.interpolators.Q_vals_interpolator(t)
  p.params.Cijk_rates_t .= p.interpolators.C_rates_interpolator(t)
   
  # Update 
  # p.p_TFs.Cijk_rates_sub_i_t[i] is replaced by p.params.Cijk_rates_t[p.p_TFs.Ci_sub_i[i]]
  # p.p_TFs.Qij_vals_sub_i_t[i]   is replaced by p.params.Qij_vals_t[p.p_TFs.Qi_eq_i[i]]
  
   
  # Populate changing "e" with time
	#terms = Vector{Float64}(undef, 4)

  @inbounds for i in 1:p.n
		p.terms .= 0.0
		
		p.terms[1], p.terms[4] = sum_Cijk_rates_Es_inbounds_simd(p.params.Cijk_rates_t[p.p_TFs.Ci_sub_i[i]], u, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=p.terms[1], term4=p.terms[4])
	
		#terms[2], terms[3] = sum_Qij_vals_inbounds_simd(p.p_TFs.Qij_vals_sub_i[i], u, p.p_TFs.Qj_sub_i[i]; term2=terms[2], term3=terms[3])
		p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals_t[p.p_TFs.Qi_eq_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		
		du[i] = p.params.mu_t_vals[i] -(p.terms[1] + p.terms[2] + p.params.mu_t_vals[i])*u[i] + p.terms[3] + p.terms[4]
  end
end # END parameterized_ClaSSE_Es_v12_simd_sums = (du,u,p,t) -> begin


# Time-varying areas & extinction rates
parameterized_ClaSSE_Ds_v12_simd_sums = (du,u,p,t) -> begin
  # Interpolate the current Q_vals_t and C_rates_t
  p.params.Qij_vals_t .= p.interpolators.Q_vals_interpolator(t)
  p.params.Cijk_rates_t .= p.interpolators.C_rates_interpolator(t)

	# Pre-calculated solution of the Es
	uE = p.sol_Es_v12(t)
	#terms = Vector{Float64}(undef, 4)
  @inbounds for i in 1:p.n
		p.terms .= 0.0

		#p.terms[1], p.terms[4] = sum_Cijk_rates_Ds_inbounds_simd(p.p_TFs.Cijk_rates_sub_i_t[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=p.terms[1], term4=p.terms[4])

# DIFFERENT:
# RIGHT: p.params.Cijk_rates_t[p.p_TFs.Ci_eq_i[i]]
# WRONG: p.params.Cijk_rates_t[p.p_TFs.Ci_sub_i[i]]
		p.terms[1], p.terms[4] = sum_Cijk_rates_Ds_inbounds_simd(p.params.Cijk_rates_t[p.p_TFs.Ci_eq_i[i]], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=p.terms[1], term4=p.terms[4])
		p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals_t[p.p_TFs.Qi_eq_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])



#		p.terms[1], p.terms[4] = sum_Cijk_rates_Ds_inbounds_simd(p.params.Cijk_rates_t[p.p_TFs.Ci_sub_i[i]], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=p.terms[1], term4=p.terms[4])
		
#		p.terms[2], p.terms[3] = sum_Qij_vals_inbounds_simd(p.params.Qij_vals_t[p.p_TFs.Qi_eq_i[i]], u, p.p_TFs.Qj_sub_i[i]; term2=p.terms[2], term3=p.terms[3])
		
		du[i] = -(p.terms[1] + p.terms[2] + p.params.mu_t_vals[i])*u[i] + p.terms[3] + p.terms[4]
  end
end # END parameterized_ClaSSE_Ds_v12_simd_sums = (du,u,p,t) -> begin








# Custom sums using inbounds and simd, but using the repeated loops
# Modeled on: https://web.eecs.umich.edu/~fessler/course/551/julia/tutor/sum-simd.html

# @inbounds version with for i=1:N range loop 
function sum_range_inbounds_simd(a::Vector)
    total = zero(eltype(a))
    @inbounds @simd for i=1:length(a)
        total += a[i]
    end
    return total
end;

function sum_Cijk_rates_Ds_inbounds_simd(Cijk_rates_sub_i, tmp_u, tmp_uE, Cj_sub_i, Ck_sub_i; term1=Float64(0.0), term4=Float64(0.0))
    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
    	term4 += Cijk_rates_sub_i[it] * (tmp_u[Ck_sub_i[it]] * tmp_uE[Cj_sub_i[it]] + tmp_u[Cj_sub_i[it]] * tmp_uE[Ck_sub_i[it]])
    end
    return term1, term4
end;

function sum_Cijk_rates_Es_inbounds_simd(Cijk_rates_sub_i, tmp_u, Cj_sub_i, Ck_sub_i; term1=Float64(0.0), term4=Float64(0.0))
    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
    	term4 += Cijk_rates_sub_i[it] * tmp_u[Cj_sub_i[it]] * tmp_u[Ck_sub_i[it]]
    end
    return term1, term4
end;



# # Tturbo turned out to not be faster...
# Simple array operation in parallel using julia
# https://stackoverflow.com/questions/68424283/simple-array-operation-in-parallel-using-julia
# using LoopVectorization
# @tturbo

function sum_Qij_vals_inbounds_simd(Qij_vals_sub_i, tmp_u, Qj_sub_i; term2=Float64(0.0), term3=Float64(0.0))
    @inbounds @simd for it=1:length(Qij_vals_sub_i)
    	term2 += Qij_vals_sub_i[it]
    	term3 += Qij_vals_sub_i[it] * tmp_u[Qj_sub_i[it]]
    end
    return term2, term3
end;



end # end of module
