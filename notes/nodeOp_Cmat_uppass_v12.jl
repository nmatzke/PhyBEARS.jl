# Uppass version of nodeOp_Cmat

# 1. For a given node, get the uppass probs at the bottom of the branch
#    (unless it's the root)
# 2. Run the Ds up from that
# 3. Run nodeOp_Cmat_get_condprobs to get uppass probs for either 
#    L or R descendant

using PhyloBits.TrUtils # for odds, Rrbind, convert_df_datatypes!



# Convert a cladogenesis table with paired events to one with non-paired events
# adjust probs / vals / rates accordingly

"""
dftxt = "Any[\"y\" 2 2 2 1 1.0 1.0 0.2222222222222222 0.2222222222222222 0.0; \"y\" 3 3 3 1 1.0 1.0 0.2222222222222222 0.2222222222222222 0.0; \"v\" 4 3 2 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"s\" 4 4 2 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"s\" 4 4 3 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"v\" 4 2 3 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"s\" 4 2 4 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0; \"s\" 4 3 4 2 2.0 0.16666666666666666 0.037037037037037035 0.037037037037037035 0.0]"

dfnames_txt = "[\"event\", \"i\", \"j\", \"k\", \"pair\", \"wt\", \"prob\", \"rate\", \"val\", \"rates_t\"]"
dfnames = Reval(dfnames_txt)

datatypes = Reval("DataType[String, Int64, Int64, Int64, Int64, Float64, Float64, Float64, Float64, Float64]")

ctable = DataFrame(Reval(dftxt), dfnames)
convert_df_datatypes!(ctable, datatypes)
ctable

"""
function make_ctable_single_events(ctable1)
	ctable2 = ctable1[ctable1.pair .== 2, :]
	tmpj = deepcopy(ctable2.j)
	tmpk = deepcopy(ctable2.k)
	ctable2.j .= tmpk
	ctable2.k .= tmpj

	ctable = Rrbind(ctable1, ctable2)
	ctable.prob .= ctable.prob ./ ctable.wt
	ctable.rate .= ctable.rate ./ ctable.wt
	ctable.val .= ctable.val ./ ctable.wt
	ctable.rates_t .= ctable.rates_t ./ ctable.wt

	return ctable
end # END function make_ctable_single_events(ctable)



#######################################################
# Run the Ds calculation from old to new (up the tree),
# reverse of the standard method.
#######################################################


# Calculate Ds UP a branch
#
# Modifies branchOp to do Ds calculation down a branch
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_ClaSSE_Ds_v7_FWD(current_nodeIndex, res; u0, tspan, p_Ds_v7, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Example slow operation
	#y = countloop(num_iterations, current_nodeIndex)
# 	print("\n")
# 	print("branchOp_ClaSSE_Ds_v5: d: ")
# 	print(p_Ds_v5.params.Qij_vals[1])
# 	print("branchOp_ClaSSE_Ds_v5: e: ")
# 	print(p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)])
# 	print("\n")
	
	prob_Ds_v7 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v7_simd_sums_FWD, deepcopy(u0), tspan, p_Ds_v7)

	#sol_Ds = solve(prob_Ds_v7, solver_options.solver, dense=false, save_start=false, save_end=true, save_everystep=false, abstol=solver_options.abstol, reltol=solver_options.reltol)
	
	sol_Ds = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, saveat=solver_options.saveat, abstol=solver_options.abstol, reltol=solver_options.reltol)
	
	
	# Error catch seems to slow it down!!
	"""
	try
		sol_Ds = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol)
	catch e
		txt = paste0(["\nError in branchOp_ClaSSE_Ds_v5(): solve(prob_Ds_v5, ", SubString(string(solver_options.solver), 1:15), "..., save_everystep=", string(solver_options.save_everystep), ", abstol=", string(solver_options.abstol), ", reltol=", string(solver_options.reltol), "): this error may indicate an impossible parameter combination (e.g. all rates are 0.0). Returning 1e-15 for all Ds."])
		print(txt)
		# Return a fake sol_Ds function
		nan_val = 1e-15
		function sol_Ds_tmp(t, tmp_n)
			res = collect(repeat([nan_val], tmp_n))
		end
		sol_Ds = t -> sol_Ds_tmp(t, length(u0)) # returns a array of nan_val, regardless of time t input
	end
	"""
# 	print(sol_Ds[length(sol_Ds)])
# 	print("\n")
	
	#nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	#nodeData_at_bottom = nodeData_at_top / 2.0
	#nodeData_at_bottom = sol_Ds.u
	
	return(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)
end



parameterized_ClaSSE_Ds_v7_simd_sums_FWD3 = (du,u,p,t) -> begin

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

#		terms[1], terms[4] = sum_Cijk_rates_Ds_inbounds_simd(p.p_TFs.Cijk_rates_sub_i[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=terms[1], term4=terms[4])
		terms[1], terms[4] = sum_Cijk_rates_Ds_inbounds_simd(p.p_TFs.Cijk_rates_sub_i[i], u, uE, p.p_TFs.Cj_sub_i[i], p.p_TFs.Ck_sub_i[i]; term1=terms[1], term4=terms[4])
	
		terms[2], terms[3] = sum_Qij_vals_inbounds_simd_FWD(p.p_TFs.Qji_vals_sub_j[i], u, p.p_TFs.Qi_sub_j[i]; term2=terms[2], term3=terms[3])
		
		du[i] = -1 * (-(terms[1] + terms[2] + mu[i])*u[i] + terms[3] + terms[4])
  end
end


# Quick
calcDs_4states2 = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_jvals		# switch these for uppass
	Qarray_jvals = p.p_indices.Qarray_ivals		# switch these for uppass
	Carray_ivals = p.p_indices.Carray_kvals   # Best so far: k, i, j; SAME as k,j,i; 0.0140907
	Carray_jvals = p.p_indices.Carray_ivals		# 2nd best: jik, jki  sum(abs.(err)) = 0.0419
	Carray_kvals = p.p_indices.Carray_jvals		# 3rd best: ijk, ikj  sum(abs.(err)) = 0.0735427
	
	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	two = 1.0
  @inbounds for i in 1:n
		# Calculation of "D" (likelihood of tip data)
#		du[i] = -(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[i] +  # case 1: no event
		du[i] = -(sum(Cijk_vals[Carray_jvals .== i]) + sum(Qij_vals[Qarray_jvals .== i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qarray_ivals .== i] .* u[(Qarray_jvals[Qarray_ivals .== i])])) #+ 	# case 2	
			#(sum(Cijk_vals[Carray_ivals .== i] .*                                               # case 34: change + eventual extinction
			#	 (u[(Carray_kvals[Carray_ivals .== i])].*uE[Carray_jvals[Carray_ivals .== i]] 
			# .+ u[(Carray_jvals[Carray_ivals .== i])].*uE[Carray_kvals[Carray_ivals .== i]]) ))
  end
end


# Quick
calcDs_4states3 = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.mu_vals
  Qij_vals = p.Qij_vals
  Cijk_vals = p.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.Qarray_jvals		# switch these for uppass
	Qarray_jvals = p.Qarray_ivals		# switch these for uppass
	Carray_ivals = p.Carray_kvals   # Best so far: k, i, j; SAME as k,j,i; 0.0140907
	Carray_jvals = p.Carray_ivals		# 2nd best: jik, jki  sum(abs.(err)) = 0.0419
	Carray_kvals = p.Carray_jvals		# 3rd best: ijk, ikj  sum(abs.(err)) = 0.0735427
	
	# Pre-calculated solution of the Es
	sol_Es = p.sol_Es_v5
	uE = p.uE
	uE = sol_Es(t)
	
	two = 1.0
  @inbounds for i in 1:n
		# Calculation of "D" (likelihood of tip data)
		#du[i] = -(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[i] +  # case 1: no event

		du[i] = -(sum(Cijk_vals[Carray_jvals .== i]) + sum(Qij_vals[Qarray_jvals .== i]) + mu[i])*u[i] +  # case 1: no event
			(sum(Qij_vals[Qarray_ivals .== i] .* u[(Qarray_jvals[Qarray_ivals .== i])])) #+ 	# case 2	
			#(sum(Cijk_vals[Carray_ivals .== i] .*                                               # case 34: change + eventual extinction
			#	 (u[(Carray_kvals[Carray_ivals .== i])].*uE[Carray_jvals[Carray_ivals .== i]] 
			# .+ u[(Carray_jvals[Carray_ivals .== i])].*uE[Carray_kvals[Carray_ivals .== i]]) ))
  end
end




# DON'T USE
# THIS ONE IS THE SAME FORMULA, JUST ADD IN "j" instead of "i" as the thing to iterate over
# (even putting in i 
function sum_Cijk_rates_Ds_inbounds_simd_FWD(Cijk_rates_sub_i, Cjik_rates_sub_j, tmp_u, tmp_uE, Cj_sub_i, Ck_sub_i, Cj_sub_j, Ck_sub_j; term1=Float64(0.0), term4=Float64(0.0))
		""" # Original:
    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
    	term4 += Cijk_rates_sub_i[it] * (tmp_u[Ck_sub_i[it]] * tmp_uE[Cj_sub_i[it]] + tmp_u[Cj_sub_i[it]] * tmp_uE[Ck_sub_i[it]])
    end
		"""

    @inbounds @simd for it=1:length(Cijk_rates_sub_i)
    	term1 += Cijk_rates_sub_i[it]
#    	term4 += Cijk_rates_sub_i[it] * (tmp_u[Ck_sub_i[it]] * tmp_uE[Cj_sub_i[it]] + tmp_u[Cj_sub_i[it]] * tmp_uE[Ck_sub_i[it]])
    	term4 += Cijk_rates_sub_i[it] * (tmp_u[Cj_sub_i[it]] * tmp_uE[Ck_sub_i[it]] + tmp_u[Cj_sub_j[it]] * tmp_uE[Ck_sub_j[it]])
    end
    return term1, term4
end;

# THIS ONE IS DIFFERENT, as i->j and j->i can have different probabilities
function sum_Qij_vals_inbounds_simd_FWD(Qji_vals_sub_j, tmp_u, Qi_sub_j; term2=Float64(0.0), term3=Float64(0.0))
		""" # Original
    @inbounds @simd for it=1:length(Qij_vals_sub_i)
    	term2 += Qij_vals_sub_i[it]
    	term3 += Qij_vals_sub_i[it] * tmp_u[Qj_sub_i[it]]
    end
		"""
		# These should sum to the same: Qij_vals_sub_i, Qji_vals_sub_j, across all i or j
#    @inbounds @simd for it=1:length(Qij_vals_sub_i)
#    	term2 += Qij_vals_sub_i[it]
#    end
    @inbounds @simd for it=1:length(Qji_vals_sub_j)
	    term2 += Qji_vals_sub_j[it]
    	term3 += Qji_vals_sub_j[it] * tmp_u[Qi_sub_j[it]] # Different on uppass; Freyman paper
    end

    return term2, term3
end;


# # Tturbo turned out to not be faster...
# Simple array operation in parallel using julia
# https://stackoverflow.com/questions/68424283/simple-array-operation-in-parallel-using-julia
# using LoopVectorization
# @tturbo
"""

function sum_Qij_vals_inbounds_simd(Qij_vals_sub_i, tmp_u, Qj_sub_i; term2=Float64(0.0), term3=Float64(0.0))
    @inbounds @simd for it=1:length(Qij_vals_sub_i)
    	term2 += Qij_vals_sub_i[it]
    	term3 += Qij_vals_sub_i[it] * tmp_u[Qj_sub_i[it]]
    end
    return term2, term3
end;
"""



function nodeOp_Cmat_uppass_v12!(res, current_nodeIndex, trdf, p_Ds_v12, solver_options)
	n = numstates = length(res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex])
	# Is this a root node?
	if (current_nodeIndex == res.root_nodeIndex)
		#uppass_probs_just_below_node = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
		uppass_probs_just_below_node = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
	
	# Tip or internal nodes require passing probs up from branch bottom
	else
		# The uppass ancestral state probs will have been previously 
		# calculated at the branch bottom
		# BGB's "relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS"
		# Multiply by saved likelihood at the node to make very small (may avoid scaling issues)
		u0 = probs_at_branch_bottom = res.uppass_probs_at_each_nodeIndex_branchBot[current_nodeIndex] .+ 0.0
		time_start = trdf.node_age[trdf.ancNodeIndex[current_nodeIndex]]
		time_end = trdf.node_age[current_nodeIndex]
		tspan = [time_start, time_end]
		
		# Seems to work 
		txt = paste0(["\nNode #", current_nodeIndex])
		print("\n")
		print(txt)
		print("\nStarting probs at branch bottom:")
		print(u0)
		
		# Uses parameterized_ClaSSE_Ds_v12_simd_sums_noNegs
		(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)= branchOp_ClaSSE_Ds_v12_noNegs(current_nodeIndex, res; u0, tspan, p_Ds_v12, solver_options=solver_options);
		
		# These are really conditional probabilities upwards, they don't 
		# have to add up to 1.0, unless normalized
		uppass_probs_just_below_node = sol_Ds.u[length(sol_Ds.u)]
		print("\nuppass_probs_just_below_node, pre-correction:")
		print(uppass_probs_just_below_node)
		
		# Correct for any values slipping below 0.0
		TF = uppass_probs_just_below_node .<= 0.0
		if (any(TF))
			# Round to zero
			#uppass_probs_just_below_node[TF] .= 0.0
			# versus add the minimum (may preserve magnitude better)
			uppass_probs_just_below_node .= uppass_probs_just_below_node .- minimum(uppass_probs_just_below_node)
		end
		
		# Normalize to sum to 1.0, *IF* sum is greater than 1
		if (sum(uppass_probs_just_below_node) > 1.0)		
			uppass_probs_just_below_node .= uppass_probs_just_below_node ./ sum(uppass_probs_just_below_node)

			print("\nuppass_probs_just_below_node, post-correction:")
			print(uppass_probs_just_below_node)
			print("\n")
		end
	end # END if (current_nodeIndex == res.root_nodeIndex)
			# END uppass from branch below
	
	# Combine info through node; Store results
	# Check if its a tip node
	if (trdf.nodeType[current_nodeIndex] == "tip")
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .* res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] = res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] ./ sum(res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex])
	
	# Internal nodes (not root)
	elseif ((trdf.nodeType[current_nodeIndex] == "intern") || (trdf.nodeType[current_nodeIndex] == "root") )
		# (Ignore direct ancestors for now)
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		
		# The root node does NOT need to be multiplied; this would produce anc_estimates.^2
		if (trdf.nodeType[current_nodeIndex] != "root")
			res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .* res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
		end
		
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] = res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] ./ sum(res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex])
		
		# Get uppass probs for Left and Right daughter branches
		node_above_Left_corner = trdf.leftNodeIndex[current_nodeIndex]
		node_above_Right_corner = trdf.rightNodeIndex[current_nodeIndex]
		# LEFT
		
		# Calculate the post-cladogenesis uppass probabilities for the Left node
		Ldownpass_likes = collect(repeat([1.0], n))
		Rdownpass_likes = res.normlikes_at_each_nodeIndex_branchBot[node_above_Right_corner]
		relprob_each_split_scenario = nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v12; use_Cijk_rates_t=true)

		for j in 1:n
			jTF = p_Ds_v12.p_indices.Carray_jvals .== j
			res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner][j] = sum(relprob_each_split_scenario[jTF])
		end
		
		# Don't normalize this
		#res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner] .= res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner] ./ sum(res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner])
		
		# RIGHT
		# Calculate the post-cladogenesis uppass probabilities for the Left node
		Rdownpass_likes = collect(repeat([1.0], n))
		Ldownpass_likes = res.normlikes_at_each_nodeIndex_branchBot[node_above_Left_corner]
		relprob_each_split_scenario = nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v12; use_Cijk_rates_t=true)
		
		for k in 1:n
			kTF = p_Ds_v12.p_indices.Carray_jvals .== k
			res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner][k] = sum(relprob_each_split_scenario[kTF])
		end
		# Don't normalize this
		#res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner] .= res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner] ./ sum(res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner])
		
		# Get ancestral range estimates for Left and Right daughter branches
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] .= Ldownpass_likes .* res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner]
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] .= res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] ./ sum(res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner])

		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] .= Rdownpass_likes .* res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner]
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] .= res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] ./ sum(res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner])

	end # END elseif internal nodes
end # END nodeOp_Cmat_uppass_v12!(res, current_nodeIndex, trdf, p_Ds_v12, solver_options)

# Get the conditional probabilities of all cladogenesis
# scenarios, conditional on ancestral range probs, Lprobs, Rprobs
#
# On an uppass, either Lprobs or Rprobs will be a vector of 1s
function nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v12; use_Cijk_rates_t=false)
	p = p_Ds_v12;
	
	# Figure out if you are solving for the left or right descendant
	# (if both are all 1.0s, assumes left)
	# (also assumes left if both are non-1st; but throws warning)
	left_or_right = ""
	if (all(Rdownpass_likes .== 1.0) == true)
		left_or_right = "right"
	elseif (all(Ldownpass_likes .== 1.0) == true)
		left_or_right = "left"
	else
		left_or_right = "left"
		txt = "WARNING in nodeOp_Cmat_get_condprobs(): Either Lprobs or Rprobs should be a vector of all 1.0. But this is not the case. Double-check your inputs."
		@warn txt
	end	
	
	
	# Go through each Ci (ancestral state index)
	# Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  if (use_Cijk_rates_t == true)
	  Cijk_vals = p.params.Cijk_rates_t # DIFFERENCE FOR V12
	else
		Cijk_vals = p.params.Cijk_vals
	end

	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	Carray_pair = p.p_indices.Carray_pair
	Carray_pair_floats = Float64.(Carray_pair)
	Carray_pair_ones = Carray_pair_floats .- 1.0
	
	# Recalculate probabilities of each cladogenesis event:
	anc_ivals = sort(unique(Carray_ivals))
	sumrates_by_i = collect(repeat([0.0], length(anc_ivals)))
	for i in (1:length(anc_ivals))
		ival = anc_ivals[i]
		iTF = Carray_ivals .== ival
		sumrates_by_i[i] = sum(Cijk_vals[iTF])
		p.params.Cijk_probs[iTF] .= Cijk_vals[iTF] ./ sumrates_by_i[i]
	end
	
	
	
	
	# See: calc_uppass_probs_v1.R
	# calc_uppass_scenario_probs_new2
	calc_uppass_scenario_probs_new2_CODE_IS_EG="""
	Lprobs = left_branch_downpass_likes[ COO_weights_columnar[[2]] + COOmat_0based_to_Qmat_1based]
	Rprobs = right_branch_downpass_likes[ COO_weights_columnar[[3]] + COOmat_0based_to_Qmat_1based]
	# Weights divided by the sum of the weights for that row
	# COO_weights_columnar[[1]] is the ancestral node, 0-based, cladogenesis-only
	# So, you will always add 1 to the index to get the rowsum conditional
	# on a particular ancestor (3 non-null states means 3 Rsp_rowsums)
	scenario_condprob = COO_weights_columnar[[4]] / Rsp_rowsums[ COO_weights_columnar[[1]] + 1 ]

	input_probs_each_split_scenario = cbind(ancprobs, Lprobs, Rprobs, scenario_condprob)
	# Multiply the probabilities through
	relprob_each_split_scenario = apply(X=input_probs_each_split_scenario, MARGIN=1, FUN=prod)
	sum_relprob_each_split_scenario = sum(relprob_each_split_scenario)
	prob_each_split_scenario = prob_each_split_scenario / sum(prob_each_split_scenario)
	prob_each_split_scenario = relprob_each_split_scenario
	"""

	ctable1 = prtCp(p)
	ctable = make_ctable_single_events(ctable1)
	
	ancprobs_by_scenario = uppass_probs_just_below_node[ctable.i] 
	lprobs_by_scenario = Ldownpass_likes[ctable.j] 
	rprobs_by_scenario = Rdownpass_likes[ctable.k] 

	relprob_each_split_scenario = ancprobs_by_scenario .* lprobs_by_scenario .* rprobs_by_scenario .* ctable.val
			
	# Because you are moving:
	# - FROM a vector of ancestral ranges
	# - TO a vector of cladogenesis scenarios
	# ...it *IS* legitimate to sum the conditional probabilities of
	#    all scenarios, then divide by the sum
	relprob_each_split_scenario = relprob_each_split_scenario ./ sum(relprob_each_split_scenario)
	
  return(relprob_each_split_scenario)
end # END nodeOp_Cmat_get_condprobs = (tmpDs; tmp1, tmp2, p_Ds_v12) -> begin



function nodeOp_Cmat_uppass_v7!(res, current_nodeIndex, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
	p = p_Ds_v7
	n = numstates = length(res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex])
	# Is this a root node?
	if (current_nodeIndex == res.root_nodeIndex)
		#uppass_probs_just_below_node = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
		# Wrong: this incorporates downpass information from both branches above the root; 
		# using information too many times!
		#uppass_probs_just_below_node = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
		
		uppass_probs_just_below_node = repeat([1/n], n)

		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
	
	# Tip or internal nodes require passing probs up from branch bottom
	else
		# The uppass ancestral state probs will have been previously 
		# calculated at the branch bottom
		# BGB's "relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS"
		# Multiply by saved likelihood at the node to make very small (may avoid scaling issues)
		u0 = probs_at_branch_bottom = res.uppass_probs_at_each_nodeIndex_branchBot[current_nodeIndex] .+ 0.0
		time_start = trdf.node_age[trdf.ancNodeIndex[current_nodeIndex]]
		time_end = trdf.node_age[current_nodeIndex]
		tspan = [time_start, time_end]
		
		# Seems to work 
		txt = paste0(["\nNode #", current_nodeIndex])
		print("\n")
		print(txt)
		print("\nStarting probs at branch bottom:")
		print(u0)
		
		# Uses parameterized_ClaSSE_Ds_v7
		# u0 = [8.322405e-13, 0.1129853, 0.677912, 0.2091026]
		(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)= branchOp_ClaSSE_Ds_v7(current_nodeIndex, res; u0, tspan, p_Ds_v7, solver_options=solver_options);
		
		"""
		u0 = [8.322405e-13, 0.1129853, 0.677912, 0.2091026]
		solver_options.solver
		solver_options.save_everystep = true
		solver_options.saveat = seq(2.0, 3.0, 0.1)
		tspan = (2.0, 3.0)
		(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)= branchOp_ClaSSE_Ds_v7(current_nodeIndex, res; u0, tspan, p_Ds_v7, solver_options=solver_options);
		
		sol_Ds(1.0)
		sol_Ds(2.1)
		sol_Ds(2.0)
		sol_Ds(2.1)
		"""
		
		# These are really conditional probabilities upwards, they don't 
		# have to add up to 1.0, unless normalized
		uppass_probs_just_below_node = sol_Ds.u[length(sol_Ds.u)]
		print("\nuppass_probs_just_below_node, pre-correction:")
		print(uppass_probs_just_below_node)
		
		# Correct for any values slipping below 0.0
		TF = uppass_probs_just_below_node .<= 0.0
		if (any(TF))
			# Round to zero
			#uppass_probs_just_below_node[TF] .= 0.0
			# versus add the minimum (may preserve magnitude better)
			uppass_probs_just_below_node .= uppass_probs_just_below_node .- minimum(uppass_probs_just_below_node)
		end
		
		# Normalize to sum to 1.0, *IF* sum is greater than 1
		if (sum(uppass_probs_just_below_node) > 1.0)		
			uppass_probs_just_below_node .= uppass_probs_just_below_node ./ sum(uppass_probs_just_below_node)

			print("\nuppass_probs_just_below_node, post-correction:")
			print(uppass_probs_just_below_node)
			print("\n")
		end
	end # END if (current_nodeIndex == res.root_nodeIndex)
			# END uppass from branch below
	
	# Combine info through node; Store results
	# Check if its a tip node
	if (trdf.nodeType[current_nodeIndex] == "tip")
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .* res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] = res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] ./ sum(res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex])
	
	# Internal nodes (not root)
	elseif ((trdf.nodeType[current_nodeIndex] == "intern") || (trdf.nodeType[current_nodeIndex] == "root") )
		# (Ignore direct ancestors for now)
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		
		# The root node does NOT need to be multiplied; this would produce anc_estimates.^2
		if (trdf.nodeType[current_nodeIndex] != "root")
			res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .* res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
		end
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] = res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] ./ sum(res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex])
		
		# Get uppass probs for Left and Right daughter branches
		node_above_Left_corner = trdf.rightNodeIndex[current_nodeIndex]
		node_above_Right_corner = trdf.leftNodeIndex[current_nodeIndex]


		# LEFT NODE UPPASS
		# Calculate the post-cladogenesis uppass probabilities for the Left node
		ctable = make_ctable_single_events(prtCp(p))
		Ldownpass_likes = collect(repeat([1.0], n))
		Rdownpass_likes = res.normlikes_at_each_nodeIndex_branchBot[node_above_Right_corner]
		relprob_each_split_scenario = nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v7; use_Cijk_rates_t=use_Cijk_rates_t)

		uppass_lprobs = repeat([0.0], n)
		uppass_rprobs = repeat([0.0], n)

		for statei in 1:n
			uppass_lprobs[statei] = sum(relprob_each_split_scenario[ctable.j .== statei])
			uppass_rprobs[statei] = sum(relprob_each_split_scenario[ctable.k .== statei]) # discarded
		end

		uppass_lprobs = uppass_lprobs ./ sum(uppass_lprobs)
		uppass_rprobs = uppass_rprobs ./ sum(uppass_rprobs)
		
		res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner] .= uppass_lprobs

		
		# RIGHT NODE UPPASS
		# Calculate the post-cladogenesis uppass probabilities for the Left node
		Rdownpass_likes = collect(repeat([1.0], n))
		Ldownpass_likes = res.normlikes_at_each_nodeIndex_branchBot[node_above_Left_corner]
		relprob_each_split_scenario = nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v7; use_Cijk_rates_t=use_Cijk_rates_t)

		uppass_lprobs = repeat([0.0], n)
		uppass_rprobs = repeat([0.0], n)

		for statei in 1:n
			uppass_lprobs[statei] = sum(relprob_each_split_scenario[ctable.j .== statei]) # discarded
			uppass_rprobs[statei] = sum(relprob_each_split_scenario[ctable.k .== statei])
		end

		uppass_lprobs = uppass_lprobs ./ sum(uppass_lprobs)
		uppass_rprobs = uppass_rprobs ./ sum(uppass_rprobs)
		res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner] .= uppass_rprobs
		
		# Get ancestral range estimates for Left and Right daughter branches
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] .= Ldownpass_likes .* res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner]
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] .= res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] ./ sum(res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner])

		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] .= Rdownpass_likes .* res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner]
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] .= res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] ./ sum(res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner])
	
	# res.uppass_probs_at_each_nodeIndex_branchBot[R_order]
	
	end # END elseif internal nodes
end # END nodeOp_Cmat_uppass_v7!(res, current_nodeIndex, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)


function uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
	uppass_edgematrix = res.uppass_edgematrix
	
	# Iterate through rows of the edge matrix, like in R
	# Do it in pairs (left branch, then right branch)
	# uppass_edgematrix: column 1 is ancestral node numbers (i.e. rows of trdf)
	#                    column 2 is descendant node numbers (i.e. rows of trdf)
	ivals_odd = odds(1:Rnrow(uppass_edgematrix))
	for i in ivals_odd
		print(i)
		j = i+1
		# Error check: Check the uppass edge matrix; 
		ancnode1 = uppass_edgematrix[i,1]
		ancnode2 = uppass_edgematrix[j,1]
		if (ancnode1 == ancnode2)
			ancnode = ancnode1
		else
			stoptxt = paste0(["STOP ERROR in uppass_ancstates_v12!(): error in res.uppass_edgematrix. This matrix should have pairs of rows, indicating pairs of branches. The node numbers in column 1 should match within this pair, but in your res.uppass_edgematrix, they do not. Error detected at res.uppass_edgematrix row i=", i, ", row j=", j, ". Printing this section of the res.uppass_edgematrix, below."])
			print("\n")
			print(stoptxt)
			print("\n")
			print(res.uppass_edgematrix[i:j,:])
			print("\n")
			error(stoptxt)
		end # END if (ancnode1 == ancnode2)
	
		Lnode = uppass_edgematrix[i,2]
		Rnode = uppass_edgematrix[j,2]
		
		# Work up through the nodes, starting from the root
		current_nodeIndex = ancnode
				
		nodeOp_Cmat_uppass_v7!(res, current_nodeIndex, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=use_Cijk_rates_t)
	
	end # END for (i in odds(1:nrow(trdf))

end # END function uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)



function uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)
	uppass_edgematrix = res.uppass_edgematrix
	
	# Iterate through rows of the edge matrix, like in R
	# Do it in pairs (left branch, then right branch)
	# uppass_edgematrix: column 1 is ancestral node numbers (i.e. rows of trdf)
	#                    column 2 is descendant node numbers (i.e. rows of trdf)
	# Go through internal nodes
	ivals_odd = odds(1:Rnrow(uppass_edgematrix))
	for i in ivals_odd
		j = i+1
		# Error check: Check the uppass edge matrix; 
		ancnode1 = uppass_edgematrix[i,1]
		ancnode2 = uppass_edgematrix[j,1]
		if (ancnode1 == ancnode2)
			ancnode = ancnode1
		else
			stoptxt = paste0(["STOP ERROR in uppass_ancstates_v12!(): error in res.uppass_edgematrix. This matrix should have pairs of rows, indicating pairs of branches. The node numbers in column 1 should match within this pair, but in your res.uppass_edgematrix, they do not. Error detected at res.uppass_edgematrix row i=", i, ", row j=", j, ". Printing this section of the res.uppass_edgematrix, below."])
			print("\n")
			print(stoptxt)
			print("\n")
			print(res.uppass_edgematrix[i:j,:])
			print("\n")
			error(stoptxt)
		end # END if (ancnode1 == ancnode2)
	
		Lnode = uppass_edgematrix[i,2]
		Rnode = uppass_edgematrix[j,2]
		
		# Work up through the nodes, starting from the root
		current_nodeIndex = ancnode
		
		current_node_time = trdf.node_age[current_nodeIndex]
		update_QC_mats_time_t!(p_Ds_v12, current_node_time)
		
		nodeOp_Cmat_uppass_v12!(res, current_nodeIndex, trdf, p_Ds_v12, solver_options)
	
	end # for (i in odds(1:nrow(trdf))
	
	# Then go through tip nodes
	rownums = (1:Rnrow(trdf))
	tipnodes = rownums[trdf.nodeType .== "tip"]
	for ancnode in tipnodes
		# Work up through the nodes, starting from the root
		current_nodeIndex = ancnode
		
		current_node_time = trdf.node_age[current_nodeIndex]
		update_QC_mats_time_t!(p_Ds_v12, current_node_time)
		
		nodeOp_Cmat_uppass_v12!(res, current_nodeIndex, trdf, p_Ds_v12, solver_options)
	end
end





function nodeOp_Cmat_uppass_v7old!(res, current_nodeIndex, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)
	n = numstates = length(res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex])
	# Is this a root node?
	if (current_nodeIndex == res.root_nodeIndex)
		# Wrong: this incorporates downpass information from both branches above the root; 
		# using information too many times!
		#uppass_probs_just_below_node = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
		
		uppass_probs_just_below_node = repeat([1/n], n)
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] .+ 0.0
	
	# Tip or internal nodes require passing probs up from branch bottom
	else
		# The uppass ancestral state probs will have been previously 
		# calculated at the branch bottom
		# BGB's "relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS"
		u0 = probs_at_branch_bottom = res.uppass_probs_at_each_nodeIndex_branchBot[current_nodeIndex] .+ 0.0
		time_start = trdf.node_age[trdf.ancNodeIndex[current_nodeIndex]]
		time_end = trdf.node_age[current_nodeIndex]
		tspan = [time_start, time_end]
		
		# Seems to work 
		txt = paste0(["Node #", current_nodeIndex])
		print("\n")
		print(txt)
		print("\nStarting probs at branch bottom:")
		print(u0)
		
		(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)= branchOp_ClaSSE_Ds_v7(current_nodeIndex, res; u0, tspan, p_Ds_v7, solver_options=solver_options);
		# These are really conditional probabilities upwards, they don't 
		# have to add up to 1.0, unless normalized
		uppass_probs_just_below_node = sol_Ds.u[length(sol_Ds.u)]
		print("\nuppass_probs_just_below_node:")
		print(uppass_probs_just_below_node)
		uppass_probs_just_below_node .= uppass_probs_just_below_node ./ sum(uppass_probs_just_below_node)
	end # END if (current_nodeIndex == res.root_nodeIndex)
			# END uppass from branch below
	
	# Combine info through node; Store results
	# Check if its a tip node
	if (trdf.nodeType[current_nodeIndex] == "tip")
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .* res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] = res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] ./ sum(res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex])
	
	# Internal nodes (not root)
	elseif ((trdf.nodeType[current_nodeIndex] == "intern") || (trdf.nodeType[current_nodeIndex] == "root") )
		# (Ignore direct ancestors for now)
		res.uppass_probs_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .+ 0.0
		# The root node does NOT need to be multiplied; this would produce anc_estimates.^2
		if (trdf.nodeType[current_nodeIndex] != "root")
			res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] .= uppass_probs_just_below_node .* res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
		end
		res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] = res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex] ./ sum(res.anc_estimates_at_each_nodeIndex_branchTop[current_nodeIndex])
		
		# Get uppass probs for Left and Right daughter branches
		node_above_Left_corner = trdf.leftNodeIndex[current_nodeIndex]
		node_above_Right_corner = trdf.rightNodeIndex[current_nodeIndex]
		# LEFT
		
		# Calculate the post-cladogenesis uppass probabilities for the Left node
		Ldownpass_likes = collect(repeat([1.0], n))
		Rdownpass_likes = res.normlikes_at_each_nodeIndex_branchBot[node_above_Right_corner]
		relprob_each_split_scenario = nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v7; use_Cijk_rates_t=use_Cijk_rates_t)

		for j in 1:n
			jTF = p_Ds_v7.p_indices.Carray_jvals .== j
			res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner][j] = sum(relprob_each_split_scenario[jTF])
		end
		res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner] .= res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner] ./ sum(res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner])
		
		# RIGHT
		# Calculate the post-cladogenesis uppass probabilities for the Left node
		Rdownpass_likes = collect(repeat([1.0], n))
		Ldownpass_likes = res.normlikes_at_each_nodeIndex_branchBot[node_above_Left_corner]
		relprob_each_split_scenario = nodeOp_Cmat_get_condprobs(uppass_probs_just_below_node, Ldownpass_likes, Rdownpass_likes, p_Ds_v7; use_Cijk_rates_t=use_Cijk_rates_t)
		
		for k in 1:n
			kTF = p_Ds_v7.p_indices.Carray_jvals .== k
			res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner][k] = sum(relprob_each_split_scenario[kTF])
		end
		res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner] .= res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner] ./ sum(res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner])
		
		# Get ancestral range estimates for Left and Right daughter branches
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] .= Ldownpass_likes .* res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Left_corner]
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] .= res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner] ./ sum(res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Left_corner])

		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] .= Rdownpass_likes .* res.uppass_probs_at_each_nodeIndex_branchBot[node_above_Right_corner]
		res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] .= res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner] ./ sum(res.anc_estimates_at_each_nodeIndex_branchBot[node_above_Right_corner])

	end # END elseif internal nodes
end # END nodeOp_Cmat_uppass_v7old!(res, current_nodeIndex, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

