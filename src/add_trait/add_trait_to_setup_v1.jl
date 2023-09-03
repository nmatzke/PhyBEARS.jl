


# Double inputs.setup for a trait
trait_states_txt = ["H", "D"]   # hermaphrodite, dioecy
trait_areanums = [3,4]
num_trait_states = length(trait_states_txt)
nt = num_trait_states

# inputs.setup: Copy the inputs.setup objects
# Edit for the new list of states
tree_height = setup.tree_height
area_names = setup.area_names
areas_list = setup.areas_list
states_list = setup.states_list[statenums_to_keep]
txt_states_list = setup.txt_states_list[statenums_to_keep]
max_range_size = setup.max_range_size
include_null_range = setup.include_null_range
root_age_mult = setup.root_age_mult
#statenums = setup.statenums[removeTF]
statenums = 1:length(statenums_to_keep)
observed_statenums = ond.(setup.observed_statenums)
numtips = setup.numtips
fossil_TF = setup.fossil_TF
direct_TF = setup.direct_TF
numstates = length(setup.states_list[statenums_to_keep])

# inputs.setup: Area-specific information (stays the same)
numareas = setup.numareas
area_of_areas = setup.area_of_areas
dmat_base = setup.dmat_base
dmat = setup.dmat
dmat_t = setup.dmat_t
jmat_base = setup.jmat_base
jmat = setup.jmat
jmat_t = setup.jmat_t
amat_base = setup.amat_base
amat = setup.amat
amat_t = setup.amat_t
elist = setup.elist
elist_base = setup.elist_base
elist_t = setup.elist_t
dispersal_multipliers_mat = setup.dispersal_multipliers_mat
distmat = setup.distmat
envdistmat = setup.envdistmat
distmat2 = setup.distmat2
distmat3 = setup.distmat3
maxent01 = setup.maxent01
bmo_rows = setup.bmo_rows




# inputs.res: subset to new space space
# node-tops (no subset needed)
regime = res.regime
node_state = res.node_state
node_method = res.node_method
node_Lparent_state = res.node_Lparent_state
node_Rparent_state = res.node_Rparent_state
root_nodeIndex = res.root_nodeIndex
numNodes = res.numNodes
uppass_edgematrix = res.uppass_edgematrix
thread_for_each_nodeOp = res.thread_for_each_nodeOp
thread_for_each_branchOp = res.thread_for_each_branchOp
calc_spawn_start = res.calc_spawn_start
calc_start_time = res.calc_start_time
calc_end_time = res.calc_end_time
calc_duration = res.calc_duration
calctime_iterations = res.calctime_iterations

# subset
sampling_f = res.sampling_f[statenums_to_keep]
tipsamp_f = res.tipsamp_f[statenums_to_keep]

# branch-tops (nodes)
sumLikes_at_node_at_branchTop = res.sumLikes_at_node_at_branchTop
lnL_at_node_at_branchTop = res.lnL_at_node_at_branchTop
lq_at_branchBot = res.lq_at_branchBot
like_at_branchBot = res.like_at_branchBot

# Subset each node/row to the new number of states / columns
Es_at_each_nodeIndex_branchTop = subvv(res.Es_at_each_nodeIndex_branchTop, statenums_to_keep)
Es_at_each_nodeIndex_branchBot = subvv(res.Es_at_each_nodeIndex_branchBot, statenums_to_keep)
fakeX0s_at_each_nodeIndex_branchTop = subvv(res.fakeX0s_at_each_nodeIndex_branchTop, statenums_to_keep)
likes_at_each_nodeIndex_branchTop = subvv(res.likes_at_each_nodeIndex_branchTop, statenums_to_keep)
normlikes_at_each_nodeIndex_branchTop = subvv(res.normlikes_at_each_nodeIndex_branchTop, statenums_to_keep)
likes_at_each_nodeIndex_branchBot = subvv(res.likes_at_each_nodeIndex_branchBot, statenums_to_keep)
normlikes_at_each_nodeIndex_branchBot = subvv(res.normlikes_at_each_nodeIndex_branchBot, statenums_to_keep)
uppass_probs_at_each_nodeIndex_branchBot = subvv(res.uppass_probs_at_each_nodeIndex_branchBot, statenums_to_keep)
anc_estimates_at_each_nodeIndex_branchBot = subvv(res.anc_estimates_at_each_nodeIndex_branchBot, statenums_to_keep)
uppass_probs_at_each_nodeIndex_branchTop = subvv(res.uppass_probs_at_each_nodeIndex_branchTop, statenums_to_keep)
anc_estimates_at_each_nodeIndex_branchTop = subvv(res.anc_estimates_at_each_nodeIndex_branchTop, statenums_to_keep)
fixNodesMult_at_each_nodeIndex_branchBot = subvv(res.fixNodesMult_at_each_nodeIndex_branchBot, statenums_to_keep)
fixNodesMult_at_each_nodeIndex_branchTop = subvv(res.fixNodesMult_at_each_nodeIndex_branchTop, statenums_to_keep)



# Subset inputs.p_Ds_v5
# inputs.p_Ds_v5.params
# inputs.p_Ds_v5.p_indices
# inputs.p_Ds_v5.p_TFs

# inputs.p_Ds_v5.uE
uE = p_Ds_v5.uE[statenums_to_keep]

# inputs.p_Ds_v5.p_indices
rn(inputs.p_Ds_v5.p_indices)

# Subset the Q matrix to only states you are keeping
Qdf = prtQp(inputs.p_Ds_v5)
keepTF = R_in(Qdf.i, statenums_to_keep);
TF = R_in(Qdf.j, statenums_to_keep);
keepTF[TF.==false] .= false;
Qdf[keepTF,:]
sum(keepTF)

# inputs.setup: Anagenetic transitions shortcuts - edit after re-doing anagenesis table
d_rows = setup.d_rows
d_froms = setup.d_froms
d_tos = setup.d_tos
d_drows = setup.d_drows
a_rows = setup.a_rows
a_froms = setup.a_froms
a_tos = setup.a_tos
a_arows = setup.a_arows
e_rows = setup.e_rows

# Areas gained and lost;
# when areas 3 and 4 are traits (H and D), no d events should add them
# when areas 3 and 4 are traits (H and D), no e events should lose them

gains = setup.gains
losses = setup.losses

dTF = Qdf.event .== "d"
gain_of_trait_TF = R_in_vv_ints(gains, trait_areanums)
prohibitTF = (dTF .+ gain_of_trait_TF) .== 2
keepTF[prohibitTF] .= false
sum(keepTF)
Qdf[keepTF,:]

eTF = Qdf.event .== "e"
loss_of_trait_TF = R_in_vv_ints(losses, trait_areanums)
prohibitTF = (eTF .+ loss_of_trait_TF) .== 2
keepTF[prohibitTF] .= false
sum(keepTF)


# The "a" events should have only gains/losses in trait states (3->4 and 4->3)
z = DataFrame(Rcbind(vv_to_v_ints(gains), vv_to_v_ints(losses)), :auto)
z = rename_df!(z, ["gains", "losses"])
Qdf_gl = Rcbind(Qdf, z)

aTF = Qdf.event .== "a"
loss_of_trait_TF = R_in_vv_ints(losses, trait_areanums)
gain_of_trait_TF = R_in_vv_ints(gains, trait_areanums)
a_events_to_keepTF = (aTF .+ loss_of_trait_TF .+ gain_of_trait_TF) .== 3
# Apply only to a events
keepTF[aTF][a_events_to_keepTF[aTF] .== false] .= false
sum(keepTF)  # no change, as the matching "a" events were all excluded based on statenums_to_keep

# Reduced Q matrix
Qdf_reduced = Qdf_gl[keepTF,:]

# Convert the i and j states to the new numbering
Qdf_reduced.i = ond.(Qdf_reduced.i);
Qdf_reduced.j = ond.(Qdf_reduced.j);
Qdf_reduced

# Add the necessary "a" events in the Q matrix (off-diagonal events)
num_nonzero_rates = length(statenums_to_keep) * (length(statenums_to_keep) - 1)

# Initialize empty arrays
Qarray_ivals = repeat(Int64[0], num_nonzero_rates)
Qarray_jvals = repeat(Int64[0], num_nonzero_rates)
Qarray_i_losses = repeat(Int64[0], num_nonzero_rates)
Qarray_j_gains = repeat(Int64[0], num_nonzero_rates)

Qij_vals =  repeat(Float64[0.0], num_nonzero_rates)  # 0-element Array{Any,1}  
Qij_vals_t =  repeat(Float64[0.0], num_nonzero_rates)  # 0-element Array{Any,1}  
										 # This is populated by calculating through the others
#base_vals = Array{Float64, num_nonzero_rates}  # base rates -- d, e, a
#mod_vals = Array{Float64, num_nonzero_rates}  # base modifiers -- e.g., "2" for AB -> ABC
#mult_vals = Array{Float64, num_nonzero_rates} # multipliers from mult_mat, e_mult 
# (ie modification by distance, area, etc.)
Qarray_event_types = repeat(String[""], num_nonzero_rates)
index = 0

numstates_new = length(states_list)
index = 0
for i in 1:(numstates_new-1)
	for j in (i+1):numstates_new
		areas_i = Vector{Int64}(states_list[i])
		areas_j = Vector{Int64}(states_list[j])
		
		# Must have same number of areas
		if (length(areas_i) != length(areas_j))
			continue
		end
		
		#shared_areas = intersect(areas_j, areas_i)
		#all_areas = union(areas_j, areas_i)
		diff_areas = sort!(symdiff(areas_j, areas_i))
		diff_areas = Vector{Int64}(diff_areas)
		
		# If there are 2 different areas, and they are both in the trait_areanums vector,
		# then it's a trait-transition "a" event:
		TF = R_in(diff_areas, trait_areanums);
		if (sum(TF) == length(TF))
			print(i)
			print(j)
			
			index = index + 1
			Qarray_event_types[index] = "a"
			Qarray_ivals[index] = i
			Qarray_jvals[index] = j
			Qarray_i_losses[index] = intersect(areas_i, diff_areas)[1]
			Qarray_j_gains[index] = intersect(areas_j, diff_areas)[1]
			Qij_vals[index] = 0.123
			Qij_vals_t[index] = 0.123
			
			# Reverse event, as well
			index = index + 1
			Qarray_event_types[index] = "a"
			Qarray_ivals[index] = j
			Qarray_jvals[index] = i
			Qarray_i_losses[index] = intersect(areas_j, diff_areas)[1]
			Qarray_j_gains[index] = intersect(areas_i, diff_areas)[1]
			Qij_vals[index] = 0.123
			Qij_vals_t[index] = 0.123

		end
	end
end

Qdf_add_As = DataFrame(event=Qarray_event_types, i=Qarray_ivals, j=Qarray_jvals, val=Qij_vals, vals_t=Qij_vals_t, gains=Qarray_j_gains, losses=Qarray_i_losses);
Qdf_add_As = Qdf_add_As[1:index,:]

Qdf_new = Rrbind(Qdf_reduced, Qdf_add_As)

# Q matrix DONE


# inputs.setup: Cladogenetic transitions shortcuts - edit after re-doing cladogenesis table
j_rows = setup.j_rows
j_froms = setup.j_froms
j_tos = setup.j_tos
j_jrows = setup.j_jrows
j_numdispersals = setup.j_numdispersals
v_rows = setup.v_rows
vicdist_base = setup.vicdist_base
vicdist = setup.vicdist
vicdist_t = setup.vicdist_t
s_rows = setup.s_rows

# inputs.setup: Other
mu_func = setup.mu_func
max_extinction_rate = setup.max_extinction_rate
multi_area_ranges_have_zero_mu = setup.multi_area_ranges_have_zero_mu
min_stepsize = setup.min_stepsize











