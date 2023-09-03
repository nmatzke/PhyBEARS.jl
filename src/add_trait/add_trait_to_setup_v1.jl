


# Double inputs.setup for a trait
trait_states_txt = ["H", "D"]   # hermaphrodite, dioecy
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
gains = setup.gains
losses = setup.losses

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











