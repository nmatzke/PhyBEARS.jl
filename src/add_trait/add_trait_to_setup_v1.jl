


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











