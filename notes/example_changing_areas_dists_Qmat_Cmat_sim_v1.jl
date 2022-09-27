#######################################################
# PhyBEARS model to simulator
#
# Design a simulation to produce something like the
# Psychotria tree (i.e. use the ML-inferred paramters
# on Psychotria to do a simulation)
# 
########################################################

using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using PhyBEARS
using DataFrames
using CSV

trfn = "/GitHub/PhyBEARS.jl/data/Psychotria_tree.newick"
tr = readTopology(trfn)

lgdata_fn = "/GitHub/PhyBEARS.jl/data/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
include_null_range = false
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr)
bmo.est[bmo.rownames .== "deathRate"] .= 0.2
bmo.est[bmo.rownames .== "d"] .= 0.034
bmo.est[bmo.rownames .== "e"] .= 0.028
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.00001
bmo.est[bmo.rownames .== "u"] .= -1.0
bmo.est[bmo.rownames .== "x"] .= -1.0

bmo.est .= bmo_updater_v1(bmo);

# Set up inputs 
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Es_v5, Es_tspan) = inputs;
numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])
root_age = maximum(trdf[!, :node_age])
p = p_Es_v5;
prtCp(p)


#######################################################
# An area-interpolator by hand:
#######################################################
times = [0.0, 4.0, 20.0, 22.0, 30.0, 80.0, 600.0]
# current NZ area = 268,021 km2
areas = [1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 10.0]

area_interpolator = interpolate((times,), areas, Gridded(Linear()));
tvals = seq(0.0, 40, 1.0);
area_interpolator(tvals)


#######################################################
# Area interpolator for a series of areas, 
# reading the areas in from a text file
#######################################################
areas_list = inputs.setup.areas_list;
states_list = inputs.setup.states_list;
area_of_areas_file = "/GitHub/PhyBEARS.jl/notes/area_of_areas_Hawaii_1s_v2.txt";
area_of_areas_table = CSV.read(area_of_areas_file, DataFrame; delim="\t");
area_of_areas_table = readtable(area_of_areas_file);
area_of_areas_vec = PhyBEARS.TimeDep.area_of_areas_df_to_vectors(area_of_areas_table);

# Check if the area of areas vector has the same number of areas as the geog_df
if (Rncol(geog_df)-1) != length(area_of_areas_vec[1])
	txt = paste0(["STOP ERROR. The area_of_areas_table, used as the input to area_of_areas_df_to_vectors() and then interpolate(), must have the same number of areas as the geographic ranges file. Fix this before proceeding.\n\nSpecifically, these should match, but don't: length(area_of_areas_vec[1])=", length(area_of_areas_vec[1]), ". Rncol(geog_df)-1=", Rncol(geog_df)-1, ""])
	print(txt)
	error(txt)
end

area_of_areas_interpolator = interpolate((area_of_areas_table.times,), area_of_areas_vec, Gridded(Linear()))

# Calculate the area at different times
tval = 19.0
state_as_areas_list = [1,2];
area_of_areas = area_of_areas_interpolator(tval);
total_area = PhyBEARS.TimeDep.get_area_of_range(tval, state_as_areas_list, area_of_areas);

get_areas_of_range_at_t = x -> PhyBEARS.TimeDep.get_area_of_range_using_interpolator(x, state_as_areas_list, area_of_areas_interpolator);
tvals = seq(18.0, 23.0, 0.25)
get_areas_of_range_at_t.(tvals)


#######################################################
# Distances interpolator for a series of distances
#######################################################
dists1 = [0.00 0.33 0.58 1.00;
0.33 0.00 0.25 0.67;
0.58 0.25 0.00 0.17;
1.00 0.67 0.17 0.00];
dists2 = [0.00 0.17 0.29 0.50;
0.17 0.00 0.13 0.33;
0.29 0.13 0.00 0.08;
0.50 0.33 0.08 0.00];
dists3 = [0.00 0.08 0.15 0.25;
0.08 0.00 0.06 0.17;
0.15 0.06 0.00 0.04;
0.25 0.17 0.04 0.00];
dists4 = [0.00 0.04 0.07 0.13;
0.04 0.00 0.03 0.08;
0.07 0.03 0.00 0.02;
0.13 0.08 0.02 0.00];
dists5 = [0.00 0.04 0.07 0.13;
0.04 0.00 0.03 0.08;
0.07 0.03 0.00 0.02;
0.13 0.08 0.02 0.00];

changing_distances = [dists1, dists2, dists3, dists4, dists5];

times = [0.0, 1.5, 3.0, 5.5, 600];

distances_interpolator = interpolate((times,), changing_distances, Gridded(Linear()));
tvals = seq(0.0, 10, 1.0);
distances_interpolator(tvals)


#######################################################
# Set up an inference
#######################################################

# Parameters vectors & model structure in "p_Es_v10"
p_Es_v10 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, area_of_areas_interpolator=area_of_areas_interpolator, distances_interpolator=distances_interpolator, use_distances=true, bmo=bmo);

p=p_Es_v10;
p.distances_interpolator.(tvals);

# Calculate the Es (extinction probabilities at any time t, for each state/range)
prob_Es_v10 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v10_simd_sums, p_Es_v10.uE, Es_tspan, p_Es_v10);
# This solution is an interpolator
sol_Es_v10 = solve(prob_Es_v10, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v10;

# Check the interpolator
Es_interpolator(1.0)
Es_interpolator(1.5)


# Calculate the Ds (likelihoods of tip data/diversity)
# parameters / model structure, with the Es interpolator added
p_Ds_v10 = (n=p_Es_v10.n, params=p_Es_v10.params, p_indices=p_Es_v10.p_indices, p_TFs=p_Es_v10.p_TFs, uE=p_Es_v10.uE, terms=p_Es_v10.terms, setup=p_Es_v10.setup, states_as_areas_lists=p_Es_v10.states_as_areas_lists, area_of_areas_interpolator=p_Es_v10.area_of_areas_interpolator, distances_interpolator=p_Es_v10.distances_interpolator, bmo=p_Es_v10.bmo, sol_Es_v10=sol_Es_v10);


# Log-likelihood calculation is slow the first time, but then speeds up
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)









#######################################################
# Note that the area of areas changes, as Hawaiian 
# islands sink as we go back in time
#######################################################
tvals = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
area_of_areas_interpolator(tvals)

#######################################################
# Get the Q matrix rates, and C matrix rates, at time t=0.0
#######################################################
t = 0.0

# At time t, get the size of areas, and distances between areas
p.setup.area_of_areas .= p.area_of_areas_interpolator(t)
p.setup.distmat .= p.distances_interpolator(t)

# Update extinction rates (mu_t)
max_extinction_rate = p.setup.max_extinction_rate
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


# Update the Qmatrix
# Get the e_vals for the Qij matrix, at time t
# elist_actual = elist_base * area_of_area_lost^u_e
update_Qij_e_vals!(p);

# Get the d_vals for the Qij matrix, at time t
# Update the dmat (dmat contains all the dispersal rates for 
# individual pairs of areas; size numareas x numareas
# It is modified by the distance_matrix^u, etc.
update_Qij_d_vals!(p);
##prtQp(p)

# Update the Cmatrix
# Using the current t's distmat, etc. update the jmat_t, then
# propagate through the C matrix
##x1 = prtCp(p).rates_t
update_Cijk_j_rates!(p);

# Save the matrices at time 0.0
orig_mus = deepcopy(p.params.mu_t_vals)
orig_Qmat = deepcopy(prtQp(p))
orig_Cmat = deepcopy(prtCp(p))




#######################################################
# Get the Q matrix rates, and C matrix rates, at time t=5.0
#######################################################
t = 5.0

# At time t, get the size of areas, and distances between areas
p.setup.area_of_areas .= p.area_of_areas_interpolator(t)
 p.setup.distmat .= p.distances_interpolator(t)

# Update extinction rates (mu_t)
max_extinction_rate = p.setup.max_extinction_rate
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


# Update the Qmatrix
# Get the e_vals for the Qij matrix, at time t
# elist_actual = elist_base * area_of_area_lost^u_e
update_Qij_e_vals!(p);

# Get the d_vals for the Qij matrix, at time t
# Update the dmat (dmat contains all the dispersal rates for 
# individual pairs of areas; size numareas x numareas
# It is modified by the distance_matrix^u, etc.
update_Qij_d_vals!(p);
##prtQp(p)

# Update the Cmatrix
# Using the current t's distmat, etc. update the jmat_t, then
# propagate through the C matrix
##x1 = prtCp(p).rates_t
update_Cijk_j_rates!(p);

# Save the matrices at time 0.0
new_mus = deepcopy(p.params.mu_t_vals)
new_Qmat = deepcopy(prtQp(p))
new_Cmat = deepcopy(prtCp(p))



#######################################################
# Compare t=0.0 and t=1.234
# 
#######################################################

# Extinction rates of each of the 16 states, at time 0.0 and 5.0
orig_mus
new_mus

# Note how this matches the areas:
# -- extinction rates increase when an area or range gets smaller
area_of_areas_interpolator(tvals)



# Q matrix -- first rows, note that vals_t gets updated
orig_Qmat[1:10,:]
new_Qmat[1:10,:]

# Q matrix -- last rows
rownums = collect((Rnrow(orig_Qmat)-10):Rnrow(orig_Qmat))
orig_Qmat[rownums,:]
new_Qmat[rownums,:]
# Note how the "e" -- range loss rate -- increases dramatically
# for losing areas that are sunk to minimum size (0.00001)

# C matrix -- first rows, note that rates_t gets updated,
#             BUT ONLY FOR THE "J" EVENTS
# (a new model could make e.g. "v"/vicariance depend on
#  distance)
orig_Cmat[1:10,:]
new_Cmat[1:10,:]

# C matrix -- rows, note that rates_t gets updated
# Note that the rates_t doesn't change for the last
# row, because these aren't jump dispersals
rownums = collect((Rnrow(orig_Cmat)-10):Rnrow(orig_Cmat))
orig_Cmat[rownums,:]
new_Cmat[rownums,:]


