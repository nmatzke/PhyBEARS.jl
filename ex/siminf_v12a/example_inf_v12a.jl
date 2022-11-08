
#######################################################
# Interpolate a Carray matrix
#
# Notes: Saving interpolated function into a separate file in Julia
# https://stackoverflow.com/questions/58499330/saving-interpolated-function-into-a-separate-file-in-julia
#######################################################

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
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.0
bmo.est[bmo.rownames .== "e"] .= 0.0
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.11
bmo.est[bmo.rownames .== "u"] .= 0.0
bmo.est[bmo.rownames .== "x"] .= 0.0
bmo.est .= bmo_updater_v1(bmo);

bmo.est[bmo.rownames .== "d"] .= 0.034
bmo.est[bmo.rownames .== "e"] .= 0.028
bmo.est[bmo.rownames .== "j"] .= 0.00001

bmo.est[bmo.rownames .== "birthRate"] .= 0.1
bmo.est[bmo.rownames .== "deathRate"] .= 0.1

bmo.est[bmo.rownames .== "d"] .= 0.034
bmo.est[bmo.rownames .== "e"] .= 0.028
bmo.est[bmo.rownames .== "j"] .= 0.11
bmo.est[bmo.rownames .== "xv"]

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;

# setup.vicdist_base  # a list of distances between ranges (e.g. minimum distance), 
# corresponding to v_rows of the Carray
Rnames(setup)
setup.vicdist_base

p = p_Ds_v5
solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-12;
solver_options.reltol = 1e-12;

p_Ds_v5 = inputs.p_Ds_v5;
inputs.setup.dmat_base
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
inputs.setup.dmat_base
prtCp(p_Ds_v5)


#######################################################
# OK let's make a function that interpolates different rates
# based on changing distances / areas
#######################################################
times = [0.0, 4.0, 20.0, 22.0, 30.0, 80.0, 600.0]
# current NZ area = 268,021 km2
areas = [1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 10.0]

area_interpolator = interpolate((times,), areas, Gridded(Linear()))
tvals = seq(0.0, 40, 1.0)
area_interpolator(tvals)


#######################################################
# Area interpolator for a series of states
#######################################################
areas_list = inputs.setup.areas_list
states_list = inputs.setup.states_list
area_of_areas_file = "/GitHub/PhyBEARS.jl/notes/area_of_areas_Hawaii_1s_v2.txt"
area_of_areas_table = CSV.read(area_of_areas_file, DataFrame; delim="\t")
area_of_areas_table = readtable(area_of_areas_file)

#area_of_areas_vec = area_of_areas_df_to_vectors(area_of_areas_table; cols=2)
area_of_areas_vec = PhyBEARS.TimeDep.area_of_areas_df_to_vectors(area_of_areas_table)

# Check if the area of areas vector has the same number of areas as the geog_df
if (Rncol(geog_df)-1) != length(area_of_areas_vec[1])
	txt = paste0(["STOP ERROR. The area_of_areas_table, used as the input to area_of_areas_df_to_vectors() and then interpolate(), must have the same number of areas as the geographic ranges file. Fix this before proceeding.\n\nSpecifically, these should match, but don't: length(area_of_areas_vec[1])=", length(area_of_areas_vec[1]), ". Rncol(geog_df)-1=", Rncol(geog_df)-1, ""])
	print(txt)
	error(txt)
end

area_of_areas_interpolator = interpolate((area_of_areas_table.times,), area_of_areas_vec, Gridded(Linear()))

# Calculate the area at different times
tval = 19.0
state_as_areas_list = [1,2]
area_of_areas = area_of_areas_interpolator(tval)
total_area = PhyBEARS.TimeDep.get_area_of_range(tval, state_as_areas_list, area_of_areas)

get_areas_of_range_at_t = x -> PhyBEARS.TimeDep.get_area_of_range_using_interpolator(x, state_as_areas_list, area_of_areas_interpolator)
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

# dists1 .= 1.0;
# dists2 .= 1.0;
# dists3 .= 1.0;
# dists4 .= 1.0;
# dists5 .= 1.0;


changing_distances = [dists1, dists2, dists3, dists4, dists5];

times = [0.0, 1.5, 3.0, 5.5, 600];

distances_interpolator = interpolate((times,), changing_distances, Gridded(Linear()));
tvals = seq(0.0, 10, 1.0);
distances_interpolator(tvals)



#######################################################
# Construct vicariance minimum-distances interpolator
# (for the v_rows of the matrix)
#######################################################
# A Vector of # times Vectors, each of length v_rows
# THIS AVOIDS THE LINKED VECTORS ISSUE
changing_mindists = [Vector{Float64}(undef, length(setup.v_rows)) for _ = 1:length(times)]

# Populate the minimum distances
for i in 1:length(times)
	distmat = distances_interpolator(times[i])
	update_min_vdist_at_time_t!(changing_mindists[i], setup.v_rows, distmat, setup.states_list, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals; mindist=1.0e9)
end

vicariance_mindists_interpolator = interpolate((times,), changing_mindists, Gridded(Linear()));
vicariance_mindists_interpolator(tvals)



# 1. create xv parameter
# 2. create function to interpolate Carray and Qarray
# 3. put these interpolators into SSEs v11






# Set up inputs 
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Es_v5, Es_tspan) = inputs;
numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])
root_age = maximum(trdf[!, :node_age])

inputs.setup.dmat_base
p_Es_v5.params.Qij_vals
p_Es_v5.params.Qij_vals_t


prtCp(p_Es_v5)

# Problem: this simulation allowed null-range speciation (null->null, null or 000->000,000)
# Solution for now: Add the null -> null, null cladogenesis event to the p_Es_v5
p_Es_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, area_of_areas_interpolator=area_of_areas_interpolator, distances_interpolator=distances_interpolator, vicariance_mindists_interpolator=vicariance_mindists_interpolator, use_distances=true, bmo=bmo)

birthRate = bmo.est[bmo.rownames .== "birthRate"];
#add_111_to_Carray!(p_Es_v5, birthRate);

prtCp(p_Es_v5) # now 1->1,1 is an allowed cladogenesis event



p_Es_v10 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, area_of_areas_interpolator=area_of_areas_interpolator, distances_interpolator=distances_interpolator, vicariance_mindists_interpolator=vicariance_mindists_interpolator, use_distances=true, bmo=bmo)
Rnames(p_Es_v10)
prtCp(p_Es_v10)

p = p_Es_v10

p.bmo.est[p.setup.bmo_rows.u_e] = -1.0
p.bmo.est[p.setup.bmo_rows.x] = -1.0
p.bmo.est[p.setup.bmo_rows.xv] = 1.0

prtQp(p).vals_t
update_QC_mats_time_t!(p, 0.0)
prtQp(p).vals_t

prtQp(p).vals_t
update_QC_mats_time_t!(p, 1.0)
prtQp(p).vals_t

prtQp(p).vals_t
update_QC_mats_time_t!(p, 3.0)
prtQp(p).vals_t


prtCp(p).rates_t
update_QC_mats_time_t!(p, 0.0)
prtCp(p).rates_t


prtCp(p).rates_t
update_QC_mats_time_t!(p, 1.0)
prtCp(p).rates_t



prtCp(p).rates_t
update_QC_mats_time_t!(p, 3.0)
prtCp(p).rates_t


p2 = construct_QC_interpolators(p, tvals);


# Interpolators
p2.Q_vals_interpolator(0.0)[1:3]
p2.Q_vals_interpolator(1.0)[1:3]
p2.Q_vals_interpolator(2.0)[1:3]
p2.Q_vals_interpolator(3.0)[1:3]

p2.C_rates_interpolator(0.0)[1:3]
p2.C_rates_interpolator(1.0)[1:3]
p2.C_rates_interpolator(2.0)[1:3]
p2.C_rates_interpolator(3.0)[1:3]

p_Es_v12 = p2


# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")

du = collect(repeat([0.0], numstates)); 
u = collect(repeat([0.0], numstates));
p = p_Es_v10;
t = 0.0;
PhyBEARS.SSEs.parameterized_ClaSSE_Es_v10_simd_sums(du, u, p, t);
du
du = collect(repeat([0.0], numstates));
u = collect(repeat([0.0], numstates)) ;
p = p_Es_v10;
t = 0.0;
PhyBEARS.SSEs.parameterized_ClaSSE_Es_v11_simd_sums(du, u, p, t);
du

prob_Es_v10 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v10_simd_sums, p_Es_v10.uE, Es_tspan, p_Es_v10);
# This solution is an interpolator
sol_Es_v10 = solve(prob_Es_v10, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v10;
p_Ds_v7 = (n=p_Es_v10.n, params=p_Es_v10.params, p_indices=p_Es_v10.p_indices, p_TFs=p_Es_v10.p_TFs, uE=p_Es_v10.uE, terms=p_Es_v10.terms, setup=p_Es_v10.setup,states_as_areas_lists=p_Es_v10.states_as_areas_lists, area_of_areas_interpolator=p_Es_v10.area_of_areas_interpolator, vicariance_mindists_interpolator=p_Es_v10.vicariance_mindists_interpolator, bmo=p_Es_v10.bmo, sol_Es_v5=sol_Es_v10);

# Check the interpolator
p_Ds_v7.sol_Es_v5(1.0)
Es_interpolator(1.0)



prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
# This solution is an interpolator
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v12;





# Calculate the Ds & total lnL via downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

p_Ds_v10 = (n=p_Es_v10.n, params=p_Es_v10.params, p_indices=p_Es_v10.p_indices, p_TFs=p_Es_v10.p_TFs, uE=p_Es_v10.uE, terms=p_Es_v10.terms, setup=p_Es_v10.setup, states_as_areas_lists=p_Es_v10.states_as_areas_lists, area_of_areas_interpolator=p_Es_v10.area_of_areas_interpolator, distances_interpolator=p_Es_v10.distances_interpolator, vicariance_mindists_interpolator=p_Es_v10.vicariance_mindists_interpolator, bmo=p_Es_v10.bmo, sol_Es_v10=sol_Es_v10);


p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, area_of_areas_interpolator=p_Es_v12.area_of_areas_interpolator, distances_interpolator=p_Es_v12.distances_interpolator, vicariance_mindists_interpolator=p_Es_v12.vicariance_mindists_interpolator, Q_vals_interpolator=p_Es_v12.Q_vals_interpolator, C_rates_interpolator=p_Es_v12.C_rates_interpolator, bmo=p_Es_v12.bmo, sol_Es_v12=sol_Es_v12);


# Use ONLY with add_111
#p_Ds_v10.params.Cijk_rates[1] = 0.0;
#p_Ds_v10.params.Cijk_rates_t[1] = 0.0;
p = p_Ds_v10;

all(abs.(prtQp(p_Ds_v10).val .- prtQp(p_Ds_v7).val) .< 1e-20)
all(abs.(prtQp(p_Ds_v10).vals_t .- prtQp(p_Ds_v7).val) .< 1e-20)
all(abs.(prtCp(p_Ds_v10).rate .- prtCp(p_Ds_v5).rate) .< 1e-6)

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

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


# WORKS
@benchmark PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)






#######################################################
# Conclusions on calculation time changes
#
# deathRate = 0.0, u=-1.0  - Fails as null range with 0.0^-1.0 gives Inf
# deathRate = 0.0, u= 1.0  - Time  (median):     30.675 ms  
# deathRate = 0.01,u= 1.0  - Time  (median):     31.224 ms
# deathRate = 0.1, u=-1.0 -  Time  (median):     39.746 ms 
# deathRate = 0.1, u=-1.0 - Time  (median):     111.056 ms 
#######################################################







p_Ds_v10.bmo.est[p_Ds_v10.bmo.rownames .== "u"] .= 1.0

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

p_Ds_v10.bmo.est[p_Ds_v10.bmo.rownames .== "u"] .= -10.0

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)




using PProf
@pprof PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v10!(res; trdf=trdf, p_Ds_v10=p_Ds_v10, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)




#######################################################
# Now make dispersal events dependent on distance
#######################################################



#######################################################
# Now make "e" events depend on "u"
#######################################################
# Update elist, i.e. the multiplier on the base e rate (fixed time)
inputs.setup.elist .= inputs.setup.elist_base .* inputs.setup.area_of_areas.^bmo.est[bmo.rownames .== "u"][1]

# Update elist_t, i.e. the multiplier on the base e rate (at time t)
p = p_Ds_v10;
tval = 5.1
PhyBEARS.TimeDep.update_Qij_e_vals!(p);
p.params.Qij_vals_t
p.setup.elist_t .= p.setup.elist_base .* area_of_areas_interpolator(tval) .^ bmo.est[p.setup.u_row]

# Update the Qmat, using elist_t
prtQp(p)
Rnames(p.p_indices)
e_rows = (1:length(p.p_indices.Qarray_event_types))[p.p_indices.Qarray_event_types .== "e"]

for i in 1:length(e_rows)
	starting_statenum = p.p_indices.Qarray_ivals[e_rows[i]]
	ending_statenum = p.p_indices.Qarray_jvals[e_rows[i]]
	area_lost = symdiff(inputs.setup.states_list[starting_statenum], inputs.setup.states_list[ending_statenum])
	# actual rate of e = base_rate_of_e * area_of_area_lost ^ u
	p.params.Qij_vals_t[e_rows[i]] = p.params.Qij_vals[e_rows[i]] * inputs.setup.elist_t[area_lost]
end


p.params.Qij_vals_t

# Now make an extinction_rate_interpolator
function get_extinction_rate_multiplier(tval, uval, state_as_areas_list, area_of_areas_interpolator)
	extinction_rate_multiplier = 1.0 * (get_area_of_range_using_interpolator(tval, state_as_areas_list, area_of_areas_interpolator) ^ uval)
	# Error check
	extinction_rate_multiplier = minimum([extinction_rate_multiplier, 10000.0])
	return(extinction_rate_multiplier)
end

get_extinction_rate_multipliers = x -> get_extinction_rate_multiplier(x, uval, state_as_areas_list, area_of_areas_interpolator)

# https://docs.juliahub.com/JLD2/O1EyT/0.1.13/
cd("/GitHub/PhyBEARS.jl/src")
@save "list_of_interpolators.jld2" list_of_interpolators

# Reloading
using JLD2
using Interpolations # for Interpolations.scale, Interpolations.interpolate
@load "list_of_interpolators.jld2" list_of_interpolators # <- It loads to "list_of_interpolators"...has to be this

# Find location of a module
pathof(PhyBEARS)
interp_jld2_path = join([pathof(PhyBEARS), ])

itp = interpolate((x,), y, Gridded(Linear()))

