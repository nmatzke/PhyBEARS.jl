
#######################################################
# Example inference on a simulated tree & geography dataset
# (from Wallis's ss8_sim_001)
# 2022-11-09
#
# We will compare inference under 
# * SSE likelihood calculator v7 (constant rates through time, very fast)
# * SSE likelihood calculator v12 (changing rates through time, very fast)
#   - Here, to start, we will only change the distances through time
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

# Change the working directory as needed
wd = "/GitHub/PhyBEARS.jl/ex/siminf_v12a/"
cd(wd)

# This simulation has 50 living species
trfn = "tree.newick"
tr = readTopology(trfn)
trdf = prt(tr)
oldest_possible_age = 100.0

lgdata_fn = "rangedata.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
include_null_range = true
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr)
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.034
bmo.est[bmo.rownames .== "e"] .= 0.028
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
bmo.est[bmo.rownames .== "u"] .= 0.0
bmo.est[bmo.rownames .== "x"] .= 0.0
bmo.est .= bmo_updater_v1(bmo);

bmo.est[bmo.rownames .== "birthRate"] .= 0.1
bmo.est[bmo.rownames .== "deathRate"] .= 0.1


# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

p = p_Ds_v5
solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-12;
solver_options.reltol = 1e-12;

#######################################################
# Load Wallis's distances
#######################################################
Times_3 = repeat([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],3)
Land1_3 = repeat(["test1","test1","test2"],inner = 11)
Land2_3 = repeat(["test2","test3","test3"],inner = 11)
Distances_km_3 = [0, 1000, 2000, 3500, 4100, 4900, 6000, 5800, 7000, 7500, 8000,0, 2000, 3000, 4000, 5100, 5900, 7000, 6800, 8000, 8500, 9000,0, 500, 1000, 2500, 3100, 3900, 5000, 4800, 6000, 6500, 7000] 
testdf = DataFrame(Times = repeat([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],3), Land1 = repeat(["test1","test1","test2"],inner = 11), Land2 = repeat(["test2","test3","test3"],inner = 11), Distances_km = [0, 1000, 2000, 3500, 4100, 4900, 6000, 5800, 7000, 7500, 8000,0, 2000, 3000, 4000, 5100, 5900, 7000, 6800, 8000, 8500, 9000,0, 500, 1000, 2500, 3100, 3900, 5000, 4800, 6000, 6500, 7000])

# Let's put in the form of distance matrices
# This will create 11 distances matrices, of size 3x3
numareas = 3
unique_times = 1.0 .* sort(unique(testdf.Times))
distmats = [Array{Float64}(undef, numareas, numareas) for _ = 1:(length(unique_times)+1)]

# Fill in the distmats
for i in 1:(length(unique_times))
	distmats[i] .= 0.0  # zero-out the distances matrix
	TF = testdf.Times .== unique_times[i]
	tmpdf = testdf[TF,:]
	for j in 1:Rnrow(tmpdf)
		land1 = tmpdf.Land1[j]
		land2 = tmpdf.Land2[j]
		area1_num = extract_first_integer_from_string(land1)
		area2_num = extract_first_integer_from_string(land2)
		distmats[i][area1_num, area2_num] = convert(Float64, testdf.Distances_km[i])
		distmats[i][area2_num, area1_num] = convert(Float64, testdf.Distances_km[i])
	end
end

# Add one last distmat on the end to cover the bottom of the tree
distmats[length(distmats)] .= distmats[length(distmats)-1]
push!(unique_times, oldest_possible_age)
sort!(unique_times)
times = unique_times

show(distmats)

# Let's divide the distances by the minimum nonzero distance
dists_vector = collect(Iterators.flatten(vec.(distmats)))
dists_vector2 = dists_vector[dists_vector .> 0.0]
maxval = maximum(dists_vector2)

for i in 1:length(distmats)
	distmats[i] .= distmats[i] ./ maxval
end
distmats

#######################################################
# Distances interpolator for a series of distances
#######################################################
distances_interpolator = interpolate((times,), distmats, Gridded(Linear()));
tvals = seq(0.0, oldest_possible_age, 1.0);
distances_interpolator(tvals)





#######################################################
# Changing (or not) area of areas
#######################################################
times = unique_times
# Set area of areas to 1.0 for now
area_of_areas = [Vector{Float64}(undef, numareas) for _ = 1:(length(unique_times))]
for i in 1:length(area_of_areas)
	area_of_areas[i] .= 1.0
end

area_interpolator = interpolate((times,), area_of_areas, Gridded(Linear()))
tvals = seq(0.0, maximum(times), 1.0)
area_interpolator(tvals)


#######################################################
# Area interpolator for a series of states
#######################################################
areas_list = inputs.setup.areas_list
states_list = inputs.setup.states_list
#area_of_areas_file = "/GitHub/PhyBEARS.jl/notes/area_of_areas_Hawaii_1s_v2.txt"
#area_of_areas_table = CSV.read(area_of_areas_file, DataFrame; delim="\t")
#area_of_areas_table = readtable(area_of_areas_file)

###area_of_areas_vec = area_of_areas_df_to_vectors(area_of_areas_table; cols=2)
#area_of_areas_vec = PhyBEARS.TimeDep.area_of_areas_df_to_vectors(area_of_areas_table)
area_of_areas_vec = area_of_areas


# Check if the area of areas vector has the same number of areas as the geog_df
if (Rncol(geog_df)-1) != length(area_of_areas_vec[1])
	txt = paste0(["STOP ERROR. The area_of_areas_table, used as the input to area_of_areas_df_to_vectors() and then interpolate(), must have the same number of areas as the geographic ranges file. Fix this before proceeding.\n\nSpecifically, these should match, but don't: length(area_of_areas_vec[1])=", length(area_of_areas_vec[1]), ". Rncol(geog_df)-1=", Rncol(geog_df)-1, ""])
	print(txt)
	error(txt)
end

#area_of_areas_interpolator = interpolate((area_of_areas_table.times,), area_of_areas_vec, Gridded(Linear()))
area_of_areas_interpolator = interpolate((times,), area_of_areas_vec, Gridded(Linear()))

# Calculate the area at different times
tval = 10.0
state_as_areas_list = states_list[8]
area_of_areas = area_of_areas_interpolator(tval)
total_area = PhyBEARS.TimeDep.get_area_of_range(tval, state_as_areas_list, area_of_areas)

get_areas_of_range_at_t = x -> PhyBEARS.TimeDep.get_area_of_range_using_interpolator(x, state_as_areas_list, area_of_areas_interpolator)
tvals = seq(1.0, 5.0, 1.0)
get_areas_of_range_at_t.(tvals)






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


# Set up inputs 
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;
numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])
root_age = maximum(trdf[!, :node_age])

p_Es_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, area_of_areas_interpolator=area_of_areas_interpolator, distances_interpolator=distances_interpolator, vicariance_mindists_interpolator=vicariance_mindists_interpolator, use_distances=true, bmo=bmo)


p_Es_v10 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, area_of_areas_interpolator=area_of_areas_interpolator, distances_interpolator=distances_interpolator, vicariance_mindists_interpolator=vicariance_mindists_interpolator, use_distances=true, bmo=bmo)


p2 = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v10, times);


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


(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)



files.times_fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/v12a_times.txt"
files.distances_fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/v12a_distances.txt"
files.area_of_areas_fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/v12a_area_of_areas.txt"

interpolators = files_to_interpolators(files, p.setup.numareas, p.setup.states_list, p.setup.v_rows, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals; oldest_possible_age=1000.0)

# include("/Users/nickm/GitHub/PhyBEARS.jl/ex/siminf_v12a/example_inf_v12a.jl")


