#######################################################
# Compare the speeds of tradSSE and Gflow versions on 
# a larger (197 species, 9 areas, 512 states) dataset
# (Podocarpaceae + Araucariaceae dataset, from Klaus & 
#  Matzke 2020, Systematic Biology)
#
# Result:
# total_calctime_in_sec_nFv6 -- traditional approach, SSE on each branch
#  seconds
# total_calctime_in_sec_GFv6 -- Louca & Pennell approach, Gflow matrices saved
#  seconds
# 
# lnLs match to <0.1 with abstol and reltol set to 1e-9 (but not 1e-6
#######################################################
#Pkg.instantiate()

"""
# Run with:
julia -t auto -p auto

cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/speedtests_Cyrtandra_wExtinction+J_v2speed.jl")
"""


using Distributed
Distributed.nprocs()

using Hwloc
Hwloc.num_physical_cores()
Hwloc.num_virtual_cores()


using Dates									# for e.g. Dates.now(), DateTime
using DataFrames
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames						# for DataFrame()
using DelimitedFiles				# for readdlm()
using NLopt									# seems to be the best gradient-free, box-constrained								
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Base.Threads  # <-- this is good for @spawn, Distributed.@spawn is BAD, it produces Futures that have to be scheduled etc]

# List each PhyBEARS code file prefix here
using PhyloBits
using PhyBEARS
using PhyloBits.TrUtils			# for e.g. numstxt_to_df()
using PhyloBits.TreeTable
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.ModelLikes # e.g. setup_DEC_SSE2
using PhyBEARS.Uppass

correct_results = (1.803, 12, -295.10996921099974, -11.470247869093797, -306.5802170800935, -119.7379445017597)
Julia_sum_lq 			= correct_results[3]
rootstates_lnL 		= correct_results[4]
Julia_total_lnLs1 = correct_results[5]
bgb_lnL 					= correct_results[6]


# List each PhyBEARS code file prefix here
doesnt_work="""
@everywhere using Pkg
@everywhere cd("/GitHub/PhyBEARS.jl")
@everywhere Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl")) # Gives error:
# ERROR: TaskFailedException
#   nested task error: TOML Parser error:
#   none:2:5 error: expected equal sign after key
# ...but the next line still works.
@everywhere using PhyBEARS
"""


#@testset "speedtests_Cyrtandra_wExtinction+J_v2speed.jl" begin

#include("/GitHub/PhyBEARS.jl/src/Gmaps.jl")
#import .Gmaps
#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes
#include("/GitHub/PhyBEARS.jl/notes/jl")
#import .Flow

#trfn = "/GitHub/PhyBEARS.jl/data/Klaus_Matzke_2020_PodoArau_197sp.newick"
#tr = readTopology(trfn)
#lgdata_fn = "/GitHub/PhyBEARS.jl/data/Podocarpaceae_197_9areas_5Araucariaceae.data"
#geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
trfn = "/GitHub/PhyBEARS.jl/data/Cyrtandra.newick"
tr = readTopology(trfn)

lgdata_fn = "/GitHub/PhyBEARS.jl/data/Cyrtandra_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.1
bmo.est[bmo.rownames .== "deathRate"] .= 0.01
bmo.est[bmo.rownames .== "d"] .= 0.02
bmo.est[bmo.rownames .== "e"] .= 0.02
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.11
bmo.est[:] = bmo_updater_v1_SLOW(bmo);
numareas = 7
n = 2^numareas          # 4 areas, 16 states

# Set up the model
root_age_mult=1.5; max_range_size=NaN; include_null_range=true; max_range_size=NaN
max_range_size = NaN # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;
vfft(res.likes_at_each_nodeIndex_branchTop)

#@everywhere trdf = trdf

solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-8;
solver_options.reltol = 1e-8;

p_Es_v5 = inputs.p_Es_v5;
p_Ds_v5_updater_v1!(p_Es_v5, inputs);

# Solve the Es
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Es_v5.uE, Es_tspan, p_Es_v5);
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, sol_Es_v5=sol_Es_v5);

res_nonFlow_v5 = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv5, iteration_number_nFv5, Julia_sum_lq_nFv5, rootstates_lnL_nFv5, Julia_total_lnLs1_nFv5, bgb_lnl_nFv5) = res_nonFlow_v5

res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6, iteration_number_nFv6, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6

res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7

res_nonFlow_v5 = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv5, iteration_number_nFv5, Julia_sum_lq_nFv5, rootstates_lnL_nFv5, Julia_total_lnLs1_nFv5, bgb_lnl_nFv5) = res_nonFlow_v5

res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6, iteration_number_nFv6, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6

res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7


# 2022-03-29 Podocarps: (324.287, 22, -1271.480048120258, -14.85776329108454, -1286.3378114113425, -576.3087517374857)
# 2022-03-29 Cyrtandra: (21.556, 12, -295.1099959548604, -11.470227634304381, -306.58022358916475, -119.7379632144717)
# 2022-03-30 Cyrtandra: (23.295, 12, -295.1099959548604, -11.470227634304381, -306.58022358916475, -119.7379632144717)
# 2022-03-30 Cyrtandra: (24.813, 12, -295.1099959548604, -11.470227634304381, -306.58022358916475, -119.7379632144717)

#######################################################
# Parallelized -- turn off for default tests
#######################################################

# res_nonFlow_v6par = iterative_downpass_parallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
# (total_calctime_in_sec_nFv6par, iteration_number_nFv6par, Julia_sum_lq_nFv6par, rootstates_lnL_nFv6par, Julia_total_lnLs1_nFv6par, bgb_lnl_nFv6par) = res_nonFlow_v6par

# res_nonFlow_v6par = iterative_downpass_parallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
# (total_calctime_in_sec_nFv6par, iteration_number_nFv6par, Julia_sum_lq_nFv6par, rootstates_lnL_nFv6par, Julia_total_lnLs1_nFv6par, bgb_lnl_nFv6par) = res_nonFlow_v6par
# 2022-03-30 Cyrtandra: (24.011, 12, -295.1099959548604, -11.470227634304381, -306.58022358916475, -119.7379632144717)
# 2022-03-30 Cyrtandra: (21.178, 12, -295.1099959548604, -11.470227634304381, -306.58022358916475, -119.7379632144717)
# After changing from Distributed.@spawn to Base.Threads.@spawn:
# 2022-03-30 Cyrtandra: (6.027, 99251, -295.1099959548604, -11.470227634304381, -306.58022358916475, -119.7379632144717)


solver_options.abstol = 1e-15;
solver_options.reltol = 1e-15;

# Version 7/2 ClaSSE Gflow calculations
G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 1.01 * maximum(trdf.node_age))
prob_Gs_v5 = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE, G0, tspan, pG);
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6
archived_Gflow_v6 = deepcopy(res);
# 2022-03-29: (7.323, 22, -1271.0625625625657, -14.712302169676475, -1285.7748647322421, -575.9248397325294)
# 2022-03-30: (total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6

nd = 3
@test round(res_nonFlow_v5[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_nonFlow_v5[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_nonFlow_v5[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_nonFlow_v5[6], digits=nd)	== round(bgb_lnL, digits=nd)

@test round(res_nonFlow_v6[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_nonFlow_v6[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_nonFlow_v6[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_nonFlow_v6[6], digits=nd)	== round(bgb_lnL, digits=nd)

@test round(res_nonFlow_v7[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_nonFlow_v7[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_nonFlow_v7[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_nonFlow_v7[6], digits=nd)	== round(bgb_lnL, digits=nd)

# Leave out parallel stuff
"""
@test round(res_nonFlow_v6par[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_nonFlow_v6par[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_nonFlow_v6par[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_nonFlow_v6par[6], digits=nd)	== round(bgb_lnL, digits=nd)

@test round(res_Gflow_v6[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v6[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v6[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v6[6], digits=nd)	== round(bgb_lnL, digits=nd)
"""


root_age = trdf[tr.root,:node_age]
num_incs = 100
Gseg_times = seq(0.0, root_age, root_age/num_incs);
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);


# Calculate array of Gflow matrices with float64 matrix multiplication
(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap = Gmaps.construct_Gmap_interpolator(pG, Gseg_times; abstol=solver_options.abstol, reltol=solver_options.reltol);

(Gseg_timesD, Gflows_arrayD, Gflows_array_totalsD, Gflows_dictD) = GmapD = construct_Gmap_interpolator_float64D(pG, Gseg_times; abstol=solver_options.abstol, reltol=solver_options.reltol);


# These should be DIFFERENT, if extinction is positive!
Gflows_dict[1](0.1)
minimum(Gflows_dict[10].t)
maximum(Gflows_dict[10].t)
Gflows_dict[10](maximum(Gflows_dict[10].t))

Gflows_dictD[1](0.1)
minimum(Gflows_dictD[10].t)
maximum(Gflows_dictD[10].t)
Gflows_dictD[10](maximum(Gflows_dictD[10].t))


# Calculate array of Gflow matrices with double64 matrix multiplication
(Gseg_timesDF, Gflows_arrayDF, Gflows_array_totalsDF, Gflows_dictDF) = Gmap_Double64 = Gmaps.construct_Gmap_interpolator_double64(pG, Gseg_times; abstol=solver_options.abstol, reltol=solver_options.reltol);

Gflows_dictDF[1](0.1)
minimum(Gflows_dictDF[10].t)
maximum(Gflows_dictDF[10].t)
Gflows_dictDF[10](maximum(Gflows_dictDF[10].t))


Gflows_array_totals[:,:,1][1,:]
Gflows_array_totalsD[:,:,1][1,:]
Gflows_array_totalsDF[:,:,1][1,:]

tval = 0.1
Gmaps.interp_from_Gmap(tval, Gmap)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, GmapD)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, Gmap_Double64)[1:5,1:5]
Gflow_to_01_GMRES(tval)[1:5,1:5]

tval = 0.2
Gmaps.interp_from_Gmap(tval, Gmap)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, GmapD)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, Gmap_Double64)[1:5,1:5]
Gflow_to_01_GMRES(tval)[1:5,1:5]

tval = 1.0
Gmaps.interp_from_Gmap(tval, Gmap)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, GmapD)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, Gmap_Double64)[1:5,1:5]
Gflow_to_01_GMRES(tval)[1:5,1:5]


Gflows_array_totals[1:5,1:5,10]
Gflows_array_totalsD[1:5,1:5,10]
Gflows_array_totalsDF[1:5,1:5,10]

tval = 5.0
Gmaps.interp_from_Gmap(tval, Gmap)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, GmapD)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, Gmap_Double64)[1:5,1:5]
Gflow_to_01_GMRES(tval)[1:5,1:5]


tval = 40.0
Gmaps.interp_from_Gmap(tval, Gmap)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, GmapD)[1:5,1:5]
Gmaps.interp_from_Gmap(tval, Gmap_Double64)[1:5,1:5]
Gflow_to_01_GMRES(tval)[1:5,1:5]


@test mean(abs.(Gmaps.interp_from_Gmap(0.1, Gmap) .- Gflow_to_01_GMRES(0.1))) < 0.00001
@test mean(abs.(Gmaps.interp_from_Gmap(0.1, Gmap_Double64) .- Gflow_to_01_GMRES(0.1))) < 0.00001


Gflows_array_totals[:,:,2]
Gmaps.interp_from_Gmap(0.2, Gmap)
Gflow_to_01_GMRES(0.2)
@test mean(abs.(Gmaps.interp_from_Gmap(0.2, Gmap) .- Gflow_to_01_GMRES(0.2))) < 0.0001
@test mean(abs.(Gmaps.interp_from_Gmap(0.2, Gmap_Double64) .- Gflow_to_01_GMRES(0.2))) < 0.0001


Gflows_array_totals[:,:,3]
Gmaps.interp_from_Gmap(0.3, Gmap)
Gflow_to_01_GMRES(0.3)
@test mean(abs.(Gmaps.interp_from_Gmap(0.3, Gmap) .- Gflow_to_01_GMRES(0.3))) < 0.0001
@test mean(abs.(Gmaps.interp_from_Gmap(0.3, Gmap_Double64) .- Gflow_to_01_GMRES(0.3))) < 0.0001

Gflows_array_totals[:,:,52]
Gmaps.interp_from_Gmap(root_age, Gmap)
Gflow_to_01_GMRES(root_age)
@test mean(abs.(Gmaps.interp_from_Gmap(root_age, Gmap) .- Gflow_to_01_GMRES(root_age))) < 0.0001
@test mean(abs.(Gmaps.interp_from_Gmap(root_age, Gmap_Double64) .- Gflow_to_01_GMRES(root_age))) < 0.0001


Gflow_via_Gmap = t -> Gmaps.interp_from_Gmap(t, Gmap)
Gflow_via_GmapD = t -> Gmaps.interp_from_Gmap(t, GmapD)

sum(Gflow_via_Gmap(root_age)[1,:] .== Gflow_to_01_GMRES(root_age)[1,:])

DataFrame(Gflow_via_Gmap(root_age) .- Gflow_to_01_GMRES(root_age), :auto)

@test all(abs.(Gflow_via_Gmap(root_age) .- Gflow_to_01_GMRES(root_age)) .< 0.0001)
@test all(abs.(Gflow_via_Gmap(root_age) .- Gflow_to_01_GMRES(root_age)) .< 0.001)



# The Gmap strategy works with Float64 or Double64
res_Gflow_v6a = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_via_Gmap, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6a, iteration_number_GFv6a, Julia_sum_lq_GFv6a, rootstates_lnL_GFv6a, Julia_total_lnLs1_GFv6a, bgb_lnl_GFv6a) = res_Gflow_v6a
archived_Gflow_v6a = deepcopy(res);

res_Gflow_v6b = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_via_GmapD, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6b, iteration_number_GFv6b, Julia_sum_lq_GFv6b, rootstates_lnL_GFv6b, Julia_total_lnLs1_GFv6b, bgb_lnl_GFv6b) = res_Gflow_v6b
archived_Gflow_v6b = deepcopy(res);

res_Gflow_v6
res_Gflow_v6a
res_Gflow_v6b

R_order = sort(trdf, :Rnodenums).nodeIndex

vfft(archived_Gflow_v6.normlikes_at_each_nodeIndex_branchBot[R_order]) .== vfft(archived_Gflow_v6a.normlikes_at_each_nodeIndex_branchBot[R_order])


rnode = 1
vfft(archived_Gflow_v6.normlikes_at_each_nodeIndex_branchBot[R_order])[rnode,:]
vfft(archived_Gflow_v6a.normlikes_at_each_nodeIndex_branchBot[R_order])[rnode,:]



# 2022-03-30: (0.589, 12, -285.1005967032069, -12.025331095939737, -297.12592779914667, -110.07183134705335)

# Identical results with Double64 (so probably unnecessary here)
Gflow_Double64 = t -> Gmaps.interp_from_Gmap(t, Gmap_Double64)
res_Gflow_v6_Double64 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_Double64, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6_Double64, iteration_number_GFv6_Double64, Julia_sum_lq_GFv6_Double64, rootstates_lnL_GFv6_Double64, Julia_total_lnLs1_GFv6_Double64, bgb_lnl_GFv6_Double64) = res_Gflow_v6_Double64
# 2022-03-30: (8.086, 12, -285.1005966444711, -12.025331096149882, -297.125927740621, -110.07183128847518)


print("\nTesting DEC+J Gflow SSE likelihood downpass v6 vs. Gflow_arrays v7, with half-matrix:\n")
@test abs(Julia_sum_lq_nFv6 - Julia_sum_lq_GFv6) < 0.1
@test abs(rootstates_lnL_nFv6 - rootstates_lnL_GFv6) < 0.1
@test abs(Julia_total_lnLs1_nFv6 - Julia_total_lnLs1_GFv6) < 0.1
@test abs(bgb_lnl_nFv6 - bgb_lnl_GFv6) < 0.1

print("\nTesting DEC+J traditional SSE likelihood downpass v6 vs. Gflow_arrays v7 using Float64, with half-matrix:\n")
@test abs(Julia_sum_lq_GFv6 - Julia_sum_lq_GFv6a) < 0.1
@test abs(rootstates_lnL_GFv6 - rootstates_lnL_GFv6a) < 0.1
@test abs(Julia_total_lnLs1_GFv6 - Julia_total_lnLs1_GFv6a) < 0.1
@test abs(bgb_lnl_GFv6 - bgb_lnl_GFv6a) < 0.1

print("\nTesting DEC+J traditional SSE likelihood downpass v6 vs. Gflow_arrays v7 using Double64, with half-matrix:\n")
@test abs(Julia_sum_lq_nFv6 - Julia_sum_lq_GFv6_Double64) < 0.1
@test abs(rootstates_lnL_nFv6 - rootstates_lnL_GFv6_Double64) < 0.1
@test abs(Julia_total_lnLs1_nFv6 - Julia_total_lnLs1_GFv6_Double64) < 0.1
@test abs(bgb_lnl_nFv6 - bgb_lnl_GFv6_Double64) < 0.1

total_calctime_in_sec_nFv6
total_calctime_in_sec_GFv6
total_calctime_in_sec_GFv6a
total_calctime_in_sec_GFv6_Double64



# Check the condition numbers of linear dynamics A and Gflow G (I think Julia does this automatically)
tvals = seq(0.0, 5.2, 0.1);
kappa_Arates_df = check_linearDynamics_of_As(tvals, p_Ds_v5; max_condition_number=1e8)
G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 2.5);
prob_Gs_v5_condnums = DifferentialEquations.ODEProblem(calc_Gs_SSE_condnums!, G0, tspan, pG)

end
