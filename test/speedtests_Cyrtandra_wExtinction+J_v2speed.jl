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

cd(expanduser("~/GitHub/PhyBEARS.jl/test/"))
include(expanduser("~/GitHub/PhyBEARS.jl/test/speedtests_Cyrtandra_wExtinction+J_v2speed.jl"))
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
using PhyBEARS.Flow
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
@everywhere cd(expanduser("~/GitHub/PhyBEARS.jl"))
@everywhere Pkg.add(PackageSpec(path=expanduser("~/GitHub/PhyBEARS.jl"))) # Gives error:
# ERROR: TaskFailedException
#   nested task error: TOML Parser error:
#   none:2:5 error: expected equal sign after key
# ...but the next line still works.
@everywhere using PhyBEARS
"""


#@testset "speedtests_Cyrtandra_wExtinction+J_v2speed.jl" begin

#include(expanduser("~/GitHub/PhyBEARS.jl/src/Gmaps.jl"))
#import .Gmaps
#include(expanduser("~/GitHub/PhyBEARS.jl/notes/ModelLikes.jl"))
#import .ModelLikes
#include(expanduser("~/GitHub/PhyBEARS.jl/notes/jl"))
#import .Flow

#trfn = expanduser("~/GitHub/PhyBEARS.jl/data/Klaus_Matzke_2020_PodoArau_197sp.newick")
#tr = readTopology(trfn)
#lgdata_fn = expanduser("~/GitHub/PhyBEARS.jl/data/Podocarpaceae_197_9areas_5Araucariaceae.data")
#geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
trfn = expanduser("~/GitHub/PhyBEARS.jl/data/Cyrtandra.newick")
tr = readTopology(trfn)

lgdata_fn = expanduser("~/GitHub/PhyBEARS.jl/data/Cyrtandra_geog.data")
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
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;
vfft(res.likes_at_each_nodeIndex_branchTop)

#@everywhere trdf = trdf

solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-8;
solver_options.reltol = 1e-8;

p_Es_v5 = inputs.p_Ds_v5;
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


#@benchmark iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true)

# BenchmarkTools.Trial: 34 samples with 1 evaluation.
#  Range (min … max):  136.110 ms … 178.505 ms  ┊ GC (min … max): 0.00% … 14.49%
#  Time  (median):     140.932 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   149.529 ms ±  15.212 ms  ┊ GC (mean ± σ):  5.07% ±  6.93%
# 
#   ▁   ▁█  ▁                                                      
#   █▇▄▁██▇▇█▄▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▄▁▄▄▇▇▁▁▁▁▁▄▁▄ ▁
#   136 ms           Histogram: frequency by time          179 ms <
# 
#  Memory estimate: 58.13 MiB, allocs estimate: 330407.
# 



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

#@benchmark Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol)

# BenchmarkTools.Trial: 1 sample with 1 evaluation.
#  Single result which took 9.507 s (3.09% GC) to evaluate,
#  with a memory estimate of 2.62 GiB, over 13001642 allocations.

res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6
archived_Gflow_v6 = deepcopy(res);

print("\nTesting DEC+J traditional SSE likelihood downpass v6 vs. Gflow:\n")
@test abs(Julia_sum_lq_nFv6 - Julia_sum_lq_GFv6) < 0.1
@test abs(rootstates_lnL_nFv6 - rootstates_lnL_GFv6) < 0.1
@test abs(Julia_total_lnLs1_nFv6 - Julia_total_lnLs1_GFv6) < 0.1
@test abs(bgb_lnl_nFv6 - bgb_lnl_GFv6) < 0.1


# 2022-03-29: (7.323, 22, -1271.0625625625657, -14.712302169676475, -1285.7748647322421, -575.9248397325294)
# 2022-03-30: (total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6

#@benchmark res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true)

# BenchmarkTools.Trial: 42 samples with 1 evaluation.
#  Range (min … max):  104.521 ms … 209.983 ms  ┊ GC (min … max): 0.00% … 44.14%
#  Time  (median):     112.160 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   120.821 ms ±  24.679 ms  ┊ GC (mean ± σ):  5.26% ± 11.54%
# 
#     ▃▃█                                                          
#   ▃▃███▅▅▁▆▅▁▁▁▃▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▁▃▁▁▁▁▁▃ ▁
#   105 ms           Histogram: frequency by time          210 ms <
# 
#  Memory estimate: 63.48 MiB, allocs estimate: 168882.




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



# Repeat for Gflow v7
# Gflow_v7 ClaSSE Gflow calculations
G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 1.01 * maximum(trdf.node_age))

dG = deepcopy(G0);
G = similar(A);
t = 1.0;
dG = calc_Gs_SSE_v7simd(dG, G, pG, t);
minimum(dG)
maximum(dG)

dG = deepcopy(G0);
G = similar(A);
dG = calc_Gs_SSE(dG, G, pG, t);
minimum(dG)
maximum(dG)


A1 = parameterized_ClaSSE_As_v7!(A, t, pG.p_Ds_v5);
A2 = parameterized_ClaSSE_As_v7_simd!(A, t, pG.p_Ds_v5);
@test all(A1 .== A2)


prob_Gs_v7 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE_v7, G0, tspan, pG);
Gflow_to_01_GMRES_v7  = solve(prob_Gs_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);

#@benchmark Gflow_to_01_GMRES_v7  = solve(prob_Gs_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol)

# BenchmarkTools.Trial: 3 samples with 1 evaluation.
#  Range (min … max):  1.731 s …   1.760 s  ┊ GC (min … max): 0.00% … 3.46%
#  Time  (median):     1.749 s              ┊ GC (median):    1.54%
#  Time  (mean ± σ):   1.747 s ± 14.661 ms  ┊ GC (mean ± σ):  1.68% ± 1.73%
# 
#   █                                 █                     █  
#   █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
#   1.73 s         Histogram: frequency by time        1.76 s <
# 
#  Memory estimate: 576.62 MiB, allocs estimate: 68893.

res_Gflow_v7 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
res_Gflow_v7 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv7, iteration_number_GFv7, Julia_sum_lq_GFv7, rootstates_lnL_GFv7, Julia_total_lnLs1_GFv7, bgb_lnl_GFv7) = res_Gflow_v7
archived_Gflow_v7 = deepcopy(res);

#@benchmark res_Gflow_v7 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true)

# BenchmarkTools.Trial: 48 samples with 1 evaluation.
#  Range (min … max):   95.578 ms … 187.513 ms  ┊ GC (min … max): 0.00% … 43.14%
#  Time  (median):     101.731 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   105.946 ms ±  19.204 ms  ┊ GC (mean ± σ):  4.58% ± 10.64%
# 
#    ▁ ▃█▆                                                         
#   ▇█▇███▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▁▃▁▁▁▁▁▁▃ ▁
#   95.6 ms          Histogram: frequency by time          188 ms <
# 
#  Memory estimate: 63.48 MiB, allocs estimate: 168882.

nd = 3
@test round(res_Gflow_v7[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v7[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v7[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v7[6], digits=nd)	== round(bgb_lnL, digits=nd)

@test round(res_Gflow_v7[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v7[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v7[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v7[6], digits=nd)	== round(bgb_lnL, digits=nd)

@test round(res_Gflow_v7[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v7[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v7[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v7[6], digits=nd)	== round(bgb_lnL, digits=nd)





G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 1.01 * maximum(trdf.node_age))

prob_Gs_v7simd = DifferentialEquations.ODEProblem(calc_Gs_SSE_v7simd, G0, tspan, pG);
Gflow_to_01_GMRES_v7simd  = solve(prob_Gs_v7simd, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);

#@benchmark Gflow_to_01_GMRES_v7simd  = solve(prob_Gs_v7simd, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol)

# BenchmarkTools.Trial: 3 samples with 1 evaluation.
#  Range (min … max):  1.689 s …   1.797 s  ┊ GC (min … max): 0.00% … 3.25%
#  Time  (median):     1.721 s              ┊ GC (median):    1.61%
#  Time  (mean ± σ):   1.735 s ± 55.165 ms  ┊ GC (mean ± σ):  1.65% ± 1.63%
# 
#   █               █                                       █  
#   █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
#   1.69 s         Histogram: frequency by time         1.8 s <
# 
#  Memory estimate: 576.62 MiB, allocs estimate: 68893.


res_Gflow_v7simd = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7simd, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
res_Gflow_v7simd = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7simd, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv7, iteration_number_GFv7, Julia_sum_lq_GFv7, rootstates_lnL_GFv7, Julia_total_lnLs1_GFv7, bgb_lnl_GFv7) = res_Gflow_v7simd
archived_Gflow_v7 = deepcopy(res);

#@benchmark res_Gflow_v7simd = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7simd, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true)

# BenchmarkTools.Trial: 47 samples with 1 evaluation.
#  Range (min … max):   92.745 ms … 177.313 ms  ┊ GC (min … max): 0.00% …  0.00%
#  Time  (median):     101.278 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   106.959 ms ±  20.788 ms  ┊ GC (mean ± σ):  4.37% ± 10.58%
# 
#    ▆  ▁▄▃█                                                       
#   ▆█▆▄████▇▄▁▄▁▄▄▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▁▁▁▁▆▁▄ ▁
#   92.7 ms          Histogram: frequency by time          177 ms <
# 
#  Memory estimate: 63.48 MiB, allocs estimate: 168882.
# 
nd = 3
@test round(res_Gflow_v7simd[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v7simd[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v7simd[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v7simd[6], digits=nd)	== round(bgb_lnL, digits=nd)

@test round(res_Gflow_v7simd[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v7simd[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v7simd[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v7simd[6], digits=nd)	== round(bgb_lnL, digits=nd)

@test round(res_Gflow_v7simd[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v7simd[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v7simd[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v7simd[6], digits=nd)	== round(bgb_lnL, digits=nd)









G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 1.01 * maximum(trdf.node_age))

prob_Gs_v7simd_B = DifferentialEquations.ODEProblem(calc_Gs_SSE_v7simd_B!, G0, tspan, pG);
Gflow_to_01_GMRES_v7simd_B  = solve(prob_Gs_v7simd_B, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);

#@benchmark Gflow_to_01_GMRES_v7simd_B  = solve(prob_Gs_v7simd_B, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol)

# BenchmarkTools.Trial: 3 samples with 1 evaluation.
#  Range (min … max):  1.780 s …   1.829 s  ┊ GC (min … max): 1.94% … 0.00%
#  Time  (median):     1.793 s              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   1.801 s ± 25.113 ms  ┊ GC (mean ± σ):  0.64% ± 1.12%
# 
#   █              █                                        █  
#   █▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
#   1.78 s         Histogram: frequency by time        1.83 s <
# 
#  Memory estimate: 576.62 MiB, allocs estimate: 68893.


res_Gflow_v7simd_B = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7simd_B, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
res_Gflow_v7simd_B = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7simd_B, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv7_B, iteration_number_GFv7_B, Julia_sum_lq_GFv7_B, rootstates_lnL_GFv7_B, Julia_total_lnLs1_GFv7_B, bgb_lnl_GFv7_B) = res_Gflow_v7simd_B
archived_Gflow_v7_B = deepcopy(res_Gflow_v7simd_B);

#@benchmark res_Gflow_v7simd_B = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7simd_B, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true)

# BenchmarkTools.Trial: 47 samples with 1 evaluation.
#  Range (min … max):   96.718 ms … 215.058 ms  ┊ GC (min … max): 0.00% … 49.98%
#  Time  (median):     101.383 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   107.950 ms ±  24.712 ms  ┊ GC (mean ± σ):  5.91% ± 12.33%
# 
#    █▅                                                            
#   ███▇▇▃▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▁▁▃ ▁
#   96.7 ms          Histogram: frequency by time          215 ms <
# 
#  Memory estimate: 63.48 MiB, allocs estimate: 168882.

nd = 3
@test round(res_Gflow_v7simd_B[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v7simd_B[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v7simd_B[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v7simd_B[6], digits=nd)	== round(bgb_lnL, digits=nd)

@test round(res_Gflow_v7simd_B[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v7simd_B[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v7simd_B[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v7simd_B[6], digits=nd)	== round(bgb_lnL, digits=nd)

@test round(res_Gflow_v7simd_B[3], digits=nd)	== round(Julia_sum_lq, digits=nd)
@test round(res_Gflow_v7simd_B[4], digits=nd)	== round(rootstates_lnL, digits=nd)
@test round(res_Gflow_v7simd_B[5], digits=nd)	== round(Julia_total_lnLs1, digits=nd)
@test round(res_Gflow_v7simd_B[6], digits=nd)	== round(bgb_lnL, digits=nd)












#######################################################
# The Gflow-matrices approach doesn't work (wrong answers if there is multiplication)
# so is being cut
# (and it's slow)
# 2023-03-09
#######################################################


# end