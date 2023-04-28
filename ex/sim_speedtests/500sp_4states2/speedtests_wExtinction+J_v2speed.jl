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

wd = "/Users/nmat471/HD/GitHub/PhyBEARS.jl/ex/sim_speedtests/500sp_4states2/"
cd(wd)
#trfn = "/GitHub/PhyBEARS.jl/data/Klaus_Matzke_2020_PodoArau_197sp.newick"
#tr = readTopology(trfn)
#lgdata_fn = "/GitHub/PhyBEARS.jl/data/Podocarpaceae_197_9areas_5Araucariaceae.data"
#geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
trfn = "living_tree_noNodeLabels.newick"
tr = readTopology(trfn)

lgdata_fn = "geog_living.data"
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
numareas = Rncol(geog_df) - 1
include_null_range=true
n = (2^numareas) - (include_null_range==false)          # 4 areas, 16 states

# Set up the model
#root_age_mult=1.5; max_range_size=NaN; include_null_range=true; max_range_size=NaN
max_range_size = NaN # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=include_null_range, bmo=bmo);
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

(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6





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


res_Gflow_v7 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
res_Gflow_v7 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv7, iteration_number_GFv7, Julia_sum_lq_GFv7, rootstates_lnL_GFv7, Julia_total_lnLs1_GFv7, bgb_lnl_GFv7) = res_Gflow_v7



# v7 simd


G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 1.01 * maximum(trdf.node_age))

prob_Gs_v7simd = DifferentialEquations.ODEProblem(calc_Gs_SSE_v7simd, G0, tspan, pG);
Gflow_to_01_GMRES_v7simd  = solve(prob_Gs_v7simd, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);

res_Gflow_v7simd = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7simd, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
res_Gflow_v7simd = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_v7simd, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv7simd, iteration_number_GFv7simd, Julia_sum_lq_GFv7simd, rootstates_lnL_GFv7simd, Julia_total_lnLs1_GFv7simd, bgb_lnl_GFv7simd) = res_Gflow_v7simd


res_nonFlow_v5
res_nonFlow_v6
res_nonFlow_v7 
res_Gflow_v6
res_Gflow_v7
res_Gflow_v7simd

total_calctime_in_sec_nFv5
total_calctime_in_sec_nFv6
total_calctime_in_sec_nFv7
total_calctime_in_sec_GFv6
total_calctime_in_sec_GFv7
total_calctime_in_sec_GFv7simd

tdf = Rrbind(res_nonFlow_v5,
res_nonFlow_v6[1:6],
res_nonFlow_v7 ,
res_Gflow_v6,
res_Gflow_v7,
res_Gflow_v7simd)

df = DataFrame(tdf)

writedlm("times_v1.txt", eachrow(df), "\t")
# end