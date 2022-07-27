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
#using PhyBEARS
#using PhyBEARS.BGExample			# default examples
#using PhyBEARS.TrUtils			# basic utility functions 
#using PhyBEARS.MaxentInterp	# preconstructed interpolator for weighting rangesize of smaller daughter
#using PhyBEARS.TreeTable			# for prt() tree tables (DFs), bd_liks(), etc.
#using PhyBEARS.StateSpace	# set up lists of areas and states (geographic ranges)
#using PhyBEARS.SSEs				# SSE calculations with various amounts of speed optimization
#using PhyBEARS.Parsers			# Parsers to read e.g. geography file
#using PhyBEARS.TreePass		# downpass and uppass through the phylogeny; prt() etc.
#using PhyBEARS.ModelLikes		# likelihood calculations
#using PhyBEARS.Flow		# downpass and uppass through the phylogeny
#using PhyBEARS.Gmaps		# Gmaps arrays etc.

using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using Statistics 			# for mean(), max()
using PhyBEARS.TreeTable # for prt()
using PhyBEARS.TrUtils # for flat2() (similar to unlist)
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.ModelLikes
using PhyBEARS.Flow
using PhyBEARS.Parsers
using PhyBEARS.Gmaps
using PhyBEARS.Optimizers

using Profile					# for @profile
using BenchmarkTools	# for @benchmark
#using PProf						# for pprof()
"""
# Run with:
cd /GitHub/PhyBEARS.jl
JULIA_NUM_THREADS=23 julia
Threads.nthreads()

# Launch multi-core Julia with e.g.: 
#julia --procs auto  # (uses as many cores as are available)
#julia -p auto       # (same)
#julia --procs 10    # (for 10 cores)
#
# Once you are inside Julia:
# length(Sys.cpu_info())
# nprocs
#julia --procs auto   # these require PhyBEARS installed (!)
Threads.nthreads()

cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/speedtests_PodoArau_wExtinction+J_downPass_v2.jl")
"""

#@testset "speedtests_Cyrtandra_wExtinction+J_v3.jl" begin

include("/GitHub/PhyBEARS.jl/notes/BranchSpeeds.jl")
import .BranchSpeeds
include("/GitHub/PhyBEARS.jl/notes/tmp_par_downpass_v1.jl")
#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes
#include("/GitHub/PhyBEARS.jl/notes/jl")
#import .Flow

# Podocarp + Araucariaceae data
#trfn = "/GitHub/PhyBEARS.jl/data/Klaus_Matzke_2020_PodoArau_197sp.newick"
#tr = readTopology(trfn)
#lgdata_fn = "/GitHub/PhyBEARS.jl/data/Podocarpaceae_197_9areas_5Araucariaceae.data"
#lgdata_fn = "/GitHub/PhyBEARS.jl/data/Podocarpaceae_197_9areas_2traits_5Araucariaceae.data"
#geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

trfn = "/GitHub/PhyBEARS.jl/data/Psychotria_tree.newick"
tr = readTopology(trfn)

lgdata_fn = "/GitHub/PhyBEARS.jl/data/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
include_null_range = false


# Cyrtandra
trfn = "/GitHub/PhyBEARS.jl/data/Cyrtandra.newick"
tr = readTopology(trfn)
lgdata_fn = "/GitHub/PhyBEARS.jl/data/Cyrtandra_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Geography data
include_null_range = false
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# Phylogeny
# Divergence times (setting tips near 0.0 mya to 0.0)
trtable = prt(tr)
trtable_node_ages = trtable.node_age
trtable_node_ages[trtable_node_ages .< 1.0e-6] .= 0.0
node_ages = sort!(unique(trtable_node_ages))[2:length(unique(trtable_node_ages))]


# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.1
bmo.est[bmo.rownames .== "deathRate"] .= 0.01
bmo.est[bmo.rownames .== "d"] .= 0.02
bmo.est[bmo.rownames .== "e"] .= 0.02
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.11
bmo.est[:] = bmo_updater_v1(bmo);

# Set up the model
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
prtCi(inputs)




(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
p_Es_v5 = deepcopy(inputs.p_Ds_v5);
p_Ds_v5_updater_v1!(p_Es_v5, inputs);

solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-6;
solver_options.reltol = 1e-3;

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol);

p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, sol_Es_v5=sol_Es_v7);

include("/GitHub/PhyBEARS.jl/notes/tmp_par_downpass_v1.jl")
res_tmppar_v6 = iterative_downpass_parallel_ClaSSE_v6tmp!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
print("\n\nwait asdf1\n\n")
res_tmppar_v7 = iterative_downpass_parallel_ClaSSE_v7tmp!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
print("\n\nwait asdf2\n\n")

@time res_tmppar_v6 = iterative_downpass_parallel_ClaSSE_v6tmp!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
@time res_tmppar_v7 = iterative_downpass_parallel_ClaSSE_v7tmp!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);

res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);

@time res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);


@time res_par_v7 = iterative_downpass_parallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);

@time res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);



(total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7
res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7




# Parallel version...not sure if it will work, @simd may not like parallelization
res_parallel_v7 = iterative_downpass_parallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_par7, iteration_number_par7, Julia_sum_lq_par7, rootstates_lnL_par7, Julia_total_lnLs1_par7, bgb_lnl_par7) = res_parallel_v7
res_parallel_v7 = iterative_downpass_parallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_par7, iteration_number_par7, Julia_sum_lq_par7, rootstates_lnL_par7, Julia_total_lnLs1_par7, bgb_lnl_par7) = res_parallel_v7
res_parallel_v7outputs = deepcopy(res)
Rnames(res_parallel_v7outputs)
unique(res_parallel_v7outputs.thread_for_each_nodeOp)
res_parallel_v7outputs.calc_duration

res_parallel_v6 = iterative_downpass_parallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_par6, iteration_number_par6, Julia_sum_lq_par6, rootstates_lnL_par6, Julia_total_lnLs1_par6, bgb_lnl_par6) = res_parallel_v6
res_parallel_v6 = iterative_downpass_parallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_par6, iteration_number_par6, Julia_sum_lq_par6, rootstates_lnL_par6, Julia_total_lnLs1_par6, bgb_lnl_par6) = res_parallel_v6



res_nonFlow_v8 = iterative_downpass_nonparallel_ClaSSE_v8!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv8, iteration_number_nFv8, Julia_sum_lq_nFv8, rootstates_lnL_nFv8, Julia_total_lnLs1_nFv8, bgb_lnl_nFv8) = res_nonFlow_v8
res_nonFlow_v8 = iterative_downpass_nonparallel_ClaSSE_v8!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv8, iteration_number_nFv8, Julia_sum_lq_nFv8, rootstates_lnL_nFv8, Julia_total_lnLs1_nFv8, bgb_lnl_nFv8) = res_nonFlow_v8


# 1 cores:
"""
julia> (total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7
(22.91, 22, -1834.6057811694277, -17.5498716247096, -1852.1556527941373, -1140.7386575347782)

julia> res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);

julia> (total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7
(20.962, 22, -1834.6057811694277, -17.5498716247096, -1852.1556527941373, -1140.7386575347782)

julia> res_nonFlow_v8 = iterative_downpass_nonparallel_ClaSSE_v8!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);

julia> (total_calctime_in_sec_nFv8, iteration_number_nFv8, Julia_sum_lq_nFv8, rootstates_lnL_nFv8, Julia_total_lnLs1_nFv8, bgb_lnl_nFv8) = res_nonFlow_v8
(22.371, 22, -1834.605781170132, -17.549871624637444, -1852.1556527947696, -1140.7386575354103)

julia> res_nonFlow_v8 = iterative_downpass_nonparallel_ClaSSE_v8!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);

julia> (total_calctime_in_sec_nFv8, iteration_number_nFv8, Julia_sum_lq_nFv8, rootstates_lnL_nFv8, Julia_total_lnLs1_nFv8, bgb_lnl_nFv8) = res_nonFlow_v8
(21.506, 22, -1834.605781170132, -17.549871624637444, -1852.1556527947696, -1140.7386575354103)
""" 

# 23 cores: slightly slower actually; try standard parallel
"""
julia> (total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7
(24.199, 22, -1834.6057811694277, -17.5498716247096, -1852.1556527941373, -1140.7386575347782)

julia> res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);

julia> (total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7
(20.78, 22, -1834.6057811694277, -17.5498716247096, -1852.1556527941373, -1140.7386575347782)

julia> res_nonFlow_v8 = iterative_downpass_nonparallel_ClaSSE_v8!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);

julia> (total_calctime_in_sec_nFv8, iteration_number_nFv8, Julia_sum_lq_nFv8, rootstates_lnL_nFv8, Julia_total_lnLs1_nFv8, bgb_lnl_nFv8) = res_nonFlow_v8
(38.728, 22, -1834.605781170132, -17.549871624637444, -1852.1556527947696, -1140.7386575354103)

julia> res_nonFlow_v8 = iterative_downpass_nonparallel_ClaSSE_v8!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);

julia> (total_calctime_in_sec_nFv8, iteration_number_nFv8, Julia_sum_lq_nFv8, rootstates_lnL_nFv8, Julia_total_lnLs1_nFv8, bgb_lnl_nFv8) = res_nonFlow_v8
(22.01, 22, -1834.605781170132, -17.549871624637444, -1852.1556527947696, -1140.7386575354103)
"""




end
