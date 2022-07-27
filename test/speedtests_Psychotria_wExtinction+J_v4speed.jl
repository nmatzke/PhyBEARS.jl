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
using PhyBEARS.Parsers
using PhyBEARS.ModelLikes
using PhyBEARS.Gmaps

"""
# Run with:
cd /GitHub/PhyBEARS.jl
JULIA_NUM_THREADS=23 julia
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/speedtests_PodoArau_wExtinction+J.jl")
"""

#@testset "speedtests_PodoArau_wExtinction+J.jl" begin

#include("/GitHub/PhyBEARS.jl/src/Gmaps.jl")
#import .Gmaps
#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes
#include("/GitHub/PhyBEARS.jl/notes/jl")
#import .Flow

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
bmo.est[bmo.rownames .== "birthRate"] .= 0.1
bmo.est[bmo.rownames .== "deathRate"] .= 0.01
bmo.est[bmo.rownames .== "d"] .= 0.02
bmo.est[bmo.rownames .== "e"] .= 0.02
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.11
bmo.est[:] = bmo_updater_v1(bmo);

# Set up the model
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-12;
solver_options.reltol = 1e-12;

p_Ds_v5 = inputs.p_Ds_v5;
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

# Solve the Es
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
maximum(sol_Es_v5.t)

p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);

#Rnames(p_Ds_v5.p_TFs)
#p_Ds_v5.p_TFs.Qij_singleNum_sub_i

res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6, iteration_number_nFv6, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6
# (5.417, 8, -66.58313029898488, -4.938958329245979, -71.52208862823086, -33.4908373922376)

res_nonFlow_v6par = iterative_downpass_parallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6par, iteration_number_nFv6par, Julia_sum_lq_nFv6par, rootstates_lnL_nFv6par, Julia_total_lnLs1_nFv6par, bgb_lnl_nFv6par) = res_nonFlow_v6par
# (0.594, 15442, -66.58313029898488, -4.938958329245979, -71.52208862823086, -33.4908373922376)

solver_options.abstol = 1e-16;
solver_options.reltol = 1e-16;



# Version 7/2 ClaSSE Gflow calculations
G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 1.01 * maximum(trdf.node_age))
prob_Gs_v5 = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE, G0, tspan, pG);
prob_Gs_v5_sub_i = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE_sub_i, G0, tspan, pG);
prob_Gs_v5_parallel = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE_parallel, G0, tspan, pG);

starttime = Dates.now();
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
endtime = Dates.now();
endtime - starttime # 0.196 sec


Gflow_to_01_GMRES0  = solve(prob_Gs_v5_sub_i, save_everystep=true);
starttime = Dates.now();
Gflow_to_01_GMRES0  = solve(prob_Gs_v5_sub_i, save_everystep=true);
endtime = Dates.now();
endtime - starttime  # 0.008 sec

Gflow_to_01_GMRESa  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-3, reltol=1e-3);
starttime = Dates.now();
Gflow_to_01_GMRESa  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-3, reltol=1e-3);
endtime = Dates.now();
endtime - starttime # 0.006 sec 

Gflow_to_01_GMRESb  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-6, reltol=1e-6);
starttime = Dates.now();
Gflow_to_01_GMRESb  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-6, reltol=1e-6);
endtime = Dates.now();
endtime - starttime # 0.008 sec

Gflow_to_01_GMRESc  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-9, reltol=1e-9);
starttime = Dates.now();
Gflow_to_01_GMRESc  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-9, reltol=1e-9);
endtime = Dates.now();
endtime - starttime # 0.012 sec

Gflow_to_01_GMRESd  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-12, reltol=1e-12);
starttime = Dates.now();
Gflow_to_01_GMRESd  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-12, reltol=1e-12);
endtime = Dates.now();
endtime - starttime # 0.16 sec

Gflow_to_01_GMRESe  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-17, reltol=1e-17);
starttime = Dates.now();
Gflow_to_01_GMRESe  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-17, reltol=1e-17);
endtime = Dates.now();
endtime - starttime # 0.332 sec

Gflow_to_01_GMRESf  = solve(prob_Gs_v5_sub_i, Tsit5(), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
starttime = Dates.now();
Gflow_to_01_GMRESf  = solve(prob_Gs_v5_sub_i, Tsit5(), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
endtime = Dates.now();
endtime - starttime # 0.279 sec

starttime = Dates.now();
Gflow_to_01_GMRESg  = solve(prob_Gs_v5_sub_i, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
endtime = Dates.now();
endtime - starttime # 1.6 sec


#@benchmark Gflow_to_01_GMRES  = solve(prob_Gs_v5_sub_i, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol)


# Very slow on large trees!
starttime = Dates.now();
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
endtime = Dates.now();
endtime - starttime
# 2022-03-30: Gmaps.calc_Gs_SSE: 0.186 sec

# "Best case" -- ie matches tradSSE
starttime = Dates.now();
Gflow_to_01_GMRES_res16 = solve(prob_Gs_v5_sub_i, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=1e-16, reltol=1e-16);
endtime = Dates.now();
endtime - starttime
# 2022-03-31: Gmaps.calc_Gs_SSE_sub_i: 0.137 sec


starttime = Dates.now();
Gflow_to_01_GMRES  = solve(prob_Gs_v5_parallel, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
endtime = Dates.now();
endtime - starttime
# 2022-03-30: Gmaps.calc_Gs_SSE_parallel: 2.158 sec


Gflow_to_01_GMRES.alg
Gflow_to_01_GMRES0.alg
Gflow_to_01_GMRESa.alg
Gflow_to_01_GMRESb.alg
Gflow_to_01_GMRESc.alg
Gflow_to_01_GMRESd.alg
Gflow_to_01_GMRESe.alg
Gflow_to_01_GMRESf.alg
Gflow_to_01_GMRESg.alg
Gflow_to_01_GMRES_res16.alg



# This is fast, once Gflow is already computed
res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6
# 2022-03-31_e16: (1.344, 8, -66.658592796525, -4.926860033675147, -71.58545283020014, -33.55928414900776)

res_Gflow_v6_res16 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES_res16, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6_res16
# (0.088, 8, -66.658592796525, -4.926860033675147, -71.58545283020014, -33.55928414900776)

res_Gflow_v60 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES0, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v60
# (0.66, 8, -66.65859302670283, -4.92686010272035, -71.58545312942319, -33.559284453004544)

res_Gflow_v6a = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRESa, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6a
# (0.006, 8, -66.65859905550221, -4.926860435086798, -71.585459490589, -33.559290854307996)

res_Gflow_v6b = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRESb, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6b
# (0.006, 8, -66.65859302670283, -4.92686010272035, -71.58545312942319, -33.559284453004544)


res_Gflow_v6c = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRESc, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6c
# (0.006, 8, -66.65859285466826, -4.926860027493808, -71.58545288216206, -33.559284210106675)


res_Gflow_v6d = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRESd, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6d
# (0.006, 8, -66.65859285064882, -4.926860025192063, -71.58545287584087, -33.559284204284715)

res_Gflow_v6e = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRESe, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6e
# (0.006, 8, -66.6585928501906, -4.926860021367958, -71.58545287155856, -33.55928420111669)

res_Gflow_v6f = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRESf, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6f
# (0.162, 8, -66.65859279652445, -4.926860033675143, -71.58545283019959, -33.559284149007205)


res_Gflow_v6g = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRESg, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6g
# (0.006, 8, -66.658592796525, -4.926860033675147, -71.58545283020014, -33.55928414900776)




root_age = trdf[tr.root,:node_age]
num_incs = 10
Gseg_times = seq(0.0, root_age, root_age/num_incs);
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);


# Calculate array of Gflow matrices with float64 matrix multiplication
starttime = Dates.now();
(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap = Gmaps.construct_Gmap_interpolator_float64(pG, Gseg_times; abstol=solver_options.abstol, reltol=solver_options.reltol);
endtime = Dates.now();
endtime - starttime
# 0.791 sec

starttime = Dates.now();
(Gseg_timesDF, Gflows_arrayDF, Gflows_array_totalsDF, Gflows_dictDF) = Gmap_Double64 = Gmaps.construct_Gmap_interpolator_double64(pG, Gseg_times; abstol=solver_options.abstol, reltol=solver_options.reltol);
endtime = Dates.now();
endtime - starttime
# 0.474 sec

# These should be DIFFERENT, if extinction is positive!
Gflows_dict[1](0.1)
Gflows_dict[10](5.2)

# Calculate array of Gflow matrices with double64 matrix multiplication


Gflows_array_totals[:,:,1]
Gmaps.interp_from_Gmap(0.1, Gmap)
Gflow_to_01_GMRES(0.1)

@test mean(abs.(Gmaps.interp_from_Gmap(0.1, Gmap) .- Gflow_to_01_GMRES(0.1))) < 0.0001
@test mean(abs.(Gmaps.interp_from_Gmap(0.1, Gmap_Double64) .- Gflow_to_01_GMRES(0.1))) < 0.0001


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

Gflows_array_totals[:,:,10]
Gmaps.interp_from_Gmap(root_age, Gmap)
Gflow_to_01_GMRES(root_age)
@test mean(abs.(Gmaps.interp_from_Gmap(root_age, Gmap) .- Gflow_to_01_GMRES(root_age))) < 0.0001
@test mean(abs.(Gmaps.interp_from_Gmap(root_age, Gmap_Double64) .- Gflow_to_01_GMRES(root_age))) < 0.0001


Gflow_via_Gmap = t -> Gmaps.interp_from_Gmap(t, Gmap)

# The Gmap strategy works with Float64 or Double64
res_Gflow_v6a = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_via_Gmap, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6a, iteration_number_GFv6a, Julia_sum_lq_GFv6a, rootstates_lnL_GFv6a, Julia_total_lnLs1_GFv6a, bgb_lnl_GFv6a) = res_Gflow_v6a

# Identical results with Double64 (so probably unnecessary here)
Gflow_Double64 = t -> Gmaps.interp_from_Gmap(t, Gmap_Double64)
res_Gflow_v6_Double64 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_Double64, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6_Double64, iteration_number_GFv6_Double64, Julia_sum_lq_GFv6_Double64, rootstates_lnL_GFv6_Double64, Julia_total_lnLs1_GFv6_Double64, bgb_lnl_GFv6_Double64) = res_Gflow_v6_Double64



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
