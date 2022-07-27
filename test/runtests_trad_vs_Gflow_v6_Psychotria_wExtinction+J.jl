#######################################################
# Test the lnLs from trad SSE v6 (uses a half matrix)
# and Gflow v7 (uses a half matrix)
# 2022-03-24: WORKS!
#######################################################

using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyBEARS.Parsers
using PhyBEARS.Gmaps

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_trad_vs_Gflow_v5_Psychotria_wExtinction+J.jl")
"""

#@testset "runtests_trad_vs_Gflow_v5_Psychotria_wExtinction+J.jl" begin

#include("/GitHub/PhyBEARS.jl/notes/Gmaps.jl")
#import .Gmaps
#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes
#include("/GitHub/PhyBEARS.jl/notes/jl")
#import .Flow

tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")

lgdata_fn = "/GitHub/PhyBEARS.jl/data/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.32881638319078066
bmo.est[bmo.rownames .== "deathRate"] .= 0.1
bmo.est[bmo.rownames .== "d"] .= 1e-10
bmo.est[bmo.rownames .== "e"] .= 1e-10
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.11
bmo.est[:] = bmo_updater_v1(bmo);
numareas = 4
n = 15          # 4 areas, 16 states

# Set up the model #
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;
solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-6;
solver_options.reltol = 1e-6;

p_Ds_v5 = inputs.p_Ds_v5;
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

# Solve the Es
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);

res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6

# Version 7/2 ClaSSE Gflow calculations
G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 1.01 * maximum(trdf.node_age))
prob_Gs_v5 = DifferentialEquations.ODEProblem(Gmaps.calc_Gs_SSE, G0, tspan, pG);
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-6, reltol = 1e-6);
res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRES, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6

root_age = trdf[tr.root,:node_age]
Gseg_times = seq(0.0, root_age, 0.1);
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);


# Calculate array of Gflow matrices with float64 matrix multiplication
(Gseg_times, Gflows_array, Gflows_array_totals, Gflows_dict) = Gmap = Gmaps.construct_Gmap_interpolator(pG, Gseg_times; abstol=1e-6, reltol=1e-6);

# These should be DIFFERENT, if extinction is positive!
Gflows_dict[1](0.1)
Gflows_dict[20](1.91)

# Calculate array of Gflow matrices with double64 matrix multiplication
(Gseg_timesF, Gflows_arrayF, Gflows_array_totalsF, Gflows_dictF) = Gmap_Double64 = Gmaps.construct_Gmap_interpolator_double64(pG, Gseg_times; abstol=1e-6, reltol=1e-6);


Gflows_array_totals[:,:,1]
Gmaps.interp_from_Gmap(0.1, Gmap)
Gflow_to_01_GMRES(0.1)

@test mean(abs.(Gmaps.interp_from_Gmap(0.1, Gmap) .- Gflow_to_01_GMRES(0.1))) < 0.00001
@test mean(abs.(Gmaps.interp_from_Gmap(0.1, Gmap_Double64) .- Gflow_to_01_GMRES(0.1))) < 0.00001


Gflows_array_totals[:,:,2]
Gmaps.interp_from_Gmap(0.2, Gmap)
Gflow_to_01_GMRES(0.2)
@test mean(abs.(Gmaps.interp_from_Gmap(0.2, Gmap) .- Gflow_to_01_GMRES(0.2))) < 0.00001
@test mean(abs.(Gmaps.interp_from_Gmap(0.2, Gmap_Double64) .- Gflow_to_01_GMRES(0.2))) < 0.00001


Gflows_array_totals[:,:,3]
Gmaps.interp_from_Gmap(0.3, Gmap)
Gflow_to_01_GMRES(0.3)
@test mean(abs.(Gmaps.interp_from_Gmap(0.3, Gmap) .- Gflow_to_01_GMRES(0.3))) < 0.00001
@test mean(abs.(Gmaps.interp_from_Gmap(0.3, Gmap_Double64) .- Gflow_to_01_GMRES(0.3))) < 0.00001

Gflows_array_totals[:,:,50]
Gmaps.interp_from_Gmap(5.0, Gmap)
Gflow_to_01_GMRES(5.0)
@test mean(abs.(Gmaps.interp_from_Gmap(5.0, Gmap) .- Gflow_to_01_GMRES(5.0))) < 0.00001
@test mean(abs.(Gmaps.interp_from_Gmap(5.0, Gmap_Double64) .- Gflow_to_01_GMRES(5.0))) < 0.00001


Gflow_via_Gmap = t -> Gmaps.interp_from_Gmap(t, Gmap)

# The Gmap strategy works with Float64 or Double64
res_Gflow_v6a = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_via_Gmap, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6a, rootstates_lnL_GFv6a, Julia_total_lnLs1_GFv6a, bgb_lnl_GFv6a) = res_Gflow_v6a

# Identical results with Double64 (so probably unnecessary here)
Gflow_Double64 = t -> Gmaps.interp_from_Gmap(t, Gmap_Double64)
res_Gflow_v6_Double64 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_Double64, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6_Double64, rootstates_lnL_GFv6_Double64, Julia_total_lnLs1_GFv6_Double64, bgb_lnl_GFv6_Double64) = res_Gflow_v6_Double64



print("\nTesting DEC+J Gflow SSE likelihood downpass v6 vs. Gflow_arrays v7, with half-matrix:\n")
@test round(Julia_sum_lq_nFv6, digits=2) == round(Julia_sum_lq_GFv6, digits=2)
@test round(rootstates_lnL_nFv6, digits=2) == round(rootstates_lnL_GFv6, digits=2)
@test round(Julia_total_lnLs1_nFv6, digits=2) == round(Julia_total_lnLs1_GFv6, digits=2)
@test round(bgb_lnl_nFv6, digits=2) == round(bgb_lnl_GFv6, digits=2)

print("\nTesting DEC+J traditional SSE likelihood downpass v6 vs. Gflow_arrays v7 using Float64, with half-matrix:\n")
@test round(Julia_sum_lq_GFv6, digits=2) == round(Julia_sum_lq_GFv6a, digits=2)
@test round(rootstates_lnL_GFv6, digits=2) == round(rootstates_lnL_GFv6a, digits=2)
@test round(Julia_total_lnLs1_GFv6, digits=2) == round(Julia_total_lnLs1_GFv6a, digits=2)
@test round(bgb_lnl_GFv6, digits=2) == round(bgb_lnl_GFv6a, digits=2)

print("\nTesting DEC+J traditional SSE likelihood downpass v6 vs. Gflow_arrays v7 using Double64, with half-matrix:\n")
@test round(Julia_sum_lq_nFv6, digits=2) == round(Julia_sum_lq_GFv6_Double64, digits=2)
@test round(rootstates_lnL_nFv6, digits=2) == round(rootstates_lnL_GFv6_Double64, digits=2)
@test round(Julia_total_lnLs1_nFv6, digits=2) == round(Julia_total_lnLs1_GFv6_Double64, digits=2)
@test round(bgb_lnl_nFv6, digits=2) == round(bgb_lnl_GFv6_Double64, digits=2)



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
#Gflow_to_25_condnums  = solve(prob_Gs_v5_condnums, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol = 1e-6, reltol = 1e-6)

end
