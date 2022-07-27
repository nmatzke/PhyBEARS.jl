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
using PProf						# for pprof()
"""
# Run with:
cd /GitHub/PhyBEARS.jl
JULIA_NUM_THREADS=23 julia
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/speedtests_PodoArau_wExtinction+J_branchTests_v1.jl")
"""

#@testset "speedtests_Cyrtandra_wExtinction+J_v3.jl" begin

include("/GitHub/PhyBEARS.jl/notes/BranchSpeeds.jl")
import .BranchSpeeds
#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes
#include("/GitHub/PhyBEARS.jl/notes/jl")
#import .Flow

# Podocarp + Araucariaceae data
trfn = "/GitHub/PhyBEARS.jl/data/Klaus_Matzke_2020_PodoArau_197sp.newick"
tr = readTopology(trfn)
lgdata_fn = "/GitHub/PhyBEARS.jl/data/Podocarpaceae_197_9areas_5Araucariaceae.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Cyrtandra
#trfn = "/GitHub/PhyBEARS.jl/data/Cyrtandra.newick"
#tr = readTopology(trfn)
#lgdata_fn = "/GitHub/PhyBEARS.jl/data/Cyrtandra_geog.data"
#geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

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
prob_Es_v5 = DifferentialEquations.ODEProblem(BranchSpeeds.parameterized_ClaSSE_Es_v5orig, p_Es_v5.uE, Es_tspan, p_Es_v5);
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol);

include("/GitHub/PhyBEARS.jl/notes/BranchSpeeds.jl")
import .BranchSpeeds


prob_Es_v7 = DifferentialEquations.ODEProblem(BranchSpeeds.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol);

@time sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol);

@time sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol);


@benchmark sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol)

@benchmark sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol)


#@profile sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, 
#abstol=solver_options.abstol, reltol=solver_options.reltol)
#pprof(;webport=58699)

#@profile sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, 
#abstol=solver_options.abstol, reltol=solver_options.reltol)
#pprof(;webport=58700)



sol_Es_v5.u[length(sol_Es_v5.u)]
sol_Es_v7.u[length(sol_Es_v7.u)]

sol_Es_v5.u[length(sol_Es_v5.u)] - sol_Es_v7.u[length(sol_Es_v7.u)]


p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);


p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v7);

#Rnames(p_Ds_v5.p_TFs)
#p_Ds_v5.p_TFs.Qij_singleNum_sub_i

Ds_tspan = maximum(trtable.node_age)
solver_options.save_everystep = false
solver_options.saveat = node_ages
tip_Ds = inputs.res.likes_at_each_nodeIndex_branchTop[1]
prob_Ds_v5 = DifferentialEquations.ODEProblem(BranchSpeeds.parameterized_ClaSSE_Ds_v5orig, tip_Ds, Ds_tspan, p_Ds_v5);
sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);

include("/GitHub/PhyBEARS.jl/notes/BranchSpeeds.jl")
import .BranchSpeeds

du = similar(tip_Ds)



prob_Ds_v7 = DifferentialEquations.ODEProblem(BranchSpeeds.parameterized_ClaSSE_Ds_v7orig, tip_Ds, Ds_tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);
@time sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);
#@code_warntype BranchSpeeds.parameterized_ClaSSE_Ds_v7orig(du, tip_Ds, p_Ds_v7, Ds_tspan)

prob_Ds_v7 = DifferentialEquations.ODEProblem(BranchSpeeds.parameterized_ClaSSE_Ds_v7_simd_sums, tip_Ds, Ds_tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);
@time sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);
#@code_warntype BranchSpeeds.parameterized_ClaSSE_Ds_v7_simd_sums(du, tip_Ds, p_Ds_v7, Ds_tspan)



du = similar(tip_Ds)
@time BranchSpeeds.parameterized_ClaSSE_Ds_v7orig(du, tip_Ds, p_Ds_v7, Ds_tspan)
#@code_warntype BranchSpeeds.parameterized_ClaSSE_Ds_v7orig(du, tip_Ds, p_Ds_v7, Ds_tspan)
@benchmark BranchSpeeds.parameterized_ClaSSE_Ds_v7orig(du, tip_Ds, p_Ds_v7, Ds_tspan)



@time sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);

@time sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);


@benchmark sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

#@profile sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=false, 
#abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);
#pprof(;webport=58699)

#@profile sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=false, 
#abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);
#pprof(;webport=58700)

# Traceur is essentially a codified version of the Julia performance tips
# https://github.com/JunoLab/Traceur.jl
using Traceur
#@trace sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=false, 
#abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);



#######################################################
# Speed tests with version 7 simd inbounds shared loops
#######################################################

# 40.883ms
@benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(), save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, Tsit5(), save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

#@benchmark sol_Ds_v7 = solve(prob_Ds_v7, QNDF(), save_everystep=false, 
#abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

#@benchmark sol_Ds_v7 = solve(prob_Ds_v7, FBDF(), save_everystep=false, 
#abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, Vern9(), save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, lsoda(), save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, VCABM(), save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:Band, jac_upper=1, jac_lower=1), save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:BCG), save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:TFQMR), save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)


# Results: CVODE_BDF(linear_solver=:GMRES) wins by a lot

"""
julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
BenchmarkTools.Trial: 164 samples with 1 evaluation.
 Range (min … max):  24.284 ms … 61.546 ms  ┊ GC (min … max): 0.00% … 22.44%
 Time  (median):     27.532 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   30.614 ms ±  7.042 ms  ┊ GC (mean ± σ):  4.61% ±  9.83%

   ▅█▂                                                         
  ▄███▅▅▄▄▄▄▅▄▆▄▃▃▄▄▃▃▃▃▁▄▁▄▃▁▁▄▃▃▅▄▃▁▁▃▁▃▃▃▁▃▁▃▁▁▁▁▁▁▁▁▁▁▃▁▃ ▃
  24.3 ms         Histogram: frequency by time        53.1 ms <

 Memory estimate: 14.75 MiB, allocs estimate: 640575.

julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
BenchmarkTools.Trial: 20 samples with 1 evaluation.
 Range (min … max):  216.064 ms … 308.667 ms  ┊ GC (min … max): 0.00% … 7.57%
 Time  (median):     255.206 ms               ┊ GC (median):    2.84%
 Time  (mean ± σ):   252.656 ms ±  22.349 ms  ┊ GC (mean ± σ):  3.37% ± 3.44%

  ▃                            ▃       █                         
  █▁▁▁▁▁▁▁▇▇▇▁▁▁▇▁▁▇▇▁▇▁▁▁▁▇▇▁▇█▁▇▁▁▁▁▇█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇ ▁
  216 ms           Histogram: frequency by time          309 ms <

 Memory estimate: 74.30 MiB, allocs estimate: 3379515.

julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, Tsit5(), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
BenchmarkTools.Trial: 76 samples with 1 evaluation.
 Range (min … max):  44.407 ms … 241.721 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     58.974 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   66.026 ms ±  26.329 ms  ┊ GC (mean ± σ):  5.35% ± 9.69%

  ▁ ▄ ▁ █   ▃                                                   
  █▆█▄█▇█▇▆▆█▄▆▇▇▄▇▄▁▇▇▄▆▁▁▁▆▄▄▄▄▁▄▁▆▁▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▁▁▁▄ ▁
  44.4 ms         Histogram: frequency by time          128 ms <

 Memory estimate: 25.60 MiB, allocs estimate: 1134391.

julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, QNDF(), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
# ERROR: MethodError: no method matching Float64(::ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEq.OrdinaryDiffEqTag, Float64}, Float64, 12})

julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, FBDF(), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
# ERROR: MethodError: no method matching Float64(::ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEq.OrdinaryDiffEqTag, Float64}, Float64, 12})


julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, Vern9(), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
BenchmarkTools.Trial: 34 samples with 1 evaluation.
 Range (min … max):  125.347 ms … 213.894 ms  ┊ GC (min … max): 0.00% … 10.31%
 Time  (median):     146.169 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   149.924 ms ±  19.716 ms  ┊ GC (mean ± σ):  6.28% ±  6.57%

      ▄▁▄ ▁▄         █        ▁                                  
  ▆▁▆▁███▁██▁▁▆▆▆▁▁▆▁█▆▆▆▆▆▁▁▆█▁▁▁▁▆▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▆▁▁▁▁▁▆ ▁
  125 ms           Histogram: frequency by time          214 ms <

 Memory estimate: 72.86 MiB, allocs estimate: 3232848.

julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, lsoda(), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
# ERROR: Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.


julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, VCABM(), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
BenchmarkTools.Trial: 95 samples with 1 evaluation.
 Range (min … max):  41.549 ms … 84.828 ms  ┊ GC (min … max): 0.00% … 32.93%
 Time  (median):     50.506 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   52.746 ms ±  9.023 ms  ┊ GC (mean ± σ):  5.98% ± 10.49%

   ▃▆ ▁ ▆ ▆█▃ ▃  ▃▁▁█           ▁                              
  ▄██▇█▄█▇███▇█▇▇████▇▇▄▁▄▄▁▇▁▁▄█▄▇▁▇▇▁▄▁▄▁▁▄▄▁▄▄▁▁▄▁▄▁▁▁▁▁▁▇ ▁
  41.5 ms         Histogram: frequency by time        77.6 ms <

 Memory estimate: 25.32 MiB, allocs estimate: 1115827.

julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:Band, jac_upper=1, jac_lower=1), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
BenchmarkTools.Trial: 14 samples with 1 evaluation.
 Range (min … max):  265.727 ms … 576.697 ms  ┊ GC (min … max): 0.00% … 7.23%
 Time  (median):     318.471 ms               ┊ GC (median):    5.38%
 Time  (mean ± σ):   374.632 ms ± 113.716 ms  ┊ GC (mean ± σ):  6.27% ± 3.49%

  █  █ █  ▁   ▁   ▁                  ▁      ▁       ▁  ▁      ▁  
  █▁▁█▁█▁▁█▁▁▁█▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁▁█▁▁▁▁▁▁█ ▁
  266 ms           Histogram: frequency by time          577 ms <

 Memory estimate: 129.06 MiB, allocs estimate: 5898736.

julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:BCG), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
BenchmarkTools.Trial: 151 samples with 1 evaluation.
 Range (min … max):  26.004 ms … 67.858 ms  ┊ GC (min … max): 0.00% … 47.35%
 Time  (median):     30.229 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   33.300 ms ±  8.605 ms  ┊ GC (mean ± σ):  6.65% ± 11.96%

   █  ▄▁▂                                                      
  ███████▆█▆▄▃█▅▃▃▃▁▅▁▃▁▅▁▃▁▃▁▁▃▃▃▁▁▁▃▁▁▃▃▁▁▃▃▃▁▁▁▁▁▃▃▁▁▁▁▁▁▃ ▃
  26 ms           Histogram: frequency by time        64.1 ms <

 Memory estimate: 16.02 MiB, allocs estimate: 698981.

julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:TFQMR), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
BenchmarkTools.Trial: 111 samples with 1 evaluation.
 Range (min … max):  36.337 ms … 77.494 ms  ┊ GC (min … max): 0.00% … 26.53%
 Time  (median):     42.163 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   45.394 ms ±  8.909 ms  ┊ GC (mean ± σ):  5.62% ± 10.19%

  ▃█   ▃ ▂▂                                                    
  █████████▇▇▆▅▆▃▅▇▅▃▃▁▃▃▅▃▅▅▁▃▅▅▃▆▁▁▁▃▁▃▁▃▃▃▁▃▁▁▁▁▁▃▁▁▁▁▁▁▁▃ ▃
  36.3 ms         Histogram: frequency by time        75.5 ms <

 Memory estimate: 22.04 MiB, allocs estimate: 975636.

julia> @benchmark sol_Ds_v7 = solve(prob_Ds_v7, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, 
       abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)
BenchmarkTools.Trial: 153 samples with 1 evaluation.
 Range (min … max):  23.905 ms … 74.067 ms  ┊ GC (min … max): 0.00% … 20.92%
 Time  (median):     30.137 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   32.696 ms ±  9.165 ms  ┊ GC (mean ± σ):  6.58% ± 11.73%

   ▁▁ ▁▁▄▁▅█▁▁▁                                                
  ▄██▆█████████▆▄▅▃▃▁▃▄▁▁▃▄▄▁▁▃▁▁▁▁▃▁▁▃▁▁▁▁▃▁▃▃▃▁▁▄▁▁▃▁▃▁▁▁▃▃ ▃
  23.9 ms         Histogram: frequency by time          64 ms <

 Memory estimate: 14.75 MiB, allocs estimate: 640575.
"""








# Automatically Accelerating ODEProblem Code
# https://github.com/SciML/ModelingToolkit.jl/blob/master/docs/src/mtkitize_tutorials/modelingtoolkitize.md
using ModelingToolkit
sys = modelingtoolkitize(prob_Ds_v7; u0=tip_Ds, tspan=Ds_tspan, p=p_Ds_v5)

# If we want to get a symbolic representation, we can simply call modelingtoolkitize on the prob, which will return an ODESystem:
prob_jac = ODEProblem(prob_Ds_v7, [], Ds_tspan, p_Ds_v5, jac=true)
prob_jac = ODEProblem(prob_Ds_v7, tip_Ds, Ds_tspan, p_Ds_v5, jac=true)

# Error:
# ERROR: No methods were found for the model function passed to the equation solver.
# The function `f` needs to have dispatches, for example, for an ODEProblem
# `f` must define either `f(u,p,t)` or `f(du,u,p,t)`. For more information
# on how the model function `f` should be defined, consult the docstring for 
# the appropriate `AbstractSciMLFunction`.



# Version 6
# Solve the Es
p_Es_v6 = deepcopy(inputs.p_Ds_v5);

prob_Es_v6 = DifferentialEquations.ODEProblem(BranchSpeeds.parameterized_ClaSSE_Es_v6orig, p_Es_v6.uE, Es_tspan, p_Es_v6);
sol_Es_v6 = solve(prob_Es_v6, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol);

p_Ds_v6 = (n=p_Es_v6.n, params=p_Es_v6.params, p_indices=p_Es_v6.p_indices, p_TFs=p_Es_v6.p_TFs, uE=p_Es_v6.uE, sol_Es_v5=sol_Es_v5);

Ds_tspan = maximum(trtable.node_age)
solver_options.save_everystep = false
solver_options.saveat = node_ages
tip_Ds = inputs.res.likes_at_each_nodeIndex_branchTop[1]
prob_Ds_v6 = DifferentialEquations.ODEProblem(BranchSpeeds.parameterized_ClaSSE_Ds_v6orig, tip_Ds, Ds_tspan, p_Ds_v6);
sol_Ds_v6 = solve(prob_Ds_v6, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);

round.(sol_Ds_v6 .- sol_Ds_v5, digits=3)




# Version 7
# Solve the Es
p_Es_v7 = deepcopy(inputs.p_Ds_v5);

prob_Es_v7 = DifferentialEquations.ODEProblem(BranchSpeeds.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, inputs.Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol);

p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, sol_Es_v5=sol_Es_v5);

Ds_tspan = maximum(trtable.node_age)
solver_options.save_everystep = false
solver_options.saveat = node_ages
tip_Ds = inputs.res.likes_at_each_nodeIndex_branchTop[1]
prob_Ds_v7 = DifferentialEquations.ODEProblem(BranchSpeeds.parameterized_ClaSSE_Ds_v7_simd_sums, tip_Ds, Ds_tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat);

round.(sol_Ds_v7 .- sol_Ds_v5, digits=3)





@benchmark sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol)

@benchmark sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)



@benchmark sol_Es_v6 = solve(prob_Es_v6, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol)

@benchmark sol_Ds_v6 = solve(prob_Ds_v6, solver_options.solver, save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)


@benchmark sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=true, 
abstol=solver_options.abstol, reltol=solver_options.reltol)

@benchmark sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=false, 
abstol=solver_options.abstol, reltol=solver_options.reltol, saveat=solver_options.saveat)


res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v7, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7
# (3.293, 22, -1253.7954418304166, -14.74717146725072, -1268.5426132976672, -558.5136417895891)


res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6, iteration_number_nFv6, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6
# 

res_nonFlow_v6par = iterative_downpass_parallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6par, iteration_number_nFv6par, Julia_sum_lq_nFv6par, rootstates_lnL_nFv6par, Julia_total_lnLs1_nFv6par, bgb_lnl_nFv6par) = res_nonFlow_v6par
# (0.594, 15442, -66.58313029898488, -4.938958329245979, -71.52208862823086, -33.4908373922376)

solver_options.abstol = 1e-9;
solver_options.reltol = 1e-6;
res_nonFlow_v6par_loRes = iterative_downpass_parallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6par, iteration_number_nFv6par, Julia_sum_lq_nFv6par, rootstates_lnL_nFv6par, Julia_total_lnLs1_nFv6par, bgb_lnl_nFv6par) = res_nonFlow_v6par_loRes


solver_options.abstol = 1e-12;
solver_options.reltol = 1e-9;
res_nonFlow_v6par_hiRes = iterative_downpass_parallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv6par, iteration_number_nFv6par, Julia_sum_lq_nFv6par, rootstates_lnL_nFv6par, Julia_total_lnLs1_nFv6par, bgb_lnl_nFv6par) = res_nonFlow_v6par_hiRes


res_nonFlow_v7 = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7par
# (0.594, 15442, -66.58313029898488, -4.938958329245979, -71.52208862823086, -33.4908373922376)

solver_options.abstol = 1e-9;
solver_options.reltol = 1e-6;
res_nonFlow_v7_loRes = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7_loRes


solver_options.abstol = 1e-12;
solver_options.reltol = 1e-9;
res_nonFlow_v7_hiRes = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7, iteration_number_nFv7, Julia_sum_lq_nFv7, rootstates_lnL_nFv7, Julia_total_lnLs1_nFv7, bgb_lnl_nFv7) = res_nonFlow_v7_hiRes



solver_options.abstol = 1e-9;
solver_options.reltol = 1e-6;
res_nonFlow_v7par = iterative_downpass_parallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7par, iteration_number_nFv7par, Julia_sum_lq_nFv7par, rootstates_lnL_nFv7par, Julia_total_lnLs1_nFv7par, bgb_lnl_nFv7par) = res_nonFlow_v7par

solver_options.abstol = 1e-6;
solver_options.reltol = 1e-3;
res_nonFlow_v7par_loRes = iterative_downpass_parallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7par_loRes, iteration_number_nFv7par_loRes, Julia_sum_lq_nFv7par_loRes, rootstates_lnL_nFv7par_loRes, Julia_total_lnLs1_nFv7par_loRes, bgb_lnl_nFv7par_loRes) = res_nonFlow_v7par_loRes



solver_options.abstol = 1e-12;
solver_options.reltol = 1e-9;
res_nonFlow_v7par_hiRes = iterative_downpass_parallel_ClaSSE_v7!(res; trdf, p_Ds_v7=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_nFv7par_hiRes, iteration_number_nFv7par_hiRes, Julia_sum_lq_nFv7par_hiRes, rootstates_lnL_nFv7par_hiRes, Julia_total_lnLs1_nFv7par_hiRes, bgb_lnl_nFv7par_hiRes) = res_nonFlow_v7par_hiRes



res_nonFlow_v6
res_nonFlow_v6par
res_nonFlow_v6par_loRes
res_nonFlow_v6par_hiRes
res_nonFlow_v6par_QNDF
res_nonFlow_v6par_FBDF

res_nonFlow_v7
res_nonFlow_v7_loRes
res_nonFlow_v7_hiRes




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

#Gflow_to_01_GMRESe  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-17, reltol=1e-17);
#starttime = Dates.now();
#Gflow_to_01_GMRESe  = solve(prob_Gs_v5_sub_i, save_everystep=true, abstol=1e-17, reltol=1e-17);
#endtime = Dates.now();
#endtime - starttime # 0.332 sec

Gflow_to_01_GMRESf  = solve(prob_Gs_v5_sub_i, Tsit5(), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
starttime = Dates.now();
Gflow_to_01_GMRESf  = solve(prob_Gs_v5_sub_i, Tsit5(), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
endtime = Dates.now();
endtime - starttime # 0.279 sec

starttime = Dates.now();
Gflow_to_01_GMRESg  = solve(prob_Gs_v5_sub_i, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
endtime = Dates.now();
endtime - starttime # 1.6 sec
# 2022-03-31 Gmaps.calc_Gs_SSE_sub_i: 10.338 sec



#Gflow_to_01_GMRESb_QNDF  = solve(prob_Gs_v5_sub_i, QNDF(), save_everystep=true, abstol=1e-6, reltol=1e-6);
#starttime = Dates.now();
#Gflow_to_01_GMRESb_QNDF  = solve(prob_Gs_v5_sub_i, QNDF(), save_everystep=true, abstol=1e-6, reltol=1e-6);
#endtime = Dates.now();
#endtime - starttime # 0.008 sec

Gflow_to_01_GMRESb_FBDF = solve(prob_Gs_v5_sub_i, FBDF(), save_everystep=true, abstol=1e-6, reltol=1e-6);
starttime = Dates.now();
Gflow_to_01_GMRESb_FBDF = solve(prob_Gs_v5_sub_i, FBDF(), save_everystep=true, abstol=1e-6, reltol=1e-6);
endtime = Dates.now();
endtime - starttime # 0.008 sec

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


#res_Gflow_v6b_QNDF = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRESb_QNDF, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
#(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6b_QNDF

res_Gflow_v6b_FBDF = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow_to_01_GMRESb_FBDF, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec_GFv6, iteration_number_GFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6b_FBDF



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
