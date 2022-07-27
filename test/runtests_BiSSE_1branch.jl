using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DataFrames  # for DataFrame
using DifferentialEquations
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA


# List each PhyBEARS code file prefix here
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreeTable
using PhyBEARS.TreePass
using PhyBEARS.TrUtils
using PhyBEARS.SSEs


"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_BiSSE_1branch.jl")
"""

@testset "runtests_SSE.jl" begin
	@test hello("runtests_BiSSE_1branch.jl") == "Hello, runtests_BiSSE_1branch.jl"
#	@test domath(2.0) ≈ 7.0
end

#######################################################
#######################################################
# Testing BiSSE models against:
# 
# * parameterized_ClaSSE         (slower, dumber implementation)
# * parameterized_ClaSSE_v5      (adds Julia speedup tricks)
# * parameterized_ClaSSE_Es_v5 + parameterized_ClaSSE_Ds_v5
#                                (more tricks + precalc Es)
# * calc_Gs_SSE (the "Gflow" approach)
#                                (precalculate Es, then Gflow)
# 
# 
# Do a bunch of tests of the SSE calculation of 
# Ds, Es, and likelihoods, on
# branches, nodes, and trees,
# under a variety of simple and more complex models
#######################################################
#######################################################
#######################################################



@testset "biSSE_1branch_n0" begin

#######################################################
# Calculation of Es and Ds on a single branch
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, birth-death, no Q transitions, branchlength=1)
#
# birthRate = 0.222222222
# deathRate = 0.1
#
# Run with:
# source("/GitHub/PhyBEARS.jl/test/biSSE_1branch_n0.R")
# Truth:
R_result_EsDs = [1 0 0 0 0.800737406147247]

#######################################################
# Note that the result for biSSE_1branch_n0 is 
exp(-0.2222222)
# 0.8007374
#
#######################################################


include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
inputs = ModelLikes.setup_MuSSE_biogeo()


# Repeat calculation in Julia
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.2222222, deathRate=0.0, d_val=0.0, e_val=0.0, a_val=0.0, j_val=0.0)
numstates = 2
n = 2
inputs = ModelLikes.setup_MuSSE_biogeo(numstates, tr; root_age_mult=1.5, in_params=in_params)
prtQi(inputs)
prtCi(inputs)

trdf = inputs.trdf
root_age = maximum(trdf[!, :node_age])
# Es_interpolator = inputs.p_Ds_v5.sol_Es_v5;

# Solve the Es
n = inputs.p_Ds_v5.n
uE = collect(repeat([0.0], n))
Es_tspan = (0.0, 1.4*root_age)
p_Es_v5 = inputs.p_Ds_v5
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, uE, Es_tspan, p_Es_v5)
#print(Es_tspan)

# This solution is a linear interpolator
solver_options=construct_SolverOpt()
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_Es_v5(1.0)

Es_interpolator = sol_Es_v5
prtQi(inputs)
prtCi(inputs)


#######################################################
# Check Es at t=1.0
#######################################################
R_Es = R_result_EsDs[((1+1):(n+1))]      # Get the saved R result
Julia_result_Es = Es_interpolator(1.0)   # Get the interpolated "E" values

# Test:
@test all(round.(Julia_result_Es; digits=4) .== round.(R_Es; digits=4))
# julia> Julia_result_Es
# 2-element Vector{Float64}:
#  0.06845321412227141
#  0.06845321412227141
# 
# julia> R_Es
# 2-element Vector{Float64}:
#  0.0860322055692215
#  0.0860322055692215

### 2021-07-20_ ABCD CHECK HERE!!!
### 2021-09-03_the solution is: ALWAYS use save_everystep=true
### if you want a generic interpolator.  Otherwise you
### only get good answers at the end, and at the saveat=
### points.


#######################################################
# Check Ds at t=1.0
#######################################################
n = inputs.p_Ds_v5.n
p_Ds_v5 = inputs.p_Ds_v5
u0 = collect(repeat([0.0], n))
u0[2] = 1.0
tspan = (0.0, 1.2*root_age)
p_Ds_v5
p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)
prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p)

ground_truth_Ds_interpolatorT = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorG = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorL = solve(prob_Ds_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)


#######################################################
# Check Ds at t=1.0
#######################################################
ground_truth_Ds_interpolatorT(1.0)
ground_truth_Ds_interpolatorG(1.0)
ground_truth_Ds_interpolatorL(1.0)

# Test Tsit5():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=4) .== round.(R_Ds; digits=4))

# Test GMRES():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=4) .== round.(R_Ds; digits=4))

# Test lsoda():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=2) .== round.(R_Ds; digits=2))



#######################################################
# Test likelihood function, v1 (parameterized_ClaSSE)
#######################################################
uEsDs = [0.0, 0.0, 0.0, 1.0]
prob_EsDs = DifferentialEquations.ODEProblem(parameterized_ClaSSE, uEsDs, Es_tspan, p_Es_v5)
sol_EsDs = solve(prob_EsDs, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_EsDs(1.0)

ground_truth_EsDs_v1_interpolatorT = solve(prob_EsDs, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v1_interpolatorG = solve(prob_EsDs, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v1_interpolatorL = solve(prob_EsDs, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Test Tsit5():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test GMRES():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test lsoda():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=2) .== round.(R_EsDs; digits=2))



#######################################################
# Test likelihood function, v5 (parameterized_ClaSSE_v5)
#######################################################

uEsDs = [0.0, 0.0, 0.0, 1.0]
prob_uEsDs_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_v5, uEsDs, Es_tspan, p_Es_v5)
sol_uEsDs_v5 = solve(prob_uEsDs_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

ground_truth_EsDs_v5_interpolatorT = solve(prob_uEsDs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v5_interpolatorG = solve(prob_uEsDs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v5_interpolatorL = solve(prob_uEsDs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

sol_uEsDs_v5 = solve(prob_uEsDs_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_uEsDs_v5(1.0)

ground_truth_EsDs_v5_interpolatorT = solve(prob_uEsDs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);
ground_truth_EsDs_v5_interpolatorG = solve(prob_uEsDs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9);
ground_truth_EsDs_v5_interpolatorL = solve(prob_uEsDs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);

# Test Tsit5():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test GMRES():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test lsoda():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=2) .== round.(R_EsDs; digits=2))


# Check if all of the results are equal, using standard Ds SSE calcs
Tsit5_results = round.(ground_truth_EsDs_v5_interpolatorT.u[length(ground_truth_EsDs_v5_interpolatorT.u)]; digits=4)
GMRES_results = round.(ground_truth_EsDs_v5_interpolatorG.u[length(ground_truth_EsDs_v5_interpolatorG.u)]; digits=4)
LSODA_results = round.(ground_truth_EsDs_v5_interpolatorL.u[length(ground_truth_EsDs_v5_interpolatorL.u)]; digits=4)

@test all(Tsit5_results .== GMRES_results)
@test all(Tsit5_results .== LSODA_results)





#######################################################
# Check if all of the results are equal, using "Flow" Ds SSE calcs
#######################################################
include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
import .Flow

# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2)
A = reshape(tmpzero, (n,n))

# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
# The "I" requires "include NumericAlgebra"
G0 = Matrix{Float64}(I, n, n) 
#p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)
pG = (n=n, p_Ds_v5=p, A=A)
tspan = (0.0, 1.1*root_age)
prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG)


Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Tsit5  = solve(prob_Gs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Lsoda  = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Check that the different interpolators match
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=5)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=5)
@test all(Gmat_GMRES .== Gmat_Tsit5)

# Lsoda requires a less precise match (!)
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=3)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=3)
Gmat_Lsoda = round.(Gflow_to_01_Lsoda(1.0); digits=3)
@test all(Gmat_GMRES .== Gmat_Lsoda)
@test all(Gmat_Tsit5 .== Gmat_Lsoda)


# Calculate the flow, on a single branch, starting from u0
# (tip values, so no fakeX0 calculation needed)
X0 = u0
Xc_GMRES = Gflow_to_01_GMRES(1.0) * X0
Xc_Tsit5 = Gflow_to_01_Tsit5(1.0) * X0
Xc_Lsoda = Gflow_to_01_Lsoda(1.0) * X0

# Compare standard to Flow
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=6) .== round.(Xc_GMRES; digits=6))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=6) .== round.(Xc_Tsit5; digits=6))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

ground_truth_Ds_interpolatorG(1.0)
ground_truth_Ds_interpolatorT(1.0)
ground_truth_Ds_interpolatorL(1.0)

Ds_indices = 1 .+ collect((n+1):(2*n));
R_Ds = R_result_EsDs[Ds_indices]

# Test the R Ds, against the Julia Ds
# Standard calc of Ds
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(R_Ds; digits=3)) # Much worse match
# Flow calc of Ds
@test all(round.(R_Ds; digits=4) .== round.(Xc_GMRES; digits=4))
@test all(round.(R_Ds; digits=4) .== round.(Xc_Tsit5; digits=4))
@test all(round.(R_Ds; digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

R_EDs = vcat(R_Es, R_Ds)
Julia_EDs_GMRES = vcat(Julia_result_Es, Xc_GMRES)
Julia_EDs_Tsit5 = vcat(Julia_result_Es, Xc_Tsit5)
Julia_EDs_Lsoda = vcat(Julia_result_Es, Xc_Lsoda)

GMRES_Rlsoda_diffs = R_EDs .- Julia_EDs_GMRES
Tsit5_Rlsoda_diffs = R_EDs .- Julia_EDs_Tsit5
Lsoda_Rlsoda_diffs = R_EDs .- Julia_EDs_Lsoda

print("\nDifferences between Julia and R biSSE_1branch_n0 calculation:\n")
print("GMRES: ")
print(GMRES_Rlsoda_diffs)
print("\n")
print("Tsit5: ")
print(Tsit5_Rlsoda_diffs)
print("\n")
print("LSODA: ")
print(Lsoda_Rlsoda_diffs)
print("\n")


#######################################################
# Compare against the analytic solution
#######################################################
length_of_result = length(ground_truth_Ds_interpolatorG(1.0))
@test all(round.(ground_truth_Ds_interpolatorG(1.0)[length_of_result]; digits=6) .== round.(exp(-0.2222222); digits=6))
@test all(round.(ground_truth_Ds_interpolatorT(1.0)[length_of_result]; digits=6) .== round.(exp(-0.2222222); digits=6))
@test all(round.(ground_truth_Ds_interpolatorL(1.0)[length_of_result]; digits=3) .== round.(exp(-0.2222222); digits=3))


end # END @testset "biSSE_1branch_n0" begin






@testset "biSSE_1branch_n1" begin

#######################################################
# Calculation of Es and Ds on a single branch
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, birth-death, no Q transitions, branchlength=1)
#
# birthRate = 0.222222222
# deathRate = 0.1
#
# Run with:
# source("/GitHub/PhyBEARS.jl/test/biSSE_1branch_n1.R")
# Truth:
R_result_EsDs = [1 0.0860322055692215 0.0860322055692215 0 0.739232931655601]
#######################################################


include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
inputs = ModelLikes.setup_MuSSE_biogeo()


# Repeat calculation in Julia
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.2222222, deathRate=0.1, d_val=0.0, e_val=0.0, a_val=0.0, j_val=0.0)
numstates = 2
n = 2
inputs = ModelLikes.setup_MuSSE_biogeo(numstates, tr; root_age_mult=1.5, in_params=in_params)
prtQi(inputs)
prtCi(inputs)

trdf = inputs.trdf
root_age = maximum(trdf[!, :node_age])
# Es_interpolator = inputs.p_Ds_v5.sol_Es_v5;

# Solve the Es
n = inputs.p_Ds_v5.n
uE = collect(repeat([0.0], n))
Es_tspan = (0.0, 1.4*root_age)
p_Es_v5 = inputs.p_Ds_v5
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, uE, Es_tspan, p_Es_v5)
print(Es_tspan)

# This solution is a linear interpolator
solver_options=construct_SolverOpt()
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_Es_v5(1.0)

Es_interpolator = sol_Es_v5
prtQi(inputs)
prtCi(inputs)


#######################################################
# Check Es at t=1.0
#######################################################
R_Es = R_result_EsDs[((1+1):(n+1))]      # Get the saved R result
Julia_result_Es = Es_interpolator(1.0)   # Get the interpolated "E" values

# Test:
@test all(round.(Julia_result_Es; digits=4) .== round.(R_Es; digits=4))
# julia> Julia_result_Es
# 2-element Vector{Float64}:
#  0.06845321412227141
#  0.06845321412227141
# 
# julia> R_Es
# 2-element Vector{Float64}:
#  0.0860322055692215
#  0.0860322055692215

### 2021-07-20_ ABCD CHECK HERE!!!
### 2021-09-03_the solution is: ALWAYS use save_everystep=true
### if you want a generic interpolator.  Otherwise you
### only get good answers at the end, and at the saveat=
### points.


#######################################################
# Check Ds at t=1.0
#######################################################
n = inputs.p_Ds_v5.n
p_Ds_v5 = inputs.p_Ds_v5
u0 = collect(repeat([0.0], n))
u0[2] = 1.0
tspan = (0.0, 1.2*root_age)
p_Ds_v5
p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)
prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p)

ground_truth_Ds_interpolatorT = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorG = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorL = solve(prob_Ds_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)


#######################################################
# Check Ds at t=1.0
#######################################################
ground_truth_Ds_interpolatorT(1.0)
ground_truth_Ds_interpolatorG(1.0)
ground_truth_Ds_interpolatorL(1.0)

# Test Tsit5():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=4) .== round.(R_Ds; digits=4))

# Test GMRES():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=4) .== round.(R_Ds; digits=4))

# Test lsoda():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=2) .== round.(R_Ds; digits=2))



#######################################################
# Test likelihood function, v1 (parameterized_ClaSSE)
#######################################################
uEsDs = [0.0, 0.0, 0.0, 1.0]
prob_EsDs = DifferentialEquations.ODEProblem(parameterized_ClaSSE, uEsDs, Es_tspan, p_Es_v5)
sol_EsDs = solve(prob_EsDs, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_EsDs(1.0)

ground_truth_EsDs_v1_interpolatorT = solve(prob_EsDs, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v1_interpolatorG = solve(prob_EsDs, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v1_interpolatorL = solve(prob_EsDs, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Test Tsit5():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test GMRES():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test lsoda():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=2) .== round.(R_EsDs; digits=2))



#######################################################
# Test likelihood function, v5 (parameterized_ClaSSE_v5)
#######################################################

uEsDs = [0.0, 0.0, 0.0, 1.0]
prob_uEsDs_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_v5, uEsDs, Es_tspan, p_Es_v5)
sol_uEsDs_v5 = solve(prob_uEsDs_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

ground_truth_EsDs_v5_interpolatorT = solve(prob_uEsDs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v5_interpolatorG = solve(prob_uEsDs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v5_interpolatorL = solve(prob_uEsDs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

sol_uEsDs_v5 = solve(prob_uEsDs_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_uEsDs_v5(1.0)

ground_truth_EsDs_v5_interpolatorT = solve(prob_uEsDs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);
ground_truth_EsDs_v5_interpolatorG = solve(prob_uEsDs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9);
ground_truth_EsDs_v5_interpolatorL = solve(prob_uEsDs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);

# Test Tsit5():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test GMRES():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test lsoda():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=2) .== round.(R_EsDs; digits=2))


# Check if all of the results are equal, using standard Ds SSE calcs
Tsit5_results = round.(ground_truth_EsDs_v5_interpolatorT.u[length(ground_truth_EsDs_v5_interpolatorT.u)]; digits=4)
GMRES_results = round.(ground_truth_EsDs_v5_interpolatorG.u[length(ground_truth_EsDs_v5_interpolatorG.u)]; digits=4)
LSODA_results = round.(ground_truth_EsDs_v5_interpolatorL.u[length(ground_truth_EsDs_v5_interpolatorL.u)]; digits=4)

@test all(Tsit5_results .== GMRES_results)
@test all(Tsit5_results .== LSODA_results)





#######################################################
# Check if all of the results are equal, using "Flow" Ds SSE calcs
#######################################################
include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
import .Flow

# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2)
A = reshape(tmpzero, (n,n))

# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
# The "I" requires "include NumericAlgebra"
G0 = Matrix{Float64}(I, n, n) 
#p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)
pG = (n=n, p_Ds_v5=p, A=A)
tspan = (0.0, 1.1*root_age)
prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG)


Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Tsit5  = solve(prob_Gs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Lsoda  = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Check that the different interpolators match
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=5)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=5)
@test all(Gmat_GMRES .== Gmat_Tsit5)

# Lsoda requires a less precise match (!)
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=3)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=3)
Gmat_Lsoda = round.(Gflow_to_01_Lsoda(1.0); digits=3)
@test all(Gmat_GMRES .== Gmat_Lsoda)
@test all(Gmat_Tsit5 .== Gmat_Lsoda)


# Calculate the flow, on a single branch, starting from u0
# (tip values, so no fakeX0 calculation needed)
X0 = u0
Xc_GMRES = Gflow_to_01_GMRES(1.0) * X0
Xc_Tsit5 = Gflow_to_01_Tsit5(1.0) * X0
Xc_Lsoda = Gflow_to_01_Lsoda(1.0) * X0

# Compare standard to Flow
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=6) .== round.(Xc_GMRES; digits=6))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=6) .== round.(Xc_Tsit5; digits=6))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

ground_truth_Ds_interpolatorG(1.0)
ground_truth_Ds_interpolatorT(1.0)
ground_truth_Ds_interpolatorL(1.0)

Ds_indices = 1 .+ collect((n+1):(2*n));
R_Ds = R_result_EsDs[Ds_indices]

# Test the R Ds, against the Julia Ds
# Standard calc of Ds
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(R_Ds; digits=3)) # Much worse match
# Flow calc of Ds
@test all(round.(R_Ds; digits=4) .== round.(Xc_GMRES; digits=4))
@test all(round.(R_Ds; digits=4) .== round.(Xc_Tsit5; digits=4))
@test all(round.(R_Ds; digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

R_EDs = vcat(R_Es, R_Ds)
Julia_EDs_GMRES = vcat(Julia_result_Es, Xc_GMRES)
Julia_EDs_Tsit5 = vcat(Julia_result_Es, Xc_Tsit5)
Julia_EDs_Lsoda = vcat(Julia_result_Es, Xc_Lsoda)

GMRES_Rlsoda_diffs = R_EDs .- Julia_EDs_GMRES
Tsit5_Rlsoda_diffs = R_EDs .- Julia_EDs_Tsit5
Lsoda_Rlsoda_diffs = R_EDs .- Julia_EDs_Lsoda

print("\nDifferences between Julia and R biSSE_1branch_n1 calculation:\n")
print("GMRES: ")
print(GMRES_Rlsoda_diffs)
print("\n")
print("Tsit5: ")
print(Tsit5_Rlsoda_diffs)
print("\n")
print("LSODA: ")
print(Lsoda_Rlsoda_diffs)
print("\n")

end # END @testset "biSSE_1branch_n1" begin






@testset "biSSE_1branch_n2" begin

#######################################################
# Calculation of Es and Ds on a single branch
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, pure birth, no Q transitions, branchlength=1)
#
# Run with:
# source("/GitHub/PhyBEARS.jl/test/biSSE_1branch_n2.R")
# Truth:
R_result_EsDs = [1 0.304473654892522 0.304473654892522 0.0523757642945074 0.525500584759885]
#######################################################

include("/GitHub/PhyBEARS.jl/src/StateSpace.jl")
import .StateSpace
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
inputs = ModelLikes.setup_MuSSE_biogeo()
prtQi(inputs)
prtCi(inputs)


# Repeat calculation in Julia
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.222222222, deathRate=0.4, d_val=0.0, e_val=0.0, a_val=0.1, j_val=0.0)
numstates = 2
n = 2
inputs = ModelLikes.setup_MuSSE_biogeo(numstates, tr; root_age_mult=1.5, in_params=in_params)
prtQi(inputs)
prtCi(inputs)

trdf = inputs.trdf
root_age = maximum(trdf[!, :node_age])
# Es_interpolator = inputs.p_Ds_v5.sol_Es_v5;

# Solve the Es
n = inputs.p_Ds_v5.n
uE = collect(repeat([0.0], n))
Es_tspan = (0.0, 1.4*root_age)
p_Es_v5 = inputs.p_Ds_v5
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, uE, Es_tspan, p_Es_v5)
print(Es_tspan)

# This solution is a linear interpolator
solver_options=construct_SolverOpt()
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_Es_v5(1.0)

Es_interpolator = sol_Es_v5
prtQi(inputs)
prtCi(inputs)


#######################################################
# Check Es at t=1.0
#######################################################
R_Es = R_result_EsDs[((1+1):(n+1))]      # Get the saved R result
Julia_result_Es = Es_interpolator(1.0)   # Get the interpolated "E" values

# Test:
@test all(round.(Julia_result_Es; digits=4) .== round.(R_Es; digits=4))
# julia> Julia_result_Es
# 2-element Vector{Float64}:
#  0.06845321412227141
#  0.06845321412227141
# 
# julia> R_Es
# 2-element Vector{Float64}:
#  0.0860322055692215
#  0.0860322055692215

### 2021-07-20_ ABCD CHECK HERE!!!
### 2021-09-03_the solution is: ALWAYS use save_everystep=true
### if you want a generic interpolator.  Otherwise you
### only get good answers at the end, and at the saveat=
### points.


#######################################################
# Check Ds at t=1.0
#######################################################
n = inputs.p_Ds_v5.n
p_Ds_v5 = inputs.p_Ds_v5
u0 = collect(repeat([0.0], n))
u0[2] = 1.0
tspan = (0.0, 1.2*root_age)
p_Ds_v5
p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)
prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p)

ground_truth_Ds_interpolatorT = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorG = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorL = solve(prob_Ds_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)


#######################################################
# Check Ds at t=1.0
#######################################################
ground_truth_Ds_interpolatorT(1.0)
ground_truth_Ds_interpolatorG(1.0)
ground_truth_Ds_interpolatorL(1.0)

# Test Tsit5():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=4) .== round.(R_Ds; digits=4))

# Test GMRES():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=4) .== round.(R_Ds; digits=4))

# Test lsoda():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=2) .== round.(R_Ds; digits=2))



#######################################################
# Test likelihood function, v1 (parameterized_ClaSSE)
#######################################################
uEsDs = [0.0, 0.0, 0.0, 1.0]
prob_EsDs = DifferentialEquations.ODEProblem(parameterized_ClaSSE, uEsDs, Es_tspan, p_Es_v5)
sol_EsDs = solve(prob_EsDs, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_EsDs(1.0)

ground_truth_EsDs_v1_interpolatorT = solve(prob_EsDs, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v1_interpolatorG = solve(prob_EsDs, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v1_interpolatorL = solve(prob_EsDs, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Test Tsit5():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test GMRES():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test lsoda():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=2) .== round.(R_EsDs; digits=2))



#######################################################
# Test likelihood function, v5 (parameterized_ClaSSE_v5)
#######################################################

uEsDs = [0.0, 0.0, 0.0, 1.0]
prob_uEsDs_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_v5, uEsDs, Es_tspan, p_Es_v5)
sol_uEsDs_v5 = solve(prob_uEsDs_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

ground_truth_EsDs_v5_interpolatorT = solve(prob_uEsDs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v5_interpolatorG = solve(prob_uEsDs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v5_interpolatorL = solve(prob_uEsDs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

sol_uEsDs_v5 = solve(prob_uEsDs_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_uEsDs_v5(1.0)

ground_truth_EsDs_v5_interpolatorT = solve(prob_uEsDs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);
ground_truth_EsDs_v5_interpolatorG = solve(prob_uEsDs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9);
ground_truth_EsDs_v5_interpolatorL = solve(prob_uEsDs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);

# Test Tsit5():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test GMRES():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test lsoda():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=2) .== round.(R_EsDs; digits=2))


# Check if all of the results are equal, using standard Ds SSE calcs
Tsit5_results = round.(ground_truth_EsDs_v5_interpolatorT.u[length(ground_truth_EsDs_v5_interpolatorT.u)]; digits=4)
GMRES_results = round.(ground_truth_EsDs_v5_interpolatorG.u[length(ground_truth_EsDs_v5_interpolatorG.u)]; digits=4)
LSODA_results = round.(ground_truth_EsDs_v5_interpolatorL.u[length(ground_truth_EsDs_v5_interpolatorL.u)]; digits=4)

@test all(Tsit5_results .== GMRES_results)
@test all(Tsit5_results .== LSODA_results)





#######################################################
# Check if all of the results are equal, using "Flow" Ds SSE calcs
#######################################################
include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
import .Flow

# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2)
A = reshape(tmpzero, (n,n))

# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
# The "I" requires "include NumericAlgebra"
G0 = Matrix{Float64}(I, n, n) 
#p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)
pG = (n=n, p_Ds_v5=p, A=A)
tspan = (0.0, 1.1*root_age)
prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG)


Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Tsit5  = solve(prob_Gs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Lsoda  = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Check that the different interpolators match
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=5)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=5)
@test all(Gmat_GMRES .== Gmat_Tsit5)

# Lsoda requires a less precise match (!)
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=3)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=3)
Gmat_Lsoda = round.(Gflow_to_01_Lsoda(1.0); digits=3)
@test all(Gmat_GMRES .== Gmat_Lsoda)
@test all(Gmat_Tsit5 .== Gmat_Lsoda)


# Calculate the flow, on a single branch, starting from u0
# (tip values, so no fakeX0 calculation needed)
X0 = u0
Xc_GMRES = Gflow_to_01_GMRES(1.0) * X0
Xc_Tsit5 = Gflow_to_01_Tsit5(1.0) * X0
Xc_Lsoda = Gflow_to_01_Lsoda(1.0) * X0

# Compare standard to Flow
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=6) .== round.(Xc_GMRES; digits=6))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=6) .== round.(Xc_Tsit5; digits=6))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

ground_truth_Ds_interpolatorG(1.0)
ground_truth_Ds_interpolatorT(1.0)
ground_truth_Ds_interpolatorL(1.0)

Ds_indices = 1 .+ collect((n+1):(2*n));
R_Ds = R_result_EsDs[Ds_indices]

# Test the R Ds, against the Julia Ds
# Standard calc of Ds
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(R_Ds; digits=3)) # Much worse match
# Flow calc of Ds
@test all(round.(R_Ds; digits=4) .== round.(Xc_GMRES; digits=4))
@test all(round.(R_Ds; digits=4) .== round.(Xc_Tsit5; digits=4))
@test all(round.(R_Ds; digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

R_EDs = vcat(R_Es, R_Ds)
Julia_EDs_GMRES = vcat(Julia_result_Es, Xc_GMRES)
Julia_EDs_Tsit5 = vcat(Julia_result_Es, Xc_Tsit5)
Julia_EDs_Lsoda = vcat(Julia_result_Es, Xc_Lsoda)

GMRES_Rlsoda_diffs = R_EDs .- Julia_EDs_GMRES
Tsit5_Rlsoda_diffs = R_EDs .- Julia_EDs_Tsit5
Lsoda_Rlsoda_diffs = R_EDs .- Julia_EDs_Lsoda

print("\nDifferences between Julia and R biSSE_1branch_n1 calculation:\n")
print("GMRES: ")
print(GMRES_Rlsoda_diffs)
print("\n")
print("Tsit5: ")
print(Tsit5_Rlsoda_diffs)
print("\n")
print("LSODA: ")
print(Lsoda_Rlsoda_diffs)
print("\n")

end # END @testset "biSSE_1branch_n2" begin






@testset "biSSE_1branch_n3" begin

#######################################################
# Calculation of Es and Ds on a single branch
# Example BiSSE calculation
# result_EsDs_biSSE_1branch_pureBirth_bl1
# (1 branch, pure birth, no Q transitions, branchlength=1)
#
# Run with:
# source("/GitHub/PhyBEARS.jl/test/biSSE_1branch_n3.R")
# Truth:
R_result_EsDs = [1 0.0860322055692216 0.0860322055692216 0.0670001652615721 0.672232766394028]
#######################################################

include("/GitHub/PhyBEARS.jl/src/StateSpace.jl")
import .StateSpace
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
inputs = ModelLikes.setup_MuSSE_biogeo()
prtQi(inputs)
prtCi(inputs)


# Repeat calculation in Julia
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.222222222, deathRate=0.1, d_val=0.0, e_val=0.0, a_val=0.1, j_val=0.0)
numstates = 2
n = 2
inputs = ModelLikes.setup_MuSSE_biogeo(numstates, tr; root_age_mult=1.5, in_params=in_params)
prtQi(inputs)
prtCi(inputs)

trdf = inputs.trdf
root_age = maximum(trdf[!, :node_age])
# Es_interpolator = inputs.p_Ds_v5.sol_Es_v5;

# Solve the Es
n = inputs.p_Ds_v5.n
uE = collect(repeat([0.0], n))
Es_tspan = (0.0, 1.4*root_age)
p_Es_v5 = inputs.p_Ds_v5
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, uE, Es_tspan, p_Es_v5)
print(Es_tspan)

# This solution is a linear interpolator
solver_options=construct_SolverOpt()
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_Es_v5(1.0)

Es_interpolator = sol_Es_v5
prtQi(inputs)
prtCi(inputs)


#######################################################
# Check Es at t=1.0
#######################################################
R_Es = R_result_EsDs[((1+1):(n+1))]      # Get the saved R result
Julia_result_Es = Es_interpolator(1.0)   # Get the interpolated "E" values

# Test:
@test all(round.(Julia_result_Es; digits=4) .== round.(R_Es; digits=4))
# julia> Julia_result_Es
# 2-element Vector{Float64}:
#  0.06845321412227141
#  0.06845321412227141
# 
# julia> R_Es
# 2-element Vector{Float64}:
#  0.0860322055692215
#  0.0860322055692215

### 2021-07-20_ ABCD CHECK HERE!!!
### 2021-09-03_the solution is: ALWAYS use save_everystep=true
### if you want a generic interpolator.  Otherwise you
### only get good answers at the end, and at the saveat=
### points.


#######################################################
# Check Ds at t=1.0
#######################################################
n = inputs.p_Ds_v5.n
p_Ds_v5 = inputs.p_Ds_v5
u0 = collect(repeat([0.0], n))
u0[2] = 1.0
tspan = (0.0, 1.2*root_age)
p_Ds_v5
p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)
prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p)

ground_truth_Ds_interpolatorT = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorG = solve(prob_Ds_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_Ds_interpolatorL = solve(prob_Ds_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)


#######################################################
# Check Ds at t=1.0
#######################################################
ground_truth_Ds_interpolatorT(1.0)
ground_truth_Ds_interpolatorG(1.0)
ground_truth_Ds_interpolatorL(1.0)

# Test Tsit5():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=4) .== round.(R_Ds; digits=4))

# Test GMRES():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=4) .== round.(R_Ds; digits=4))

# Test lsoda():
R_Ds = R_result_EsDs[((1+n+1):(1+n+n))]      # Get the saved R result
Julia_result_Ds = ground_truth_Ds_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_Ds; digits=2) .== round.(R_Ds; digits=2))



#######################################################
# Test likelihood function, v1 (parameterized_ClaSSE)
#######################################################
uEsDs = [0.0, 0.0, 0.0, 1.0]
prob_EsDs = DifferentialEquations.ODEProblem(parameterized_ClaSSE, uEsDs, Es_tspan, p_Es_v5)
sol_EsDs = solve(prob_EsDs, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_EsDs(1.0)

ground_truth_EsDs_v1_interpolatorT = solve(prob_EsDs, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v1_interpolatorG = solve(prob_EsDs, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v1_interpolatorL = solve(prob_EsDs, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Test Tsit5():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test GMRES():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test lsoda():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v1_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=2) .== round.(R_EsDs; digits=2))



#######################################################
# Test likelihood function, v5 (parameterized_ClaSSE_v5)
#######################################################

uEsDs = [0.0, 0.0, 0.0, 1.0]
prob_uEsDs_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_v5, uEsDs, Es_tspan, p_Es_v5)
sol_uEsDs_v5 = solve(prob_uEsDs_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

ground_truth_EsDs_v5_interpolatorT = solve(prob_uEsDs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v5_interpolatorG = solve(prob_uEsDs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
ground_truth_EsDs_v5_interpolatorL = solve(prob_uEsDs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

sol_uEsDs_v5 = solve(prob_uEsDs_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_uEsDs_v5(1.0)

ground_truth_EsDs_v5_interpolatorT = solve(prob_uEsDs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);
ground_truth_EsDs_v5_interpolatorG = solve(prob_uEsDs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9);
ground_truth_EsDs_v5_interpolatorL = solve(prob_uEsDs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);

# Test Tsit5():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorT(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test GMRES():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorG(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=4) .== round.(R_EsDs; digits=4))

# Test lsoda():
R_EsDs = R_result_EsDs[(1+1):(1+n+n)]      # Get the saved R result
Julia_result_EsDs = ground_truth_EsDs_v5_interpolatorL(1.0)   # Get the interpolated "E" values
@test all(round.(Julia_result_EsDs; digits=2) .== round.(R_EsDs; digits=2))


# Check if all of the results are equal, using standard Ds SSE calcs
Tsit5_results = round.(ground_truth_EsDs_v5_interpolatorT.u[length(ground_truth_EsDs_v5_interpolatorT.u)]; digits=4)
GMRES_results = round.(ground_truth_EsDs_v5_interpolatorG.u[length(ground_truth_EsDs_v5_interpolatorG.u)]; digits=4)
LSODA_results = round.(ground_truth_EsDs_v5_interpolatorL.u[length(ground_truth_EsDs_v5_interpolatorL.u)]; digits=4)

@test all(Tsit5_results .== GMRES_results)
@test all(Tsit5_results .== LSODA_results)





#######################################################
# Check if all of the results are equal, using "Flow" Ds SSE calcs
#######################################################
include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
import .Flow

# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2)
A = reshape(tmpzero, (n,n))

# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
# The "I" requires "include NumericAlgebra"
G0 = Matrix{Float64}(I, n, n) 
#p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)

p = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5)
pG = (n=n, p_Ds_v5=p, A=A)
tspan = (0.0, 1.1*root_age)
prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG)


Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Tsit5  = solve(prob_Gs_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
Gflow_to_01_Lsoda  = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)

# Check that the different interpolators match
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=5)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=5)
@test all(Gmat_GMRES .== Gmat_Tsit5)

# Lsoda requires a less precise match (!)
Gmat_GMRES = round.(Gflow_to_01_GMRES(1.0); digits=3)
Gmat_Tsit5 = round.(Gflow_to_01_Tsit5(1.0); digits=3)
Gmat_Lsoda = round.(Gflow_to_01_Lsoda(1.0); digits=3)
@test all(Gmat_GMRES .== Gmat_Lsoda)
@test all(Gmat_Tsit5 .== Gmat_Lsoda)


# Calculate the flow, on a single branch, starting from u0
# (tip values, so no fakeX0 calculation needed)
X0 = u0
Xc_GMRES = Gflow_to_01_GMRES(1.0) * X0
Xc_Tsit5 = Gflow_to_01_Tsit5(1.0) * X0
Xc_Lsoda = Gflow_to_01_Lsoda(1.0) * X0

# Compare standard to Flow
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=6) .== round.(Xc_GMRES; digits=6))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=6) .== round.(Xc_Tsit5; digits=6))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

ground_truth_Ds_interpolatorG(1.0)
ground_truth_Ds_interpolatorT(1.0)
ground_truth_Ds_interpolatorL(1.0)

Ds_indices = 1 .+ collect((n+1):(2*n));
R_Ds = R_result_EsDs[Ds_indices]

# Test the R Ds, against the Julia Ds
# Standard calc of Ds
@test all(round.(ground_truth_Ds_interpolatorG(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorT(1.0); digits=4) .== round.(R_Ds; digits=4))
@test all(round.(ground_truth_Ds_interpolatorL(1.0); digits=3) .== round.(R_Ds; digits=3)) # Much worse match
# Flow calc of Ds
@test all(round.(R_Ds; digits=4) .== round.(Xc_GMRES; digits=4))
@test all(round.(R_Ds; digits=4) .== round.(Xc_Tsit5; digits=4))
@test all(round.(R_Ds; digits=3) .== round.(Xc_Lsoda; digits=3)) # Much worse match

R_EDs = vcat(R_Es, R_Ds)
Julia_EDs_GMRES = vcat(Julia_result_Es, Xc_GMRES)
Julia_EDs_Tsit5 = vcat(Julia_result_Es, Xc_Tsit5)
Julia_EDs_Lsoda = vcat(Julia_result_Es, Xc_Lsoda)

GMRES_Rlsoda_diffs = R_EDs .- Julia_EDs_GMRES
Tsit5_Rlsoda_diffs = R_EDs .- Julia_EDs_Tsit5
Lsoda_Rlsoda_diffs = R_EDs .- Julia_EDs_Lsoda

print("\nDifferences between Julia and R biSSE_1branch_n1 calculation:\n")
print("GMRES: ")
print(GMRES_Rlsoda_diffs)
print("\n")
print("Tsit5: ")
print(Tsit5_Rlsoda_diffs)
print("\n")
print("LSODA: ")
print(Lsoda_Rlsoda_diffs)
print("\n")

end # END @testset "biSSE_1branch_n3" begin


