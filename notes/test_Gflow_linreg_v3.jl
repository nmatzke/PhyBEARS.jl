#######################################################
# Shows the equivalence of traditional and Gflow
# likelihood calculations.
#
# Speeding up the Gflow calculation
# - using half-transition matrix
# - avoiding loops where possible 
# 
# Notes:
# * using the [m] loop should be fine; 
# * using Cijk_vals[Ci_eq_i] instead of Cijk_rates_sub_i, as the latter are another thing to update.
#######################################################

using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Base.Threads  # <-- this is good for @spawn, Distributed.@spawn is BAD, it produces Futures that have to be scheduled etc]
using Combinatorics					# for e.g. combinations()
using DataFrames

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DataFrames  # for DataFrame
using DifferentialEquations
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
using CSV					# for CSV.write (like write.table)
using RCall				# to call e.g. plot()

# List each PhyBEARS code file prefix here
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreeTable
using PhyBEARS.TreePass
using PhyBEARS.TrUtils
using PhyBEARS.SSEs
using PhyBEARS.ModelLikes
using PhyBEARS.Flow

cd("/GitHub/PhyBEARS.jl/notes/")
"""
# Run with:
cd("/GitHub/PhyBEARS.jl/notes/")
include("/GitHub/PhyBEARS.jl/notes/test_Gflow_linreg_v3.jl")
"""
# 
# """
# # Run with:
# include("/GitHub/PhyBEARS.jl/notes/test_Gflow_linreg_v3.jl")
# """
# 
# @testset "Example" begin
# 	@test hello("runtests_BiSSE_tree_n3.jl") == "Hello, runtests_BiSSE_tree_n3.jl"
# #	@test domath(2.0) ≈ 7.0
# end
# 
# 
# #######################################################
# # Do a bunch of tests of the SSE calculation of 
# # Ds, Es, and likelihoods, on
# # branches, nodes, and trees,
# # under a variety of simple and more complex models
# #######################################################
# 
#@testset "test_Gflow_linreg_v1.jl" begin

# Repeat calculation in Julia
#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes
#include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
#import .Flow

tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.2, deathRate=1.0, d_val=0.5, e_val=0.4, a_val=0.0, j_val=1.5)
numareas = 2;
n = 3;

bmo = construct_BioGeoBEARS_model_object();
bmo.type[bmo.rownames .== "j"] .= "free";
bmo.est[bmo.rownames .== "birthRate"] .= in_params.birthRate;
bmo.est[bmo.rownames .== "deathRate"] .= in_params.deathRate;
bmo.est[bmo.rownames .== "d"] .= in_params.d_val;
bmo.est[bmo.rownames .== "e"] .= in_params.e_val;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= in_params.j_val;
bmo.est[:] = bmo_updater_v1(bmo);

# Set up the model
geog_df = DataFrame(tipnames=["chimp","human","gorilla"],A=[1,1,1],B=[0,1,0]); # humans widespread, apes not
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;

#######################################################
# Parameter arrays
#######################################################
d_vals = seq(0.0, 1.0, 0.5);
e_vals = seq(0.0, 1.0, 0.5);
j_vals = [0.0, 0.1, 1.0];
lambda_vals = [0.3];
mu_vals1 = seq(0.0, 1.0, 0.5);
mu_vals2 = seq(0.0, 1.0, 0.5);

length_results = length(d_vals) * length(e_vals) * length(j_vals) * length(lambda_vals) * length(mu_vals1) * length(mu_vals2)
res_table = Array{Float64}(undef, length_results, 11);
counter = 0

trdf = inputs.trdf;
p_Ds_v5 = inputs.p_Ds_v5;
root_age = maximum(trdf[!, :node_age]);

print("\n\nCounting through loops, max= ")
print(length_results)
print("\n")


# Loop through all parameters
# example:
d_val=d_vals[2]; e_val=e_vals[1]; j_val=j_vals[1]; lambda_val=lambda_vals[1]; mu_val1=mu_vals1[1]; mu_val2=mu_vals2[1];
counter = 0
"""
for d_val in d_vals
for e_val in e_vals
for j_val in j_vals
for lambda_val in lambda_vals
for mu_val1 in mu_vals1
"""

counter = counter + 1
print(counter)
print(" ")

# Change parameter inputs 
inputs.bmo.est[inputs.bmo.rownames .== "birthRate"] .= lambda_val;
inputs.bmo.est[inputs.bmo.rownames .== "deathRate"] .= mu_val1;
inputs.bmo.est[inputs.bmo.rownames .== "d"] .= d_val;
inputs.bmo.est[inputs.bmo.rownames .== "e"] .= e_val;
inputs.bmo.est[inputs.bmo.rownames .== "a"] .= 0.0;
inputs.bmo.est[inputs.bmo.rownames .== "j"] .= j_val;
inputs.bmo.est[:] = bmo_updater_v1(inputs.bmo);

#p_Ds_v5 = inputs.p_Ds_v5;
p_Ds_v5_updater_v1!(p_Ds_v5, inputs)
#@time p_Ds_v5_updater_v1!(p_Ds_v5, inputs)

print("\nsol_Es_v6...")

# Solve the Es
#print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v6 = DifferentialEquations.ODEProblem(Flow.parameterized_ClaSSE_Es_v6, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v6 = solve(prob_Es_v6, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
#@time sol_Es_v6 = solve(prob_Es_v6, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

Es_interpolator = sol_Es_v6;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v6);
print("\n...done")


# Solving the Ds once as well

print("\nres_nonFlow_v6...")
# Version 6 ClaSSE standard calculation
res_nonFlow_v6 = Flow.iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6) = res_nonFlow_v6
#@time res_nonFlow_v6 = Flow.iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
print("\n...done")

# Version 6 ClaSSE Gflow calculation
# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
# The "I" requires "include NumericAlgebra"


G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 2.5);

print("\nGflow_to_01_GMRES...")

prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE, G0, tspan, pG);
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-6, reltol = 1e-6);
#@time Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-6, reltol = 1e-6);
Gflow = Gflow_to_01_GMRES;
print("\n...done")

print("\nres_Gflow_v6...")
res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6) = res_Gflow_v6

#@time res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
print("\n...done")

tmprow = [d_val, e_val, j_val, lambda_val, mu_val1, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6]

res_table[counter,:] .= tmprow

end
end
end
end
end

tmprow_names = ["d_val", "e_val", "j_val", "lambda_val", "mu_val1", "Julia_sum_lq_nFv5", "rootstates_lnL_nFv5", "Julia_total_lnLs1_nFv5", "Julia_sum_lq_nFv6", "rootstates_lnL_nFv6", "Julia_total_lnLs1_nFv6", "Julia_sum_lq_GFv6", "rootstates_lnL_GFv6", "Julia_total_lnLs1_GFv6"]

res_table_df = DataFrame(res_table, :auto);

# Add names
rename!(res_table_df, Symbol.(tmprow_names) ) 

resfn = "/GitHub/PhyBEARS.jl/notes/res_table_df_tradSSEv6_vs_GflowV7.txt"
CSV.write(resfn, res_table_df; delim="\t")

#######################################################
# Call some R code to plot the results
#######################################################
resfn2 = getfn(resfn)
@rput resfn2
R"""
#R CODE HERE
# resfn = "/GitHub/PhyBEARS.jl/notes/res_table_df_tradSSEv6_vs_GflowV7.txt"
res_table_df = read.table(resfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(res_table_df)

tmpcolors = "grey70"
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6,  col=tmpcolors, xlab="trad SSE", ylab="Gflow SSE", main=resfn2, xlim=c(-16,0), ylim=c(-16,0))
segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")

# Histogram
diffs2 = res_table_df$Julia_sum_lq_GFv6 - res_table_df$Julia_sum_lq_nFv6
titletxt = paste0("Errors betw. trad & Gflow SSE: ",resfn2)
hist(diffs2, breaks=100, main=titletxt)

lmres = lm(res_table_df$Julia_sum_lq_GFv6~res_table_df$Julia_sum_lq_nFv6)
intercept = lmres$coefficients[1]
slope = lmres$coefficients[2]
rsq = summary(lmres)$r.squared

"""
@rget diffs2;
@rget slope;
@rget intercept;
@rget rsq;


# end


