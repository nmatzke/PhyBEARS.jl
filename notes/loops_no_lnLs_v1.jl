#######################################################
# Test the lnLs from trad SSE v6 (uses a half matrix)
# and Gflow v7 (uses a half matrix)
# 2022-03-24: WORKS!
#######################################################

using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
#using RCall
print("\n...finished loading packages.")

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/notes/")
include("/GitHub/PhyBEARS.jl/notes/loops_no_lnLs_v1.jl")
"""

@testset "loops_no_lnLs_v1.jl" begin


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

solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
#solver_options.solver = Tsit5()
solver_options.save_everystep = true;
solver_options.abstol = 1e-6;
solver_options.reltol = 1e-6;


#######################################################
# Parameter arrays
#######################################################
d_vals = seq(0.0, 1.0, 0.5)
e_vals = seq(0.0, 1.0, 0.5)
j_vals = [0.0, 0.1, 1.0]
lambda_vals = [0.3]
mu_vals1 = seq(0.0, 1.0, 0.5)

length_results = length(d_vals) * length(e_vals) * length(j_vals) * length(lambda_vals) * length(mu_vals1)
res_table = Array{Float64}(undef, length_results, 13);
counter = 0

# Loop through all parameters
# example:
d_val=d_vals[2]; e_val=e_vals[1]; j_val=j_vals[1]; lambda_val=lambda_vals[1]; mu_val1=mu_vals1[1]

for d_val in d_vals
for e_val in e_vals
for j_val in j_vals
for lambda_val in lambda_vals
for mu_val1 in mu_vals1

counter = counter + 1
print(counter)
print(" ")


inputs.bmo.est[inputs.bmo.rownames .== "birthRate"] .= lambda_val;
inputs.bmo.est[inputs.bmo.rownames .== "deathRate"] .= mu_val1;
inputs.bmo.est[inputs.bmo.rownames .== "d"] .= d_val;
inputs.bmo.est[inputs.bmo.rownames .== "e"] .= e_val;
inputs.bmo.est[inputs.bmo.rownames .== "a"] .= 0.0;
inputs.bmo.est[inputs.bmo.rownames .== "j"] .= j_val;
inputs.bmo.est[:] = bmo_updater_v1(inputs.bmo);
p_Ds_v5_updater_v1!(p_Ds_v5, inputs)


# Solve the Es
prob_Es_v6 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v6, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
sol_Es_v6 = solve(prob_Es_v6, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
p_Ds_v5_new = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v6);

# Version 6 ClaSSE standard calculation (triangle of cladogenetic transition matrix instead of full matrix)
#res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
#(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6
# (0.023, 2, -7.780033846857172, -3.8046659173817328, -11.584699764238906, -6.269812932269745)



# Version 7/2 ClaSSE Gflow calculations
#G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
#tmpzero = repeat([0.0], n^2);
#A = reshape(tmpzero, (n,n));
#pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
#tspan = (0.0, 2.5);

#prob_Gs_v5 = DifferentialEquations.ODEProblem(calc_Gs_SSE_v7, G0, tspan, pG);
#Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-6, reltol = 1e-6);
#Gflow = Gflow_to_01_GMRES;

#res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
#(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6
# (0.0, 2, -7.779971880045283, -3.8046557418245324, -11.584627621869815, -6.269741097248855)

end
end
end
end
end

end
