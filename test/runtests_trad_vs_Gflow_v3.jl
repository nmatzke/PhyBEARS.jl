#######################################################
# Manual comparison on lnLs, on full-size transition matrix (includes Cijk and Cikj) of:
# trad SSE v5
# trad SSE v6
# Gflow SSE v7
#######################################################


using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_trad_vs_Gflow_v3.jl")
"""

@testset "runtests_trad_vs_Gflow_v3.jl" begin


#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes
#include("/GitHub/PhyBEARS.jl/notes/jl")
#import .Flow

tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.2, deathRate=1.0, d_val=0.5, e_val=0.4, a_val=0.0, j_val=1.5)
numareas = 2;
n = 3;

# Set up the model
inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, in_params=in_params);
(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs;

Rnames(p_Ds_v5.p_TFs)
p_Ds_v5.params
prtCp(p_Ds_v5)

# Set the tip likelihoods so they aren't identical
res.likes_at_each_nodeIndex_branchTop[1] = [1.0, 0.0, 0.0]	;		# state 1 for tip #1 (chimps in Africa)
res.normlikes_at_each_nodeIndex_branchTop[1] = [1.0, 0.0, 0.0];	# state 1 for tip #1 (chimps in Africa)
res.likes_at_each_nodeIndex_branchTop[2] = [0.0, 0.0, 1.0]	;		# state 3 for tip #2 (humans are everywhere)
res.normlikes_at_each_nodeIndex_branchTop[2] = [0.0, 0.0, 1.0];	# state 3 for tip #2 (humans are everywhere)
res.likes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0, 0.0];			# state 1 for tip #4 (gorillas in Africa)
res.normlikes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0, 0.0];	# state 1 for tip #4 (gorillas in Africa)

# Solve the Es
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

# Solve the Ds
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);


# Do downpass
solver_options=construct_SolverOpt();
solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
#solver_options.solver = Tsit5()
solver_options.save_everystep = true;
solver_options.abstol = 1e-6;
solver_options.reltol = 1e-6;

# Version 5 ClaSSE standard calculation
res_nonFlow_v5 = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv5, rootstates_lnL_nFv5, Julia_total_lnLs1_nFv5, bgb_lnl_nFv5) = res_nonFlow_v5


# Version 6:
inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, in_params=in_params);
(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs;

# Set the tip likelihoods so they aren't identical
res.likes_at_each_nodeIndex_branchTop[1] = [1.0, 0.0, 0.0];			# state 1 for tip #1 (chimps in Africa)
res.normlikes_at_each_nodeIndex_branchTop[1] = [1.0, 0.0, 0.0];	# state 1 for tip #1 (chimps in Africa)
res.likes_at_each_nodeIndex_branchTop[2] = [0.0, 0.0, 1.0];			# state 3 for tip #2 (humans are everywhere)
res.normlikes_at_each_nodeIndex_branchTop[2] = [0.0, 0.0, 1.0];	# state 3 for tip #2 (humans are everywhere)
res.likes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0, 0.0];			# state 1 for tip #4 (gorillas in Africa)
res.normlikes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0, 0.0];	# state 1 for tip #4 (gorillas in Africa)

solver_options=construct_SolverOpt();
solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
#solver_options.solver = Tsit5()
solver_options.save_everystep = true;
solver_options.abstol = 1e-6;
solver_options.reltol = 1e-6;


# Solve the Es
prob_Es_v6 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v6, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
sol_Es_v6 = solve(prob_Es_v6, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

# Solve the Ds
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v6);

prtCp(p_Ds_v5)

# Version 6 ClaSSE standard calculation (triangle of cladogenetic transition matrix instead of full matrix)
res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6

# Version 7/2 ClaSSE Gflow calculations
G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 2.5);

prob_Gs_v5 = DifferentialEquations.ODEProblem(calc_Gs_SSE!, G0, tspan, pG);
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-6, reltol = 1e-6);
Gflow = Gflow_to_01_GMRES;

res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6


print("\nTesting traditional SSE likelihood downpass v5 vs. v6, with full matrix:\n")
@test round(Julia_sum_lq_nFv5, digits=3) == round(Julia_sum_lq_nFv6, digits=3)
@test round(rootstates_lnL_nFv5, digits=3) == round(rootstates_lnL_nFv6, digits=3)
@test round(Julia_total_lnLs1_nFv5, digits=3) == round(Julia_total_lnLs1_nFv6, digits=3)
@test round(bgb_lnl_nFv5, digits=3) == round(bgb_lnl_nFv6, digits=3)

print("\nTesting traditional SSE likelihood downpass v6 vs. Gflow v7, with full matrix:\n")
@test round(Julia_sum_lq_nFv6, digits=3) == round(Julia_sum_lq_GFv6, digits=3)
@test round(rootstates_lnL_nFv6, digits=3) == round(rootstates_lnL_GFv6, digits=3)
@test round(Julia_total_lnLs1_nFv6, digits=2) == round(Julia_total_lnLs1_GFv6, digits=2)
@test round(bgb_lnl_nFv6, digits=2) == round(bgb_lnl_GFv6, digits=2)


# Check the condition numbers of linear dynamics A and Gflow G (I think Julia does this automatically)
tvals = seq(0.0, 2.5, 0.01);
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
