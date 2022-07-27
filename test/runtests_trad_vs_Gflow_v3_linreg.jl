using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using CSV
using RCall

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_trad_vs_Gflow_v3_linreg.jl")
"""

@testset "runtests_trad_vs_Gflow_v3_linreg.jl" begin


#include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
#import .ModelLikes
#include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
#import .Flow

tr = readTopology("((chimp:1,human:1):1,gorilla:2);")

in_params = (birthRate=0.2, deathRate=1.0, d_val=0.5, e_val=0.4, a_val=0.0, j_val=1.5)
numareas = 2;
n = 3;

# Set up the model
inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, in_params=in_params);
(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs;

# Set the tip likelihoods so they aren't identical
res.likes_at_each_nodeIndex_branchTop[1] = [1.0, 0.0, 0.0];			# state 1 for tip #1 (chimps in Africa)
res.normlikes_at_each_nodeIndex_branchTop[1] = [1.0, 0.0, 0.0];	# state 1 for tip #1 (chimps in Africa)
res.likes_at_each_nodeIndex_branchTop[2] = [0.0, 0.0, 1.0];			# state 3 for tip #2 (humans are everywhere)
res.normlikes_at_each_nodeIndex_branchTop[2] = [0.0, 0.0, 1.0];	# state 3 for tip #2 (humans are everywhere)
res.likes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0, 0.0];			# state 1 for tip #4 (gorillas in Africa)
res.normlikes_at_each_nodeIndex_branchTop[4] = [1.0, 0.0, 0.0];	# state 1 for tip #4 (gorillas in Africa)


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

# Change parameter inputs manually
inputs.p_Ds_v5.params.Qij_vals[1:2] .= d_val
inputs.p_Ds_v5.params.Qij_vals[3:4] .= e_val
area1_weight = (3.0-j_val)/3.0 + j_val + j_val
area1_nonj_prob = ((3.0-j_val)/3.0) / area1_weight
area1_j_prob = (j_val) / area1_weight
area2_prob = 1/6
inputs.p_Ds_v5.params.Cijk_weights[1] = (3.0-j_val)/3.0
inputs.p_Ds_v5.params.Cijk_weights[2:5] .= j_val
inputs.p_Ds_v5.params.Cijk_weights[6:12] .= (3.0-j_val)/3.0
inputs.p_Ds_v5.params.Cijk_vals[1] = area1_nonj_prob * lambda_val
inputs.p_Ds_v5.params.Cijk_vals[2:5] .= area1_j_prob * lambda_val
inputs.p_Ds_v5.params.Cijk_vals[6] = area1_nonj_prob * lambda_val
inputs.p_Ds_v5.params.Cijk_vals[7:12] .= area2_prob * lambda_val
inputs.p_Ds_v5.params.mu_vals[:] .= mu_val1

p_Ds_v5 = inputs.p_Ds_v5;




solver_options=construct_SolverOpt();
solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
#solver_options.solver = Tsit5()
solver_options.save_everystep = true;
solver_options.abstol = 1e-6;
solver_options.reltol = 1e-6;


# Solve the Es
prob_Es_v6 = DifferentialEquations.ODEProblem(Flow.parameterized_ClaSSE_Es_v6, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
sol_Es_v6 = solve(prob_Es_v6, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

# Solve the Ds
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v6);

prtCp(p_Ds_v5)

# Version 6 ClaSSE standard calculation (triangle of cladogenetic transition matrix instead of full matrix)
res_nonFlow_v6 = Flow.iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6) = res_nonFlow_v6

# Version 7/2 ClaSSE Gflow calculations
G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_Ds_v5, A=A);
tspan = (0.0, 2.5);

prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE!, G0, tspan, pG);
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-6, reltol = 1e-6);
Gflow = Gflow_to_01_GMRES;

res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6) = res_Gflow_v6


#print("\nTesting with setup_DEC_SSE, full Cijk and Cikj transition matrix, traditional SSE likelihood downpass v6 vs. Gflow v7:\n")
#@test round(Julia_sum_lq_nFv6, digits=3) == round(Julia_sum_lq_GFv6, digits=3)
#@test round(rootstates_lnL_nFv6, digits=3) == round(rootstates_lnL_GFv6, digits=3)
#@test round(Julia_total_lnLs1_nFv6, digits=2) == round(Julia_total_lnLs1_GFv6, digits=2)
#@test round(bgb_lnl_nFv6, digits=2) == round(bgb_lnl_GFv6, digits=2)


tmprow = [d_val, e_val, j_val, lambda_val, mu_val1, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6, bgb_lnl_GFv6]

res_table[counter,:] .= tmprow

end
end
end
end
end

tmprow_names = ["d_val", "e_val", "j_val", "lambda_val", "mu_val1", "Julia_sum_lq_nFv6", "rootstates_lnL_nFv6", "Julia_total_lnLs1_nFv6", "bgb_lnl_nFv6", "Julia_sum_lq_GFv6", "rootstates_lnL_GFv6", "Julia_total_lnLs1_GFv6", "bgb_lnl_GFv6"]

res_table_df = DataFrame(res_table, :auto);

# Add names
rename!(res_table_df, Symbol.(tmprow_names) ) 

resfn = "/GitHub/PhyBEARS.jl/notes/res_table_df_fullMatrix_tradSSEv6_vs_GflowV7.txt"
CSV.write(resfn, res_table_df; delim="\t")

#######################################################
# Call some R code to plot the results
#######################################################
resfn2 = getfn(resfn)
@rput resfn
@rput resfn2
R"""
#R CODE HERE
# resfn = "/GitHub/PhyBEARS.jl/notes/res_table_df_fullMatrix_tradSSEv6_vs_GflowV7.txt"
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

print("\nTesting if lnLs matrix across a range of parameters between trad_SSE_v6 and Gflow_v7, with full matrix:\n")

@test mean(diffs2) < 0.01
@test round(slope, digits=2) == 1.0
@test round(intercept, digits=2) == 0.0
@test round(rsq, digits=2) == 1.0


end
