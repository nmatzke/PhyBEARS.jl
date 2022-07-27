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
using CSV					# for CSV.write (like write.table)


# List each PhyBEARS code file prefix here
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreeTable
using PhyBEARS.TreePass
using PhyBEARS.TrUtils
using PhyBEARS.SSEs

cd("/GitHub/PhyBEARS.jl/notes/")
"""
# Run with:
cd("/GitHub/PhyBEARS.jl/notes/")
include("/GitHub/PhyBEARS.jl/notes/test_Gflow_linreg_v1.jl")
"""
# 
# """
# # Run with:
# include("/GitHub/PhyBEARS.jl/test/runtests_BiSSE_tree_n3.jl")
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
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes
include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
import .Flow

tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
in_params = (birthRate=0.2, deathRate=1.0, d_val=0.5, e_val=0.4, a_val=0.0, j_val=1.5)
numareas = 2
n = 3

#######################################################
# Parameter arrays
#######################################################
d_vals = seq(0.0, 1.0, 0.5)
e_vals = seq(0.0, 1.0, 0.5)
j_vals = [0.0, 0.1, 1.0]
lambda_vals = [0.3]
mu_vals1 = seq(0.0, 1.0, 0.5)
mu_vals2 = seq(0.0, 1.0, 0.5)

length_results = length(d_vals) * length(e_vals) * length(j_vals) * length(lambda_vals) * length(mu_vals1) * length(mu_vals2)
res_table = Array{Float64}(undef, length_results, 15);
counter = 0

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, in_params=in_params)
(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs

res.likes_at_each_nodeIndex_branchTop[2] = [1.0, 0.0, 0.0]			# state 1 for tip #2
res.normlikes_at_each_nodeIndex_branchTop[2] = [1.0, 0.0, 0.0]	# state 1 for tip #2
inputs.res.likes_at_each_nodeIndex_branchTop
inputs.res.normlikes_at_each_nodeIndex_branchTop

trdf = inputs.trdf
p_Ds_v5 = inputs.p_Ds_v5
root_age = maximum(trdf[!, :node_age])

print("\n\nCounting through loops, max= ")
print(length_results)
print("\n")


# Loop through all parameters
# example:
d_val=d_vals[2]; e_val=e_vals[1]; j_val=j_vals[1]; lambda_val=lambda_vals[1]; mu_val1=mu_vals1[1]; mu_val2=mu_vals2[1];

for d_val in d_vals
for e_val in e_vals
for j_val in j_vals
for lambda_val in lambda_vals
for mu_val1 in mu_vals1
for mu_val2 in mu_vals2

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
inputs.p_Ds_v5.params.mu_vals[1] = mu_val1
inputs.p_Ds_v5.params.mu_vals[2] = mu_val2
inputs.p_Ds_v5.params.mu_vals[3] = mu_val1


p_Ds_v5 = inputs.p_Ds_v5

# Solve the Es
#print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);


# Solving the Ds once as well
uD = collect(repeat([1.0], n))
prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, uD, Es_tspan, p_Ds_v5)

# This solution is a linear interpolator
sol_Ds_v5 = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_Ds_v5(1.5)
sol_Ds_v5(2.5)

p_inputsG = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uD=uD, uE=p_Ds_v5.uE,  sol_Ds_v5=sol_Ds_v5, sol_Es_v5=sol_Es_v5)



# Do downpass
solver_options=construct_SolverOpt()
solver_options.solver = CVODE_BDF(linear_solver=:GMRES)
#solver_options.solver = Tsit5()
solver_options.save_everystep = true
solver_options.abstol = 1e-6
solver_options.reltol = 1e-6


# Version 5 ClaSSE standard calculation
res_nonFlow_v5 = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true)
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv5, rootstates_lnL_nFv5, Julia_total_lnLs1_nFv5) = res_nonFlow_v5

# Version 6 ClaSSE standard calculation
res_nonFlow_v6 = iterative_downpass_nonparallel_ClaSSE_v6!(res; trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10, return_lnLs=true);
(total_calctime_in_sec, iteration_number, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6) = res_nonFlow_v6

# Version 6 ClaSSE Gflow calculation
# Map the likelihood "flow" of Ds, G (or Gmap or Psi).
# Start with an identity matrix
# The "I" requires "include NumericAlgebra"


G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_inputsG, A=A);
tspan = (0.0, 2.5);

prob_Gs_v5 = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE_v7!, G0, tspan, pG);

saveats = seq(tspan[1], tspan[2], tspan[2]/10000);
Gflow_to_01_GMRES  = solve(prob_Gs_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol = 1e-6, reltol = 1e-6);#, saveat=saveats);
#Gflow_to_01_GMRES  = solve(prob_Gs_v5, lsoda(), save_everystep=true, abstol = 1e-6, reltol = 1e-6, saveat=saveats);

# Downpass with the Gflow system
Gflow = Gflow_to_01_GMRES;



# Gflow calculated in steps
#include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
#import .Flow
Gseg_times = seq(0.0, 2.5, 0.1)  # Times for Gsegments
pG = (n=n, p_Ds_v5=p_inputsG, A=A);
#Gflows_array = Flow.calc_Gmaps_SSE(pG, Gseg_times)

#res_Gflow_v6 = Flow.iterative_downpass_Gflow_nonparallel_v3!(res; trdf, p_Ds_v5, Gflow, Gseg_times, Gflows_array, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true, include_null_range=true)

# error < 1 lnL
res_Gflow_v6 = iterative_downpass_Gflow_nonparallel_v2!(res; trdf, p_Ds_v5=p_Ds_v5, Gflow=Gflow, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=true);

(total_calctime_in_sec, iteration_number, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6) = res_Gflow_v6


lnl_diff = res_nonFlow_v6[3] - res_Gflow_v6[3]

tmprow = [d_val, e_val, j_val, lambda_val, mu_val1, mu_val2, Julia_sum_lq_nFv5, rootstates_lnL_nFv5, Julia_total_lnLs1_nFv5, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, Julia_sum_lq_GFv6, rootstates_lnL_GFv6, Julia_total_lnLs1_GFv6]

res_table[counter,:] .= tmprow

end
end
end
end
end
end

tmprow_names = ["d_val", "e_val", "j_val", "lambda_val", "mu_val1", "mu_val2", "Julia_sum_lq_nFv5", "rootstates_lnL_nFv5", "Julia_total_lnLs1_nFv5", "Julia_sum_lq_nFv6", "rootstates_lnL_nFv6", "Julia_total_lnLs1_nFv6", "Julia_sum_lq_GFv6", "rootstates_lnL_GFv6", "Julia_total_lnLs1_GFv6"]

res_table_df = DataFrame(res_table, :auto);

# Add names
rename!(res_table_df, Symbol.(tmprow_names) ) 

resfn = "/GitHub/PhyBEARS.jl/notes/res_table_df_Ci_eq_i_GflowMap_v9.txt"
CSV.write(resfn, res_table_df; delim="\t")



include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
import .Flow

tvals = seq(0.0, 2.5, 0.01)
kappa_Arates_df = Flow.check_linearDynamics_of_As(tvals, p_Ds_v5; max_condition_number=1e8)
calc_Gs_SSE_condnums! = (dG, G, pG, t) 



include("/GitHub/PhyBEARS.jl/notes/Flow.jl")
import .Flow
G0 = Matrix{Float64}(I, n, n) ;
# build an A (the linear dynamics, i.e. Q and C matrices combined into a square matrix)
tmpzero = repeat([0.0], n^2);
A = reshape(tmpzero, (n,n));
pG = (n=n, p_Ds_v5=p_inputsG, A=A);
tspan = (0.0, 2.5);
prob_Gs_v5_condnums = DifferentialEquations.ODEProblem(Flow.calc_Gs_SSE_condnums!, G0, tspan, pG)
Gflow_to_25_condnums  = solve(prob_Gs_v5_condnums, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol = 1e-6, reltol = 1e-6)


# Call some R code to plot the results
@rput resfn
R"""
#R CODE HERE
res_table_df = read.table(resfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(res_table_df)

tmpcolors = "grey70"
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6,  col=tmpcolors, xlab="trad SSE", ylab="Gflow SSE", main=resfn, xlim=c(-16,0), ylim=c(-16,0))
segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")

# Histogram
diffs = res_table_df$Julia_sum_lq_GFv6 - res_table_df$Julia_sum_lq_nFv6
hist(diffs, breaks=100)

lmres = lm(res_table_df$Julia_sum_lq_GFv6~res_table_df$Julia_sum_lq_nFv6)
summary(lmres)
"""
@rget diffs;
a

"""
# R code
wd = "/GitHub/PhyBEARS.jl/notes/"
setwd(wd)

# ci_eq_i TFs, j and k correct, including Ds
resfn = "res_table_df_Ci_eq_i_GflowMap_v9.txt" # BAD
mean(diffs)




TF = abs(lmres$residuals) > 1
sum(TF)
TF = abs(lmres$residuals) > 0.75
printall(res_table_df[TF,])
hist(lmres$residuals)


diffs = res_table_df$Julia_sum_lq_GFv6 - res_table_df$Julia_sum_lq_nFv6
hist(diffs, breaks=100)

# Almost there except for Ci_sub_i -> Ci_eq_i
resfn = "res_table_df_Ci_eq_i_Gflow_v6a.txt" # LESS BAD
res_table_df = read.table(resfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(res_table_df)

tmpcolors = "grey"
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6, col=tmpcolors, xlab="trad SSE", ylab="Gflow SSE", main=resfn, xlim=c(-16,0), ylim=c(-16,0))
segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")
lmres = lm(res_table_df$Julia_sum_lq_GFv6~res_table_df$Julia_sum_lq_nFv6)
summary(lmres)
diffs = res_table_df$Julia_sum_lq_GFv6 - res_table_df$Julia_sum_lq_nFv6
hist(diffs, breaks=100)




# ci_eq_i TFs, j and k correct, accuracy 1e-12 
resfn = "res_table_df_Ci_eq_i_highres.txt" # BAD
res_table_df = read.table(resfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(res_table_df)

tmpcolors = "grey"
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6, pch=".", col=tmpcolors, xlab="trad SSE", ylab="Gflow SSE")

segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")
lmres = lm(res_table_df$Julia_sum_lq_GFv6~res_table_df$Julia_sum_lq_nFv6)
summary(lmres)
diffs = res_table_df$Julia_sum_lq_GFv6 - res_table_df$Julia_sum_lq_nFv6
hist(diffs, breaks=100)


TF = abs(lmres$residuals) > 1
sum(TF)
TF = abs(lmres$residuals) > 0.75
printall(res_table_df[TF,])

# switching j and k in cladogenesis SSE calc
# BAD
resfn = "res_table_df_Ci_eq_i2.txt" # BAD
res_table_df = read.table(resfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(res_table_df)

tmpcolors = "grey"
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6, pch=".", col=tmpcolors, xlab="trad SSE", ylab="Gflow SSE")

segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")
lmres = lm(res_table_df$Julia_sum_lq_GFv6~res_table_df$Julia_sum_lq_nFv6)
summary(lmres)
TF = abs(lmres$residuals) > 1
sum(TF)



# Bad - uses i ancstates instead of TF on i ancstates
resfn = "res_table_df_Ci_sub_i.txt" # BAD
res_table_df = read.table(resfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(res_table_df)

tmpcolors = "grey"
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6, pch=".", col=tmpcolors, xlab="trad SSE", ylab="Gflow SSE")

segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")


# BEST YET
# pretty good; errors ~1 lnL...
resfn = "res_table_df_Ci_eq_i.txt"
res_table_df = read.table(resfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(res_table_df)

tmpcolors = "grey"
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6, pch=".", col=tmpcolors, xlab="trad SSE", ylab="Gflow SSE")

segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")

lmres = lm(res_table_df$Julia_sum_lq_GFv6~res_table_df$Julia_sum_lq_nFv6)
summary(lmres)
TF = abs(lmres$residuals) > 1
sum(TF)




resfn = "res_table_df.txt"
res_table_df = read.table(resfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(res_table_df)

map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

map2color(0.5, pal=rainbow(nrow(res_table_df)), limits=c(0,1))

max_mu = max(c(res_table_df$mu_val1, res_table_df$mu_val2))
mean_mus = (2*res_table_df$mu_val1+res_table_df$mu_val2) / 3


levels = 1-(max_mu-mean_mus)/max_mu
tmpcolors = "grey"
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6, pch=".", col=tmpcolors, xlab="trad SSE", ylab="Gflow SSE")

segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")


levels = 1-(max_mu-mean_mus)/max_mu
tmpcolors = map2color(levels, pal=heat.colors(nrow(res_table_df)), limits=c(0,1))
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6, pch=".", col=tmpcolors)

segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")

max_j = max(c(res_table_df$j_val))
j_scaled = (max_j-res_table_df$j_val) / max_j

levels = 1-j_scaled
tmpcolors = map2color(levels, pal=rainbow(nrow(res_table_df)), limits=c(0,1))
plot(res_table_df$Julia_sum_lq_nFv6, res_table_df$Julia_sum_lq_GFv6, pch=".", col=tmpcolors)

segments(x0=-14, y0=-14, x1=0, y1=0, lty="dashed", col="black")



lmres = lm(res_table_df$Julia_sum_lq_GFv6~res_table_df$Julia_sum_lq_nFv6)
TF = abs(lmres$residuals) > 1
sum(TF)

res_table_df[TF,]

"""


# end


