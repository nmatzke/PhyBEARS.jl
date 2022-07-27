#######################################################
# Test the lnLs from trad SSE v6 (uses a half matrix)
# and Gflow v7 (uses a half matrix)
# 2022-03-24: WORKS!
#######################################################

using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using CSV
using RCall
using PhyBEARS.Parsers
using PhyBEARS.Gmaps

print("\n...finished loading packages.")

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests_trad_vs_Gflow_v5_Psychotria_linreg.jl")
"""

@testset "runtests_trad_vs_Gflow_v5_Psychotria_linreg.jl" begin



trfn = "/GitHub/PhyBEARS.jl/notes/amph_shl_dates_frogs2.tre"
tr = readTopology(trfn)
trdf = prt(tr);

lgdata_fn = "/GitHub/PhyBEARS.jl/notes/frog_ecoregions_reduced.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.est[bmo.rownames .== "birthRate"] .= 0.32881638319078066
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
bmo.est[:] = bmo_updater_v1(bmo);
numareas = 8
n = 2^8            # 4 areas, 16 states

# Set up the model
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;

p_Ds_v5 = inputs.p_Ds_v5
p_Ds_v5_updater_v1!(p_Ds_v5, inputs)

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
d_val=in_params.d_val; e_val=in_params.e_val; j_val=in_params.j_val; lambda_val=in_params.birthRate; mu_val1=in_params.deathRate
#for d_val in d_vals
#for e_val in e_vals
#for j_val in j_vals
#for lambda_val in lambda_vals
#for mu_val1 in mu_vals1

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
p_Ds_v5 = inputs.p_Ds_v5  # 2022-03-24_THIS LINE IS CRUCIAL -- SUPER-SLOW WITHOUT IT (scoping issue)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs)


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


tmprow = [d_val, e_val, j_val, lambda_val, mu_val1, Julia_sum_lq_nFv6, rootstates_lnL_nFv6, Julia_total_lnLs1_nFv6, bgb_lnl_nFv6, Julia_sum_lq_GFv6a, rootstates_lnL_GFv6a, Julia_total_lnLs1_GFv6a, bgb_lnl_GFv6a]

res_table[counter,:] .= tmprow


end
end
end
end
end

tmprow_names = ["d_val", "e_val", "j_val", "lambda_val", "mu_val1", "Julia_sum_lq_nFv6", "rootstates_lnL_nFv6", "Julia_total_lnLs1_nFv6", "bgb_lnl_nFv6", "Julia_sum_lq_GFv6a", "rootstates_lnL_GFv6a", "Julia_total_lnLs1_GFv6a", "bgb_lnl_GFv6a"]

res_table_df = DataFrame(res_table, :auto);

# Add names
rename!(res_table_df, Symbol.(tmprow_names) ) 

# Filter for NaNs
res_table_df = (filter(row->!any(isnan(x) for x in row), res_table_df))

resfn = "/GitHub/PhyBEARS.jl/notes/res_table_df_halfMatrix_tradSSEv6_vs_GflowArrays_v7.txt"
CSV.write(resfn, res_table_df; delim="\t")

#######################################################
# Call some R code to plot the results
#######################################################
resfn2 = getfn(resfn)
@rput resfn
@rput resfn2
R"""
#R CODE HERE
# resfn = "/GitHub/PhyBEARS.jl/notes/res_table_df_halfMatrix_tradSSEv6_vs_GflowV7.txt"
res_table_df = read.table(resfn, header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(res_table_df)

tmpcolors = "grey70"
xvals = res_table_df$Julia_total_lnLs1_nFv6
yvals = res_table_df$Julia_total_lnLs1_GFv6
xmin = min(c(min(xvals), min(yvals)))
xmax = max(c(max(xvals), max(yvals)))
ymin = xmin
ymax = xmax


plot(xvals, yvals, col=tmpcolors, xlab="trad SSE", ylab="GflowArrays SSE", main=paste0("TradSSE_julia_sumlq_vs_BGB\n", resfn2), xlim=c(xmin,xmax), ylim=c(ymin,ymax))
segments(x0=xmin, y0=ymin, x1=xmax, y1=ymax, lty="dashed", col="black")

# Histogram
diffs2 = xvals - yvals
titletxt = paste0("Diff betw. trad SSE sumLq & GflowArrays_v7:\n",resfn2)
hist(diffs2, breaks=100, main=titletxt)

lmres = lm(yvals~xvals)
summary(lmres)

intercept = lmres$coefficients[1]
slope = lmres$coefficients[2]
rsq = summary(lmres)$r.squared

"""
@rget diffs2;
@rget slope;
@rget intercept;
@rget rsq;

print("\nTesting if lnLs matrix across a range of parameters between trad_SSE_v6 and GflowArrays_v7, with half matrix:\n")

@test abs(mean(diffs2)) > 0.01
@test round(slope, digits=2) == 1.0
#@test round(intercept, digits=2) != 0.0
@test round(rsq, digits=2) > 0.95



end
