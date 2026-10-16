
#######################################################
# Example inference on a simulated tree & geography dataset
# (from Wallis's ss8_sim_001)
# 2022-11-09
#
# We will compare inference under 
# * SSE likelihood calculator v7 (constant rates through time, very fast)
# * SSE likelihood calculator v12 (changing rates through time, very fast)
#   - Here, to start, we will only change the distances through time
#######################################################

using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using PhyBEARS
using DataFrames
using CSV

# Change the working directory as needed
wd = "~/GitHub/PhyBEARS.jl/ex/siminf_v12a/"
cd(expanduser(wd))

# This simulation has 50 living species
trfn = "tree.newick"
tr = readTopology(trfn)
trdf = prt(tr)
oldest_possible_age = 100.0

lgdata_fn = "rangedata.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = true
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
#bmo.type[bmo.rownames .== "j"] .= "free";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 0.034;
bmo.est[bmo.rownames .== "e"] .= 0.028;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.est[bmo.rownames .== "x"] .= 0.0;
bmo.est[bmo.rownames .== "xv"] .= 0.1;
bmo.max[bmo.rownames .== "xv"] .= 10.0;

bmo.est .= bmo_updater_v1(bmo);

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

#######################################################
# Read in and parse distances and area-of-areas
#######################################################
files.times_fn = "v12a_times.txt"
files.distances_fn = "v12a_distances.txt"
files.area_of_areas_fn = "v12a_area_of_areas.txt"

# Construct interpolators, times-at-which-to-interpolate QC
p = p_Ds_v5
interpolators = files_to_interpolators(files, setup.numareas, setup.states_list, setup.v_rows, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals, trdf; oldest_possible_age=100.0);

p_Es_v12 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo, interpolators=interpolators);

# Add Q, C interpolators
p_Es_v12 = p = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v12, p_Es_v12.interpolators.times_for_SSE_interpolators);

# Solve the Es
prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


#######################################################
# Maximum likelihood inference
#######################################################
bmo.type[bmo.rownames .== "x"] .= "free"
bmo.type[bmo.rownames .== "xv"] .= "free"
pars = bmo.est[bmo.type .== "free"]
parnames = bmo.rownames[bmo.type .== "free"]
func = x -> func_to_optimize_v12(x, parnames, inputs, p_Ds_v12; returnval="bgb_lnL", printlevel=1)
pars = [0.4, 0.001, -0.001, 4.0]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


using NLopt
opt = NLopt.Opt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2;
lower = bmo.min[bmo.type .== "free"];
upper = bmo.max[bmo.type .== "free"];
opt.lower_bounds = lower::Union{AbstractVector,Real};
opt.upper_bounds = upper::Union{AbstractVector,Real};
#opt.ftol_abs = 0.01 # tolerance on log-likelihood
#opt.ftol_rel = 0.01 # tolerance on log-likelihood
#opt.xtol_abs = 0.0001 # tolerance on parameters
#opt.xtol_rel = 0.001 # tolerance on parameters
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################




# Get the inputs & res:
pars = optx;
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars;
bmo_updater_v1!(inputs.bmo);
res = inputs.res;

# Solution, under best ML parameters
p_Ds_v5_updater_v1!(p_Ds_v12, inputs);
p_Es_v12 = TimeDep.construct_QC_interpolators(p_Ds_v12, p_Ds_v12.interpolators.times_for_SSE_interpolators);

# Solve the Es
prob_Es_v12 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, inputs.Es_tspan, p_Es_v12)
# This solution is an interpolator
sol_Es_v12 = solve(prob_Es_v12, inputs.solver_options.solver, save_everystep=inputs.solver_options.save_everystep, abstol=inputs.solver_options.abstol, reltol=inputs.solver_options.reltol);
p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

# Calculate the Ds, and final log-likelihood etc.
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)









#######################################################
# Calculate ancestral range probabilities under the PhyBEARS DEC model,
# ML parameters
#######################################################

# (updates the "res" results object)
PhyBEARS.Uppass.uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=false)

# Output the updated res to a series of .Rdata objects that can be read in R
#######################################################
# Ancestral states estimation and plotting
#######################################################
tmpres = deepcopy(inputs.res);
geogfn = lgdata_fn
lnLs_tuple = (total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, total_loglikelihood=bgb_lnL)
optim_result = build_optim_result(opt, optf, optx, ret)
juliaRes_to_Rdata(inputs.res, trdf, inputs, lnLs_tuple, optim_result, geogfn, trfn; outwd=getwd(), outfns=NaN)
tmpres_archive = deepcopy((inputs.res, inputs, lnLs_tuple, optim_result, geogfn, trfn));
"""
# Then run, in R:
"""
# (NOTE: NO COMMENTS are allowed *after* R commands, in the same line)
# $ -- This changes the Julia window to an R window
using RCall

# Load the R string WITHOUT dollar-sign symbols (using pipes instead), then replace
rstring = """
library(ape)
library(cladoRcpp)
library(diversitree)
library(BioGeoBEARS)
wd = "~/GitHub/PhyBEARS.jl/ex/siminf_v12a/"  # CHANGE THIS
setwd(wd)
sourceall("~/GitHub/PhyBEARS.jl/Rsrc/")
res = PhyBEARS_res_to_BGB_res(outfns=NaN)
resDECj = res  # CHANGE THIS
results_object = res

trfn = res|inputs|trfn				# The pipe characters will be converted to dollar signs
geogfn = res|inputs|geogfn		# The pipe characters will be converted to dollar signs
tr = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

max_range_size = res|inputs|max_range_size
include_null_range = res|inputs|include_null_range

pdffn = "example_inf_v12a_v4.jl.pdf"  # CHANGE THIS
pdf(pdffn, height=24, width=9)
analysis_titletxt ="example_inf_v12a_v4.jl"  # CHANGE THIS
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)
"""
rstring = replace(rstring, "|"=>"\$") # replacing the pipes
reval(rstring)
