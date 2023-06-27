
using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames						# for DataFrame()
using DelimitedFiles				# for readdlm()
using NLopt									# seems to be the best gradient-free, box-constrained								

# List each PhyBEARS code file prefix here
using PhyloBits
using PhyloBits.TrUtils			# for e.g. numstxt_to_df()
using PhyloBits.TreeTable
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.ModelLikes # e.g. setup_DEC_SSE2
using PhyBEARS.Uppass

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/simulator/examples/ex3/")
include("/Users/nickm/GitHub/PhyBEARS.jl/simulator/examples/ex3/01_DECjBD_inf_on_Psychotria_v1.jl")
#
Check against:
cd /GitHub/PhyBEARS.jl/ex/cinthy/
Rscript /GitHub/PhyBEARS.jl/ex/cinthy/epic_M0_ML_BSM_v1.R
"""

wd = "/GitHub/PhyBEARS.jl/simulator/examples/ex3/"
cd(wd)
getwd()

# Load files
geogfn = lgdata_fn = "geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Epacridoideae tree
trfn = "tree.newick"
tr = PhyloBits.PNreadwrite.readTopology(trfn)
trdf = prt(tr)
root_age = maximum(trdf[!, :node_age])

ML_yule_birthRate(tr)
ML_yule_birthRate_wRoot(tr)


# DEC model on Hawaiian Epacridoideae
bmo = construct_BioGeoBEARS_model_object()
birthRate = ML_yule_birthRate(tr)
bmo.est[bmo.rownames .== "birthRate"] .= birthRate
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.0001
bmo.est[bmo.rownames .== "e"] .= 1.000000e-12 
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.11
bmo.type[bmo.rownames .== "j"] .= "free"

numareas = Rncol(geog_df)-1
max_range_size = numareas
oldest_possible_age = 100.0

inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;

numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])


files.times_fn = "times.txt"
files.distances_fn = "distances.txt"
files.area_of_areas_fn = "area_of_areas.txt"

# Construct interpolators, times-at-which-to-interpolate QC
p = p_Es_v5;
interpolators = files_to_interpolators(files, setup.numareas, setup.states_list, setup.v_rows, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals, trdf; oldest_possible_age=oldest_possible_age);

ts = [0.0, 1.0, 3.0]
interpolators.distances_interpolator(ts)
interpolators.area_of_areas_interpolator(ts)

p_Es_v12 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo, interpolators=interpolators);

# Add Q, C interpolators
p_Es_v12 = p = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v12, p_Es_v12.interpolators.times_for_SSE_interpolators);


# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")

# v7 likelihood (no interpolation / changing distances etc.) 
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v7);
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


# v12 likelihood (allows interpolation / changing distances etc.) 
prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)



#######################################################
# ML inference on DEC+J
#######################################################
parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]
pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v12; returnval="lnL", printlevel=1)
#pars = [0.9, 0.9, 0.9]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)

using NLopt
pars = [0.9, 0.9, 0.9]
func(pars)
opt = NLopt.Opt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2
opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
opt.lower_bounds
opt.upper_bounds
#opt.ftol_abs = 0.001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################
# For saving to dataframe
optim_result = build_optim_result(opt, optf, optx, ret)

# Get the inputs & res:
pars = optx
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
bmo_updater_v1!(inputs.bmo)
p_Ds_v5_updater_v1!(p_Es_v5, inputs);
inputs.solver_options.save_everystep=true	# WORKS!! Can make a difference EVEN ON THE Ds!!

# Solve the Es
p_Es_v12 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo, interpolators=interpolators);

# Add Q, C interpolators
p_Es_v12 = p = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v12, p_Es_v12.interpolators.times_for_SSE_interpolators);

# v12 likelihood (allows interpolation / changing distances etc.) 
prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

#######################################################
# Ancestral states estimation and plotting
#######################################################
uppass_ancstates_v12!(inputs.res, trdf, p_Ds_v12, solver_options);
resDECj_Yule = deepcopy(inputs.res);
geogfn = lgdata_fn
lnLs_tuple = (total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, total_loglikelihood=bgb_lnL)
optim_result = build_optim_result(opt, optf, optx, ret)
juliaRes_to_Rdata(inputs.res, trdf, inputs, lnLs_tuple, optim_result, geogfn, trfn; outfns=NaN)
resDECj_Yule_archive = deepcopy((inputs.res, inputs, lnLs_tuple, optim_result, geogfn, trfn));
"""
# Then run, in R:
"""
# (NOTE: NO COMMENTS are allowed *after* R commands, in the same line)
# $ -- This changes the Julia window to an R window
using RCall

# Plot with BioGeoBEARS by evaluating this R string with reval()
# NOTE: Julia turns $wd into the full pre-specificed working directory
# NOTE: The replace() function changes the "|" characters "$" for evaluation in R
# NOTE: reval() runs the modified string
rstring = """
library(ape)
library(cladoRcpp)
library(BioGeoBEARS)
wd = "$wd"
setwd(wd)
sourceall("/GitHub/PhyBEARS.jl/Rsrc/")
res = PhyBEARS_res_to_BGB_res(outfns=NaN)
resDECj = res
results_object = res

trfn = res|inputs|trfn
geogfn = res|inputs|geogfn
tr = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

max_range_size = res|inputs|max_range_size
include_null_range = res|inputs|include_null_range

pdffn = "Psychotria_DEC+J+Yule_M0_unconstrained_v1.pdf"
pdf(pdffn, height=6, width=6)
analysis_titletxt ="PhyBEARS DEC+J+Yule on Psychotria M0_unconstrained"
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
rstring = replace(rstring, "|"=>"\$")
reval(rstring)







#######################################################
# Output this model setup to files for a castor simulator
#######################################################
# Parsers.jl
timepoints = [0.0, 5.0, 100.0]
model_to_text_v12(p_Ds_v12, timepoints; prefix="")

# Created these files:
simfns = ["setup_df.txt",
"timepoints.txt",
"mu_vals_by_t.txt",
"Qvals_by_t.txt",
"Crates_by_t.txt",
"Qarray.txt",
"Carray.txt",
"area_names.txt",
"states_list.R"]

#######################################################
# Run castor simulator in R
#######################################################
using RCall

# NOTE: Julia turns $wd into the full pre-specificed working directory
# NOTE: The replace() function changes the "|" characters "$" for evaluation in R
# NOTE: reval() runs the modified string
rstring = """
library(ape)
library(cladoRcpp)
library(BioGeoBEARS)
library(castor)
sourceall("/GitHub/PhyBEARS.jl/Rsrc/")
wd = "$wd"
setwd(wd)

start_state = 2 # number of the starting state
max_simulation_time = 15.0 # Set to 0 if the user doesn't want to set a max simulation time
min_tips = 50
max_tips=100
simfns=default_simfns()
seedval = 54321
max_rate=10.0
numtries = 1000

start_state=2; max_simulation_time=15; min_tips=50; max_tips=100; simfns=default_simfns(); seedval=54321; max_rate=10.0; numtries=1000

simulation2 = simulate_tdsse2_for_timeperiod(wd, start_state, max_simulation_time, min_tips, max_tips, simfns, seedval, max_rate, numtries)

simulation = remove_last_tip_from_simulation(simulation2)
write_out_original_castor_simfiles(simulation, wd, fn="rawsim")
area_names = readLines("area_names.txt")
source("states_list.R")
write_out_reordered_castor_simfiles(simulation, wd, area_names, states_list, tip_age_tolerance=1e-6

"""
rstring = replace(rstring, "|"=>"\$")
reval(rstring)



#######################################################
# Look at all the files!
#######################################################
opd()


