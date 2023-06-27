
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

# Load files
geogfn = lgdata_fn = "geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Epacridoideae tree
trfn = "tree.newick"
tr = PhyloBits.PNreadwrite.readTopology(trfn)

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

inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;

numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, sol_Es_v5=sol_Es_v5);

# Do downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

#######################################################
# ML inference on DEC+J
#######################################################
parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]
pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL", printlevel=1)
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
# SEE runtests_ClaSSE_tree_n13_DECj_WORKS.jl
# save_everystep_EQ_false_CAN_MATTER_EVEN_ON_THE_Ds
#inputs.solver_options.save_everystep=false # CAN PRODUCE A -20.9 vs. -20.6 difference!
inputs.solver_options.save_everystep=true	# WORKS!! Can make a difference EVEN ON THE Ds!!

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


#######################################################
# Ancestral states estimation and plotting
#######################################################
uppass_ancstates_v7!(inputs.res, trdf, p_Ds_v7, solver_options);
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

# Load the R string WITHOUT "S" symbols, then replace
rstring = """
library(ape)
library(cladoRcpp)
library(BioGeoBEARS)
wd = "/GitHub/PhyBEARS.jl/simulator/examples/ex3/"
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
pdf(pdffn, height=18, width=6)
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
# ML inference on DEC+J + birth-death
#######################################################
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= 0.1
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.1
bmo.est[:] = bmo_updater_v1(bmo, inputs.setup.bmo_rows) # works


root_age_mult=1.5; max_range_size=NaN; include_null_range=false; max_range_size=NaN
max_range_size = NaN # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;



parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo); # works
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);  # WORKS 2022-03-10

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, sol_Es_v5=sol_Es_v5);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL")
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)

pars = [0.9, 0.9, 0.1, 0.3, 0.2]
func(pars)
func2(pars, [])

#######################################################
# Best optimizer so far - 2022-03-15
#######################################################
using NLopt
opt = NLopt.Opt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2
opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
opt.lower_bounds
opt.upper_bounds
#opt.ftol_abs = 0.00001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################

# Get the inputs & res:
pars = optx

# 2023-02-08_archived pars
optx = [1.0e-12, 1.0e-12, 0.10999899813382318, 0.34816245870282203, 0.0]
pars = optx

inputs.bmo.est[inputs.bmo.type .== "free"] .= optx
bmo_updater_v1!(inputs.bmo);
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

func(pars)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);


# inputs.solver_options.save_everystep = false
#	inputs.solver_options.saveat = nodetimes(trdf)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v6!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(inputs.res; trdf=trdf, p_Ds_v5=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


# Ancestral states estimation
uppass_ancstates_v7!(inputs.res, trdf, p_Ds_v7, solver_options);
resDECj_BD = deepcopy(inputs.res);

Rnames(resDECj_BD)
vfft(resDECj_BD.likes_at_each_nodeIndex_branchTop)
vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchTop)
vfft(resDECj_BD.likes_at_each_nodeIndex_branchBot)
vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchBot)
vfft(resDECj_BD.anc_estimates_at_each_nodeIndex_branchBot)
vfft(resDECj_BD.anc_estimates_at_each_nodeIndex_branchTop)

rowSums(vfft(resDECj_BD.likes_at_each_nodeIndex_branchTop))
rowSums(vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchTop))
rowSums(vfft(resDECj_BD.likes_at_each_nodeIndex_branchBot))
rowSums(vfft(resDECj_BD.normlikes_at_each_nodeIndex_branchBot))
rowSums(vfft(resDECj_BD.anc_estimates_at_each_nodeIndex_branchBot))
rowSums(vfft(resDECj_BD.anc_estimates_at_each_nodeIndex_branchTop))


#######################################################
# Ancestral states estimation and plotting
#######################################################
uppass_ancstates_v7!(inputs.res, trdf, p_Ds_v7, solver_options);
resDECj_BD = deepcopy(inputs.res);
geogfn = lgdata_fn
lnLs_tuple = (total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, total_loglikelihood=bgb_lnL)
optim_result = build_optim_result(opt, optf, optx, ret)
juliaRes_to_Rdata(inputs.res, trdf, inputs, lnLs_tuple, optim_result, geogfn, trfn; outfns=NaN)
resDECj_BD_archive = deepcopy((inputs.res, inputs, lnLs_tuple, optim_result, geogfn, trfn))
"""
# Then run, in R:
"""
# (NOTE: NO COMMENTS are allowed *after* R commands, in the same line)
using RCall
# Load the R string WITHOUT "S" symbols, then replace
rstring = """
library(ape)
library(cladoRcpp)
library(BioGeoBEARS)
wd = "/GitHub/PhyBEARS.jl/ex/cinthy/"
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

pdffn = "Epic_DEC+J+BD_M0_unconstrained_v1.pdf"
pdf(pdffn, height=18, width=6)
analysis_titletxt ="PhyBEARS DEC+J+BD on Epic M0_unconstrained"
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

