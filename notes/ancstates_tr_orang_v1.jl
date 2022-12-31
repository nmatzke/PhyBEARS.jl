
#######################################################
# Example inference: Imagine NZ sunk, recolonized from elsewhere
#######################################################

using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using DataFrames
using CSV

using PhyBEARS
#using PhyBEARS.Parsers


# Change the working directory as needed
wd = "/GitHub/PhyBEARS.jl/notes/"
cd(wd)

# Simple tree
tr = readTopology("(((chimp:1,human:1):1,gorilla:2):1,orang:3);")

trdf = prt(tr);
oldest_possible_age = 4.0

lgdata_fn = "ancstates_tr_orang.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = false
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
#bmo.type[bmo.rownames .== "j"] .= "free";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 0.0;
bmo.est[bmo.rownames .== "e"] .= 0.0;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.min[bmo.rownames .== "u"] .= 0.0;
bmo.max[bmo.rownames .== "u"] .= 0.0;

bmo.type[bmo.rownames .== "j"] .= "fixed";
bmo.type[bmo.rownames .== "u"] .= "fixed";
bmo.type[bmo.rownames .== "birthRate"] .= "free";
bmo.type[bmo.rownames .== "deathRate"] .= "fixed";



# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=include_null_range, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

bmo.est[:] = bmo_updater_v2(bmo, inputs.setup.bmo_rows);

inputs.setup.txt_states_list


# Solve the Es
p_Es_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo);

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Es_v7(0.0)
sol_Es_v7(1.0)
sol_Es_v7(1.5)

p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


#######################################################
# Maximum likelihood inference
#######################################################
inputs.bmo.type[inputs.bmo.rownames .== "j"] .= "free"
inputs.bmo.type[inputs.bmo.rownames .== "birthRate"] .= "free"
inputs.bmo.type[inputs.bmo.rownames .== "deathRate"] .= "birthRate"

inputs.bmo.type[inputs.bmo.rownames .== "u"] .= "fixed"
inputs.bmo.est[inputs.bmo.rownames .== "u"] .= -1.0
inputs.bmo.init[inputs.bmo.rownames .== "u"] .= -1.0

pars = deepcopy(inputs.bmo.est[inputs.bmo.type .== "free"])
parnames = inputs.bmo.rownames[inputs.bmo.type .== "free"]
func = x -> func_to_optimize_v12(x, parnames, inputs, p_Ds_v12; returnval="lnL", printlevel=1)
#pars = [0.04, 0.001, 0.0001, 0.1, inputs.bmo.est[bmo.rownames .== "birthRate"][1], 0.0]



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
#opt.ftol_abs = 0.0001 # tolerance on log-likelihood
#opt.ftol_rel = 0.01 # tolerance on log-likelihood
#opt.xtol_abs = 0.00001 # tolerance on parameters
#opt.xtol_rel = 0.001 # tolerance on parameters
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################


# Get the inputs & res:
pars = optx;

# Give the simulation a substantial death rate
func(pars)
#pars[parnames .== "deathRate"] .= 0.5*pars[parnames .== "birthRate"]
#pars[parnames .== "u"] .= -1.0
func(pars)

inputs.bmo.est[inputs.bmo.type .== "free"] .= pars;
inputs.bmo.est[bmo.rownames .== "birthRate"] = inputs.bmo.est[bmo.rownames .== "birthRate"] / 5
inputs.bmo.est[:] = bmo_updater_v2(bmo, inputs.setup.bmo_rows)
res = inputs.res;

# Solution, under best ML parameters
p_Ds_v5_updater_v1!(p_Ds_v12, inputs);
p_Es_v12 = TimeDep.construct_QC_interpolators(p_Ds_v12, p_Ds_v12.interpolators.times_for_SSE_interpolators);

# Solve the Es
prob_Es_v12 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, inputs.Es_tspan, p_Es_v12)
# This solution is an interpolator
sol_Es_v12 = solve(prob_Es_v12, inputs.solver_options.solver, save_everystep=inputs.solver_options.save_everystep, abstol=inputs.solver_options.abstol, reltol=inputs.solver_options.reltol);
p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

Rnames(p_Ds_v12.interpolators)

p_Ds_v12.interpolators.area_of_areas_interpolator(20.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(21.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(22.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(23.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(24.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(25.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(26.0)



p_Ds_v12.interpolators.mu_vals_interpolator(0.0)
p_Ds_v12.interpolators.mu_vals_interpolator(1.0)
p_Ds_v12.interpolators.mu_vals_interpolator(20.0)
p_Ds_v12.interpolators.mu_vals_interpolator(21.0)
p_Ds_v12.interpolators.mu_vals_interpolator(22.0)
p_Ds_v12.interpolators.mu_vals_interpolator(23.0)
p_Ds_v12.interpolators.mu_vals_interpolator(23.5)
p_Ds_v12.interpolators.mu_vals_interpolator(24.0)
p_Ds_v12.interpolators.mu_vals_interpolator(60.0)


# Calculate the Ds, and final log-likelihood etc.
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

Rnames(res)
round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]; digits=3)

# 0.06
# 0.613
# 0.327


# ancestral_range_estimation
# This term is preferable to e.g. "ancestral area reconstruction"

Rnames(res)

rootnode = inputs.res.root_nodeIndex

lnode = trdf[rootnode,"leftNodeIndex"]
rnode = trdf[rootnode,"rightNodeIndex"]

# ACE for left descendant
nodenum = rootnode
nodelikes = res.normlikes_at_each_nodeIndex_branchTop[nodenum]



uppass_edgematrix = res.uppass_edgematrix

include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
current_nodeIndex = 6
x = nodeOp_Cmat_uppass_v12!(res, current_nodeIndex, trdf, p_Ds_v12, solver_options)

solver_options.abstol = 1.0e-9
solver_options.reltol = 1.0e-9
uppass_ancstates_v12(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=true)

res.uppass_probs_at_each_nodeIndex_branchBot
res.anc_estimates_at_each_nodeIndex_branchBot
res.uppass_probs_at_each_nodeIndex_branchTop
res.anc_estimates_at_each_nodeIndex_branchTop
res.fixNodesMult_at_each_nodeIndex_branchTop
res.fixNodesMult_at_each_nodeIndex_branchBot


tspan

uppass_Ds_v12 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v12_simd_sums, deepcopy(u0), tspan, p_Ds_v12)

	sol_Ds = solve(prob_Ds_v12, solver_options.solver, dense=false, save_start=false, save_end=true, save_everystep=false, abstol=solver_options.abstol, reltol=solver_options.reltol)


# Install modified "castor" package in R
# install.packages(pkgs="/GitHub/PhyBEARS.jl/simulator/castor_1.7.2.000004.tar.gz", lib="/Library/Frameworks/R.framework/Resources/library/", repos=NULL, type="source")

# Write model out to text files that can be read in to simulator
geog_interpolator_times = parse_times_fn(files.times_fn)
timepoints = sort(unique(vcat(seq(0.0, maximum(geog_interpolator_times), 1.0), geog_interpolator_times)))
# (the best way to do this is to do simulations for a fixed period of time; the number of taxa
#  will vary, but have an average)
outfns = model_to_text_v12(p_Ds_v12, timepoints; prefix="")


Rcode = """
library(cladoRcpp)
library(BioGeoBEARS)
library(ape)
library(castor)

# for: reorder_castor_sim_to_default_ape_node_order(simulation)
source("/GitHub/PhyBEARS.jl/Rsrc/castor_helpers.R")

wd = "/GitHub/PhyBEARS.jl/sims/sunkNZ_v1/"
setwd(wd)
simfns = c("setup_df.txt",
"timepoints.txt", 
"mu_vals_by_t.txt", 
"Qvals_by_t.txt",
"Crates_by_t.txt",
"Qarray.txt",
"Carray.txt",
"area_names.txt",
"states_list.R")


simulation2 = simulate_tdsse2_for_timeperiod(wd, start_state=2, max_simulation_time=100.0, min_tips=50, max_tips=500, simfns=default_simfns(), seedval=543221, max_rate=10.0, numtries=250)
get_root_age(simulation2$tree)
get_root_age(simulation2$living_tree)
"""
