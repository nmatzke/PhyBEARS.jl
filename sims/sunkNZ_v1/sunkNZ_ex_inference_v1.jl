
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
using DataFrames
using CSV

using PhyBEARS
#using PhyBEARS.Parsers


# Change the working directory as needed
wd = "/GitHub/PhyBEARS.jl/sims/sunkNZ_v1/"
cd(wd)

# This simulation has 148 living species
trfn = "living_tree_noNodeLabels.newick"
tr = readTopology(trfn)
trdf = prt(tr);
oldest_possible_age = 125.0

lgdata_fn = "geog_living.data"
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
bmo.est[bmo.rownames .== "deathRate"] .= 0.9*ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "d"] .= 0.034;
bmo.est[bmo.rownames .== "e"] .= 0.028;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.1;
bmo.est[bmo.rownames .== "u"] .= -1.0;
bmo.min[bmo.rownames .== "u"] .= -2.5;
bmo.max[bmo.rownames .== "u"] .= 0.0;

bmo.type[bmo.rownames .== "j"] .= "free";
bmo.type[bmo.rownames .== "u"] .= "free";
bmo.type[bmo.rownames .== "birthRate"] .= "free";
bmo.type[bmo.rownames .== "deathRate"] .= "free";
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "free"


bmo.est .= bmo_updater_v1(bmo);

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=include_null_range, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

#######################################################
# Read in and parse distances and area-of-areas
#######################################################
files.times_fn = "sunkNZ_times.txt"
files.distances_fn = "sunkNZ_distances.txt"
files.area_of_areas_fn = "sunkNZ_area_of_areas.txt"

# Construct interpolators, times-at-which-to-interpolate QC
p = p_Ds_v5;
interpolators = files_to_interpolators(files, setup.numareas, setup.states_list, setup.v_rows, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals, trdf; oldest_possible_age=oldest_possible_age);

interpolators.area_of_areas_interpolator

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
inputs.bmo.type[inputs.bmo.rownames .== "j"] .= "free"
inputs.bmo.type[inputs.bmo.rownames .== "birthRate"] .= "free"
inputs.bmo.type[inputs.bmo.rownames .== "deathRate"] .= "free"
pars = deepcopy(inputs.bmo.est[inputs.bmo.type .== "free"])
parnames = inputs.bmo.rownames[inputs.bmo.type .== "free"]
func = x -> func_to_optimize_v12(x, parnames, inputs, p_Ds_v12; returnval="lnL", printlevel=1)
#pars = [0.04, 0.001, 0.0001, 0.1, inputs.bmo.estinputs.[bmo.rownames .== "birthRate"][1], 0.0]

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
# d=0.09035,	e=0.00116,	x=-0.53656,	xv=7.41251,	birthRate=0.25956,	deathRate=0.21085,	Julia_sum_lq=-175.6323, rootstates_lnL=13.6898,	Julia_total_lnLs1=-161.9425, bgb_lnL=-44.7349
# (-44.73494901514184, [0.09035419093152777, 0.001159414738073889, -0.5365638714798167, 7.41251181819682, 0.2595586392951748, 0.21085253647485827], :ROUNDOFF_LIMITED)
# d=0.09035,	e=0.00116,	x=-0.53656,	xv=7.41251,	birthRate=0.25956,	deathRate=0.21085,	Julia_sum_lq=-175.6323, rootstates_lnL=13.6898,	Julia_total_lnLs1=-161.9425, bgb_lnL=-44.7349
# (-44.73494901514184, [0.09035419093152777, 0.001159414738073889, -0.5365638714798167, 7.41251181819682, 0.2595586392951748, 0.21085253647485827], :ROUNDOFF_LIMITED)



# Get the inputs & res:
pars = optx;

# Give the simulation a substantial death rate
func(pars)
pars[parnames .== "deathRate"] .= 0.5*pars[parnames .== "birthRate"]
pars[parnames .== "u"] .= -1.0
func(pars)

inputs.bmo.est[inputs.bmo.type .== "free"] .= pars;
inputs.bmo.est[bmo.rownames .== "birthRate"] = inputs.bmo.est[bmo.rownames .== "birthRate"] / 5
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

Rnames(p_Ds_v12.interpolators)
p_Ds_v12.interpolators.mu_vals_interpolator(0.0)
p_Ds_v12.interpolators.mu_vals_interpolator(1.0)
p_Ds_v12.interpolators.mu_vals_interpolator(20.0)
p_Ds_v12.interpolators.mu_vals_interpolator(21.0)
p_Ds_v12.interpolators.mu_vals_interpolator(22.0)
p_Ds_v12.interpolators.mu_vals_interpolator(23.0)
p_Ds_v12.interpolators.mu_vals_interpolator(23.5)
p_Ds_v12.interpolators.mu_vals_interpolator(24.0)
p_Ds_v12.interpolators.mu_vals_interpolator(60.0)

p_Ds_v12.interpolators.area_of_areas_interpolator(20.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(21.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(22.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(23.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(24.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(25.0)
p_Ds_v12.interpolators.area_of_areas_interpolator(26.0)


# Calculate the Ds, and final log-likelihood etc.
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

Rnames(res)
round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]; digits=3)

# 0.0
#  0.0
#  0.0
#  0.0
#  0.296
#  0.028
#  0.292
#  0.384



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
