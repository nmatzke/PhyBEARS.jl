
#######################################################
# Example inference on a example tree & geography dataset
#
# We will compare inference under 
# * SSE likelihood calculator v7 (constant rates through time, very fast)
# * SSE likelihood calculator v12 (changing rates through time, fast)
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
wd = "/GitHub/PhyBEARS.jl/data/"
cd(wd)

# Psychotria tree from Ree & Smith 2008
trfn = "Psychotria_tree.newick"
tr = readTopology(trfn)
trdf = prt(tr)
oldest_possible_age = 100.0

lgdata_fn = "Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = true
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0000001, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
#bmo.type[bmo.rownames .== "j"] .= "free";
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 0.034;
bmo.est[bmo.rownames .== "e"] .= 0.028;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0000001;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.est[bmo.rownames .== "x"] .= 0.0;
#bmo.est[bmo.rownames .== "xv"] .= 0.1;
bmo.max[bmo.rownames .== "xv"] .= 10.0;

bmo.est .= bmo_updater_v1(bmo);

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

#######################################################
# Read in and parse distances and area-of-areas
#######################################################
files.times_fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/v12a_times.txt"
files.distances_fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/v12a_distances.txt"
files.area_of_areas_fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/v12a_area_of_areas.txt"

files.times_fn = "Hawaii_KOMH_times.txt"
files.distances_fn = "Hawaii_KOMH_distances_changing_fractions_PhyBEARS.txt"
files.area_of_areas_fn = ""

# Construct interpolators, times-at-which-to-interpolate QC
p = p_Ds_v5;
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
bmo.type[bmo.rownames .== "xv"] .= "free"
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "free"
bmo.type[bmo.rownames .== "x"] .= "free"
pars = bmo.est[bmo.type .== "free"]
parnames = bmo.rownames[bmo.type .== "free"]
func = x -> func_to_optimize_v12(x, parnames, inputs, p_Ds_v12; returnval="bgb_lnL", printlevel=1)
pars = [0.04, 0.01, -0.1, 0.1, 0.34, 0.0]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


using NLopt
opt = NLopt.Opt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]
opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
#opt.ftol_abs = 0.0001 # tolerance on log-likelihood
#opt.ftol_rel = 0.01 # tolerance on log-likelihood
#opt.xtol_abs = 0.00001 # tolerance on parameters
#opt.xtol_rel = 0.001 # tolerance on parameters
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################
#opt.ftol_rel = 0.01 # tolerance on log-likelihood
#opt.xtol_abs = 0.0001 # tolerance on parameters
# (-16.91509708891156, [0.04786084799206467, 1.0e-12, 5.467412494208904], :SUCCESS)
# Unconstrained:
# d=0.04791,	e=0.0,	xv=5.47159,	Julia_sum_lq=-56.8738, rootstates_lnL=3.8986,	Julia_total_lnLs1=-52.9752, bgb_lnL=-15.0759
# (-15.075897880145341, [0.04791423273879585, 1.0e-12, 5.4715887419793905], :ROUNDOFF_LIMITED)
# 
# opt.ftol_abs = 0.001 # tolerance on log-likelihood
# opt.xtol_abs = 0.0001 # tolerance on parameters
# d=0.04793,	e=0.0,	xv=5.46762,	Julia_sum_lq=-77.7007, rootstates_lnL=3.7721,	Julia_total_lnLs1=-73.9286, bgb_lnL=-36.7417
# (-16.91509708891156, [0.04786084799206467, 1.0e-12, 5.467412494208904], :SUCCESS)
# d=0.07304,	e=1.0e-5,	xv=5.43854,	Julia_sum_lq=-60.3596, rootstates_lnL=3.5499,	Julia_total_lnLs1=-56.8098, bgb_lnL=-19.2144
# 16.077195154750136, [0.07303955699440073, 1.0e-12, 5.438560655327038], :SUCCESS)

# d=0.07303,	e=0.0,	xv=5.43862,	Julia_sum_lq=-57.9226, rootstates_lnL=3.4612,	Julia_total_lnLs1=-54.4613, bgb_lnL=-16.3902
# (-16.390191978494094, [0.07303402776871623, 1.0e-12, 5.438618806149209], :ROUNDOFF_LIMITED)

# d=0.07303,	e=0.0,	xv=5.43862,	Julia_sum_lq=-57.9226, rootstates_lnL=3.4612,	Julia_total_lnLs1=-54.4613, bgb_lnL=-16.3902
# (-16.390191978494094, [0.07303402776871623, 1.0e-12, 5.438618806149209], :ROUNDOFF_LIMITED)

# Unconstrained, 5 parameters:
# d=0.05993,	e=0.0108,	x=-0.23754,	xv=5.10425,	birthRate=0.3646,	deathRate=0.19304,	Julia_sum_lq=-56.6835, rootstates_lnL=6.7275,	Julia_total_lnLs1=-49.9559, bgb_lnL=-12.8284
# (-12.828352684984715, [0.059932699102433235, 0.010797111791669806, -0.23754106458148194, 5.104251203844556, 0.3645992824829854, 0.19303903497496325], :ROUNDOFF_LIMITED)



# Get the inputs & res:
pars = optx
#pars = [0.9747407112459348, 0.8, 0.11]
#pars = [100.0, 1.8, 0.11]
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
bmo_updater_v1!(inputs.bmo)
p_Ds_v5_updater_v1!(p_Ds_v12, inputs);
p_Es_v12 = TimeDep.construct_QC_interpolators(p_Ds_v12, p_Ds_v12.interpolators.times_for_SSE_interpolators);

# Solve the Es
prob_Es_v12 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, inputs.Es_tspan, p_Es_v12)
# This solution is an interpolator
sol_Es_v12 = solve(prob_Es_v12, inputs.solver_options.solver, save_everystep=inputs.solver_options.save_everystep, abstol=inputs.solver_options.abstol, reltol=inputs.solver_options.reltol);
p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

# Calculate the Ds, and final log-likelihood etc.
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

Rnames(res)
round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]; digits=3)
#  0.0  0.0  0.0  0.0  0.0  0.002  0.003  0.947  0.0  0.003  0.0  0.001  0.024  0.02  0.0  0.0
