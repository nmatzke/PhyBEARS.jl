#######################################################
# Tutorial: On an example simulated dataset
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
using PhyBEARS.Uppass


# Change the working directory as needed
wd = "/Users/nmat471/HD/GitHub/PhyBEARS.jl/ex/tutorials/4sims/ss8_sim_001/"
cd(wd)

# This simulation has 148 living species
#trfn = "tree.newick"
trfn = "tree2.newick"
tr = readTopology(trfn)
trdf = prt(tr);
oldest_possible_age = 100.0

#lgdata_fn = "rangedata.data"
#lgdata_fn = "rangedata_noNull.data"
lgdata_fn = "rangedata_noSp50.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = true
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
bmo.est[bmo.rownames .== "d"] .= 0.01;
bmo.est[bmo.rownames .== "e"] .= 0.01;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.min[bmo.rownames .== "u"] .= -2.5;
bmo.max[bmo.rownames .== "u"] .= 0.0;

bmo.type[bmo.rownames .== "j"] .= "fixed";
bmo.type[bmo.rownames .== "u"] .= "fixed";
bmo.type[bmo.rownames .== "birthRate"] .= "free";
bmo.type[bmo.rownames .== "deathRate"] .= "fixed";



# Set up the model
testing_can_ignore="""
using PhyBEARS.ModelLikes
root_age_mult=1.5
max_range_size=NaN
include_null_range=include_null_range
bmo=bmo
allow_null_cladogenesis=true
manual_states_list=NaN
area_names=LETTERS(1:numareas)
fossils_older_than=1e-5

#predeclare_array_length=Integer(min(length(states_list)*length(states_list)*round((length(states_list)/2)), 10000000))
min_precision=1.0e-9
"""
# ADDING OPTION: allow_null_cladogenesis=true
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=include_null_range, bmo=bmo, allow_null_cladogenesis=true);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

prtCi(inputs)

bmo.est[:] = bmo_updater_v2(bmo, inputs.setup.bmo_rows);

inputs.setup.txt_states_list


p_Es_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo);

prtQp(p_Es_v7)
prtCp(p_Es_v7)



# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

t = [0.0, 1.0, 2.0]
sol_Es_v7.(t)

# Construct parameter object for calculating the Ds
p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

#######################################################
# NJM gets:
#######################################################
# (0.135, 13, -221.9734049941743, -5.551995962419084, -227.5254009565934, -118.49721930590728)


#######################################################
# Run optimization
#######################################################
pars = bmo.est[bmo.type .== "free"]
parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL", printlevel=1)
#pars = [0.01, 0.01, 0.01]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


#######################################################
# Best optimizer so far - 2022-03-15
#######################################################
using NLopt
#pars = [0.9, 0.9, 0.9]
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




# Get the inputs & res:
pars = optx
#pars = [0.9747407112459348, 0.8, 0.11]
#pars = [100.0, 1.8, 0.11]
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
bmo_updater_v1!(inputs.bmo)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);
res = inputs.res

prtQp(p_Ds_v5)
prtCp(p_Ds_v5)

# results object before uppass
rn(res)
res.uppass_probs_at_each_nodeIndex_branchBot
res.uppass_probs_at_each_nodeIndex_branchTop

res.anc_estimates_at_each_nodeIndex_branchBot
res.anc_estimates_at_each_nodeIndex_branchTop

uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

# results object after uppass
rn(res)
res.uppass_probs_at_each_nodeIndex_branchBot
res.uppass_probs_at_each_nodeIndex_branchTop

res.anc_estimates_at_each_nodeIndex_branchBot
res.anc_estimates_at_each_nodeIndex_branchTop


# The res tables are hard to read, the vfft function helps
# by putting into DataFrame
vfft(res.anc_estimates_at_each_nodeIndex_branchBot)
vfft(res.anc_estimates_at_each_nodeIndex_branchTop)

# View the ancestral range probabilities in R's default node order
R_node_order = sort(trdf, :Rnodenums).nodeIndex

# Tree table is reordered
trdf[R_node_order,:]

# Ancestral state probs are reorderd
vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_node_order])
vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_node_order])

