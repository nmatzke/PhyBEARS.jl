#######################################################
# This script shows:
#
# * Inference on a simulated tree / geographic data file
#
#######################################################

"""
cd /GitHub/PhyBEARS.jl/ex/siminf1/ss8_simulation_001
julia
"""

using Test
using DataFrames
using Dates									# for e.g. Dates.now(), DateTime
#using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames
using StatsBase
#using Optim                 # for e.g. LBFGS Maximum Likelihood optimization,optimize
# Optim really sucks, try LBFGSB: https://github.com/JuliaNLSolvers/Optim.jl/issues/953
# LBFGSB also sucks: https://github.com/Gnimuc/LBFGSB.jl
using NLopt									# seems to be the best gradient-free, box-constrained								

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DataFrames  # for DataFrame
using DifferentialEquations
using OrdinaryDiffEq, Sundials, DiffEqDevTools, ODEInterfaceDiffEq, ODE, LSODA

# List each PhyBEARS code file prefix here
# using PhyBEARS.BGExample
# using PhyloBits.PNreadwrite  # for readTopology etc.
# using PhyloBits.TrUtils
# using PhyBEARS.StateSpace
# using PhyloBits.TreeTable		# for prt, nodetimes
# using PhyBEARS.TreePass
# using PhyBEARS.Parsers
# using PhyBEARS.SSEs
# using PhyBEARS.ModelLikes
# using PhyBEARS.Optimizers
using PhyloBits
using PhyBEARS

#####################################################
# Check multiple cores/workers
#####################################################
using Distributed
Distributed.nprocs()

using Hwloc
Hwloc.num_physical_cores()
Hwloc.num_virtual_cores()


#t = @async Distributed.addprocs(7) # Adds workers
Distributed.nprocs()
# Check that you have same number of processors and threads
Distributed.nprocs()
numthreads = Base.Threads.nthreads()


#####################################################
# The newick file has the species as "6", "7" etc;
# should be "sp6", "sp7"
#####################################################
trfn = "newtree.newick"
tr = readTopology(trfn)

# Change the tip names from e.g. "6" to "sp6"
newnames = deepcopy(tr.names)
for i in 1:length(newnames)
	newnames[i] = Rpaste0(["sp", tr.names[i]])
end
tr.names .= newnames

# Change the names on the nodes
for i in 1:length(tr.node)
	if tr.node[i].leaf == true
		tr.node[i].name = Rpaste0(["sp", tr.node[i].name])
	end
end

writeTopology(tr, "tree.newick")

#######################################
# Load the tree and geog file
#######################################
trfn = "tree.newick"
tr = readTopology(trfn)
# Tree table
trdf = prt(tr)

lgdata_fn = "rangedata.data"
geog_df = PhyBEARS.Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Do they all match?
all(sort(tr.names) .== sort(geog_df.tipnames))
# Yay!


#######################################
# Create a standard DEC model on the simulated dataset
#######################################
# Yule birthrate
# The Yule-process ML birthrate is just (# internal nodes - 1)/total_tree_length
numTips = sum(trdf.nodeType .== "tip")
numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")
ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
yuleBirthRate = (numInternal-1) / ttl_tree_length

bmo = PhyBEARS.StateSpace.construct_BioGeoBEARS_model_object()
bmo.est[bmo.rownames .== "birthRate"] .= 0.24
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.001
bmo.est[bmo.rownames .== "e"] .= 0.001
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
numareas = Rncol(geog_df)-1

# Set up inputs 
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Es_v5, Es_tspan) = inputs;
numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])
root_age = maximum(trdf[!, :node_age])

prtCp(p_Es_v5)

# Problem: this simulation allowed null-range speciation (null->null, null or 000->000,000)
# Solution for now: Add the null -> null, null cladogenesis event to the p_Es_v5
birthRate = bmo.est[bmo.rownames .== "birthRate"]
add_111_to_Carray!(p_Es_v5, birthRate)

prtCp(p_Es_v5) # now 1->1,1 is an allowed cladogenesis event

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);

# Check the interpolator
p_Ds_v5.sol_Es_v5(1.0)
Es_interpolator(1.0)

# Calculate the Ds & total lnL via downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)





#######################################################
#######################################################
#######################################################
# Optimize the DEC model (fixed birthRate and deathRate=0.0)
#######################################################
#######################################################
#######################################################

pars = bmo.est[bmo.type .== "free"]
parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL", printlevel=1)
pars = [0.9, 0.9]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


#######################################################
# Best optimizer so far - 2022-03-15
#######################################################
using NLopt
pars = [0.9, 0.9]
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
opt.ftol_abs = 0.00001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################






# Get the inputs & "res" results object:
pars = optx
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
bmo_updater_v1!(inputs.bmo)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
prtCp(p_Ds_v5)
p_Ds_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
prtCp(p_Ds_v5)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

# The approximate ancestral state probabilities are in the "res" results object here:
# (I say "approximate" because the normalized downpass state likelihoods are not quite
#  identical with Felsenstein's ancestral state probabilities, which require an uppass as
#  well as a downpass; but they are typically the dominant influence.)
#
Rnames(res)

# E.g. root state probabilities:
round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root], digits=3)
# The rows are the nodes as listed in the tree table at prt(tr)
# The columns are the 8 ranges: null, A, B, C, AB, AC, BC, ABC









#######################################################
#######################################################
#######################################################
# ML inference on DEC + birth-death
#######################################################
#######################################################
#######################################################
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "fixed"
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.3288164
bmo.est[bmo.rownames .== "deathRate"] .= 0.1
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
bmo.est[:] = bmo_updater_v1(bmo) # works


inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Es_v5, Es_tspan) = inputs;

prtCp(p_Es_v5)

# Problem: this simulation allowed null-range speciation (null->null, null or 000->000,000)
# Solution for now: Add the null -> null, null cladogenesis event to the p_Es_v5
birthRate = bmo.est[bmo.rownames .== "birthRate"]
add_111_to_Carray!(p_Es_v5, birthRate)

prtCp(p_Es_v5) # now 1->1,1 is an allowed cladogenesis event


inputs.bmo.type[bmo.rownames .== "j"] .= "fixed"
parnames = bmo.rownames[bmo.type .== "free"]
lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]

pars = bmo.est[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works
inputs.bmo

prtCp(p_Es_v5)
p_Ds_v5_updater_v1!(p_Es_v5, inputs);  # WORKS 2022-03-10
prtCp(p_Es_v5)

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);

# Check the interpolator
p_Ds_v5.sol_Es_v5(1.0)
Es_interpolator(1.0)

# Calculate the Ds & total lnL via downpass
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

pars = bmo.est[bmo.type .== "free"]
#func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL")
func = x -> func_to_optimize(x, parnames, inputs, p_Ds_v5; returnval="lnL")

pars = [0.9, 0.9, 0.3, 0.2]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


#######################################################
# Best optimizer so far - 2022-03-15
#######################################################
using NLopt
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
opt.ftol_abs = 0.00001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################





# Get the inputs & res:
pars = optx
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
inputs_updater_v1!(inputs)
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

printlevel=1
returnval="bgb_lnL"
func(pars);

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Ds_v5.uE, Es_tspan, p_Ds_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
prtCp(p_Ds_v5)
p_Ds_v7 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, sol_Es_v5=sol_Es_v5);
prtCp(p_Ds_v5)


(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(inputs.res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)

# The approximate ancestral state probabilities are in the "res" results object here:
# (I say "approximate" because the normalized downpass state likelihoods are not quite
#  identical with Felsenstein's ancestral state probabilities, which require an uppass as
#  well as a downpass; but they are typically the dominant influence.)
#
Rnames(res)

# E.g. root state probabilities:
round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root], digits=3)
# The rows are the nodes as listed in the tree table at prt(tr)
# The columns are the 8 ranges: null, A, B, C, AB, AC, BC, ABC









