
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
wd = "/GitHub/PhyBEARS.jl/test/apes_SSE/"
cd(wd)

# Simple tree
tr = readTopology("(((chimp:1,human:1):1,gorilla:2):1,orang:3);")

trdf = prt(tr);
oldest_possible_age = 4.0

lgdata_fn = "ancstates_tr_orang.data"
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
bmo.est[bmo.rownames .== "d"] .= 1.010557e-01;
bmo.est[bmo.rownames .== "e"] .= 1.000000e-12;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.min[bmo.rownames .== "u"] .= 0.0;
bmo.max[bmo.rownames .== "u"] .= 0.0;

bmo.type[bmo.rownames .== "j"] .= "fixed";
bmo.type[bmo.rownames .== "u"] .= "fixed";
bmo.type[bmo.rownames .== "birthRate"] .= "fixed";
bmo.type[bmo.rownames .== "deathRate"] .= "fixed";

birthRate = 0.2222222222
deathRate = 0.0
bd_liks(tr, birthRate, deathRate)

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

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);

# Solve the Ds
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


# Solve the Ds, single branch
u0 = res.likes_at_each_nodeIndex_branchTop[1]
prob_Ds_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Ds_v7_simd_sums, u0, Es_tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_Ds_v7(0.0)
sol_Ds_v7(1.0)
sol_Ds_v7(1.5)

branch_bottom_time = 1.5
branch_top_time = 0.0
reverse_tspan = (branch_bottom_time, branch_top_time)
u0rev = sol_Ds_v7(branch_bottom_time)
prob_Ds_v7rev = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Ds_v7_simd_sums, u0rev, reverse_tspan, p_Ds_v7);
sol_Ds_v7rev = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_Ds_v7rev(1.5)
sol_Ds_v7rev(1.0)
sol_Ds_v7rev(0.0)

all(sol_Ds_v7(0.0) .== sol_Ds_v7rev(0.0))
all(sol_Ds_v7(1.0) .== sol_Ds_v7rev(1.0))
all(sol_Ds_v7(1.5) .== sol_Ds_v7rev(1.5))






#######################################################
# Estimate ancestral states
#######################################################


#solver_options.solver = Tsit5()
#solver_options.solver = Vern9()
solver_options.abstol = 1e-12
solver_options.reltol = 1e-12
solver_options.save_everystep = false
# ancestral_range_estimation
# This term is preferable to e.g. "ancestral area reconstruction"


# Check if solver is functional
u0 = [8.322405e-13, 0.1129853, 0.677912, 0.2091026]
solver_options.solver
solver_options.save_everystep = true
solver_options.saveat = seq(2.0, 3.0, 0.1)
tspan = (2.0, 3.0)
(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)= branchOp_ClaSSE_Ds_v7(current_nodeIndex, res; u0, tspan, p_Ds_v7, solver_options=solver_options);

sol_Ds(1.0)
sol_Ds(2.1)
sol_Ds(2.0)
sol_Ds(2.1)



Rnames(res)

rootnode = inputs.res.root_nodeIndex

lnode = trdf[rootnode,"leftNodeIndex"]
rnode = trdf[rootnode,"rightNodeIndex"]

# ACE for left descendant
nodenum = rootnode
nodelikes = res.normlikes_at_each_nodeIndex_branchTop[nodenum]





R_order =  sort(trdf, :Rnodenums).nodeIndex

uppass_edgematrix = res.uppass_edgematrix

include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
current_nodeIndex = 7
x = nodeOp_Cmat_uppass_v7!(res, current_nodeIndex, trdf, p_Ds_v7, solver_options)

solver_options.abstol = 1.0e-13
solver_options.reltol = 1.0e-13
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

res.uppass_probs_at_each_nodeIndex_branchBot[R_order,:]
res.anc_estimates_at_each_nodeIndex_branchBot[R_order,:]
res.uppass_probs_at_each_nodeIndex_branchTop[R_order,:]
res.anc_estimates_at_each_nodeIndex_branchTop[R_order,:]

res.normlikes_at_each_nodeIndex_branchTop[R_order,:]

p = p_Ds_v7;
ctable = prtCp(p)
ctable2 = ctable[ctable.pair .== 2, :]
tmpj = deepcopy(ctable2.j)
tmpk = deepcopy(ctable2.k)
ctable2.j .= tmpk
ctable2.k .= tmpj

ctable = Rrbind(ctable1, ctable2)
ctable.prob .= ctable.prob ./ ctable.wt
ctable.rate .= ctable.rate ./ ctable.wt
ctable.val .= ctable.val ./ ctable.wt
ctable.rates_t .= ctable.rates_t ./ ctable.wt
ctable

# Calculating uppass LEFT branch for Julia node 5 (R 6), sister Julia 6 (R 4)
ancnode = 7
lnode = 5
rnode = 6

res.normlikes_at_each_nodeIndex_branchTop[rnode]

ancprobs = repeat([1.0/n], n)
lprobs = repeat([1.0], n)
rprobs = res.normlikes_at_each_nodeIndex_branchBot[rnode]

ancprobs_by_scenario = ancprobs[ctable.i] 
lprobs_by_scenario = lprobs[ctable.j] 
rprobs_by_scenario = rprobs[ctable.k] 

relprob_each_split_scenario = ancprobs_by_scenario .* lprobs_by_scenario .* rprobs_by_scenario .* ctable.val

relprob_each_split_scenario = relprob_each_split_scenario ./ sum(relprob_each_split_scenario)

uppass_lprobs = repeat([0.0], n)
uppass_rprobs = repeat([0.0], n)

for statei in 1:n
	# Left
	uppass_lprobs[statei] = sum(relprob_each_split_scenario[ctable.j .== statei])
	uppass_rprobs[statei] = sum(relprob_each_split_scenario[ctable.k .== statei])
end

uppass_lprobs
uppass_rprobs
uppass_lprobs ./ sum(uppass_lprobs)
uppass_rprobs ./ sum(uppass_rprobs)

