
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
wd = "/GitHub/PhyBEARS.jl/test/simtree_5taxa_SSE/"
cd(wd)

# Simple tree
trfn = "simtree_5taxa.newick"
tr = readTopology("(sp5:9.201374725,(sp6:3.944281182,(sp7:1.688329525,(sp8:0.07452322119,sp9:0.07452322119):1.613806304):2.255951656):5.257093543);")
tr = readTopology(trfn)


trdf = prt(tr);
oldest_possible_age = 4.0

lgdata_fn = "simtree_tipstates.txt"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = false
numareas = Rncol(geog_df)-1
max_range_size = 1
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
#bmo.type[bmo.rownames .== "j"] .= "free";
bmo.est[bmo.rownames .== "birthRate"] .= 0.22222222;
bmo.est[bmo.rownames .== "deathRate"] .= 0.11111111; # (for one of them)
bmo.est[bmo.rownames .== "d"] .= 0.0;
bmo.est[bmo.rownames .== "e"] .= 0.0;
bmo.est[bmo.rownames .== "a"] .= 0.03; # (for one of them
bmo.est[bmo.rownames .== "j"] .= 0.0;
bmo.est[bmo.rownames .== "u"] .= 0.0;
bmo.min[bmo.rownames .== "u"] .= 0.0;
bmo.max[bmo.rownames .== "u"] .= 0.0;

bmo.type[bmo.rownames .== "j"] .= "fixed";
bmo.type[bmo.rownames .== "u"] .= "fixed";
bmo.type[bmo.rownames .== "birthRate"] .= "fixed";
bmo.type[bmo.rownames .== "deathRate"] .= "fixed";

birthRate = 0.2222222222
deathRate = 0.1111111111
bd_liks(tr, birthRate, deathRate)

# Update the bmo, BEFORE setup_DEC_SSE2
bmo_updater_v1!(bmo)

# Set up the model
manual_states_list=NaN; area_names=LETTERS(1:numareas)
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=include_null_range, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

prtQi(inputs)
prtCi(inputs)

bmo.est[:] = bmo_updater_v2(bmo, inputs.setup.bmo_rows);

inputs.setup.txt_states_list



#######################################################
# Put in the parameters straight from diversitree
#######################################################
p_Ds_v5.params.Cijk_vals[1] = 0.22222222
p_Ds_v5.params.Cijk_vals[2] = 0.22222222
p_Ds_v5.params.mu_vals[1] = 0.11111111
p_Ds_v5.params.mu_vals[2] = 0.01111111
p_Ds_v5.params.Qij_vals[1] = 0.06000000
p_Ds_v5.params.Qij_vals[2] = 0.05000000

prtQp(p_Ds_v5)
prtCp(p_Ds_v5)
p_Ds_v5.params.mu_vals

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

prtQp(p_Ds_v7)
prtCp(p_Ds_v7)
p_Ds_v7.params.mu_vals

# R:
# lik(pars)
# [1] -13.08989


prtQp(p)
p.p_TFs.Qj_sub_i
p.p_TFs.Qij_vals_sub_i


prtCp(p)
p.p_TFs.Ci_sub_i
p.p_TFs.Cj_sub_i
p.p_TFs.Ck_sub_i
p.p_TFs.Qij_vals_sub_i


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

solver_options.solver = CVODE_BDF{:Newton, :GMRES, Nothing, Nothing}(0, 0, 0, false, 10, 5, 7, 3, 10, nothing, nothing, 0)
#solver_options.solver = Tsit5()
solver_options.solver = Vern9()
solver_options.abstol = 1e-12
solver_options.reltol = 1e-12
solver_options.save_everystep = false
# ancestral_range_estimation
# This term is preferable to e.g. "ancestral area reconstruction"


# Check if solver is functional
truth = [8.322405e-13, 0.1129853, 0.677912, 0.2091026]
#u0 = [0, 0.1129853, 0.677912, 0.2091026]
u0 = [0, 0.125, 0.75, 0.125]
solver_options.solver = CVODE_BDF(linear_solver=:GMRES)
#solver_options.solver = Vern9()
solver_options.abstol = 1.0e-13
solver_options.reltol = 1.0e-13

solver_options.save_everystep = true
solver_options.saveat = seq(2.0, 3.0, 0.1)
tspan = (3.0, 2.0)
#tspan = (2.0, 3.0)
current_nodeIndex = 5
(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)= branchOp_ClaSSE_Ds_v7(current_nodeIndex, res; u0, tspan, p_Ds_v7, solver_options=solver_options);
sol_Ds

sol_Ds(2.0)
sol_Ds(2.1)
sol_Ds(3.0)


prob_Ds_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Ds_v7_simd_sums, u0, tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);
sol_Ds_v7(2.0)
sol_Ds_v7(2.1)
sol_Ds_v7(3.0)

include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
tspan = (3.0, 2.0)

prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4states2, u0, tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
sol_Ds_v7(2.0)
sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))

err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)

tspan = (2.0, 3.0)

prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4states2A, u0, tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
sol_Ds_v7(3.0)
sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))

best = (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)

Rcbind(truth, best, err)

# CLOSEST

d_val = 0.1010557
e_val = 1.0e-12


Qmat = [0.0 0.0 0.0 0.0;
e_val 0.0 0.0 d_val;
e_val 0.0 0.0 d_val;
0.0 e_val e_val 0.0]


Qmat[make_diag_TF(dim(Qmat)[1])] .= -1 .* flat2(sum(Qmat; dims=2))
unscaled = exp(transpose(Qmat)*1.0)*u0
scaled = unscaled ./ sum(unscaled)

exp_err = truth .- scaled
sum(exp_err)
sum(abs.(exp_err))







#######################################################
# Try with full clado calculation, positive extinction
#######################################################
include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")

tspan = (2.0, 3.0)

birthRate = p_Ds_v7.bmo.est[p_Ds_v7.setup.bmo_rows.birthRate]
p_Ds_v7.bmo.est[p_Ds_v7.setup.bmo_rows.deathRate] = 0.5 * birthRate
p_Ds_v5_updater_v1!(p_Ds_v7, inputs);
p_Ds_v5_updater_v1!(p_Es_v7, inputs);
p_Ds_v7.params.mu_vals
p_Es_v7.params.mu_vals

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Es_v7(2.0)
sol_Es_v7(2.5)
sol_Es_v7(3.0)

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);



# Extinction equiprobable, clado left out, no change in ancestral states
prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4states2A, u0, tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
sol_Ds_v7(3.0)
sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))

best = (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)

ncA = Rcbind(truth, best, err)


# Extinction equiprobable, clado included; no change in ancestral states
prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4states2B, u0, tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
sol_Ds_v7(3.0)
sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))

best = (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)

ncB = Rcbind(truth, best, err)


# Extinction equiprobable, clado included; no change in ancestral states
prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4states2C, u0, tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
sol_Ds_v7(3.0)
sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))

best = (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)

ncC = Rcbind(truth, best, err)



#######################################################
# Try with unequal extinction
#######################################################
include("/Users/nickm/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")

p_Es_v7.params.mu_vals[1] = 10.0	# null range has high extinction rate
p_Es_v7.params.mu_vals[2] = 2 * p_Ds_v7.params.mu_vals[3]		#	range AB has 0.0 extinction rate
p_Es_v7.params.mu_vals[4] = 0.0		#	range AB has 0.0 extinction rate
p_Es_v7.params.mu_vals

# Solve the Es
prob_Es_v7 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v7.uE, Es_tspan, p_Es_v7);
sol_Es_v7 = solve(prob_Es_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Es_v7(2.0)
sol_Es_v7(2.5)
sol_Es_v7(3.0)

p = p_Ds_v7 = (n=p_Es_v7.n, params=p_Es_v7.params, p_indices=p_Es_v7.p_indices, p_TFs=p_Es_v7.p_TFs, uE=p_Es_v7.uE, terms=p_Es_v7.terms, setup=p_Es_v7.setup, states_as_areas_lists=p_Es_v7.states_as_areas_lists, use_distances=p_Es_v7.use_distances, bmo=p_Es_v7.bmo, sol_Es_v5=sol_Es_v7);


# Extinction equiprobable, clado left out, no change in ancestral states
prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4statesA, u0, tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
solA = sol_Ds_v7(3.0)
sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))

best = (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)

wcA = Rcbind(truth, best, err)


# Extinction equiprobable, clado included; no change in ancestral states
prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4states2B, u0, tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
solB = sol_Ds_v7(3.0)
sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))

best = (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)

wcB = Rcbind(truth, best, err)


# Extinction equiprobable, clado included; no change in ancestral states
prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4states2C, u0, tspan, p_Ds_v7);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
solC = sol_Ds_v7(3.0)
sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))

best = (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)

wcC = Rcbind(truth, best, err)

ncA .== ncB
ncA .== ncC

ncA .== wcA
ncA .== wcB
ncA .== wcC

wcA .== wcB
wcA .== wcC


solA
solB
solC











# Try with written-out matrix
include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")
Cmat = make_ctable_single_events(prtCp(p))

ptmp = (n=n, mu_vals=p.params.mu_vals, Qij_vals = p.params.Qij_vals, Cijk_vals=Cmat.val, Qarray_ivals=p.p_indices.Qarray_ivals, Qarray_jvals=p.p_indices.Qarray_jvals, Carray_ivals=Cmat.i, 	Carray_jvals=Cmat.j, Carray_kvals=Cmat.k, sol_Es_v5=p.sol_Es_v5, uE=p.uE);


tspan = (3.0, 2.0)

prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4states3, u0, tspan, ptmp);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
sol_Ds_v7(3.0)
sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))

best = (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)



include("/GitHub/PhyBEARS.jl/notes/nodeOp_Cmat_uppass_v12.jl")

solver_options.solver = AutoTsit5(Rosenbrock23())
#solver_options.solver = CVODE_BDF(linear_solver=:GMRES)
#solver_options.solver = Vern9()
solver_options.abstol = 1.0e-6
solver_options.reltol = 1.0e-6


tspan = (2.0, 3.0)

prob_Ds_v7 = DifferentialEquations.ODEProblem(calcDs_4states3, u0, tspan, ptmp);
sol_Ds_v7 = solve(prob_Ds_v7, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol, dtmax=1.0/100000);

sol_Ds_v7(2.0) ./ sum(sol_Ds_v7(2.0))
sol_Ds_v7(2.1) ./ sum(sol_Ds_v7(2.1))
sol_Ds_v7(3.0)
sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))

best = (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
err = truth .- (sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0)))
sum(abs.(err))
sum(err)



# julia> sol_Ds_v7(3.0) ./ sum(sol_Ds_v7(3.0))
# 4-element Vector{Float64}:
#  8.039148447906605e-13
#  0.11484497782742717
#  0.6890698669637856
#  0.19608515520798325
# 
# julia> sol_Ds_v7(3.0)
# 4-element Vector{Float64}:
#  8.749999999997373e-13
#  0.12500000000016923
#  0.7500000000001692
#  0.21342373749997345




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

