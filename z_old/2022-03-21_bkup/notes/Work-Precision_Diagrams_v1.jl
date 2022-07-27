
#######################################################
# Work-Precision Diagrams on ClaSSE Es & Ds calculations
#######################################################

# Modified from:
# 100 Independent Linear Work-Precision Diagrams
# Chris Rackauckas
# https://benchmarks.juliadiffeq.org/html/NonStiffODE/linear_wpd.html

using Pkg
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
Pkg.add(PackageSpec(url="https://github.com/JuliaDiffEq/deSolveDiffEq.jl"))
using deSolveDiffEq # https://docs.juliadiffeq.org/stable/solvers/ode_solve/index.html


using Random
Random.seed!(123)
gr()

# 2D Linear ODE
function f(du,u,p,t)
  @inbounds for i in eachindex(u)
    du[i] = 1.01*u[i]
  end
end
function f_analytic(u₀,p,t)
  u₀*exp(1.01*t)
end
tspan = (0.0,10.0)
prob = ODEProblem(ODEFunction(f,analytic=f_analytic),rand(100,100),tspan)

abstols = 1.0 ./ 10.0 .^ (3:13)
reltols = 1.0 ./ 10.0 .^ (0:10);



setups = [Dict(:alg=>DP5())
          Dict(:alg=>ode45())
          Dict(:alg=>ARKODE(Sundials.Explicit(),etable=Sundials.DORMAND_PRINCE_7_4_5))
          Dict(:alg=>Tsit5())]
solnames = ["OrdinaryDiffEq";"ODE";"Sundials ARKODE";"OrdinaryDiffEq Tsit5"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;names=solnames,save_everystep=false,numruns=10)
plot(wp)


# Set up a 6-area, 64-state problem

cd("/GitHub/BioGeoJulia.jl/notes/")
include("ModelLikes.jl")
import .ModelLikes
#using .Tmp

using Profile     # for @profile
using DataFrames  # for DataFrame
using PhyloNetworks
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs


	numareas = 3
	areas_list = collect(1:numareas)
	states_list = areas_list_to_states_list(areas_list, numareas, false)
	n = length(states_list)
	max_numareas = length(areas_list)

	
	res = construct_Res(tr, n)
	rootnodenum = tr.root
	trdf = prt(tr, rootnodenum)
	tipnodes = trdf[!,1][trdf[!,10].=="tip"]
	
	birthRate = 0.222222
	deathRate = 0.1
	
	d_val = 0.01
	e_val = 0.001
	j_val = 0.0
	
	dmat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list)))
	amat=reshape(repeat([1.0], (length(areas_list)^2)), (length(areas_list),length(areas_list)))
	elist = repeat([1.0], length(areas_list))
	
	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	print("\n")
	print("elist")
	print(elist)
	print("\n")
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	event_type_vals = Qmat.event_type_vals


	Cparams = default_Cparams()
	maxent_constraint_01 = 0.0
	maxent01symp = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01sub = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01jump = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent_constraint_01 = 0.0
	maxent01vic = relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01)
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
	Carray = setup_DEC_Cmat(areas_list, states_list, maxent01, Cparams)

	# Possibly varying parameters
	mu_vals = repeat([deathRate], n)

	params = (mu_vals=mu_vals, Qij_vals=Qmat.Qij_vals, Cijk_vals=birthRate.*Carray.Cijk_vals)

	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	p_indices = (Qarray_ivals=Qmat.Qarray_ivals, Qarray_jvals=Qmat.Qarray_jvals, Carray_ivals=Carray.Carray_ivals, Carray_jvals=Carray.Carray_jvals, Carray_kvals=Carray.Carray_kvals)

	# True/False statements by index
	# The calculation of dEi and dDi for state i involves many
	# ==i and !=i operations across Q and C. These only need to be 
	# done once per problem (may or may not save time to 
	# pre-calculate).
	# 
	# Pre-allocating the Carray_ivals .== i, Qarray_jvals[Qarray_ivals .== i
	# Reduces GC (Garbage Collection) from 40% to ~5%
	# 10+ times speed improvement (!)
	Qi_eq_i = Any[]
	Ci_eq_i = Any[]

	Qi_sub_i = Any[]
	Qj_sub_i = Any[]

	# These are the (e.g.) j state-indices (left descendant) when the ancestor==state i
	Ci_sub_i = Any[]
	Cj_sub_i = Any[]
	Ck_sub_i = Any[]

	# The push! operation may get slow at huge n
	# This will have to change for non-Mk models
	for i in 1:n
		push!(Qi_eq_i, Qmat.Qarray_ivals .== i)
		push!(Qi_sub_i, Qmat.Qarray_ivals[Qarray_ivals .== i])
		push!(Qj_sub_i, Qmat.Qarray_jvals[Qarray_ivals .== i])

		push!(Ci_eq_i, Carray.Carray_ivals .== i)
		push!(Ci_sub_i, Carray.Carray_ivals[Carray.Carray_ivals .== i])
		push!(Cj_sub_i, Carray.Carray_jvals[Carray.Carray_ivals .== i])
		push!(Ck_sub_i, Carray.Carray_kvals[Carray.Carray_ivals .== i])
	end

	# Inputs to the Es calculation
	p_TFs = (Qi_eq_i=Qi_eq_i, Ci_eq_i=Ci_eq_i, Qi_sub_i=Qi_sub_i, Qj_sub_i=Qj_sub_i, Ci_sub_i=Ci_sub_i, Cj_sub_i=Cj_sub_i, Ck_sub_i=Ck_sub_i)
	p_orig = (n=n, params=params, p_indices=p_indices)
	p = p_orig
	p_Es_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs)
	
	# Solutions to the E vector
	u0_Es = repeat([0.0], 1*n)
	uE = repeat([0.0], n)
	Es_tspan = (0.0, 1.2*trdf[tr.root,:node_age]) # 110% of tree root age


	#######################################################
	# Downpass with ClaSSE
	#######################################################
	# Solve for the Ds
	du = repeat([0.0], n)
	u0 = repeat([0.0], n)
	u0[2] = 1.0  # Starting likelihood
	#tspan = (0.0, 2.0*trdf[tr.root,:node_age]) # 110% of tree root age
	current_nodeIndex = 1
	res = construct_Res(tr, n)
	
	res.likes_at_each_nodeIndex_branchTop
	for i in 1:tr.numTaxa
		res.likes_at_each_nodeIndex_branchTop[tipnodes[i]] = u0;
	end
	#res.likes_at_each_nodeIndex_branchTop[6] = u0;
	res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]

	# Updates res
	res_orig = res
	res_orig.likes_at_each_nodeIndex_branchTop

	solver_options = construct_SolverOpt()
	solver_options.solver=Tsit5()
	solver_options.save_everystep = false
	solver_options.abstol = 1.0e-6
	solver_options.reltol = 1.0e-6
	
	print("\nSolving the Es once, for the whole tree timespan...")
	
	# Solve the Es
	prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, u0_Es, Es_tspan, p_Es_v5)

	sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol=1e-9);

setups = [Dict(:alg=>DP5())
          Dict(:alg=>ode45())
          Dict(:alg=>ARKODE(Sundials.Explicit(),etable=Sundials.DORMAND_PRINCE_7_4_5))
          Dict(:alg=>Tsit5())]
solnames = ["OrdinaryDiffEq";"ODE";"Sundials ARKODE";"OrdinaryDiffEq Tsit5"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;names=solnames,save_everystep=false,numruns=10)
plot(wp)


setups = [Dict(:alg=>ode45())]
# deSolveDiffEq.lsoda

abstols = 1.0 ./ 10.0 .^ (3:13)
reltols = 1.0 ./ 10.0 .^ (0:10);

abstols = [0.1, 0.001, 0.0000001];
reltols = [0.1, 0.001, 0.0000001];



solnames = ["OrdinaryDiffEq Tsit5"]
wp = WorkPrecisionSet(prob_Es_v5,abstols,reltols,setups;names=solnames,save_everystep=false,numruns=100)
plot(wp)

wp = WorkPrecisionSet(prob_Es_v5,abstols,reltols,setups)