
print("\n")
print("""Running...\ninclude("/GitHub/BioGeoJulia.jl/notes/startup_basic_tests_v1.jl")\n...to run some basic BioGeoJulia code.\n""")
print("\n")


using Pkg
using RCall       # for df_to_Rdata, reval, g = globalEnv -- it loads a bunch of R, but just onces
using DataFrames			# for e.g. DataFrame()
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using LSODA          # for lsoda()
using BenchmarkTools # for @time
using InvertedIndices # for Not
using Statistics			# for e.g. mean()
	
# For SSE stuff
using OrdinaryDiffEq
using DiffEqDevTools
using ODEInterfaceDiffEq
using ODE

# For slow packages
print("using Plots...\n")
using Plots						# basic plots
print("using DifferentialEquations...\n")
using DifferentialEquations # for ODEProblem
print("using PhyloNetworks...\n")
using PhyloNetworks



print("Loading BioGeoJulia...\n")

# DON'T DO THIS:
#Pkg.rm("BioGeoJulia")
#Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
#Pkg.add(PackageSpec(url="https://github.com/nmatzke/BioGeoJulia.jl"))

# DON'T ADD BioGeoJulia, INSTEAD DO:
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

using BioGeoJulia

using BioGeoJulia.TrUtils
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs



#######################################################
# Auto-testing code
#######################################################

#great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
great_ape_newick_string = "((chimp:1,human:1):1,gorilla:2);"
tr = readTopology(great_ape_newick_string)
tr

#res2 = construct_Res(tr)
res2 = construct_Res(tr, 5)


rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
trdf


n=2

# 4 states, no Q
# 3 tips in state 2, branch is 1 Mya long
birthRate = 0.222222
deathRate = 0.0
d_val = 0.0
e_val = 0.0
j_val = 0.0

# Define Qarray - zeros
Qarray_ivals = collect(1:(n-1))
Qarray_jvals = collect(2:n)
Qarray_ivals = vcat(Qarray_ivals, collect(2:n))
Qarray_jvals = vcat(Qarray_jvals, collect(1:(n-1)))
Qij_vals = vcat(repeat([0.0],(n-1)), repeat([0.0],(n-1)))
hcat(Qarray_ivals, Qarray_jvals, Qij_vals)

# A 2-state DEC matrix
Qarray_ivals = [2,1]
Qarray_jvals = [1,2]
Qij_vals = [e_val, d_val]
hcat(Qarray_ivals, Qarray_jvals, Qij_vals)


Qij_vals[((Qarray_ivals .== 1) .+ (Qarray_jvals .!= 1)) .== 2]


# Carray: Cladogenetic parameters
# column 1: state i
# column 2: state j
# column 3: state k
# column 4: lambda_ijk (a parameter that is at least possibly nonzero, under the given model)
Carray_ivals = collect(1:n)
Carray_jvals = collect(1:n)
Carray_kvals = collect(1:n)
Cijk_vals = repeat([birthRate], n)

mu_vals = repeat([deathRate], n)

hcat(Carray_ivals, Carray_jvals, Carray_kvals, Cijk_vals)
mu_vals

# Possibly varying parameters
params = (mu_vals=mu_vals, Qij_vals=Qij_vals, Cijk_vals=Cijk_vals)

# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
p_indices = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals)

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
	push!(Qi_eq_i, Qarray_ivals .== i)
	push!(Qi_sub_i, Qarray_ivals[Qarray_ivals .== i])
	push!(Qj_sub_i, Qarray_jvals[Qarray_ivals .== i])

	push!(Ci_eq_i, Carray_ivals .== i)
	push!(Ci_sub_i, Carray_ivals[Carray_ivals .== i])
	push!(Cj_sub_i, Carray_jvals[Carray_ivals .== i])
	push!(Ck_sub_i, Carray_kvals[Carray_ivals .== i])
end

# Inputs to the Es calculation
p_TFs = (Qi_eq_i=Qi_eq_i, Ci_eq_i=Ci_eq_i, Qi_sub_i=Qi_sub_i, Qj_sub_i=Qj_sub_i, Ci_sub_i=Ci_sub_i, Cj_sub_i=Cj_sub_i, Ck_sub_i=Ck_sub_i)
p_orig = (n=n, params=params, p_indices=p_indices)
p = p_orig
p_Es_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs)


params = (mu_vals=mu_vals, Qij_vals=Qij_vals, Cijk_vals=Cijk_vals)

# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
p_indices = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals)

# Solutions to the E vector
u0_Es = repeat([0.0], 1*n)
uE = repeat([0.0], n)
tspan = (0.0, 2.0*trdf[tr.root,:node_age]) # 110% of tree root age



prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, u0_Es, tspan, p_Es_v5)
sol_Es_v5 = solve(prob_Es_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9);

#@benchmark sol_Es_v5 = solve(prob_Es_v5, lsoda(), save_everystep=true, abstol = 1e-9, reltol = 1e-9) # slower
#@benchmark sol_Es_v5 = solve(prob_Es_v5, lsoda(), save_everystep=false, abstol = 1e-9, reltol = 1e-9)


#p_Ds = (n=n, params=params, p_indices=p_indices, sol_Es=sol_Es, uE=uE)
p_Ds_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs, sol_Es_v5=sol_Es_v5, uE=uE)


#######################################################
# Downpass with ClaSSE
#######################################################
# Solve for the Ds
du = repeat([0.0], n)
u0 = repeat([0.0], n)
u0[2] = 1.0  # Starting likelihood
tspan = (0.0, 1.0) # Shorter
current_nodeIndex = 1
res = construct_Res(tr, n)
#u0, tspan, p_Ds_v5

#inputs = setup_inputs_branchOp_ClaSSE_Ds_v5(u0, tspan, p_Ds_v5; solver="Tsit5()", 
#				 save_everystep="false", abstol="1e-9", reltol="1e-9")

res.likes_at_each_nodeIndex_branchTop

res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
res.likes_at_each_nodeIndex_branchTop[1] = u0;
res.likes_at_each_nodeIndex_branchTop[2] = u0;
res.likes_at_each_nodeIndex_branchTop[4] = u0;
#res.likes_at_each_nodeIndex_branchTop[6] = u0;
res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]

# Updates res
res_orig = res
res_orig.likes_at_each_nodeIndex_branchTop


solver_options = construct_SolverOpt()
(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res, trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10);

res.likes_at_each_nodeIndex_branchTop
res.sumLikes_at_node_at_branchTop
log.(res.sumLikes_at_node_at_branchTop[res.sumLikes_at_node_at_branchTop.!=0.0])
sum(log.(res.sumLikes_at_node_at_branchTop[res.sumLikes_at_node_at_branchTop.!=0.0]))


print("\nres.likes_at_each_nodeIndex_branchTop:\n")
print(res.likes_at_each_nodeIndex_branchTop)

print("\nres.sumLikes_at_node_at_branchTop:\n")
print(res.sumLikes_at_node_at_branchTop)

print("\nlog.(res.sumLikes_at_node_at_branchTop[res.sumLikes_at_node_at_branchTop.!=0.0]):\n")
print(log.(res.sumLikes_at_node_at_branchTop[res.sumLikes_at_node_at_branchTop.!=0.0]))

print("\nsum(log.(res.sumLikes_at_node_at_branchTop[res.sumLikes_at_node_at_branchTop.!=0.0])):\n")
print(sum(log.(res.sumLikes_at_node_at_branchTop[res.sumLikes_at_node_at_branchTop.!=0.0])))


total_calctime_in_sec
iteration_number





