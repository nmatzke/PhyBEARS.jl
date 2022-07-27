
#######################################################
# Add dependencies to PhyBEARS package development
# https://discourse.julialang.org/t/how-to-manage-dependencies-of-developed-packages/25481/2
#######################################################
# To get your package to have correct-ish [deps] fields
# etc. in Project.toml and Manifest.toml:

cd("/GitHub/PhyBEARS.jl")
Pkg.activate(".")
Pkg.add("Combinatorics")
Pkg.add("DataFrames")		# for DataFrame
Pkg.add("Dates"	)				# for e.g. Dates.now(), DateTime
Pkg.add("Distributed")  # for e.g. @spawn
Pkg.add("DifferentialEquations")  # for e.g. ODEProblem
Pkg.add("Random")  # for MersenneTwister()

# BioSequences before PhyloNetworks
# https://github.com/BioJulia/BioSequences.jl
]
registry add https://github.com/BioJulia/BioJuliaRegistry.git
add BioSequences
add PhyloNetworks

# de-activate with blank "activate"
Pkg.activate()
# ^C

#######################################################
# NOTE: You may also have to then manually go into
# Project.toml to lower the minimum version numbers of
# some packages.
#
# Commit to GitHub, then continue with below
#######################################################



#######################################################
# Load dependencies of PhyBEARS (if needed)
#######################################################
import Pkg
using Pkg

# Add these packages
#Pkg.add("Combinatorics")
Pkg.add(Pkg.PackageSpec(;name="Combinatorics", version="1.0.0"))

# Then import or using as needed
import Combinatorics.combinations

Pkg.resolve()






#######################################################
# FROM FRESH JULIA: Load the PhyBEARS package
#######################################################
import Pkg
using Pkg
Pkg.rm("PhyBEARS")
Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl"))
using PhyBEARS

# Run the tests directory
Pkg.test("PhyBEARS")


#######################################################
# Re-run the tests directory
#######################################################
# First, commit to Master (quick), then:
Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl"))
using PhyBEARS
Pkg.test("PhyBEARS")



#######################################################
# FROM FRESH JULIA: Load the PhyBEARS package
#######################################################
import Pkg
using Pkg


using DifferentialEquations # for ODEProblem
using LSODA          # for lsoda()
using BenchmarkTools # for @time
using InvertedIndices # for Not
using Distributed
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloNetworks
#using PhyloPlots
# https://github.com/JuliaLang/julia/issues/28276
#using Sundials # for CVODE_BDF() e.g.





#######################################################
# Try some functions!
#######################################################
#include("/drives/Dropbox/_njm/__julia/julia4Rppl_v4.jl")

# Try some functions
Pkg.rm("PhyBEARS")
Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl"))
using PhyBEARS

using PhyBEARS.TrUtils
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs





great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

#res2 = construct_Res(tr)
res2 = construct_Res(tr, 5)


rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
trdf


n=50000

# 4 states, no Q
# 3 tips in state 2, branch is 1 Mya long
birthRate = 0.222222
deathRate = 0.1
d_val = 0.1
e_val = 0.01
j_val = 0.0

# Define Qarray - zeros
Qarray_ivals = collect(1:(n-1))
Qarray_jvals = collect(2:n)
Qarray_ivals = vcat(Qarray_ivals, collect(2:n))
Qarray_jvals = vcat(Qarray_jvals, collect(1:(n-1)))
Qij_vals = vcat(repeat([0.0],(n-1)), repeat([0.0],(n-1)))
hcat(Qarray_ivals, Qarray_jvals, Qij_vals)

# A 4-state DEC matrix
Qarray_ivals = [2,3,2,3,4,4]
Qarray_jvals = [1,1,4,4,2,3]
Qij_vals = [e_val, e_val, d_val, d_val, e_val, e_val]
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
Ci_eq_i = Any[]
Qi_eq_i = Any[]
Cj_sub_i = Any[]
Ck_sub_i = Any[]
Qj_sub_i = Any[]

for i in 1:n
	push!(Ci_eq_i, Carray_ivals .== i)
	push!(Qi_eq_i, Qarray_ivals .== i)
	push!(Cj_sub_i, Carray_jvals[Carray_ivals .== i])
	push!(Ck_sub_i, Carray_kvals[Carray_ivals .== i])
	push!(Qj_sub_i, Qarray_jvals[Qarray_ivals .== i])
end


p_TFs = (Ci_eq_i=Ci_eq_i, Qi_eq_i=Qi_eq_i, Cj_sub_i=Cj_sub_i, Ck_sub_i=Ck_sub_i, Qj_sub_i=Qj_sub_i)

p_orig = (n=n, params=params, p_indices=p_indices)
p = p_orig
p_Es_v5 = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs)

params = (mu_vals=mu_vals, Qij_vals=Qij_vals, Cijk_vals=Cijk_vals)

# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
p_indices = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals)

# Solutions to the E vector
u0_Es = repeat([0.0], 1*n)
uE = repeat([0.0], n)
tspan = (0.0, 1.1*trdf[tr.root,:node_age]) # 110% of tree root age



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
res.likes_at_each_nodeIndex_branchTop[6] = u0;
res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]

# Updates res
res_orig = res
res_orig.likes_at_each_nodeIndex_branchTop
(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res, trdf=trdf, p_Ds_v5=p_Ds_v5, max_iterations=10^10);
res.likes_at_each_nodeIndex_branchTop
res.sum_likes_at_nodes
log.(res.sum_likes_at_nodes[res.sum_likes_at_nodes.!=0.0])
sum(log.(res.sum_likes_at_nodes[res.sum_likes_at_nodes.!=0.0]))

total_calctime_in_sec
iteration_number
