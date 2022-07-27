
#######################################################
# Add dependencies to BioGeoJulia package development
# https://discourse.julialang.org/t/how-to-manage-dependencies-of-developed-packages/25481/2
#######################################################
# To get your package to have correct-ish [deps] fields
# etc. in Project.toml and Manifest.toml:

cd("/GitHub/BioGeoJulia.jl")
Pkg.activate(".")
Pkg.add("Combinatorics")
Pkg.add("DataFrames")		# for DataFrame
Pkg.add("Dates"	)				# for e.g. Dates.now(), DateTime
Pkg.add("Distributed")  # for e.g. @spawn
Pkg.add("DifferentialEquations")  # for e.g. ODEProblem
Pkg.add("Random")  # for MersenneTwister()

# BioSequences before PhyloNetworks
# 
Pkg.add(Pkg.PackageSpec(url="https://github.com/BioJulia/BioSymbols.jl"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/BioJulia/BioGenerics.jl"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/BioJulia/BioSequences.jl"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/crsl4/PhyloNetworks.jl"))
using BioSymbols
using BioGenerics
using BioSequences
using PhyloNetworks

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
# Load dependencies of BioGeoJulia (if needed)
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
# FROM FRESH JULIA: Load the BioGeoJulia package
#######################################################
# Running startup here:
# '~/.julia/config/startup.jl
import Pkg
using Pkg
Pkg.rm("BioGeoJulia")
Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
using BioGeoJulia

# Run the tests directory
Pkg.test("BioGeoJulia")


#######################################################
# Re-run the tests directory
#######################################################
# First, commit to Master (quick), then:
Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
using BioGeoJulia
Pkg.test("BioGeoJulia")



#######################################################
# FROM FRESH JULIA: Load the BioGeoJulia package
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
Pkg.rm("BioGeoJulia")
Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
using BioGeoJulia

using BioGeoJulia.TrUtils
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs





#great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
great_ape_newick_string = "((chimp:1,human:1):1,gorilla:2);"
tr = readTopology(great_ape_newick_string)
tr

#res2 = construct_Res(tr)
res2 = construct_Res(tr, 5)


rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
trdf


n=10

# 4 states, no Q
# 3 tips in state 2, branch is 1 Mya long
birthRate = 0.222222
deathRate = 0.1
d_val = 0.01
e_val = 0.001
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
#Cijk_vals = repeat([birthRate], n)

#######################################################
# Add some j events
#######################################################
# Narrow sympatry
Cijk_vals = repeat([birthRate*(1.0/(1.0+2*j_val))], n)

# j events
Carray_ivals_leftJump = collect(1:n)
Carray_jvals_leftJump = reverse(collect(1:n))
Carray_kvals_leftJump = collect(1:n)
Cijk_vals_leftJump = repeat([birthRate*(j_val/(1.0+2*j_val))], n)

Carray_ivals_rightJump = collect(1:n)
Carray_jvals_rightJump = collect(1:n)
Carray_kvals_rightJump = reverse(collect(1:n))
Cijk_vals_rightJump = repeat([birthRate*(j_val/(1.0+2*j_val))], n)


Carray_ivals = vcat(Carray_ivals, Carray_ivals_leftJump, Carray_ivals_rightJump)
Carray_jvals = vcat(Carray_jvals, Carray_jvals_leftJump, Carray_jvals_rightJump)
Carray_kvals = vcat(Carray_kvals, Carray_kvals_leftJump, Carray_kvals_rightJump)
Cijk_vals = vcat(Cijk_vals, Cijk_vals_leftJump, Cijk_vals_rightJump)





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
tspan = (0.0, 1.2*trdf[tr.root,:node_age]) # 110% of tree root age

p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
p_Es_v5 = setup_MuSSE(3; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
p_Es_v5 = setup_MuSSE(4, birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)

# Anagenetic transition matrix
hcat(p_Es_v5.p_indices.Qarray_ivals, p_Es_v5.p_indices.Qarray_jvals, p_Es_v5.params.Qij_vals)
DataFrame(p_Es_v5.p_TFs.Qi_eq_i)
p_Es_v5.p_TFs.Qi_sub_i
p_Es_v5.p_TFs.Qj_sub_i

# Cladogenetic transition matrix
hcat(p_Es_v5.p_indices.Carray_ivals, p_Es_v5.p_indices.Carray_jvals, p_Es_v5.p_indices.Carray_kvals, p_Es_v5.params.Cijk_vals)
DataFrame(p_Es_v5.p_TFs.Ci_eq_i)
p_Es_v5.p_TFs.Ci_sub_i
p_Es_v5.p_TFs.Cj_sub_i
p_Es_v5.p_TFs.Ck_sub_i

p_Es_v5 = setup_MuSSE(2, birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)

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
#tspan = (0.0, 2.0*trdf[tr.root,:node_age]) # 110% of tree root age
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
solver_options.solver=Tsit5()
solver_options.abstol = 1.0e-6
solver_options.reltol = 1.0e-6

(total_calctime_in_sec, iteration_number) = iterative_downpass_nonparallel_ClaSSE_v5!(res, trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=solver_options, max_iterations=10^10);

res.likes_at_each_nodeIndex_branchTop
res.sum_likes_at_nodes
res.logsum_likes_at_nodes
log.(res.sum_likes_at_nodes[res.sum_likes_at_nodes.!=0.0])
sum(log.(res.sum_likes_at_nodes[res.sum_likes_at_nodes.!=0.0]))

total_calctime_in_sec
iteration_number


# Log sum of each branch bottom (cumulative)
log.(sum.(res.likes_at_each_nodeIndex_branchBot))
# 5-element Array{Float64,1}:
#  -0.22222200170689896
#  -0.22222200170689896
#  -2.170744399487735  
#  -0.4444440025007777 
#  -4.1192667987652865 

# The last number (total root lnL) equals
# LnLs1 -4.119266 -2.615189
# LnLs2 -4.812413 -2.615189
# LnLs3 -4.119266 -2.615189




#######################################################
# When deathRate = 0.0
#######################################################
# Sum of branch bottom likelihoods
sum(log.(sum.(res.likes_at_each_nodeIndex_branchBot))[1:3])
# This script: -2.615188402901533
# 
# Matches sum(lq)
# -2.615189
# ...from 
# /GitHub/phyBEARS/ex/groking_ClaSSE/BiSSE_branchlikes_w_pureYule_v6_WORKING.R

# Log sum of each branch bottom (cumulative)
log.(sum.(res.likes_at_each_nodeIndex_branchBot))
# 5-element Array{Float64,1}:
#  -0.22222200170689896
#  -0.22222200170689896
#  -2.170744399487735  
#  -0.4444440025007777 
#  -4.1192667987652865 

# The last number (total root lnL) equals
# LnLs1 -4.119266 -2.615189
# LnLs2 -4.812413 -2.615189
# LnLs3 -4.119266 -2.615189

# ...i.e., flat root probabilities, condition.surv=FALSE
# ...from 
# /GitHub/phyBEARS/ex/groking_ClaSSE/BiSSE_branchlikes_w_pureYule_v6_WORKING.R
res.likes_at_each_nodeIndex_branchBot
res.likes_at_each_nodeIndex_branchTop

res.normlikes_at_each_nodeIndex_branchBot
res.normlikes_at_each_nodeIndex_branchTop
sum(log.(sum.(res.likes_at_each_nodeIndex_branchBot))[1:3])
log.(sum.(res.likes_at_each_nodeIndex_branchBot))


res.sum_likes_at_nodes
res.logsum_likes_at_nodes



#######################################################
# When deathRate = 0.0
# After putting in normalization
#######################################################
res.likes_at_each_nodeIndex_branchBot
# 5-element Array{Array{Float64,1},1}:
#  [0.0, 0.8007375794916948]
#  [0.0, 0.8007375794916948]
#  [0.0, 0.8007375796021832]
#  [0.0, 0.6411806717956291]
#  [0.0, 0.0]               
# 
res.likes_at_each_nodeIndex_branchTop
# 5-element Array{Array{Float64,1},1}:
#  [0.0, 1.0]                
#  [0.0, 1.0]                
#  [0.0, 0.14248445111767713]
#  [0.0, 1.0]                
#  [0.0, 0.11409265462308325]
# 
res.normlikes_at_each_nodeIndex_branchBot
# 5-element Array{Array{Float64,1},1}:
#  [0.0, 0.0]
#  [0.0, 0.0]
#  [0.0, 0.0]
#  [0.0, 0.0]
#  [0.0, 0.0]
# 
res.normlikes_at_each_nodeIndex_branchTop
# 5-element Array{Array{Float64,1},1}:
#  [0.0, 0.0]
#  [0.0, 0.0]
#  [0.0, 1.0]
#  [0.0, 0.0]
#  [0.0, 1.0]
# 
sum(log.(sum.(res.likes_at_each_nodeIndex_branchBot))[1:3])
# -0.6666660049827136
# 
log.(sum.(res.likes_at_each_nodeIndex_branchBot))
# 5-element Array{Float64,1}:
#    -0.22222200170689896
#    -0.22222200170689896
#    -0.22222200156891567
#    -0.4444440025007777 
#  -Inf                  
# 
res.sum_likes_at_nodes
# 5-element Array{Float64,1}:
#  0.0                
#  0.0                
#  0.14248445111767713
#  0.0                
#  0.11409265462308325
# 
res.logsum_likes_at_nodes
# 5-element Array{Float64,1}:
#   0.0               
#   0.0               
#  -1.9485224001905719
#   0.0               
#  -2.1707444008464676
# 
res.likes_at_each_nodeIndex_branchBot
# 5-element Array{Array{Float64,1},1}:
#  [0.0, 0.8007375794916948]
#  [0.0, 0.8007375794916948]
#  [0.0, 0.8007375796021832]
#  [0.0, 0.6411806717956291]
#  [0.0, 0.0]               
# 
res.sum_likes_at_nodes
# 5-element Array{Float64,1}:
#  0.0                
#  0.0                
#  0.14248445111767713
#  0.0                
#  0.11409265462308325
# 
# 
-1.9485224001905719+-2.1707444008464676
# -4.119266801037039

sum(res.logsum_likes_at_nodes)
# -4.119266801037039

# Branch likes
-2.1707444008464676 + -0.4444440025007777
-2.6151884033472452
# Matches
# LnLs1 -4.119266 -2.615189
# LnLs2 -4.812413 -2.615189
# LnLs3 -4.119266 -2.615189



res.likes_at_each_nodeIndex_branchTop
log(res.likes_at_each_nodeIndex_branchTop[3][2])
log(res.likes_at_each_nodeIndex_branchTop[5][2])
# julia> log(res.likes_at_each_nodeIndex_branchTop[3][2])
# -1.9485224001905719
# 
# julia> log(res.likes_at_each_nodeIndex_branchTop[5][2])
# -2.1707444008464676

-2.1707444008464676- -1.9485224001905719
# -0.2222220006558957


-1.9485224001905719 - -0.2222220006558957
# -1.7263003995346762

# Matches:
# lq
# [1] -0.2222222 -0.2222222 -0.4444444  0.0000000 -1.7262996





#######################################################
# When deathRate = 0.1
#######################################################
# Sum of branch bottom likelihoods
sum(log.(sum.(res.likes_at_each_nodeIndex_branchBot))[1:3])
# This script: -2.981668534925785
# 
# Matches sum(lq)
# -2.948449
# ...from 
# /GitHub/phyBEARS/ex/groking_ClaSSE/BiSSE_branchlikes_w_BD_v6_WORKING.R

# Log sum of each branch bottom (cumulative)
log.(sum.(res.likes_at_each_nodeIndex_branchBot))
# 5-element Array{Float64,1}:
#    -0.3021505419023683
#    -0.3021505419023683
#    -2.3773674511210485
#    -0.5711384700773261
#  -Inf     

sum(log.(sum.(res.likes_at_each_nodeIndex_branchBot))[1:3])
# -2.981668534925785

sum(log.(sum.(res.likes_at_each_nodeIndex_branchBot))[1:4])
# -2.981668534925785

# The last number (total root lnL) is close to
# LnLs1 -4.452527 -2.948449
# LnLs2 -5.145674 -2.948449
# LnLs3 -4.452527 -2.948449

res.logsum_likes_at_nodes
# 5-element Array{Float64,1}:
#   0.0              
#   0.0              
#  -2.108379480581511
#   0.0              
#  -4.452584317975149

# The last number equals -4.452527

# ...i.e., flat root probabilities, condition.surv=FALSE
# ...from 
# /GitHub/phyBEARS/ex/groking_ClaSSE/BiSSE_branchlikes_w_BD_v6_WORKING.R
res.likes_at_each_nodeIndex_branchBot
res.likes_at_each_nodeIndex_branchTop

res.normlikes_at_each_nodeIndex_branchBot
res.normlikes_at_each_nodeIndex_branchTop
sum(log.(sum.(res.likes_at_each_nodeIndex_branchBot))[1:3])
log.(sum.(res.likes_at_each_nodeIndex_branchBot))


res.sum_likes_at_nodes
res.logsum_likes_at_nodes









#######################################################
# When deathRate = 0.1, q01 = 0.01, q10 = 0.01
#######################################################
# Matching:
# /GitHub/phyBEARS/ex/groking_ClaSSE/BiSSE_branchlikes_w_BDq01_v6_WORKING.R

sum(log.(res.sum_likes_at_nodes[res.sum_likes_at_nodes.!=0.0]))
# -4.457350176563324

# Matches
#            [,1]      [,2]
# LnLs1 -4.457691 -2.923784
# LnLs2 -5.150440 -2.923784
# LnLs3 -4.457492 -2.923784



res.likes_at_each_nodeIndex_branchBot
# 5-element Array{Array{Float64,1},1}:
#  [0.007351758755615916, 0.7384915796335971]
#  [0.007351758755615916, 0.7384915796335971]
#  [0.007673866446696122, 0.7633169235133778]
#  [0.011174271120372939, 0.5637645381772777]
#  [0.0, 0.0]  

log.(sum.(res.likes_at_each_nodeIndex_branchBot))
# 5-element Array{Float64,1}:
#    -0.2932397029911491
#    -0.2932397029911491
#    -0.2600788510672332
#    -0.5534916624604305
#  -Inf   

sum(res.likes_at_each_nodeIndex_branchBot[4])
# 0.5749388092976506


log(res.likes_at_each_nodeIndex_branchTop[3][2])
# Likelihoods at first internal node (Humans,Chimps)
# -2.1103695548994863

# Likelihoods at the bottom of the branches above the root
log.(sum.(res.likes_at_each_nodeIndex_branchBot))
log.(sum.(res.likes_at_each_nodeIndex_branchBot))[3:4]

# julia> log.(sum.(res.likes_at_each_nodeIndex_branchBot))
# 5-element Array{Float64,1}:
#    -0.2932397029911491
#    -0.2932397029911491
#    -0.2600788510672332
#    -0.5534916624604305
#  -Inf                 
# 
# julia> log.(sum.(res.likes_at_each_nodeIndex_branchBot))[3:4]
# 2-element Array{Float64,1}:
#  -0.2600788510672332
#  -0.5534916624604305

log(res.likes_at_each_nodeIndex_branchTop[3][2]) + sum(log.(sum.(res.likes_at_each_nodeIndex_branchBot))[3:4])
# -2.92394006842715
# ...all likelihoods above the root
# Matches!
#            [,1]      [,2]
# LnLs1 -4.457691 -2.923784
# LnLs2 -5.150440 -2.923784
# LnLs3 -4.457492 -2.923784

