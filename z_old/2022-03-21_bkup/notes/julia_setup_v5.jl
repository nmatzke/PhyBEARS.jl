#######################################################
# Julia hints
#######################################################
#######################################################
# Julia setup
# (e.g. on a new computer)
# Nick Matzke
# 2020-01-08
#######################################################

# After installing a new version:
txt = '
cd /usr/local/bin
ls -la j*
ln -s /Applications/Julia-1.3.1.app/Contents/Resources/julia/bin/julia julia
rm julia
ln -s /Applications/Julia-1.3.1.app/Contents/Resources/julia/bin/julia julia
julia --version
'


# Install libraries:
import Pkg
using Pkg  # for Pkg.add, Pkg.PackageSpec
#Pkg.resolve()  # reconciles packages (speeds up?)
using Random   # for e.g.: using Random; rng = MersenneTwister(1234); randn(rng, 10)
using Dates # for e.g. Dates.now(), DateTime


# https://timholy.github.io/Revise.jl/stable/
Pkg.add("Revise")
using Revise

Pkg.add("Combinatorics") # for combinations()


Pkg.add("FFMPEG") # for Plots
Pkg.add("Plots")
using FFMPEG
using Plots

Pkg.add("Fontconfig")
#Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaGraphics/Fontconfig.jl"))
Pkg.add("FiniteDiff") # Replaces DiffEqDiffTools
# https://discourse.julialang.org/t/finitediff-jl-fast-sparse-gradients-jacobians-hessians/32976
using FiniteDiff

Pkg.add("Rmath")
using Rmath
Pkg.build("Rmath")

Pkg.add("DiffEqBase")
Pkg.add("OrdinaryDiffEq")
Pkg.add("StochasticDiffEq")
Pkg.add("DifferentialEquations") # for ODEProblem
Pkg.add("Sundials")              # for CVODE_BDF
Pkg.add("ParameterizedFunctions") # for @ode_def (macro definition of ODE)
Pkg.add("ExponentialUtilities")
Pkg.add("EllipsisNotation")
Pkg.add("Sundials")
Pkg.add("LSODA")
Pkg.add("LinearAlgebra") # for "mul!" (a replacement for A_mul_B!)
Pkg.add("DiffEqBase")
Pkg.add("OrdinaryDiffEq")
Pkg.add("StochasticDiffEq")
# DIVAnd performs an n-dimensional variational analysis of arbitrarily located observations
Pkg.add("DIVAnd") # for ind2sub, deprecated after 0.6
Pkg.add("Compat") # for @test
Pkg.add("Random") # for Random
Pkg.add("LinearAlgebra") # for Diagonal
Pkg.add("ExponentialUtilities")  # for many, many DiffEq packages (?)
Pkg.build("FFTW")
Pkg.add("GPUArrays")
Pkg.add("BenchmarkTools")  # for @benchmark

Pkg.add("DataInterpolations") # for LinearInterpolation

using Statistics     # for mean()
using DifferentialEquations
using BenchmarkTools # for @benchmark
using LSODA          # for lsoda()
using Sundials       # for CVODE_BDF()

using GPUArrays      # for GPUArray
using DataInterpolations # for LinearInterpolation


# using SharedArrays   # for SharedArray
# using Distributed
# WARNING: using Distributed.@spawn in module Main conflicts with an existing identifier.


#Pkg.add(Pkg.PackageSpec(url = "https://github.com/BenJWard/Phylogenetics.jl"))
#Pkg.add("Phylo")        # Richard Reeve, demphasize
#Pkg.add("Phylogenies")  # BioJulia, emphasize

Pkg.build("Rmath")
Pkg.add("PhyloNetworks")        # most maintained, emphasize; for HybridNetwork
Pkg.add("PhyloPlots")						# plotting phylogenies
using PhyloNetworks
using PhyloPlots


great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr


Pkg.add("RCall")      # packaage to call R from within julia
Pkg.add("CSV")        # to read from / write to text files, e.g. csv files
Pkg.add("DataFrames") # to create & manipulate data frames
Pkg.add("StatsModels")# for regression formulas

using RCall
using CSV
using DataFrames
using StatsModels

Pkg.add("ODE")			# standard slow ODE solution
Pkg.add("DiffEqProblemLibrary") # library of problems
Pkg.add("ODEInterfaceDiffEq")

Pkg.add("DiffEqDevTools") # for TestSolution


Pkg.add("ODEInterface")  # for loadODESolvers, radau()

Pkg.add("StatsModels")     # for Missing
# Pkg.add("DataArrays")    # for NA # DEPRECATED
# https://github.com/JuliaStats/DataArrays.jl/issues/306
Pkg.add("DelimitedFiles")
Pkg.add("DataFrames")
Pkg.add("RDatasets")
Pkg.add("InvertedIndices")


# Adding packages from GitHub
# This WORKS:
# https://github.com/JuliaLang/Pkg.jl/issues/629 
#Pkg.add("SpecialMatrices") # for Strang() etc.; replaces "TeachingMatrices"
Pkg.add(Pkg.PackageSpec(url = "https://github.com/ChrisRackauckas/EllipsisNotation.jl")) # dependency for SpecialMatrices
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/SpecialMatrices.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl"))

# Needed for: CUDAdrv
# For Fontconfig ─ v0.3.0
# ERROR: MethodError: no method matching getindex(::Nothing, ::String)
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaGraphics/Fontconfig.jl"))


Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/LowRankApprox.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/MatrixFactorizations.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/MatrixDepot.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/ToeplitzMatrices.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/HierarchicalMatrices.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/LAPACK.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/InfiniteBandedMatrices.jl"))
#Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/BlockBandedMatrices.jl.jl"))
#Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/BandedMatrices.jl.jl"))

# https://gitter.im/JuliaDiffEq/Lobby?at=5c088a5d4808192b03f058a8
# Christopher Rackauckas @ChrisRackauckas Dec 07 2018 04:18
# @nmatzke hey, there's no tutorial yet, but you can take a look at how we use it
# for the in-progress benchmarks for IMEX and exponential integrators: 
# http://nbviewer.jupyter.org/github/MSeeker1340/DiffEqBenchmarks.jl/blob/expRK/StiffODE/Burgers.ipynb
#Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaMatrices/DiffEqBenchmarks.jl"))
#Pkg.add(Pkg.PackageSpec(path = "/drives/GDrive/z_code/julia_code/DiffEqBenchmarks.jl-master.zip"))
#Pkg.add(Pkg.PackageSpec(name="DiffEqBenchmarks", rev="master"))
#Pkg.add(Pkg.PackageSpec(path = "/drives/GDrive/z_code/julia_code/DiffEqBenchmarks"))

#Pkg.add(Pkg.PackageSpec(url = "file:///drives/GDrive/z_code/julia_code/DiffEqBenchmarks/"))

# Pkg.add(Pkg.PackageSpec(url = "file:///drives/GDrive/z_code/julia_code/DiffEqBenchmarks.jl-master.zip"))
# Pkg.add(Pkg.PackageSpec("file:///drives/GDrive/z_code/julia_code/DiffEqBenchmarks/src/"))
# Pkg.add(Pkg.PackageSpec("file:///drives/GDrive/z_code/julia_code/DiffEqBenchmarks"))
# Pkg.add(Pkg.PackageSpec("file:///drives/GDrive/z_code/julia_code/DiffEqBenchmarks.jl-master"))
# Pkg.add(Pkg.PackageSpec("file:///drives/GDrive/z_code/julia_code/DiffEqBenchmarks.jl.zip"))
# Pkg.add(Pkg.PackageSpec("file:///drives/GDrive/z_code/julia_code/DiffEqBenchmarks.jl-master.zip"))


# https://www.juliadiff.org/
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaDiff/ForwardDiff.jl"))
#Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaDiff/BackwardDiff.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaDiff/TaylorSeries.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaDiff/DualNumbers.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaDiff/HyperDualNumbers.jl"))

using ForwardDiff
using BackwardDiff
using TaylorSeries
using DualNumbers
using HyperDualNumbers


# For GPUs
# http://www.stochasticlifestyle.com/solving-systems-stochastic-pdes-using-gpus-julia/
Pkg.add("CLArrays") #FAILS
Pkg.add("GPUArrays")

Pkg.add("OpenCL")

Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaGPU/CLBLAS.jl"))
Pkg.add(Pkg.PackageSpec(url = "https://github.com/JuliaGPU/CLArrays.jl"))


Needs clBLAS
Which needs OpenCL, which might be replaced on mac by 
pocl

which needs hwloc

hwloc
https://www.open-mpi.org/software/hwloc/v2.0/
cd /Users/nickm/Downloads/hwloc-2.0.3
./configure --prefix=/usr
make
sudo make install


cd /Users/nickm/Downloads/pocl-1.3
cmake .


#Pkg.add("RecursiveArrayTools") # for ArrayPartition
#A = zeros(M,N); B  = zeros(M,N); C = zeros(M,N)
#u0 = ArrayPartition((A,B,C))


# GPU
Pkg.add("CUDAnative")
using CUDAnative



#######################################################
# Syntax notes
#######################################################
# eye() is replaced with:
# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
# 
# I()
# one()
# Matrix{Float64}(I, 2, 2)
# Diagonal(ones(3,3))
# 
# 4. https://discourse.julialang.org/t/sparse-to-dense-matrix/8380/10
# There is no "full()" in Julia 1.0. Use Matrix(...) instead to do the conversion.
# 


# Common and Baffling Errors
# Julia is an unusual language in many respects. It is also very novel both for the 
# (rather creative) developers and for its users. This is a collection of common 
# errors and/or weird and unintuitive outcomes.
# The Julia FAQ covers many other such issues.
# For now, the following collection is unordered. They will be ordered as a final step.
# http://julia.cookbook.tips/doku.php?id=baffling
#
# julia> 1.+2
# ERROR: syntax: invalid syntax "1.+"; add space(s) to clarify


#############################################
# Differential Equations Tutorial:
# Viewing the Notebooks Locally
#############################################
# https://github.com/JuliaDiffEq/DiffEqTutorials.jl
using Pkg
Pkg.add("CUDAdrv")  # dependency of IJulia, DiffEqTutorials
Pkg.add("Weave")    # dependency of IJulia, DiffEqTutorials
Pkg.add("IJulia") # Need to do this the first time to install IJulia!
Pkg.add(Pkg.PackageSpec(url="https://github.com/JuliaDiffEq/DiffEqTutorials.jl"))
using IJulia, DiffEqTutorials
pathof(DiffEqTutorials)

# Brings up a Notebook of Planetary Orbital Motions, etc.
# /Users/nmat471/.julia/packages/DiffEqTutorials/d6xDE/src/DiffEqTutorials.jl
notebook(dir = joinpath(dirname(pathof(DiffEqTutorials)),".."))


# Local link:
# http://localhost:8888/view/html/introduction/ode_introduction.html



#######################################################
# Online tutorials
#######################################################
# Julia differential equations
#
# http://docs.juliadiffeq.org/latest/
#
# Open Julia-1.0.app (this is actually version 1.02, November 2018)
# (some older matrix functions require Julia-0.6.app)

#
#
# This tutorial:
# https://nextjournal.com/sdanisch/an-intro-to-differentialequations
#
# Used in this YouTube video:
# https://www.youtube.com/watch?v=KPEqYtEd-zY



# Introducing the methods:

# ODE Solvers
# solve(prob::ODEProblem,alg;kwargs)
# https://github.com/JuliaDiffEq/DiffEqDocs.jl/blob/master/docs/src/solvers/ode_solve.md





# From:
# 7 Julia Gotchas and How to Handle Them
# http://www.stochasticlifestyle.com/7-julia-gotchas-handle/


#Pkg.resolve
using InvertedIndices # for Not
# https://github.com/JuliaLang/julia/issues/28276

using DifferentialEquations  # for ODEProblem
using DiffEqDevTools         # for TestSolution
using Sundials               # for CVODE_BDF
using EllipsisNotation
using DiffEqBase, OrdinaryDiffEq, StochasticDiffEq
using Plots; plotly()
using ParameterizedFunctions # for @ode_def (macro definition of ODE)
using Random # for Random
using LinearAlgebra # for Diagonal
using ExponentialUtilities  # Utility functions for exponential integrators 


# https://github.com/JuliaDiffEq/ODEInterfaceDiffEq.jl/issues/12
using ODEInterface # for loadODESolvers, radau(), rodas()

# DIVAnd performs an n-dimensional variational analysis of arbitrarily located observations
using DIVAnd; # for ind2sub

using Compat
using Compat.Test # for @test


# To avoid slow recompilation
# i.e. "Recompiling stale cache file" message
using Pkg
Pkg.resolve()




# Source some handy functions
include("/drives/Dropbox/_njm/__julia/julia4Rppl_v2.jl")

# Set working directory:
prob






#######################################################
# Speed tricks
#######################################################

# Optimizing differential equations code:
#https://juliabox.com/notebook/notebooks/tutorials/introductory-tutorials/broader-topics-and-ecosystem/intro-to-solving-diffeq-in-julia/3.OptimizingDiffEqCode.ipynb

# @inbounds (if you are sure of no index errors)
# @simd (using SIMD chip stuff)
# devectorize the stencil (???)
# .= where possible

# See what @simd is doing:
function mysum(a::Vector)
		 total = zero(eltype(a))
		 @simd for x in a
				 total += x
		 end
		 return total
 end
@code_native mysum(rand(Float64 , 100000))





#######################################################
# Write data to a table
#######################################################
#using DataArrays  # for NA
using DataFrames
using RDatasets
using DelimitedFiles # writedlm() writes the contents of an object to a text file, 
										 # readdlm() reads the data from a file into an array:

# Input/Output:
# https://syl1.gitbook.io/julia-language-a-concise-tutorial/language-core/input-output
open("afile.txt", "w") do f  # "w" for writing
  for i in 1:10
  	write(f, i, "\n")   # \n for newline
  end
end

open("afile.txt", "r") do f
  for ln in eachline(f)
    println(ln)
  end
end


# For-loops
niter = 10
for i in 1:niter

end





# Output to screen

using BenchmarkTools # for @benchmark
using Statistics     # for mean

function example_add(a,b)
	c = a+b
end

x = @benchmark example_add(1, 2)

x

# Screen output -- how do I capture to a string or other variable?
#
# BenchmarkTools.Trial: 
#   memory estimate:  0 bytes
#   allocs estimate:  0
#   --------------
#   minimum time:     0.025 ns (0.00% GC)
#   median time:      0.029 ns (0.00% GC)
#   mean time:        0.029 ns (0.00% GC)
#   maximum time:     6.074 ns (0.00% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1000

minimum(x.times)
mean(x.times)
maximum(x.times)




#
# Landis & Quintero 2019
# https://landislab.github.io/assets/research/pdf/Quintero_Landis_2019_bioRxiv_biotic_interactions.pdf
#
# Software.— We denote this model as “TRIBE” (which stands for “Trait and Range
# 414 Interspecific Biogeographic Evolution”) and implement it in a new open source package
# 415 named “Tapestree” (https://github.com/ignacioq/Tapestree.jl) that we wrote in
# 416 Julia (Bezanson et al. 2017). This software makes available the tribe() function for
# 417 inference and the simulate tribe() for simulations given a fixed tree. We note that, in
# 418 the software, we allow the user to fix to 0 any or all of the parameters governing the effect
# 419 of biotic interactions (i.e., ωx, ω0, & ω1)


Pkg.add(Pkg.PackageSpec(url = "https://github.com/ignacioq/Tapestree.jl")) # installing TapesTree

using Tapestree
finches_tree_file = "/drives/Dropbox/_njm/__julia/Tapestree.jl-master/data/finches_rescaled.tre"
finches_data_file = "/drives/Dropbox/_njm/__julia/Tapestree.jl-master/data/finches_pca1.txt"
out_file  = *(homedir(),"/TapesTree_ex_inf_v1.txt")

# Run the tribe() (TRIBE: Trait and Range Interspecific Biogeographic Evolution) model:
tribe(finches_tree_file, finches_data_file, out_file)


min_dt  = 0.01                       # a float describing the percentage of tree height allowed for discretization (lower values are more precise but take longer).
niter   = 10_000                     # an integer for the number of iterations.
nburn   = 5_000                      # an integer for the number of iterations in the adaptive burn-in phase.
nthin   = 100                        # an integer for the iteration sampling frequency.
saveXY  = (true, 1_000)              # a tuple of length 2: first is a boolean to save (or not) data augmented histories, second an integer for sampling frequency.
saveDM  = (true, 1_000)              # a tuple of length 2: first is a boolean to save (or not) data augmented deterministic effects, second an integer for sampling frequency.
ωxprior = (0.,10.)                   # a tuple of length 2 for the normal prior of ωx, first the mean, second the variance.
ω1prior = (0.,10.)                   # a tuple of length 2 for the normal prior of ω1, first the mean, second the variance.
ω0prior = (0.,10.)                   # a tuple of length 2 for the normal prior of ω0, first the mean, second the variance.
σ²prior = 1e-1                       # a float for the mean of the exponential prior for σ².
λprior  = 1e-1                       # a float for the mean of the exponential prior for both λs.
weight  = (0.15,0.05,0.02,0.02,5e-3) # a tuple of length 5 specifying the probabilities to update σ², ωx, ω1 & ω0, and λ1 & λ0 respectively.
λ1i     = 1.0                        # a float for the starting value for λ1.
λ0i     = 0.5                        # a float for the starting value for λ0.
ωxi     = 0.0                        # a float for the starting value for ωx.
ω1i     = 0.0                        # a float for the starting value for ω1.
ω0i     = 0.0                        # a float for the starting value for ω0.
fix_ωx  = false                      # a boolean to make inference without ωx.
fix_ω1  = false                      # a boolean to make inference without ω1.
fix_ω0  = false                      # a boolean to make inference without ω0.


# Simulation
# Specify the path to the phylogenetic tree (in a format that ape can read):
finches_tree_file = "/drives/Dropbox/_njm/__julia/Tapestree.jl-master/data/finches_rescaled.tre"


x_init  = 0.0
n_areas = 6
tip_values, tip_areas, tree, bts = simulate_tribe(x_init, n_areas, finches_tree_file)


# Further options for simulate_tribe() are
# ωx       = 0.0   # a float for simulated value of ωx.
# σ²       = 0.5   # a float for simulated value of σ².
# ω1       = 0.0   # a float for simulated value of ω1.
# ω0       = 0.0   # a float for simulated value of ω0.
# λ1       = 0.5   # a float for simulated value of λ1.
# λ0       = 0.2   # a float for simulated value of λ0.
# const_δt = 1e-4  # a float for the delta t used to approximate the simulation (lower values are more accurate but a

out_file  = *(homedir(),"/TapesTree_ex_sim_v1.txt")
tribe(tip_values, tip_areas, tree, bts, out_file)










# From TapesTree.jl/data_initializer

"""
    read_tree(tree_file::String)

Function to read a tree using `RCall`
to call **ape** tree reading capabilities. 
"""
function read_tree(tree_file::String)

  str = reval("""
              library(\"ape\")
              tree     <- read.tree('$tree_file') 
              tree     <- reorder(tree)
              edge     <- .subset2(tree,'edge')
              Nnode    <- .subset2(tree,'Nnode')
              tiplabel <- .subset2(tree,'tip.label')
              edlength <- .subset2(tree,'edge.length')
              list(edge,Nnode,tiplabel,edlength)
              """)

  edge     = rcopy(str[1])
  edge     = convert(Array{Int64},edge)
  Nnode    = rcopy(str[2])
  Nnode    = convert(Int64,Nnode)
  tiplabel = rcopy(str[3])
  edlength = rcopy(str[4])
  edlength = convert(Array{Float64},edlength)

  tree = rtree(edge, edlength, tiplabel, Nnode)

  brtimes = reval("""
                  brtimes <- branching.times(tree)
                  """)

  brtimes = rcopy(brtimes)

  return tree, brtimes
end








"""
Immutable type of an R tree `phylo` object type.
"""
struct rtree
  ed  ::Array{Int64,2}
  el  ::Array{Float64,1}
  tlab::Array{String,1}
  nnod::Int64
end



"""
    read_data(tree_file::String, data_file::String; delim::Char = '\t', eol::Char = '\r')

Read a phylogenetic tree using **ape** package in R through 
`RCall` and the data file with the trait and biogeographic information.
"""
function read_data(tree_file::String,
                   data_file::String)

  tree, bts = read_tree(tree_file)

  tip_labels = Dict(i => val for (val,i) = enumerate(tree.tlab))

  data = DelimitedFiles.readdlm(data_file)

  if size(data,1) != (tree.nnod + 1)
    data = DelimitedFiles.readdlm(data_file, '\t', '\r')
  end

  if size(data,1) != (tree.nnod + 1)
    data = DelimitedFiles.readdlm(data_file, '\t', '\n')
  end

  if size(data,1) != (tree.nnod + 1) 
    error("Data file cannot be made of the right dimensions.\n Make sure the data file has the same number of rows as tips in the tree")
  end

  data_tlab   = convert(Array{String,1}, data[:,1])
  data_values = convert(Array{Float64,1},data[:,2])
  data_areas  = convert(Array{Int64,2},  data[:,3:end])

  # create dictionaries
  tip_areas = Dict(tip_labels[val] => data_areas[i,:] 
                   for (i,val) = enumerate(data_tlab))

  tip_values = Dict(tip_labels[val] => data_values[i] 
                    for (i,val) = enumerate(data_tlab))

  return tip_values, tip_areas, tree, bts
end

