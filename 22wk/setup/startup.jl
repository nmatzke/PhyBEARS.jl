










# 
print("\n")
print("Your '~/.julia/config/startup.jl' file is running some commands...")
print("\n")


# Perhaps useful for Revise, but couldn't get Revise() to work on my Mac
print("Adding package developing directories to LOAD_PATH (JULIA_LOAD_PATH externall)...")
print("\n")
push!(LOAD_PATH, "/GitHub/PhyBEARS.jl")
print("Printing current LOAD_PATH:\n")
print(LOAD_PATH)
print("\n\n")


# print("Setting up Revise...")
# 
# atreplinit() do repl
#     try
#         @eval using Revise          # for recompiling new code in development
#         @async Revise.wait_steal_repl_backend()
#     catch e
#         @warn(e.msg)
#     end
# end
# print("done.")
# print("\n")
# print("\n")




# Loading basic (fast, don't need compilation) default packages
print("Loading basic (fast, don't need compilation) default packages...\n")

using Pkg
using DelimitedFiles
using CSV
using DataFrames			# for e.g. DataFrame()
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using LSODA          # for lsoda()
using BenchmarkTools # for @time
using InvertedIndices # for Not
using SpecialFunctions	# for logfactorial
using Optim                 # for e.g. L-BFGS-B Maximum Likelihood optimization

print("...done loading basic default packages.\n\n")



# print("\n")
# print("Showing installed packages...")
# show(stdout, "text/plain", Pkg.installed())
# print("\n")


print("Showing environment, and installed packages, with Pkg.status()...")
Pkg.status()
#show(stdout, "text/plain", Pkg.status())
print("\n")

 

#print("NOTE: These packages, in Nick's default setup, are stored via Fezzik.brute_build_julia()")
#print("\n")

print("Your '~/.julia/config/startup.jl' file has finished running.")
print("\n")



# Loading slow packages
print("\n")
print("Loading slow packages...\n")
print("\n")
print("DifferentialEquations...\n")
using DifferentialEquations # for ODEProblem
print("PhyloBits...\n")
Pkg.add(PackageSpec(path="/GitHub/PhyloBits.jl"))

using PhyloBits

using PhyloBits.PNtypes						# Types: Tree, Node, Edge etc.
using PhyloBits.PNmanipulateNet
using PhyloBits.PNreadwrite				# Helper functions
using PhyloBits.PNreadwrite 			# Reading and writing trees; readTopology
using PhyloBits.PNdescriptive			# show() commands for trees etc.
using PhyloBits.TrUtils						# basic utility functions, R-like functions
using PhyloBits.TreeTable					# for prt() tree tables (tree dataframes, trdfs), bd_liks(), etc.

print("Unloading and re-loading PhyBEARS...\n")

# Pkg.rm("PhyloBits")
# Pkg.rm(PackageSpec(name="PhyloBits", uuid="d8844a15-86e7-4f2e-bb1d-495fc7c0502a"))
# Pkg.add(PackageSpec(path="/GitHub/PhyloBits.jl"))
# Pkg.add(url="https://github.com/nmatzke/PhyloBits.jl")

#Pkg.rm(PackageSpec(name="PhyBEARS", uuid="7876af07-990d-54b4-ab0e-23690620f79a"))
Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl"))
#Pkg.add(url="https://github.com/nmatzke/PhyBEARS.jl")

# Activate PhyBEARS
#if isfile("Project.toml") && isfile("Manifest.toml")
#    Pkg.activate(".")
#end

using PhyBEARS

#using PhyloBits.TrUtils
#using PhyloBits.TreeTable
using PhyBEARS.BGExample
using PhyBEARS.MaxentInterp
using PhyBEARS.StateSpace
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.TreePass
using PhyBEARS.Flow
using PhyBEARS.Gmaps
using PhyBEARS.Optimizers
using PhyBEARS.TimeDep


print("\n")
print("...done loading slow packages.\n")
print("\n")

print("Checking if PhyloBits is working:\n")
print("\n")

# Check if updates are working
PhyloBits.TrUtils.hello_world_TrUtils()


print("\nChecking if PhyloBits.PNreadwrite.readTopology is working:\n")
print("\n")

# Psychotria tree
Psychotria_tree = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")

show(Psychotria_tree)




