






# Nick Matzke startup, 2020-08-30
# No active Fezzik, since compiled sysimg is already done!
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
using DataFrames			# for e.g. DataFrame()
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using LSODA          # for lsoda()
using BenchmarkTools # for @time
using InvertedIndices # for Not
using SpecialFunctions	# for logfactorial

print("...done loading basic default packages.\n\n")



# print("\n")
# print("Showing installed packages...")
# show(stdout, "text/plain", Pkg.installed())
# print("\n")


print("Showing environment, and installed packages...")
Pkg.status()
#show(stdout, "text/plain", Pkg.status())
print("\n")

 

print("NOTE: These packages, in Nick's default setup, are stored via Fezzik.brute_build_julia()")
print("\n")

print("Your '~/.julia/config/startup.jl' file has finished running.")
print("\n")



# Loading slow packages
print("\n")
print("Loading slow packages...\n")
print("\n")
print("DifferentialEquations...\n")
using DifferentialEquations # for ODEProblem
print("PhyloNetworks...\n")
using PhyloNetworks

# Check if PhyloNetworks is working
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
testtr = readTopology(great_ape_newick_string)
testtr
print("\n")
print("...done loading slow packages.\n")
print("\n")


print("Unloading and re-loading PhyBEARS...\n")

#Pkg.rm("PhyBEARS")
#Pkg.rm(PackageSpec(name="PhyBEARS", uuid="7876af07-990d-54b4-ab0e-23690620f79a"))
#Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl"))

# Activate PhyBEARS
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

using PhyBEARS
using PhyBEARS.MaxentInterp
using PhyBEARS.TrUtils
using PhyBEARS.StateSpace
using PhyBEARS.TreeTable
using PhyBEARS.TreePass
using PhyBEARS.SSEs

# Check if updates are working
PhyBEARS.hello_PhyBEARS("~/.julia/startup.jl has finished.")






