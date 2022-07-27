
#######################################################
# SKIP the Fezzik block, if you have already used it
#######################################################
# FezzikAutoGenStart
# to remove remove entire block
try
    using Fezzik
    try
        Fezzik.trace()
    catch err
        @info "Something went wrong" err
    end
catch e
    try
        using Pkg
        Pkg.add("Fezzik")
        import Fezzik
        try
            Fezzik.trace()
        catch err
            @info "Something went wrong" err
        end
    catch err
        @info "Something went wrong" err
    end
end

# FezzikAutoGenEnd




#######################################################
# PASTE THE BELOW INTO YOUR startup.jl
#######################################################


# Nick Matzke startup, 2020-08-30
# No active Fezzik, since compiled sysimg is already done!
print("\n")
print("Your '~/.julia/config/startup.jl' file is running some commands...")
print("\n")


print("If you want to run some basic BioGeoJulia commands, try:\n")
print("""include("/GitHub/BioGeoJulia.jl/notes/startup_basic_tests_v1.jl")""")
print("\n\n")

# Perhaps useful for Revise, but couldn't get Revise() to work on my Mac
print("Adding package developing directories to LOAD_PATH (JULIA_LOAD_PATH externall)...")
print("\n")
push!(LOAD_PATH, "/GitHub/BioGeoJulia.jl")
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

# Check if updates are working
TrUtils.hello_world_TrUtils()


print("\n")
print("...done loading slow packages.\n")
print("\n")

print("If you want to run some basic BioGeoJulia commands, try:\n")
print("""include("/GitHub/BioGeoJulia.jl/notes/startup_basic_tests_v1.jl")""")
print("\n\n")

