print("\n\n'~/.julia/startup.jl' is loading packages precompiled with Fezzik.\n")


print("Loading Fezzik...\n")

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


# startup.jl
#
print("Loading the fast ones...\n")

using Pkg
#using Revise
using LSODA          # for lsoda()
using BenchmarkTools # for @time
using InvertedIndices # for Not
using Distributed     # for @spawn
using Distributions  # for quantile
using Convex				 # for Convex.entropy(), maximize()
using SCS						 # for SCSSolve, solve (maximize(entropy()))
using Combinatorics  # for e.g. combinations()
using DataFrames     # for e.g. DataFrame()

using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()

print("Loading the slow ones...\n")
print("Loading DifferentialEquations...\n")
using DifferentialEquations # for ODEProblem (THE SLOWEST ONE)
print("Loading PhyloNetworks...\n")
using PhyloNetworks
print("...done.\n")


