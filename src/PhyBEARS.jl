module PhyBEARS
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

print("\nPhyBEARS: loading PhyBEARS.jl.\n")

using DelimitedFiles	# for writedlm
using Distributed 	# for workers, spawnat :any, etc.
using Hwloc					# for Hwloc.num_physical_cores(), Hwloc.num_virtual_cores()
using PhyloBits			# for prt(), Node/HybridNetwork, etc.

# List each PhyBEARS code file here
# NOTE: LOAD THE DEPENDENCY .jl FILES *FIRST*, or you get "not recognized" errors
include("BGExample.jl")			# default examples
#include("TrUtils.jl")			# basic utility functions 
include("TimeDep.jl")				# time-dependency in Qmat and Cmat
include("MaxentInterp.jl")	# preconstructed interpolator for weighting rangesize of smaller daughter
#include("TreeTable.jl")			# for prt() tree tables (DFs), bd_liks(), etc.
include("StateSpace.jl")	# set up lists of areas and states (geographic ranges)
include("SSEs.jl")				# SSE calculations with various amounts of speed optimization
include("Parsers.jl")			# Parsers to read e.g. geography file
include("TreePass.jl")		# downpass and uppass through the phylogeny; prt() etc.
include("ModelLikes.jl")		# likelihood calculations
include("Flow.jl")		# downpass and uppass through the phylogeny
include("Gmaps.jl")		# Gmaps arrays etc.
include("Optimizers.jl")


export hello_PhyBEARS, add_one_PhyBEARS


print("For multithreading purposes, Threads.nthreads() = ")
print(Threads.nthreads())
print("\n")
print("For multiprocessor purposes\n:")
print("Hwloc.num_physical_cores() = ")
print(Hwloc.num_physical_cores())
print("\n")
print("Hwloc.num_virtual_cores() = ")
print(Hwloc.num_virtual_cores())
print("\n")
print("Number of workers turned on in this session: Distributed.nworkers() = ")
print(Distributed.nworkers())
print("\n")
print("List of workers turned on in this session: Distributed.workers() = ")
print(Distributed.workers())
print("\n")








"""
# Local nstallation:
using Pkg
Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl")); using PhyBEARS

# GitHub Installation
using Pkg
Pkg.add(Pkg.PackageSpec(url="https://github.com/nmatzke/PhyBEARS.jl")); using PhyBEARS

# Add to another package:
cd("/GitHub/PhyBEARS.jl")
]
activate .
add https://github.com/BioJulia/BioJuliaRegistry.git

# Loading:
Pkg.instantiate()
using PhyloBits

Pkg.instantiate()
using PhyBEARS
using PhyBEARS.BGExample
using PhyBEARS.MaxentInterp
using PhyBEARS.StateSpace
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.TreePass
using PhyBEARS.Flow
using PhyBEARS.Gmaps
using PhyBEARS.Optimizers
"""


#######################################################
# Put functions here
#######################################################
"""
    hello(who::String)

Return "PhyBEARS says, hi `who`".
"""
hello_PhyBEARS(who::String) = "hello_PhyBEARS() says '$who'"

"""
    add_one_PhyBEARS(x::Number)

Return `x + 1`.
"""
add_one_PhyBEARS(x::Number) = x + 1

end # Ends the module command
