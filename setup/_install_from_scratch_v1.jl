
# Reinstall from scratch:

import Pkg
using Pkg

Pkg.add("Compat")
Pkg.add("Preferences")

# PhyloBits dependencies
Pkg.add("DelimitedFiles")
Pkg.add("DataFrames")
Pkg.add("UUIDs")
Pkg.add("CSV")
Pkg.add("FilePathsBase")
Pkg.add("BenchmarkTools")

Pkg.add("SciMLBase")

Pkg.add("Distributed")
Pkg.add("Hwloc")
Pkg.add("Printf")
Pkg.add("SpecialFunctions")
Pkg.add("StaticArrays")
Pkg.add("StatsBase")

# PhyBEARS.jl dependencies
Pkg.add("Interpolations")
Pkg.add("RCall")
Pkg.add("Combinatorics")
Pkg.add("Dates")
Pkg.add("DiffEqNoiseProcess")
Pkg.add("Distributions")
Pkg.add("DoubleFloats")
Pkg.add("Interpolations")
Pkg.add("InvertedIndices")
Pkg.add("JLD2")
Pkg.add("LSODA")
Pkg.add("LinearAlgebra")
Pkg.add("NLopt")
Pkg.add("ODE")
Pkg.add("Optim")
Pkg.add("Polynomials")
Pkg.add("Random")
Pkg.add("Sundials")
Pkg.add("Test")
Pkg.add("DifferentialEquations")


# Loading them
using Compat
using Preferences

# PhyloBits dependencies
using DelimitedFiles
using DataFrames
using UUIDs
using CSV
using FilePathsBase
using BenchmarkTools

using SciMLBase

using Distributed
using Hwloc
using Printf
using SpecialFunctions
using StaticArrays
using StatsBase

# PhyBEARS.jl dependencies
using Combinatorics
using Dates
using DiffEqNoiseProcess
using Distributions
using DoubleFloats
using Interpolations
using InvertedIndices
using JLD2
using LSODA
using LinearAlgebra
using NLopt
using ODE
using Optim
using Polynomials
using Random
using Sundials
using Test
using DifferentialEquations

# Update and resolve
Pkg.update()
Pkg.resolve()
Pkg.gc()





using Pkg; Pkg.add(PackageSpec(path="/Users/nmat471/GitHub/PhyloBits.jl"))
# OR: 
# Pkg.add(url="https://github.com/nmatzke/PhyloBits.jl")

using PhyloBits

using PhyloBits.PNtypes						# Types: Tree, Node, Edge etc.
using PhyloBits.PNmanipulateNet
using PhyloBits.PNreadwrite				# Helper functions
using PhyloBits.PNreadwrite 			# Reading and writing trees; readTopology
using PhyloBits.PNdescriptive			# show() commands for trees etc.
using PhyloBits.TrUtils						# basic utility functions, R-like functions
using PhyloBits.TreeTable					# for prt() tree tables (tree dataframes, trdfs), bd_liks(), etc.

print("Unloading and re-loading PhyBEARS...\n")


# Activate PhyBEARS
#if isfile("Project.toml") && isfile("Manifest.toml")
#    Pkg.activate(".")
#end
#Pkg.rm(PackageSpec(name="PhyBEARS", uuid="7876af07-990d-54b4-ab0e-23690620f79a"))
using Pkg; Pkg.add(PackageSpec(path="/Users/nmat471/GitHub/PhyBEARS.jl"))
# OR: 
# Pkg.add(url="https://github.com/nmatzke/PhyBEARS.jl")
using PhyBEARS

#using PhyloBits.TrUtils
#using PhyloBits.TreeTable
using PhyBEARS.BGExample
using PhyBEARS.TimeDep
using PhyBEARS.MaxentInterp
using PhyBEARS.StateSpace
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.TreePass
using PhyBEARS.Flow
using PhyBEARS.Gmaps
using PhyBEARS.Optimizers

Pkg.update()
Pkg.resolve()
Pkg.gc()
