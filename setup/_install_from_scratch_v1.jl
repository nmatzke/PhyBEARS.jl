
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
using BenchmarkTools")

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

