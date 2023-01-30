#######################################################
# Installing Julia
#######################################################

# Download and install the latest version of Julia
# from here: 
# https://julialang.org/downloads/

# Run this at the command line to see that julia is installed:
julia --version

# If you get an error, you may have to symbolic-link to your 
# julia executable.  E.g., on a Mac:
cd /usr/local/bin/
ls -la j*
ln -s /Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia julia
ls -la j*

# Restart your command-line terminal, then:
julia --version


#######################################################
# Installing PhyBEARS
#######################################################
# This will take awhile the first time!
using Pkg
Pkg.add("BenchmarkTools") # for @time
Pkg.add("CSV")
Pkg.add("Combinatorics") # for combinations()
Pkg.add("DataFrames")			# for e.g. DataFrame()
Pkg.add("Dates")						# for e.g. DateTime, Dates.now()
Pkg.add("DelimitedFiles")
Pkg.add("DiffEqNoiseProcess")

# This is a huge package
Pkg.add("DifferentialEquations")

Pkg.add("Distributed")
Pkg.add("Distributions")
Pkg.add("DoubleFloats")
Pkg.add("Hwloc")
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
Pkg.add("SciMLBase")
Pkg.add("SpecialFunctions")
Pkg.add("StatsBase")
Pkg.add("Sundials")
Pkg.add("Test")

# Pre-compile the huge DifferentialEquations package
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

# This one can be slow
using DifferentialEquations # for ODEProblem


# Install "PhyloBits" from GitHub. This contains the
# tree-reading functionality of PhyloNetworks, but
# without all the dependencies.  It also contains
# some accessory functions by Matzke, like "prt()"
print("PhyloBits...\n")
#Pkg.add(PackageSpec(path="/GitHub/PhyloBits.jl"))
Pkg.add(url="https://github.com/nmatzke/PhyloBits.jl")

# Install "PhyBEARS" from GitHub
# Some example uninstall commands are commented out;
# useful if we change & reinstall the package

# Pkg.rm("PhyloBits")
# Pkg.rm(PackageSpec(name="PhyloBits", uuid="d8844a15-86e7-4f2e-bb1d-495fc7c0502a"))

#Pkg.rm(PackageSpec(name="PhyBEARS", uuid="7876af07-990d-54b4-ab0e-23690620f79a"))
#Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl"))
Pkg.add(url="https://github.com/nmatzke/PhyBEARS.jl")

# The "resolve" step checks everything
Pkg.resolve()

# The "status" command lists all installed packages
print("Showing environment, and installed packages, with Pkg.status()...")
Pkg.status()


# See if you can load a phylogeny
using PhyloBits
using PhyloBits.PNtypes						# Types: Tree, Node, Edge etc.
using PhyloBits.PNmanipulateNet
using PhyloBits.PNreadwrite				# Helper functions
using PhyloBits.PNreadwrite 			# Reading and writing trees; readTopology
using PhyloBits.PNdescriptive			# show() commands for trees etc.
using PhyloBits.TrUtils						# basic utility functions, R-like functions
using PhyloBits.TreeTable					# for prt() tree tables (tree dataframes, trdfs), bd_liks(), etc.

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

# See if a PhyBEARS command works
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
using PhyBEARS.TimeDep

include_null_range = true;
numstates_from_numareas(4, 4, include_null_range)
numstates_from_numareas(4, 3, include_null_range)
numstates_from_numareas(4, 2, include_null_range)
numstates_from_numareas(4, 1, include_null_range)

include_null_range = false;
numstates_from_numareas(4, 4, include_null_range)
numstates_from_numareas(4, 3, include_null_range)
numstates_from_numareas(4, 2, include_null_range)
numstates_from_numareas(4, 1, include_null_range)


#######################################################
# FINALLY: To avoid having to do all of the "using" commands,
# I recommend you copy my "startup.jl" file into your
# ~/.julia directory (~ = your default home directory)
#
# This ~/.julia/startup.jl file will run every time
# you run julia, unless you say
# julia --startup-file=no
#######################################################




