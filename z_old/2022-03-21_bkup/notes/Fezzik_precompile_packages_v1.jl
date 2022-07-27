

# Fezzik_setup
# https://github.com/TsurHerman/Fezzik

# FIRST TIME:
# "enables automatic tracing of compiler activity through adding itself to the startup.jl file."
using Pkg
#Pkg.add(PackageSpec(url="https://github.com/TsurHerman/Fezzik"))
using Fezzik
Fezzik.auto_trace()
# NOTE: To clear old Fezziks, use Fezzik.revert()


# CLOSE, THEN OPEN JULIA:
using Fezzik
using Pkg

# The last 2 seem to need blacklisting on Julia 1.4 (kept the same on Julia 1.5)
Fezzik.blacklist("BioGeoJulia")
Fezzik.blacklist("AbstractPlotting")
Fezzik.blacklist("Interact")

# Faster/smaller packages
#using Revise
using Optim       # for ML optimization
using GR          # dependency of Plots
using LSODA          # for lsoda()
using BenchmarkTools # for @time
using InvertedIndices # for Not
using DataFrames      # for e.g. DataFrame()
using Distributed     # for @spawn
using Distributions  # for quantile
using Convex				 # for Convex.entropy(), maximize()
using SCS						 # for SCSSolve, solve (maximize(entropy()))
using Combinatorics  # for e.g. combinations()
using DataFrames     # for e.g. DataFrame()
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()

using FixedPointNumbers  # for Statistics
using GeometryTypes  # for Plots
using GeometryBasics  # for Plots

using Statistics			# for e.g. mean()
using RCall       # for df_to_Rdata, reval, g = globalEnv

# For SSE stuff
using Sundials        # for CVODE_BDF(linear_solver=:GMRES)
using OrdinaryDiffEq
using DiffEqDevTools
using ODEInterfaceDiffEq
using ODE

# Large and slow...
using GR							# dependency for plots
using Plots						# basic plots
using DifferentialEquations # for ODEProblem (THE SLOWEST ONE)
using PhyloNetworks		# for phylogenetic trees
Pkg.resolve()					# Resolve package dependencies

Fezzik.brute_build_julia()

# Restart:
# remove itself from the startup.jl and deletes previous traces
using Fezzik
Fezzik.auto_trace(false)
exit()

# Then go edit ~/.julia/config/startup.jl

# Activating environment at `~/.julia/environments/v1.3/Project.toml`
# [ArrayInterface] already loaded
# [RecursiveArrayTools] already loaded
# [DiffEqBase] already loaded
# [FiniteDiff] already loaded
# [Compat] already loaded
# [DiffEqNoiseProcess] already loaded
# [Parsers] already loaded
# [RandomNumbers] already loaded
# [Base] already loaded
# [OpenBLAS_jll] already loaded
# [Pkg] already loaded
# [NLopt] already loaded
# [Setfield] already loaded
# [Arpack_jll] already loaded
# [CSV] already loaded
# [Libdl] already loaded
# [REPL] already loaded
# [SparseDiffTools] already loaded
# [LSODA] already loaded
# [FilePathsBase] already loaded
# [DelimitedFiles] already loaded
# [Logging] already loaded
# [PackageCompiler] already loaded
# [DocStringExtensions] already loaded
# [Sundials] already loaded
# [Core] already loaded
# [SpecialFunctions] already loaded
# [Requires] already loaded
# [Random] already loaded
# [Fezzik] already loaded
# [Rmath] already loaded
# Activating environment at `~/.julia/environments/v1.3/Project.toml`
# [ Info: used 443 precompile statements
# out_file = "/Users/nmat471/.julia/packages/Fezzik/BRQC8/precomp.jl"
# 
# 
# 
#  Compiling...
# Activating environment at `~/.julia/environments/v1.3/Project.toml`
# trying to import ArrayInterface
# failed to import ArrayInterface deffering
# trying to import RecursiveArrayTools
# failed to import RecursiveArrayTools deffering
# trying to import DiffEqBase
# failed to import DiffEqBase deffering
# trying to import FiniteDiff
# failed to import FiniteDiff deffering
# [Compat] already loaded
# trying to import DiffEqNoiseProcess
# failed to import DiffEqNoiseProcess deffering
# trying to import Parsers
# failed to import Parsers deffering
# trying to import RandomNumbers
# failed to import RandomNumbers deffering
# [Base] already loaded
# trying to import OpenBLAS_jll
# failed to import OpenBLAS_jll deffering
# [Pkg] already loaded
# trying to import NLopt
# failed to import NLopt deffering
# [Setfield] already loaded
# trying to import Arpack_jll
# failed to import Arpack_jll deffering
# trying to import CSV
# failed to import CSV deffering
# [Libdl] already loaded
# [REPL] already loaded
# trying to import SparseDiffTools
# failed to import SparseDiffTools deffering
# trying to import LSODA
# failed to import LSODA deffering
# [FilePathsBase] already loaded
# [DelimitedFiles] already loaded
# [Logging] already loaded
# [PackageCompiler] already loaded
# [DocStringExtensions] already loaded
# trying to import Sundials
# failed to import Sundials deffering
# [Core] already loaded
# trying to import SpecialFunctions
# failed to import SpecialFunctions deffering
# [Requires] already loaded
# [Random] already loaded
# [Fezzik] already loaded
# trying to import Rmath
# failed to import Rmath deffering
# [ArrayInterface] already loaded
# [RecursiveArrayTools] already loaded
# [DiffEqBase] already loaded
# [FiniteDiff] already loaded
# [Parsers] already loaded
# [CSV] already loaded
# [LSODA] already loaded
# [Sundials] already loaded
# [Rmath] already loaded
#   Updating registry at `~/.julia/registries/BioJuliaRegistry`
#   Updating git-repo `https://github.com/BioJulia/BioJuliaRegistry.git`
#   Updating registry at `~/.julia/registries/General`
#   Updating git-repo `https://github.com/JuliaRegistries/General.git`
#  Resolving package versions...
#   Updating `~/.julia/environments/v1.3/Project.toml`
#   [77a26b50] + DiffEqNoiseProcess v3.8.0
#   Updating `~/.julia/environments/v1.3/Manifest.toml`
#  [no changes]
# [RandomNumbers] already loaded
#  Resolving package versions...
#   Updating `~/.julia/environments/v1.3/Project.toml`
#   [4536629a] + OpenBLAS_jll v0.3.7+5
#   Updating `~/.julia/environments/v1.3/Manifest.toml`
#  [no changes]
#  Resolving package versions...
#   Updating `~/.julia/environments/v1.3/Project.toml`
#   [76087f3c] + NLopt v0.5.1
#   Updating `~/.julia/environments/v1.3/Manifest.toml`
#  [no changes]
#  Resolving package versions...
#  Installed Arpack_jll ─ v3.7.0+0
#  Installed Arpack ───── v0.3.2
#   Updating `~/.julia/environments/v1.3/Project.toml`
#   [68821587] + Arpack_jll v3.7.0+0
#   Updating `~/.julia/environments/v1.3/Manifest.toml`
#   [7d9fca2a] ↓ Arpack v0.4.0 ⇒ v0.3.2
#   [68821587] ↑ Arpack_jll v3.5.0+2 ⇒ v3.7.0+0
#   Building Arpack → `~/.julia/packages/Arpack/zCmTA/deps/build.log`
#  Resolving package versions...
#   Updating `~/.julia/environments/v1.3/Project.toml`
#   [47a9eef4] + SparseDiffTools v1.3.3
#   Updating `~/.julia/environments/v1.3/Manifest.toml`
#  [no changes]
# [SpecialFunctions] already loaded
# compiling line 201
# ┌ Info: UndefVarError(Symbol("#7#9"))
# │   LINE = 206
# └ @ Main.##anon_module#422 /Users/nmat471/.julia/packages/Fezzik/BRQC8/precomp.jl:209
# compiling line 546
# ┌ Info: UndefVarError(Symbol("#8#10"))
# │   LINE = 551
# └ @ Main.##anon_module#422 /Users/nmat471/.julia/packages/Fezzik/BRQC8/precomp.jl:554
# compiling line 2231
# done.
# .
# .
# Creating sysimg
# [ Info: DONE!!


# Julia automatically closes

# The next time you open Julia, the packages are pre-compiled!!


# To undo:
# remove itself from the startup.jl and deletes previous traces
# Fezzik.auto_trace(false)

# every now and then you'll want to start over , just call
# 
# Fezzik.revert()
# to go back to the beginning