#######################################################
# A test module to run Flow.jl, not in scope Main
# 
# Development workflow recommended by:
#
# https://docs.julialang.org/en/v1/manual/workflow-tips/
#
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst_ExampleModule.jl")

"""

module tst_ExampleModule
	# SETUP	
	cd("/GitHub/BioGeoJulia.jl/notes/")
	include("ModelLikes.jl")
	import .ModelLikes
	#using .Tmp

	include("ExampleModule.jl")
	import .ExampleModule

	using Profile     # for @profile
	using DataFrames  # for DataFrame
	using PhyloNetworks
	using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
	using BioGeoJulia.StateSpace
	using BioGeoJulia.TreePass
	using BioGeoJulia.SSEs
	
	using DifferentialEquations
	using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
	#Pkg.add(PackageSpec(url="https://github.com/JuliaDiffEq/deSolveDiffEq.jl"))
	#using deSolveDiffEq 
	# https://docs.juliadiffeq.org/stable/solvers/ode_solve/index.html

	using Profile     # for @profile
	using DataFrames  # for DataFrame
	using PhyloNetworks
	using RCall       # for df_to_Rdata, reval, g = globalEnv
	
	
	# RUN THE FUNCTIONS YOU ARE WORKING ON
	ExampleModule.example1()
	
	ExampleModule.example2()
	
	
end # End module
