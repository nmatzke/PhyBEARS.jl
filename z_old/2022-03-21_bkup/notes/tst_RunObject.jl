#######################################################
# A test module to run RunObject.jl, not in scope Main
# 
# Development workRunObject recommended by:
#
# https://docs.julialang.org/en/v1/manual/workRunObject-tips/
#
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst_RunObject.jl")

"""

module Tst_RunObject
	using BenchmarkTools # for @time
	using InvertedIndices # for Not
	using LSODA           # for lsoda()
	using Sundials        # for CVODE_BDF(linear_solver=:GMRES)
	using DifferentialEquations
	using Distributed
	using Random					# for MersenneTwister()
	using Dates						# for e.g. DateTime, Dates.now()
	using PhyloNetworks
	#using Plots						# for plot
	using DataFrames          # for DataFrame()
	using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
	using BioGeoJulia.StateSpace
	using BioGeoJulia.TreePass
	using BioGeoJulia.SSEs
	
	cd("/GitHub/BioGeoJulia.jl/notes/")
	include("Parsers.jl")
	import .Parsers

	include("RunObject.jl")
	import .RunObject

	using BioGeoJulia.StateSpace 
	using DataFrames  # for DataFrame

	RunObject.say_hello4()
	
	
	# construct the sub-objects, then the full one
	lgdata_fn = "/GitHub/BioGeoJulia.jl/ex/homs1/homs_geog.txt"
	geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

	tr_bkup = readTopology("(((human:0.5):0.5,(chimp:0.5):0.5):1.0,(gorilla:1.5):0.5);")
	trfn = "/GitHub/BioGeoJulia.jl/ex/homs1/homs_w_breakpoint.newick"
	tr = readTopology(trfn)
	trdf = prt(tr)

	ro = RunObject.construct_RunObj(trfn, lgdata_fn)
	Rnames(ro)
	Rnames(ro.meta)
	Rnames(ro.settings)
	Rnames(ro.gl)
	Rnames(ro.regs)
	
	
	#######################################################
	# Now, calculate a likelihood (assuming 1 global regime)
	#######################################################
	# 1. Set up a model (Qmat and Cmat tables)
	#
	# 2. Update the Qmat and Cmat rates with update function
	#
	# 3. Calculate the likelihood
	#    a. Es interpolator
	#    b. As interpolator (records instantaneous flow of Ds)
	#    c. Gflow interpolator
	#    d. Calc downpass
	#    e. Root state and final lnLs
	########################################################
	numareas = ro.gl.numareas
	numstates = ro.gl.numstates
	areas_list = ro.gl.areas_list
	states_list = ro.gl.states_list
	
	amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
	dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
	elist = repeat([0.123], numstates)
	allowed_event_types=["d","e"]

	Qarray = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	Qdf = prtQ(Qarray)
	Carray = setup_DEC_Cmat(areas_list, states_list)
	Cdf = prtC(Carray)
	C_row_weightvals = Carray.row_weightvals # has length numstates
	
	
	
	
	
end # End of module