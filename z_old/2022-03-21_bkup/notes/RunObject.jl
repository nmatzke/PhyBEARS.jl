#######################################################
# RunObject
#######################################################

module RunObject
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

print("\n\nStarting module 'RunObject'...loading dependencies...\n")
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

include("/GitHub/BioGeoJulia.jl/notes/Parsers.jl")
import .Parsers


# (1) List all function names here:
export say_hello4, default_BGB_params, Metadata, Settings, AllRegimes, Regime, Regimes, RunObj, construct_RunObj

#######################################################
# Goal: Set up the run object for BioGeoJulia.
#######################################################
# 
# This will include all necessary inputs and outputs
# for all regimes.
#
# The RunObject (ro) is described here:
# https://docs.google.com/document/d/1p9RXf4UGgu9NjXKKMhVBrObcVMjysqBRNsXb55xogUE/edit
# 

# Instructions to run this module:
"""
# Run instructions:
cd("/GitHub/BioGeoJulia.jl/notes/")
include("Parsers.jl")
import .Parsers


include("RunObject.jl")
ro = RunObject.setup_ro()

cd("/GitHub/BioGeoJulia.jl/notes/")
include("RunObject.jl")
include("tst_RunObject_n1.jl")

cd("/GitHub/BioGeoJulia.jl/notes/")
include("RunObject.jl")
import .RunObject
ro = RunObject.say_hello4()
ro = RunObject.setup_ro()
"""
#######################################################


#######################################################
# (2) write the functions here
#######################################################

"""
# Dummy function as an example

# To run:
say_hello4()
"""
say_hello4() = println("Function say_hello4() says hello!")
say_hello4()




function default_BGB_params()
	par = ["d", "e", "a", "b", "x", "n", "w", "u", "j", "ysv", "ys", "y", "s", "v", "mx01", "mx01j", "mx01y", "mx01s", "mx01v", "mx01r", "mf", "dp", "fdp"]
	
	type = collect(repeat(["fixed"], length(par)))
	init = collect(repeat([0.0], length(par)))
	min = collect(repeat([0.0], length(par)))
	max = collect(repeat([5.0], length(par)))
	est = collect(repeat([0.0], length(par)))
	desc = ["anagenesis: rate of 'dispersal' (range expansion)",
"anagenesis: rate of 'extinction' (range contraction)",
"anagenesis: rate of range-switching (i.e. for a standard char.)",
"anagenesis: exponent on branch lengths",
"exponent on distance (modifies d, j, a)",
"exponent on environmental distance (modifies d, j, a)",
"exponent on manual dispersal multipliers (modifies d, j, a)",
"anagenesis: exponent on extinction risk with area (modifies e)",
"cladogenesis: relative per-event weight of jump dispersal",
"cladogenesis: y+s+v",
"cladogenesis: y+s",
"cladogenesis: relative per-event weight of sympatry (range-copying)",
"cladogenesis: relative per-event weight of subset speciation",
"cladogenesis: relative per-event weight of vicariant speciation",
"cladogenesis: controls range size of smaller daughter",
"cladogenesis: controls range size of smaller daughter",
"cladogenesis: controls range size of smaller daughter",
"cladogenesis: controls range size of smaller daughter",
"cladogenesis: controls range size of smaller daughter",
"root: controls range size probabilities of root",
"mean frequency of truly sampling OTU of interest",
"detection probability per true sample of OTU of interest",
"false detection of OTU probability per true taphonomic control sample"]
	
	pardf = DataFrame(par=par, type=type, init=init, min=min, max=max, est=est, desc=desc)
	pardf
	
	return(pardf)	
end # function default_BGB_params()


# Structure for ro.meta
#
# (Metadata and input files for the run)
# 
# This structure is mutable so that the inputs can be 
# easily changed.
# 
# This is the setup for ONE tree - items that take arrays
#   do so for e.g. multiple regimes on this one tree
# 
# There is always a minimum of 1 regime
# meta: description of the run, input files, etc.
mutable struct Metadata
	# Meta: description of the run, input files, etc.
	run_num::Any
	run_name::String
	run_desc::String
	run_dir::String
	trfn::String
	
	# Geography file containing ranges for all
	# tips at all times.
	# 
	# If some areas appear/disappear/are unoccupied/merge/split etc,
	# include them all here, then set up a...
	#
	# list_of_states_lists
	#
	# ...to describe the states list in each regime.
	geogfn::String
	
	# Other files -- assuming just time-stratification here.
	# For something more complex (e.g. different regimes for
	# different clade, or both clade+time stratification), 
	# set up inputs with another function.
	
	# One file, with blank lines, specifying time-stratification
	timesfn::String
	distsfn::String
	envdistsfn::String
	dispersal_multipliers_fn::String
	area_of_areas_fn::String
	detects_fn::String
	controls_fn::String
	
	# A bunch of additional distance/dispersal multiplier matrices (for kicks)
	dists2_fn::String
	dists3_fn::String
	dists4_fn::String
	dists5_fn::String
	dists6_fn::String

	# fs: functions -- structures can't contain functions, but they can
	#                  contain userfuncs
	# calc_lnl - overall likelihood calculator
	# calc_likes_down_branch(regime) - calculate likelihoods down a branch segment during downpass
	# calc_tiplikes(regime) - calculate tip likelihoods for all tips, given trdf, pardf
	# calc_likes_down_regimes(reg1, reg2) -- convert likelihoods between regimes during downpass
	#add_two::userfunc

end # END struct Metadata

# Structure for ro.settings 
mutable struct Settings
	# Settings: user-modifiable inputs
	time_tops::Array{Float64}
	time_bot_Es::Float64  # Time before present at which to end calculation of Es (e.g. 1.5 * root)
	time_bot_Ds::Float64  # Time before present at which to end calculation of Ds (e.g. 1.25 * root)

	# Maximum maximum number of areas per range (before any manual modifications to the states_list)
	max_range_size::Int64
	include_null_range::Bool
	
	# Solver settings
	solver::Any  # Might not get used, as solver might be set in calc_likes_down_branch
	save_everystep::Bool
	abstol::Float64
	reltol::Float64

end # END struct Settings




# Structure for ro.gl : AllRegimes (items used/calculated across all regimes)
mutable struct AllRegimes
	# Areas and states can be Int numbers, letters, or other strings
	# Mostly we use Int numbers for processing; text for display
	numareas::Int64			# Number of areas (total across all regimes)
	numstates::Int64		# Number of areas (total across all regimes)
	areas_list::Array{Int64,1}           # 1D array; a list of all possible areas across all regimes
	states_list::Array{Array{Any,1},1}   # 1D array of 1D arrays; a list of all possible states across all regimes
	
	# Text versions of the areas_list and states_list
	area_names::Array{String,1}
	ranges_list::Array{String,1} # Ranges are list of areas collapse to a string, e.g. "ABC"
	geog_df::DataFrame
	
	# Tree details etc.
	tr::HybridNetwork # Tree with possible ith included direct ancestors)
	trdf::DataFrame   # tree table (with regimes in column "reg")
	pardf::DataFrame  # a df with all parameters, free/fixed/fn, regime, init, est, min, max, desc
	
	#tiplikes::Array{Float64,2} # tip likelihoods on all possible states # Might be updated
	
	# Results:
	res::Res    # likelihoods results tables (on all possible states)
	branch_lnl::Float64
	root_lnl::Float64
	ttl_lnl::Float64


end # END struct AllRegimes


# Structure for a single Regime
mutable struct Regime
	name::String
	top::Float64
	bot::Float64
	desc::String
	model::String
	areas::Array{Int64,1}    # 1D array of areas in this regime
	states::Array{Array{Int64,1},1}   # 1D array of 1D arrays; a list of all possible states in this regime
	Qdf::DataFrame
	Cdf::DataFrame
	Crow_weightvals::Array{Float64,1}    # 1D array of areas in this regime
	mu::Array{Float64,1}     # 1D array of floats with per-area extinction rates for this regime
end

function construct_Regime()
	name = ""
	top = 0.0
	bot = 0.0
	desc = ""
	model = ""
	areas = [0]
	states = [[0]]
	Qdf = DataFrame([])
	Cdf = DataFrame([])
	C_row_weightvals = [0.0]
	mu = [0.0]
	
	return(Regime(name, top, bot, desc, model, areas, states, Qdf, Cdf, C_row_weightvals, mu))
end

# Structure for ref: Regimes
mutable struct Regimes
	ra::Array{Regime,1}  # ra = Regimes Array
	num_regimes::Float64
end


# Structure for ro: the RunObject
struct RunObj
	# Meta: description of the run, input files, etc.
	meta::Metadata
	settings::Settings
	gl::AllRegimes
	regs::Regimes
end # END struct RunObject


"""
# construct the sub-objects, then the full one
lgdata_fn = "/GitHub/BioGeoJulia.jl/ex/homs1/homs_geog.txt"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

tr_bkup = readTopology("(((human:0.5):0.5,(chimp:0.5):0.5):1.0,(gorilla:1.5):0.5);")
trfn = "/GitHub/BioGeoJulia.jl/ex/homs1/homs_w_breakpoint.newick"
tr = readTopology(trfn)
trdf = prt(tr)

ro = construct_RunObject(trfn, lgdata_fn)
"""
function construct_RunObj(trfn, lgdata_fn)
	# Read the files
	geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
	tr = readTopology(trfn)
	trdf = prt(tr)

	# Text versions of the areas_list
	tmpnames = string.(names(geog_df))
	area_names = tmpnames[tmpnames .!= "tipnames"]
	numareas = length(area_names)
	root_age = trdf[tr.root,:node_age]
	
	#############################################
	# Set up the Metadata
	#############################################
	run_num = 0
	run_name = ""
	run_desc = ""
	run_dir = Rgetwd()
	trfn = trfn
	geogfn = lgdata_fn
	
	# One file, with blank lines, specifying time-stratification
	timesfn = ""
	distsfn = ""
	envdistsfn = ""
	dispersal_multipliers_fn = ""
	area_of_areas_fn = ""
	detects_fn = ""
	controls_fn = ""
	
	# A bunch of additional distance/dispersal multiplier matrices (for kicks)
	dists2_fn = ""
	dists3_fn = ""
	dists4_fn = ""
	dists5_fn = ""
	dists6_fn = ""
	
	meta = RunObject.Metadata(run_num, run_name, run_desc, run_dir, trfn, geogfn, timesfn, distsfn, envdistsfn, dispersal_multipliers_fn, area_of_areas_fn, detects_fn, controls_fn, dists2_fn, dists3_fn, dists4_fn, dists5_fn, dists6_fn)
	

	#############################################
	# Settings
	#############################################
	time_tops = [0.0]
	time_bot_Es = 1.5 * root_age
	time_bot_Ds = 1.25 * root_age

	# Maximum maximum number of areas per range (before any manual modifications to the states_list)
	max_range_size = numareas
	include_null_range = true

	# Solver settings
	solver = Tsit5()
	save_everystep = false
	abstol = 1e-9
	reltol = 1e-9

	settings = RunObject.Settings(time_tops, time_bot_Es, time_bot_Ds, max_range_size, include_null_range, solver, save_everystep, abstol, reltol)
	
	
	
	#############################################
	# AllRegimes
	# inputs for all regimes (gl: global)
	#############################################
	areas_list = collect(1:numareas)
	states_list = areas_list_to_states_list(areas_list, max_range_size, include_null_range)
	# Text versions of the areas_list and states_list
	tmpnames = string.(names(geog_df))
	area_names = tmpnames[tmpnames .!= "tipnames"]
	ranges_list = states_list_to_txt(states_list, area_names; delim="")
	
	# Tree details etc.
	# tr = tr
	# trdf = trdf
	# tiplikes::Array{Float64,2} # tip likelihoods on all possible states # Might be updated
	
	# Parameters df
	pardf = RunObject.default_BGB_params()
	
	# Results:
	numareas = length(areas_list)
	numstates = n = length(states_list)
	res = construct_Res(tr, n)
	branch_lnl = 0.0
	root_lnl = 0.0
	ttl_lnl = 0.0
	
	gl = RunObject.AllRegimes(numareas, numstates, areas_list, states_list, area_names, ranges_list, geog_df, tr, trdf, pardf, res, branch_lnl, root_lnl, ttl_lnl)
	
	#######################################################
	# Regimes
	#######################################################
	regime_nums = sort(unique(trdf[!,:reg]))
	num_regimes = length(regime_nums)
	tmp_regimes = [RunObject.construct_Regime()]
	
	# Build up the list of regimes
	if (num_regimes > 1)
		for i in 2:num_regimes
			tmp_regime = [RunObject.construct_Regime()]
			append!(tmp_regimes, tmp_regime)
		end
	end
	num_regimes = length(tmp_regimes)
	regs = Regimes(tmp_regimes, num_regimes)
	
	#######################################################
	# OK, now, finally, output the RunObject
	#######################################################
	ro = RunObj(meta, settings, gl, regs)
	
	return(ro)	
end # END function construct_RunObj(trfn, lgdata_fn)




"""
include("/GitHub/BioGeoJulia.jl/notes/RunObject.jl")
import .RunObject

ro = RunObject.setup_ro()

lgdata_fn = "/GitHub/BioGeoJulia.jl/Rsrc/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)


"""





end # ENDING RunObject