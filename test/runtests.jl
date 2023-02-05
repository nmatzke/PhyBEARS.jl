using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
#using PhyloBits.PNtypes					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames						# for DataFrame()

# List each PhyBEARS code file prefix here
using PhyloBits.TrUtils
using PhyloBits.TreeTable
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
#using PhyBEARS.TreeTable
using PhyBEARS.TreePass
#using PhyBEARS.TrUtils
using PhyBEARS.SSEs
using PhyBEARS.Parsers
#include(paste0([dd, "/notes/Parsers.jl"]))

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/test/")
include("/GitHub/PhyBEARS.jl/test/runtests.jl")
"""

# default directory for data etc.
pathof_result = pathof(PhyBEARS)
dd = pp(pathof_result, "PhyBEARS")
dd = "/GitHub/PhyBEARS.jl"


@testset "All tests" begin

@testset "Example" begin
	@test hello("Julia") == "Hello, Julia"
	@test domath(2.0) ≈ 7.0
end


@testset "PhyBEARS" begin
	@test hello_PhyBEARS("Julia") == "hello_PhyBEARS() says 'Julia'"
	@test add_one_PhyBEARS(2.0) ≈ 3.0   # note the approximate equals (compare ≈=)
end

@testset "StateSpace" begin
	@test numstates_from_numareas(3,3,false) == 7
	@test numstates_from_numareas(3,3,true) == 8
	@test numstates_from_numareas(10,1,false) == 10
	@test numstates_from_numareas(10,2,false) == 55
	@test numstates_from_numareas(10,3,false) == 175
	@test numstates_from_numareas(10,10,false) == 1023
	@test numstates_from_numareas(10,10,true) == 1024
	@test numstates_from_numareas(20,20,true) == 1048576
	
	# Set up list of areas
	area_nums = collect(1:3)
	
	tmpstr = "Array{Any,1}[[1], [2], [3]]"
	states_list_answer = eval(Meta.parse(tmpstr))
	states_list = areas_list_to_states_list(area_nums, 1, false)
	@test states_list == states_list_answer

	tmpstr = "Array{Any,1}[[], [1], [2], [3]]"
	states_list_answer = eval(Meta.parse(tmpstr))
	states_list = areas_list_to_states_list(area_nums, 1, true)
	@test states_list == states_list_answer

	tmpstr = "Array{Any,1}[[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]"
	states_list_answer = eval(Meta.parse(tmpstr))
	states_list = areas_list_to_states_list(area_nums, 3, false)
	@test states_list == states_list_answer
	
	tmpstr = "Array{Any,1}[[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]"
	states_list_answer = eval(Meta.parse(tmpstr))
	states_list = areas_list_to_states_list(area_nums, 3, true)
	@test states_list == states_list_answer

	# still need to be done:

	Cparams = default_Cparams()
	@test Cparams.y == 1
	@test Cparams.s == 1
	@test Cparams.v == 1
	@test Cparams.j == 0
	"""

	# sumy
	# sums
	# sumv
	# sumj
	unsure what these do
	"""
	# get_default_inputs prt not defined?
	# run_model
	
	tmpstr = "(mu_vals = [0.1, 0.1], Qij_vals = [0.01, 0.001], Cijk_vals = [0.222222, 0.222222])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.001)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp
	
	tmpstr = "(mu_vals = [0.1, 0.1], Qij_vals = [0.01, 0.001], Cijk_vals = [0.333333, 0.333333])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.333333, deathRate=0.1, q01=0.01, q10=0.001)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp

	tmpstr = "(mu_vals = [0.2, 0.2], Qij_vals = [0.01, 0.001], Cijk_vals = [0.222222, 0.222222])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.2, q01=0.01, q10=0.001)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp

	tmpstr = "(mu_vals = [0.1, 0.1], Qij_vals = [0.2, 0.001], Cijk_vals = [0.222222, 0.222222])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.1, q01=0.2, q10=0.001)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp

	tmpstr = "(mu_vals = [0.1, 0.1], Qij_vals = [0.01, 0.2], Cijk_vals = [0.222222, 0.222222])"
	p_Es_v5 = setup_MuSSE(2; birthRate=0.222222, deathRate=0.1, q01=0.01, q10=0.2)
	p_Es_v5_tmp =  eval(Meta.parse(tmpstr))
	@test p_Es_v5.params == p_Es_v5_tmp

	# setup_DEC_DEmat
	numareas = 3
	areas_list = collect(1:numareas)
	states_list = areas_list_to_states_list(areas_list, 3, true)
	numstates = length(states_list)
	amat = reshape(collect(1:(numareas^2)), (numareas,numareas))
	dmat = reshape(collect(1:(numareas^2)), (numareas,numareas)) ./ 100
	elist = repeat([0.123], numstates)
	allowed_event_types=["d","e"]

	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat)
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	Qarray_event_types = Qmat.Qarray_event_types
	@test Qmat.Qarray_jvals[4] == 7

	states_list = areas_list_to_states_list(areas_list, 3, false)
	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["d","e"])
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	Qarray_event_types = Qmat.Qarray_event_types
	@test Qmat.Qarray_jvals[4] == 6

	states_list = areas_list_to_states_list(areas_list, 3, false)
	Qmat = setup_DEC_DEmat(areas_list, states_list, dmat, elist, amat; allowed_event_types=["a"])
	Qarray_ivals = Qmat.Qarray_ivals
	Qarray_jvals = Qmat.Qarray_jvals
	Qij_vals = Qmat.Qij_vals
	Qarray_event_types = Qmat.Qarray_event_types
	@test Qmat.Qarray_jvals[4] == 1

	# update_Qij_vals
	Qmat1_df = hcat(Qarray_ivals, Qarray_jvals, Qij_vals, Qarray_event_types)
	dmat = reshape(repeat([0.5], numareas^2), (numareas,numareas))
	"""
	Qmat2 = update_Qij_vals(Qmat, areas_list, states_list, dmat, elist, amat )
	ERROR: MethodError: Cannot `convert` an object of type Array{Int64,2} to an object of type Float64
	I SEE NO CONVERT-LIKE MOVES WITHIN THIS FUNCTION, I AM UNSURE WHERE IT IS TRYING TO CONVERT ANYTHING
	"""
	# relative_probabilities_of_subsets
	# relative_probabilities_of_vicariants()
	# discrete_maxent_distrib_of_smaller_daughter_ranges()
	"""
	relative_probabilities_of_subsets()
	ERROR: ArgumentError: MathProgBase solvers like `solve!(problem, SCSSolver())` are no longer supported. Use instead e.g. `solve!(problem, SCS.Optimizer)`.
	
	Same issue for 
	relative_probabilities_of_vicariants()
	discrete_maxent_distrib_of_smaller_daughter_ranges()
	
	"""

	# array_in_array
	

	# is_event_vicariance
	ancstate = [1, 2, 3,4];
	lstate = [1, 2];
	rstate = [4];
	@test is_event_vicariance(ancstate, lstate, rstate) == false

	ancstate = [1, 2, 3,4];
	lstate = [1, 2];
	rstate = [2, 4];
	@test is_event_vicariance(ancstate, lstate, rstate) == false

	ancstate = [1, 2, 3,4];
	lstate = [1, 2];
	rstate = [3, 4];
	@test is_event_vicariance(ancstate, lstate, rstate) == true

	# setup_DEC_Cmat
	"""
	areas_list = [1,2,3]
	states_list = areas_list_to_states_list(areas_list, 3, true)
	Cparams=(y=1.0,s=1.0,v=1.0,j=0.0)
	max_numareas = length(areas_list)
	maxent_constraint_01 = 0.0
	maxent01symp = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01sub = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01jump = relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent_constraint_01 = 0.5
	maxent01vic = relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01)
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
	predeclare_array_length=10000000
	Carray = setup_DEC_Cmat(areas_list, states_list, Cparams)

	Cannot be run because relative_probabilities_of_subsets cannot be run!
	"""

	# update_Cijk_vals
end

@testset "TrUtils" begin
	tmpstr = "HybridNetwork"
	answer = eval(Meta.parse(tmpstr))
	great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
	tr = readTopology(great_ape_newick_string)
	@test type(tr) == answer
	
	setwd("/Users/")
	@test setwd("/Users/") == cd("/Users/")
	@test getwd() == pwd()
	@test getwd() == "/Users"
	@test Rgetwd() == pwd()
	@test Rgetwd() == "/Users"

	tmparray = recursive_find(dd)
	@test tmparray[1] == paste0([dd, "/src/BGExample.jl"])

	#tmparray = include_jls(dd)
	#@test tmparray[1] == paste0([dd, "/src/PhyBEARS.jl"])

	A = ones(3,3)
	B = ones(3,3)
	@test dim(A) == size(A)
	@test Rdim(A) == size(A)

	C = Int64[1,2,3,4,5,6,7,8,9,10]
	@test seq(1, 10, 1) == C
	@test Rchoose(10,5) == 252

	D = ones(3,6)
	@test Rcbind(A, B) == hcat(A,B)
	@test Rcbind(A, B) == D
	E = ones(6,3)
	@test Rrbind(A, B) == vcat(A,B)
	@test Rrbind(A, B) == E

	@test paste(1:12, delim="") == "123456789101112"
	@test paste0(1:12) == "123456789101112"

	F = "tester"
	@test type(F) == String
	@test class(F) == "String"
	@test Rclass(F) == "String"
	
	# Q: is there a reason slashslash() has the internal code repeated several times?
	# A: to remove any cases of ///, ////, etc.
	@test slashslash("//GitHub/PhyBEARS.jl//src//PhyBEARS.jl") == "/GitHub/PhyBEARS.jl/src/PhyBEARS.jl"
	@test slashslash(paste0([addslash("/GitHub/PhyBEARS.jl"), "/src/PhyBEARS.jl"])) == "/GitHub/PhyBEARS.jl/src/PhyBEARS.jl"

	tmpmatrix = [3 1; 3 2; 5 3; 5 4; 7 5; 7 6]
	tmpstr = repr(tmpmatrix)
	tmpstr2 = eval(Meta.parse(tmpstr))
	@test Reval(tmpstr) == tmpstr2

	tmpstr = "[[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]"
	tmpstr2 = tmpstr
	states_list = Reval(tmpstr)
	@test Rdput(states_list) == tmpstr2

	great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
	tr = readTopology(great_ape_newick_string)
	@test Rnames(tr)[1] == :numTaxa
	@test Rtypes(tr)[1] == Int64
	@test ont(tr) == Any[:numTaxa Int64; :numNodes Int64; :numEdges Int64; :node Array{PhyloBits.PNtypes.Node,1}; :edge Array{PhyloBits.PNtypes.Edge,1}; :root Int64; :names Array{String,1}; :hybrid Array{PhyloBits.PNtypes.Node,1}; :numHybrids Int64; :cladewiseorder_nodeIndex Array{Int64,1}; :visited Array{Bool,1}; :edges_changed Array{PhyloBits.PNtypes.Edge,1}; :nodes_changed Array{PhyloBits.PNtypes.Node,1}; :leaf Array{PhyloBits.PNtypes.Node,1}; :ht Array{Float64,1}; :numht Array{Int64,1}; :numBad Int64; :hasVeryBadTriangle Bool; :index Array{Int64,1}; :loglik Float64; :blacklist Array{Int64,1}; :partition Array{PhyloBits.PNtypes.Partition,1}; :cleaned Bool; :isRooted Bool]


	@test Rnrow(A) == 3
	@test Rncol(A) == 3
	@test Rsize(A) == (3,3)

	tmpDF = DataFrame(A = 1:4, B = ["M", "F", "F", "M"], C = 5:8)
	tmparray = [1,2,3,4]
	@test Rorder(tmpDF) == tmparray
	
	# tmpDF2 = DataFrame(A = 1:4, C = 5:8)
	# headLR(tmpDF, 1, 1) == tmpDF2
	# false?
	
"""
	A note! Even when the dataframes are identical, 
	pulling the left and right collumns seems to twist them?
	In this case there are ONLY 2 collumns? so the outcome should be identical to itself?

	Something to do with it being flattened? 
	Unsure on how to test atm

	tmptmp = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
	tmptmp2 = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
	headLR(tmptmp, 1, 1) == tmptmp2

	FALSE
"""
	
	tmparray = (1)
	@test single_element_array_to_scalar(tmparray) == 1

	tmpline = print("1 module PhyBEARS")
	@test headf(mpf([dd, "/src/PhyBEARS.jl"]); numlines=1) == tmpline

	# How to test these?
	# @test df_to_Rdata
	# @test source(paste0([dd, "/src/PhyBEARS.jl"])) == include(paste0([dd, "/src/PhyBEARS.jl"]))
	# @test headLR(df, num_startcols=4, num_endcols=4) ==
	# @test flat2(arr) ==
	# @test moref(fn) ==
	# @test scr2str ==
	
	#######################################################
	# Test the birth-death likelihood calculation in bd_liks()
	#######################################################
	# These numbers come from: bd_liks() in ClaSSE_mods_v2.R
	R_dev = -0.9741936
	R_lnl_topology = 36.39545
	R_lnl_numBirths = -18.90835
	R_lnl_Births_above_root = 13.58031
	R_lnl_numtips_wOneMinusDeathRate = 0.0
	R_lnl_branching_times = -30.58031
	R_lnL = 0.4870968

	tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")
	birthRate=0.3288164   # ML birthRate for Psychotria tree
	deathRate=0.0					# ML deathRate for Psychotria tree
	bd = bd_liks(tr, birthRate, deathRate)
	
	
	@test round(R_dev; digits=5) == round(bd.deviance; digits=5)
	@test round(R_lnl_topology; digits=5) == round(bd.lnl_topology; digits=5)
	@test round(R_lnl_numBirths; digits=5) == round(bd.lnl_numBirths; digits=5)
	@test round(R_lnl_Births_above_root; digits=5) == round(bd.lnl_Births_above_root; digits=5)
	@test round(R_lnl_numtips_wOneMinusDeathRate; digits=5) == round(bd.lnl_numtips_wOneMinusDeathRate; digits=5)
	@test round(R_lnl_branching_times; digits=5) == round(bd.lnl_branching_times; digits=5)
	@test round(R_lnL; digits=5) == round(bd.lnL; digits=5)


	# These numbers come from: bd_liks() in ClaSSE_mods_v2.R
	R_dev = 8.613706
	R_lnl_topology = 0.6931472
	R_lnl_numBirths = 0.0
	R_lnl_Births_above_root = 1.0
	R_lnl_numtips_wOneMinusDeathRate = 0.0
	R_lnl_branching_times = -6.0
	R_lnL = -4.306853

	tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
	birthRate=1.0
	deathRate=0.0
	bd = bd_liks(tr, birthRate, deathRate)
	
	@test round(R_dev; digits=5) == round(bd.deviance; digits=5)
	@test round(R_lnl_topology; digits=5) == round(bd.lnl_topology; digits=5)
	@test round(R_lnl_numBirths; digits=5) == round(bd.lnl_numBirths; digits=5)
	@test round(R_lnl_Births_above_root; digits=5) == round(bd.lnl_Births_above_root; digits=5)
	@test round(R_lnl_numtips_wOneMinusDeathRate; digits=5) == round(bd.lnl_numtips_wOneMinusDeathRate; digits=5)
	@test round(R_lnl_branching_times; digits=5) == round(bd.lnl_branching_times; digits=5)
	@test round(R_lnL; digits=5) == round(bd.lnL; digits=5)


	# These numbers come from: bd_liks() in ClaSSE_mods_v2.R
	R_dev = 5.782411
	R_lnl_topology = 0.6931472
	R_lnl_numBirths = -6.907755
	R_lnl_Births_above_root = 0.001
	R_lnl_numtips_wOneMinusDeathRate = -20.72327
	R_lnl_branching_times = 24.04567
	R_lnL = -2.891206

	tr = readTopology("((chimp:1,human:1):1,gorilla:2);")
	birthRate=1.0
	deathRate=0.999
	bd = bd_liks(tr, birthRate, deathRate)
	
	@test round(R_dev; digits=5) == round(bd.deviance; digits=5)
	@test round(R_lnl_topology; digits=5) == round(bd.lnl_topology; digits=5)
	@test round(R_lnl_numBirths; digits=5) == round(bd.lnl_numBirths; digits=5)
	@test round(R_lnl_Births_above_root; digits=5) == round(bd.lnl_Births_above_root; digits=5)
	@test round(R_lnl_numtips_wOneMinusDeathRate; digits=5) == round(bd.lnl_numtips_wOneMinusDeathRate; digits=5)
	@test round(R_lnl_branching_times; digits=5) == round(bd.lnl_branching_times; digits=5)
	@test round(R_lnL; digits=5) == round(bd.lnL; digits=5)
	#######################################################
	# END: Test the birth-death likelihood calculation in bd_liks()
	#######################################################

	
end

@testset "TreePass" begin
	# Test if the printed tree table from prt()
	# gets the node ages correct
	tmpstr = "[0.0, 0.0, 6.0, 0.0, 7.0, 0.0, 12.0]"
	node_ages_in_prt = eval(Meta.parse(tmpstr))
	
	great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
	tr = readTopology(great_ape_newick_string)
	rootnodenum = tr.root
	trdf = prt(tr, rootnodenum)
	@test trdf[!, :node_age] == node_ages_in_prt
end

@testset "SSEs" begin
	# parameterized_ClaSSE
	# parameterized_ClaSSE_Es
	# parameterized_ClaSSE_Ds
	# parameterized_ClaSSE_v5
	# parameterized_ClaSSE_Es_v5
	# parameterized_ClaSSE_Ds_v5

end

@testset "Parsers" begin
	lgdata_fn = paste0([dd, "/data/Psychotria/Psychotria_geog.data"])
	geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
end


@testset "Checks against BioGeoBEARS models" begin
	include("/GitHub/PhyBEARS.jl/test/apes_SSE/apes_M0_DEC_v1.jl")
	include("/GitHub/PhyBEARS.jl/test/apes_SSE/apes_M0_DEC+J_v1.jl")
	include("/GitHub/PhyBEARS.jl/test/apes_SSE/apes_M0_DIVALIKE_v1.jl")
	include("/GitHub/PhyBEARS.jl/test/apes_SSE/apes_M0_DIVALIKE+J_v1.jl")
	include("/GitHub/PhyBEARS.jl/test/apes_SSE/apes_M0_BAYAREALIKE_v1.jl")
	include("/GitHub/PhyBEARS.jl/test/apes_SSE/apes_M0_BAYAREALIKE+J_v1.jl")

end


@testset "Checks against diversitree's BiSSE" begin
	include("/GitHub/PhyBEARS.jl/test/bisse_ASR_v1.jl")
	include("/GitHub/PhyBEARS.jl/test/bisse_ASR+J_v1.jl")
	include("/GitHub/PhyBEARS.jl/test/bisse_ASR+Jv12_v1.jl")

end










end # END: @testset "All tests"
