#######################################################
# A test module to run Tmp.jl, not in scope Main
# 
# Development workflow recommended by:
#
# https://docs.julialang.org/en/v1/manual/workflow-tips/
#
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst.jl")

"""

# 
#######################################################


#######################################################
# Run all of the functions here
#######################################################

module Tst
	include("Tmp.jl")
	import .Tmp
	#using .Tmp

	include("Flow.jl")
	import .Flow

	using BioGeoJulia.StateSpace 
	using DataFrames  # for DataFrame

	Tmp.say_hello()
	# say_hello()



	# Distribution of smaller daughters:
	numareas=6
	max_numareas=6
	maxent_constraint_01=0.0
	maxent_result = Tmp.discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)
	
	max_numareas=6
	maxent_constraint_01=0.5
	maxent_result = Tmp.discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)
	
	max_numareas=6
	maxent_constraint_01=1.0
	maxent_result = Tmp.discrete_maxent_distrib_of_smaller_daughter_ranges(max_numareas, maxent_constraint_01)

	max_numareas=6
	maxent_constraint_01 = 0.0001    # ranges from 0.0001 (all weight on ranges of size 1)
																	 # to 0.9999 (all weight on ranges of max size)
	NA_val=NaN
	maxent_constraint_01=0.0
	relprob_subsets_matrix = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01, NA_val)

	maxent_constraint_01=0.5
	relprob_subsets_matrix = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01, NA_val)

	maxent_constraint_01=1.0
	relprob_subsets_matrix = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01, NA_val)


	maxent_constraint_01=0.5
	relprob_subsets_matrix = Tmp.relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01, NA_val)

	areas_list = [1,2,3]
	states_list = areas_list_to_states_list(areas_list, 3, true)
	Cparams = Tmp.default_Cparams()
	max_numareas = length(areas_list)
	maxent_constraint_01 = 0.0
	maxent01symp = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01sub = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent01jump = Tmp.relative_probabilities_of_subsets(max_numareas, maxent_constraint_01)
	maxent_constraint_01 = 0.5
	maxent01vic = Tmp.relative_probabilities_of_vicariants(max_numareas, maxent_constraint_01)
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)

	Carray = Tmp.setup_DEC_Cmat(areas_list, states_list, maxent01, Cparams)

	# Extract the values
	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals;
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	row_weightvals
	
	df = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);

	showall(df, true)
	by(df, :event, nrow)

	# DEC+J
	Cparams.j = 0.1
	Cparams.y = (3.0-Cparams.j) / 3.0
	Cparams.s = (3.0-Cparams.j) / 3.0
	Cparams.v = (3.0-Cparams.j) / 3.0
	Cparams
	Carray = Tmp.setup_DEC_Cmat(areas_list, states_list, maxent01, Cparams)

	Carray_event_types = Carray.Carray_event_types;
	Carray_ivals = Carray.Carray_ivals;
	Carray_jvals = Carray.Carray_jvals;
	Carray_kvals = Carray.Carray_kvals;
	Cijk_weights = Carray.Cijk_weights;
	Cijk_vals = Carray.Cijk_vals;
	row_weightvals = Carray.row_weightvals;
	DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	row_weightvals
	
	df = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals);
	showall(df, true)
	
	
	by(df, :event, nrow)
		
	by(df, :i, ysum = :event => Tmp.sumy, ssum = :event => Tmp.sums, vsum = :event => Tmp.sumv, jsum = :event => Tmp.sumj) 

	by(df, :i, :prob => sum)
	
	
	
	# Update the parameter values
	Cparams.j = 0.3;
	Cparams.y = (3.0-Cparams.j) / 3.0;
	Cparams.s = (3.0-Cparams.j) / 3.0;
	Cparams.v = (3.0-Cparams.j) / 3.0;
	
	Carray2 = Tmp.update_Cijk_vals(Carray, areas_list, states_list, maxent01, Cparams);
	Carray_event_types = Carray2.Carray_event_types;
	Carray_ivals = Carray2.Carray_ivals;
	Carray_jvals = Carray2.Carray_jvals;
	Carray_kvals = Carray2.Carray_kvals;
	Cijk_weights = Carray2.Cijk_weights;
	Cijk_vals = Carray2.Cijk_vals;
	row_weightvals = Carray2.row_weightvals;
	df2 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	by(df2, :i, :prob => sum)
	

	Cparams.j = 1.0;
	Cparams.y = (3.0-Cparams.j) / 3.0;
	Cparams.s = (3.0-Cparams.j) / 3.0;
	Cparams.v = (3.0-Cparams.j) / 3.0;
	
	Carray3 = Tmp.update_Cijk_vals(Carray, areas_list, states_list, maxent01, Cparams);
	Carray_event_types = Carray3.Carray_event_types;
	Carray_ivals = Carray3.Carray_ivals;
	Carray_jvals = Carray3.Carray_jvals;
	Carray_kvals = Carray3.Carray_kvals;
	Cijk_weights = Carray3.Cijk_weights;
	Cijk_vals = Carray3.Cijk_vals;
	row_weightvals = Carray3.row_weightvals;
	df3 = DataFrame(event=Carray_event_types, i=Carray_ivals, j=Carray_jvals, k=Carray_kvals, weight=Cijk_weights, prob=Cijk_vals)
	by(df3, :i, :prob => sum)
	
	
	
	#######################################################
	# NEXT:
	# 1. multiply conditional probabilities by lambda
	# 2. load some geographic states, and a phylogeny
	# 3. input into likelihood
	# 4. maximum likelihood (or speed tests)
	#######################################################
	
	
	
	
	
	
	
end # End of module