module MaxentInterp
# Use a pre-constructed interpolator for determining the weighting of 
# rangesizes of the smaller daughter.
#
# Creates maxent_interp, a predetermined list of interpolators, to 
# interpolate vectors from length 1 area-20 areas, assigning each
# rangesize a weight according to a Maximum Entropy optimization of
# a discrete multinomial variable (e.g. like a dice roll).
# 
# This way, 1 number, maxent01, can determine if the smaller daughter
# will 
# 
# * always be rangesize 1 (maxent01 = 0.0)
# * always be rangesize max (maxent01 = 1.0)
# * any daughter rangesize will have equal probability (maxent01 = 0.5)
# ...or something else.
# 
# This avoids the problem of using the SCS and Convex packages, which
# keep changing their syntax and also get quite slow as they use
# numerical optimization. For MaxentInterp, these optimizations were
# run one time, for daughter rangesize=1-20, and maxent01 = seq(0,1,0.01).
#
__precompile__(false)  # false will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/
print("PhyBEARS: loading MaxentInterp.jl dependencies...")

using JLD2					 # for @save and @load
using Interpolations # for Interpolations.scale, Interpolations.interpolate
using Distributions  # for quantile
using PhyBEARS
using PhyloBits.TrUtils # for flat2() (similar to unlist)

export maxent_interp, maxent01_defaults, discrete_maxent_distrib_of_smaller_daughter_ranges, relative_probabilities_of_subsets, relative_probabilities_of_vicariants

print("...done.\n")

# Setup
"""
include("/GitHub/PhyBEARS.jl/src/MaxentInterp.jl")
"""



#######################################################
# Load the interpolator into the module's environment
#######################################################
"""
# Deploy the maxent_interp array of interpolators
daughter_rangesize = 3
maxent01 = 0.0
maxent_interp[daughter_rangesize](maxent01)
maxent01 = 0.5
maxent_interp[daughter_rangesize](maxent01)
maxent01 = 1.0
maxent_interp[daughter_rangesize](maxent01)
"""

# Find the path where maxent_interp.jld2 is saved
path0 = replace(pathof(PhyBEARS), "PhyBEARS.jl/"=>"tmpabc.jl/")
path1 = replace(path0, "PhyBEARS.jl"=>"")
path2 = replace(path1, "tmpabc.jl/"=>"PhyBEARS.jl/")

interp_jld2_path = join([path2, "maxent_interp.jld2"])
@load interp_jld2_path maxent_interp # 
length(maxent_interp)




"""
# Default maxent multipliers for different kinds of events
include("/GitHub/PhyBEARS.jl/src/StateSpace.jl")
import .StateSpace

total_numareas = 3
maxent01 = StateSpace.maxent01_defaults(total_numareas);
maxent01symp = maxent01.maxent01symp
maxent01sub = maxent01.maxent01sub
maxent01vic = maxent01.maxent01vic
maxent01jump = maxent01.maxent01jump


total_numareas = 4
maxent01 = StateSpace.maxent01_defaults(total_numareas);
maxent01symp = maxent01.maxent01symp
maxent01sub = maxent01.maxent01sub
maxent01vic = maxent01.maxent01vic
maxent01jump = maxent01.maxent01jump

total_numareas = 4
maxent01 = StateSpace.maxent01_defaults(total_numareas; maxent_constraint_01v=0.5);
maxent01symp = maxent01.maxent01symp
maxent01sub = maxent01.maxent01sub
maxent01vic = maxent01.maxent01vic
maxent01jump = maxent01.maxent01jump

"""
function maxent01_defaults(total_numareas=4; maxent_constraint_01=0.0, maxent_constraint_01v=0.0)
	#maxent_constraint_01 = 0.0
	maxent01symp = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
	maxent01sub = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
	maxent01jump = relative_probabilities_of_subsets(total_numareas, maxent_constraint_01)
	#maxent_constraint_01v = 0.0
	maxent01vic = relative_probabilities_of_vicariants(total_numareas, maxent_constraint_01v)
	maxent01 = (maxent01symp=maxent01symp, maxent01sub=maxent01sub, maxent01vic=maxent01vic, maxent01jump=maxent01jump)
	return maxent01
end




"""
numareas = 6
total_numareas=6
maxent_constraint_01=0.0
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

total_numareas=6
maxent_constraint_01=0.5
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

total_numareas=6
maxent_constraint_01=1.0
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

numareas = 2
total_numareas=6
maxent_constraint_01=0.0
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

total_numareas=6
maxent_constraint_01=0.5
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

total_numareas=6
maxent_constraint_01=1.0
discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas, maxent_constraint_01)

"""
function discrete_maxent_distrib_of_smaller_daughter_ranges(total_numareas=6, maxent_constraint_01=0.5)
	n = total_numareas
	x = 1:n
	
	if (maxent_constraint_01 < 0.0)
		maxent_constraint_01 = 0.0
	end
	if (maxent_constraint_01 > 1.0)
		maxent_constraint_01 = 1.0
	end
	
	#discrete_values_padded = cat(0, collect(x)[:], n+1; dims=1)
	#discrete_values_padded = collect(x)
	#maxent_constraint = quantile(discrete_values_padded, maxent_constraint_01)
	#maxent_constraint = 1
	
	maxent_result = round.(maxent_interp[n](maxent_constraint_01), digits=4)
	return maxent_result
end # END function discrete_maxent_distrib_of_smaller_daughter_ranges







"""
######################################################
Weight the possible range-sizes of the smaller daughter range

E.g., if you have subset sympatry from ancestor ABC,
the smaller daugher could be size 1, 2, or 3.

DEC says the weights of these would be 1 0 0
BayArea says the weights would be 0 0 1

With this function, input 
maxent_constraint_01=0.0001 = 1 0 0 
maxent_constraint_01=0.5    = 1 1 1 
maxent_constraint_01=0.9999 = 0 0 1 
######################################################

total_numareas=6
maxent_constraint_01 = 0.0001    # ranges from 0.0001 (all weight on ranges of size 1)
						                     # to 0.9999 (all weight on ranges of max size)
NA_val=NaN
relative_probabilities_of_subsets(total_numareas, maxent_constraint_01, NA_val)
"""

function relative_probabilities_of_subsets(total_numareas=6, maxent_constraint_01=0.5, NA_val=NaN)
	# Set up a matrix to hold the maxent distributions of relative prob of 
	# smaller daughter range sizes
	relprob_subsets_matrix = reshape(repeat([NA_val], total_numareas^2), (total_numareas,total_numareas))
	for i in 1:total_numareas
		ancestor_range_size = i
		maxent_result = discrete_maxent_distrib_of_smaller_daughter_ranges(ancestor_range_size, maxent_constraint_01)

		if (i == 1)
			relprob_subsets_matrix[i,1:i] .= maxent_result
		else
			relprob_subsets_matrix[i,1:i] = maxent_result
		end
	end
	
	return relprob_subsets_matrix
end


"""
# R version:
library(BioGeoBEARS)
relative_probabilities_of_vicariants(max_numareas=3, maxent_constraint_01v=0.0001)
#      [,1] [,2] [,3]
# [1,]   NA   NA   NA
# [2,]    1   NA   NA
# [3,]    1   NA   NA
relative_probabilities_of_vicariants(max_numareas=3, maxent_constraint_01v=0.5)
#      [,1] [,2] [,3]
# [1,]   NA   NA   NA
# [2,]    1   NA   NA
# [3,]    1   NA   NA
relative_probabilities_of_vicariants(max_numareas=4, maxent_constraint_01v=0.5)
#      [,1] [,2] [,3] [,4]
# [1,]   NA   NA   NA   NA
# [2,]  1.0   NA   NA   NA
# [3,]  1.0   NA   NA   NA
# [4,]  0.5  0.5   NA   NA
relative_probabilities_of_vicariants(max_numareas=5, maxent_constraint_01v=0.5)
#      [,1] [,2] [,3] [,4] [,5]
# [1,]   NA   NA   NA   NA   NA
# [2,]  1.0   NA   NA   NA   NA
# [3,]  1.0   NA   NA   NA   NA
# [4,]  0.5  0.5   NA   NA   NA
# [5,]  0.5  0.5   NA   NA   NA
relative_probabilities_of_vicariants(max_numareas=6, maxent_constraint_01v=0.5)
#       [,1]  [,2]  [,3] [,4] [,5] [,6]
# [1,]    NA    NA    NA   NA   NA   NA
# [2,] 1.000    NA    NA   NA   NA   NA
# [3,] 1.000    NA    NA   NA   NA   NA
# [4,] 0.500 0.500    NA   NA   NA   NA
# [5,] 0.500 0.500    NA   NA   NA   NA
# [6,] 0.333 0.333 0.333   NA   NA   NA

# Julia
relative_probabilities_of_vicariants(3, 0.0001)
relative_probabilities_of_vicariants(3, 0.5)
relative_probabilities_of_vicariants(4, 0.0001)
relative_probabilities_of_vicariants(4, 0.5)
relative_probabilities_of_vicariants(5, 0.0001)
relative_probabilities_of_vicariants(5, 0.5)
relative_probabilities_of_vicariants(6, 0.0001)
relative_probabilities_of_vicariants(6, 0.5)

"""
function relative_probabilities_of_vicariants(total_numareas=4, maxent_constraint_01=0.5, NA_val=NaN)
	# Set up a matrix to hold the maxent distributions of relative prob of 
	# smaller daughter range sizes
	relprob_subsets_matrix = reshape(repeat([NA_val], total_numareas^2), (total_numareas,total_numareas))
	relprob_subsets_matrix[1,:] .= NA_val
	
	for i in 2:total_numareas
		ancestor_range_size = i
		tmpstates = collect(1:i)# .+ 0.0
		tmpstates_Floats = collect(1:i) .+ 0.0
		max_smaller_rangesize = median(tmpstates_Floats)
		possible_vicariance_smaller_rangesizes = tmpstates[tmpstates_Floats .< max_smaller_rangesize]
		vic_total_numareas = maximum(possible_vicariance_smaller_rangesizes)

		maxent_result = flat2(discrete_maxent_distrib_of_smaller_daughter_ranges(vic_total_numareas, maxent_constraint_01))

		if (i <= 3)
# 			print("\n")
# 			print(i)
# 
# 			print("\n")
# 			print(vic_total_numareas)
# 			
# 			print("\n")
# 			print(relprob_subsets_matrix[i,1:vic_total_numareas])
# 
# 			print("\n")
# 			print(maxent_result)
# 
# 			print("\n")
			relprob_subsets_matrix[i,1:vic_total_numareas] .= maxent_result
		else
			relprob_subsets_matrix[i,1:vic_total_numareas] = maxent_result
		end
	end
	
	return relprob_subsets_matrix
end









end # END module MaxentInterp
