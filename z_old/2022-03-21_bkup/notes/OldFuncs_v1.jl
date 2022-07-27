module OldFuncs
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

using BioGeoJulia.TrUtils # for e.g. flat2
using Combinatorics  # for e.g. combinations()
using DataFrames     # for e.g. DataFrame()
using PhyloNetworks

#######################################################
# Old functions that have been superceded
#######################################################
using Distributions  # for quantile
using Convex				 # for Convex.entropy(), maximize()
using SCS						 # for SCSSolve, solve (maximize(entropy()))

# This worked on 2022-03-12 with old versions:
	# ...assumes Pkg.add(name="SCS", version="0.9.0")
	# ...assumes Pkg.add(name="Convex", version="0.14.18")

# Check versions with
Pkg.status("SCS")


export discrete_maxent_distrib_of_smaller_daughter_ranges_OLD


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
function discrete_maxent_distrib_of_smaller_daughter_ranges_OLD(total_numareas=6, maxent_constraint_01=0.5)
	n = total_numareas
	x = 1:n

	#discrete_values_padded = cat(0, collect(x)[:], n+1; dims=1)
	discrete_values_padded = collect(x)
	maxent_constraint = quantile(discrete_values_padded, maxent_constraint_01)
	#maxent_constraint = 1

	# Vector of probabilities that must sum to 1
	p = Variable(length(x))
	probability_contraints = [0 <= p, p <= 1, sum(p) == 1];
	feature_constraints = sum(p'*x) == maxent_constraint

	# base or prior (e.g. Uniform)
	#h = pdf.(Uniform(1, n), x)
	
	# This solution updates p (look in p.values)
	problem = Convex.maximize(Convex.entropy(p), 0 <= p, p <= 1, sum(p) == 1, feature_constraints)
	
	# worked 2020-09-29 at work, failed 2020-09-30 at home
	#sol = Convex.solve!(problem, SCS.SCSSolver(verbose=0))
	
	# worked e.g. July 2020 (version issue I guess)
	# worked at home 2020-09-30
	# worked at home & laptop, 2022-03-11
	sol = Convex.solve!(problem, SCS.Optimizer(verbose=0))
	# ...assumes Pkg.add(name="SCS", version="0.9.0")
	# ...assumes Pkg.add(name="Convex", version="0.14.18")

	maxent_result = abs.(round.(p.value; digits=4))
	return maxent_result
end # END function discrete_maxent_distrib_of_smaller_daughter_ranges
