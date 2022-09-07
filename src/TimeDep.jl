module TimeDep
#############################################################
# Functions for time-dependency in the Qmat and Cmat
#############################################################

__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

print("PhyBEARS: loading TimeDep.jl dependencies...")

using PhyloBits.PNtypes	# for e.g. HybridNetwork
using PhyloBits.TrUtils # for e.g. flat2
using PhyloBits.TreeTable	# for e.g. get_nonrootnodes_trdf
using PhyBEARS.SSEs 
#using PhyBEARS.Flow 
using Hwloc						# for e.g. getinfo()[:Core] # the current core, I think
#using Distributed			# for e.g. Distributed.@spawn

print("...done.\n")


export area_of_areas_df_to_vectors, get_area_of_range




"""
# Area-of-areas table to a vector of vectors
# area_of_areas_file = "/GitHub/PhyBEARS.jl/notes/area_of_areas_NZ_Oz_v1.txt"
# area_of_areas_table = CSV.read(area_of_areas_file, DataFrame; delim="\t")
# repr(Matrix(area_of_areas_table))
tmpstr = "[0.0 1.0 28.6; 4.0 1.0 28.6; 20.0 0.0 28.6; 22.0 0.0 28.6; 30.0 1.0 28.6; 80.0 1.0 28.6; 600.0 1.0 28.6]"

# area_of_areas_table = Reval(tmpstr)
area_of_areas_matrix = eval(Meta.parse(tmpstr))
area_of_areas_table = DataFrame(area_of_areas_matrix, ["times","A","B"])
area_of_areas_vec = area_of_areas_df_to_vectors(area_of_areas_table)
area_of_areas_vec = area_of_areas_df_to_vectors(area_of_areas_table; cols=2)

#area_of_areas_vec = area_of_areas_df_to_vectors(area_of_areas_table; cols=2)
area_of_areas_vec = area_of_areas_df_to_vectors(area_of_areas_table)
area_of_areas_interpolator = interpolate((area_of_areas_table.times,), area_of_areas_vec, Gridded(Linear()))

"""
function area_of_areas_df_to_vectors(area_of_areas_table; cols=NaN)
	numtimes = Rnrow(area_of_areas_table)

	if isnan(cols)
		numcols = Rncol(area_of_areas_table)
		num_areas = numcols - 1
		cols = 2:numcols
	else
		numcols = length(cols)
		num_areas = numcols
	end

	area_of_areas_vec = [Vector{Float64}(undef, num_areas) for _ = 1:numtimes]
	for i in 1:numtimes
		area_of_areas_vec[i] .= flat2(area_of_areas_table[i, cols])
	end
	return(area_of_areas_vec)
end



"""
# Get the area of a range at time t
# actual_extinction_rate = base_rate * range_size^u

# Load area_of_areas_table (DataFrame)
tmpstr = "[0.0 1.0 28.6; 4.0 1.0 28.6; 20.0 0.0 28.6; 22.0 0.0 28.6; 30.0 1.0 28.6; 80.0 1.0 28.6; 600.0 1.0 28.6]"

# area_of_areas_table = Reval(tmpstr)
area_of_areas_matrix = eval(Meta.parse(tmpstr))
area_of_areas_table = DataFrame(area_of_areas_matrix, ["times","A","B"])
#area_of_areas_vec = area_of_areas_df_to_vectors(area_of_areas_table; cols=2)
area_of_areas_vec = area_of_areas_df_to_vectors(area_of_areas_table)
area_of_areas_interpolator = interpolate((area_of_areas_table.times,), area_of_areas_vec, Gridded(Linear()))

# Calculate the area at different times
tval = 19.0
state_as_areas_list = [1,2]
total_area = get_area_of_range(tval, state_as_areas_list, area_of_areas_interpolator)

get_areas_of_range = x -> get_area_of_range(x, state_as_areas_list, area_of_areas_interpolator)
tvals = seq(18.0, 23.0, 0.25)
get_areas_of_range.(tvals)

# Calculate extinction rate multipliers
state_as_areas_list = [1]
get_areas_of_range = x -> get_area_of_range(x, state_as_areas_list, area_of_areas_interpolator)
tvals = seq(18.0, 23.0, 0.25)
uval = -1.0
get_extinction_rate_multipliers = x -> get_extinction_rate_multiplier(x, uval, state_as_areas_list, area_of_areas_interpolator)
get_extinction_rate_multipliers.(tvals)
"""
function get_area_of_range(tval, state_as_areas_list, area_of_areas_interpolator)
	num_areas = length(state_as_areas_list)
	total_area = 0.0
	for i in 1:num_areas
		total_area += area_of_areas_interpolator(tval)[state_as_areas_list[i]]
	end
	return total_area
end




end # END TimeDep