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
#using Distributed			# for e.g. Distributed.@spawn

print("...done.\n")


export area_of_areas_df_to_vectors, get_area_of_range, get_area_of_range_using_interpolator, update_Qij_e_vals!, update_Qij_d_vals!, update_elist_t!, Qij_e_vals_t!(p), Qij_d_vals_t!(p), update_dmat!


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

get_areas_of_range = x -> get_area_of_range_using_interpolator(x, state_as_areas_list, area_of_areas_interpolator)
tvals = seq(18.0, 23.0, 0.25)
get_areas_of_range.(tvals)

# Calculate extinction rate multipliers
state_as_areas_list = [1]
get_areas_of_range = x -> get_area_of_range_using_interpolator(x, state_as_areas_list, area_of_areas_interpolator)
tvals = seq(18.0, 23.0, 0.25)
uval = -1.0
get_extinction_rate_multipliers = x -> get_extinction_rate_multiplier(x, uval, state_as_areas_list, area_of_areas_interpolator)
get_extinction_rate_multipliers.(tvals)
"""
function get_area_of_range(tval, state_as_areas_list, area_of_areas)
	num_areas = length(state_as_areas_list)
	total_area = 0.0
	@simd @inbounds for i in 1:num_areas
		total_area += area_of_areas[state_as_areas_list[i]]
	end
	return total_area
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

get_areas_of_range = x -> get_area_of_range_using_interpolator(x, state_as_areas_list, area_of_areas_interpolator)
tvals = seq(18.0, 23.0, 0.25)
get_areas_of_range.(tvals)

# Calculate extinction rate multipliers
state_as_areas_list = [1]
get_areas_of_range = x -> get_area_of_range_using_interpolator(x, state_as_areas_list, area_of_areas_interpolator)
tvals = seq(18.0, 23.0, 0.25)
uval = -1.0
get_extinction_rate_multipliers = x -> get_extinction_rate_multiplier(x, uval, state_as_areas_list, area_of_areas_interpolator)
get_extinction_rate_multipliers.(tvals)
"""
function get_area_of_range_using_interpolator(tval, state_as_areas_list, area_of_areas_interpolator)
	num_areas = length(state_as_areas_list)
	total_area = 0.0
	#@inbounds @simd for i in 1:num_areas
	@inbounds for i in 1:num_areas
		total_area += area_of_areas_interpolator(tval)[state_as_areas_list[i]]
	end
	return total_area
end


"""
# 1. Update parameter 'e' based on 'u' and time 't'
# 2. Propagate that through the Q matrix, in the form of updating the Qij_vals

# Adding:
# Update the 'd' rates based on the distance between areas
"""
function update_Qij_e_vals!(p)
	# Update elist_t, i.e. the base_e_rate * multiplier on the base e rate (at time t)
	#p = p_Ds_v10
	#tval = 5.1
	#@inbounds @simd for i in 1:p.setup.numareas
	#	p.setup.elist_t[i] = p.setup.elist_base[i] * p.setup.area_of_areas[i]^p.bmo.est[p.setup.u_e_row[i]]
	#end
	
	#p.setup.elist_t .= p.setup.elist_base .* p.setup.area_of_areas.^p.bmo.est[p.setup.bmo_rows.u_e]
#	@simd @inbounds for i in 1:length(p.setup.elist_t)
#		p.setup.elist_t[i] = p.setup.elist_base[i] * p.setup.area_of_areas[i]^p.bmo.est[p.setup.bmo_rows.u_e]
#	end
	update_elist_t!(p)
	
	# Update the Qmat, using elist_t
	#prtQp(p)
	#Rnames(p.p_indices)
	# PRE-ALLOCATE THIS FOR SPEED
	#e_rows = (1:length(p.p_indices.Qarray_event_types))[p.p_indices.Qarray_event_types .== "e"]
#	@simd @inbounds for i in 1:length(p.setup.e_rows)
		#starting_statenum = p.p_indices.Qarray_ivals[p.setup.e_rows[i]]
		#ending_statenum = p.p_indices.Qarray_jvals[p.setup.e_rows[i]]
		#area_lost = symdiff(p.setup.states_list[starting_statenum], p.setup.states_list[ending_statenum])
		
		# PRECALCULATE THE AREA THAT WAS LOST (AND GAINED)
		# area_lost = symdiff(p.setup.states_list[p.p_indices.Qarray_ivals[p.setup.e_rows[i]]], p.setup.states_list[p.p_indices.Qarray_jvals[p.setup.e_rows[i]]])
		#area_lost = p.setup.losses[p.setup.e_rows[i]] 
		
		# actual rate of e = base_rate_of_e * area_of_area_lost ^ u
		#p.params.Qij_vals_t[p.setup.e_rows[i]] = p.params.Qij_vals[p.setup.e_rows[i]] * p.setup.elist_t[area_lost][]
		#p.params.Qij_vals[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i][1]] ][1]
#		p.params.Qij_vals_t[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i][1]] ][1]
#	end
	Qij_e_vals_t!(p)
	
#	@inbounds @simd for i in 1:p.setup.num_e_rows
		#starting_statenum = p.p_indices.Qarray_ivals[p.setup.e_rows[i]]
		#ending_statenum = p.p_indices.Qarray_jvals[p.setup.e_rows[i]]
		#area_lost = symdiff(p.setup.states_list[starting_statenum], p.setup.states_list[ending_statenum])
		
		# PRECALCULATE THE AREA THAT WAS LOST (AND GAINED)
		# area_lost = symdiff(p.setup.states_list[p.p_indices.Qarray_ivals[p.setup.e_rows[i]]], p.setup.states_list[p.p_indices.Qarray_jvals[p.setup.e_rows[i]]])
		#area_lost = p.setup.losses[p.setup.e_rows[i]] 
		
		# actual rate of e = base_rate_of_e * area_of_area_lost ^ u
		#p.params.Qij_vals_t[p.setup.e_rows[i]] = p.params.Qij_vals[p.setup.e_rows[i]] * p.setup.elist_t[area_lost][]
		# Slower
		#p.params.Qij_vals[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i][1]] ][1]
		#p.params.Qij_vals_t[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i][1]] ][1]
		# v2:
		#p.params.Qij_vals[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i]] ][1]
#		p.params.Qij_vals_t[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i]] ][1]
#	end
	
#	p.params.Qij_vals_t[p.setup.e_rows] .= p.setup.elist_t[p.setup.losses_by_e_rows]
	
	return(p)
end # END function update_Qij_e_vals!(p)


function update_elist_t!(p)
	@simd @inbounds for i in 1:length(p.setup.elist_t)
		p.setup.elist_t[i] = p.setup.elist_base[i] * p.setup.area_of_areas[i]^p.bmo.est[p.setup.bmo_rows.u_e]
	end
end

function Qij_e_vals_t!(p)
	@simd @inbounds for i in 1:length(p.setup.e_rows)
		p.params.Qij_vals_t[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i][1]] ][1]
	end
end

function Qij_d_vals_t!(p)
	@simd @inbounds for i in 1:length(p.setup.d_drows)
		p.params.Qij_vals_t[p.setup.d_drows[i]] += p.setup.dmat[p.setup.d_froms[i], p.setup.d_tos[i]]
	end
end


function update_dmat!(p)
	@simd @inbounds for i in 1:length(p.setup.dmat)
		p.setup.dmat[i] = p.setup.dmat_base[i] * p.setup.dispersal_multipliers_mat[i]^p.bmo.est[p.setup.bmo_rows.w] * p.setup.distmat[i]^p.bmo.est[p.setup.bmo_rows.x] * p.setup.envdistmat[i]^p.bmo.est[p.setup.bmo_rows.n] * p.setup.distmat2[i]^p.bmo.est[p.setup.bmo_rows.x2] * p.setup.distmat3[i]^p.bmo.est[p.setup.bmo_rows.x3]
	end
end



"""
# 1. Update parameter 'd' based on a distance matrix (already determined for time 't')
# 2. Propagate that through the Q matrix, in the form of updating the Qij_vals
# (save multiple distance matrices for later)
"""
function update_Qij_d_vals!(p)
	
	# dmat: a numareas x numarea matrix, containing the total product of all the 
	#       dispersal multipliers (time-dependent or not)
	# initially: just based on distance
	#p.setup.dmat_t .= p.setup.dmat_base .* p.setup.distmat.^p.bmo.est[p.setup.bmo_rows.x]
	
	# dmat_base contains the original d values **between pairs of areas***
	#p.setup.dmat .= p.setup.dmat_base .* p.setup.dispersal_multipliers_mat.^p.bmo.est[p.setup.bmo_rows.w] .* p.setup.distmat.^p.bmo.est[p.setup.bmo_rows.x] .* p.setup.envdistmat.^p.bmo.est[p.setup.bmo_rows.n] .* p.setup.distmat2.^p.bmo.est[p.setup.bmo_rows.x2] .* p.setup.distmat3.^p.bmo.est[p.setup.bmo_rows.x3]
	
	#p.setup.dmat .= 0.0
	update_dmat!(p)
	
	#p.setup.amat_t .= p.setup.amat_base .* p.setup.dispersal_multipliers_mat.^p.bmo.est[p.setup.bmo_rows.w] .* p.setup.distmat.^p.bmo.est[p.setup.bmo_rows.x] .* p.setup.envdistmat.^p.bmo.est[p.setup.bmo_rows.n] .* p.setup.distmat2.^p.bmo.est[p.setup.bmo_rows.x2] .* p.setup.distmat3.^p.bmo.est[p.setup.bmo_rows.x3]
	
	# Update the Qmat, using elist_t
	#prtQp(p)
	#Rnames(p.p_indices)
	# PRE-ALLOCATE THIS FOR SPEED
	#e_rows = (1:length(p.p_indices.Qarray_event_types))[p.p_indices.Qarray_event_types .== "e"]
	
	#p.params.Qij_vals[p.setup.d_rows] .= 0.0
	p.params.Qij_vals_t[p.setup.d_rows] .= 0.0
	#@simd @inbounds for i in 1:length(p.setup.d_drows)
		#starting_statenum = p.p_indices.Qarray_ivals[p.setup.e_rows[i]]
		#ending_statenum = p.p_indices.Qarray_jvals[p.setup.e_rows[i]]
		#area_lost = symdiff(p.setup.states_list[starting_statenum], p.setup.states_list[ending_statenum])
		
		# PRECALCULATE THE AREA THAT WAS LOST (AND GAINED)
		# area_lost = symdiff(p.setup.states_list[p.p_indices.Qarray_ivals[p.setup.e_rows[i]]], p.setup.states_list[p.p_indices.Qarray_jvals[p.setup.e_rows[i]]])
		#area_lost = p.setup.losses[p.setup.e_rows[i]] 
		
		# actual rate of e = base_rate_of_e * area_of_area_lost ^ u
		#p.params.Qij_vals_t[p.setup.e_rows[i]] = p.params.Qij_vals[p.setup.e_rows[i]] * p.setup.elist_t[area_lost][]
		
		# area moving from: p.setup.d_froms[i]
		# area moving to: p.setup.d_tos[i]
		
		#p.params.Qij_vals[p.setup.d_drows[i]] += p.setup.dmat[p.setup.d_froms[i], p.setup.d_tos[i]]
	#	p.params.Qij_vals_t[p.setup.d_drows[i]] += p.setup.dmat[p.setup.d_froms[i], p.setup.d_tos[i]]
	#end
	Qij_d_vals_t!(p)
	
#	@inbounds @simd for i in 1:p.setup.num_e_rows
		#starting_statenum = p.p_indices.Qarray_ivals[p.setup.e_rows[i]]
		#ending_statenum = p.p_indices.Qarray_jvals[p.setup.e_rows[i]]
		#area_lost = symdiff(p.setup.states_list[starting_statenum], p.setup.states_list[ending_statenum])
		
		# PRECALCULATE THE AREA THAT WAS LOST (AND GAINED)
		# area_lost = symdiff(p.setup.states_list[p.p_indices.Qarray_ivals[p.setup.e_rows[i]]], p.setup.states_list[p.p_indices.Qarray_jvals[p.setup.e_rows[i]]])
		#area_lost = p.setup.losses[p.setup.e_rows[i]] 
		
		# actual rate of e = base_rate_of_e * area_of_area_lost ^ u
		#p.params.Qij_vals_t[p.setup.e_rows[i]] = p.params.Qij_vals[p.setup.e_rows[i]] * p.setup.elist_t[area_lost][]
		# Slower
		#p.params.Qij_vals[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i][1]] ][1]
		#p.params.Qij_vals_t[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i][1]] ][1]
		# v2:
		#p.params.Qij_vals[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i]] ][1]
#		p.params.Qij_vals_t[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i]] ][1]
#	end
	
#	p.params.Qij_vals_t[p.setup.e_rows] .= p.setup.elist_t[p.setup.losses_by_e_rows]
	
	return(p)
end # END function update_Qij_d_vals!(p)



end # END TimeDep