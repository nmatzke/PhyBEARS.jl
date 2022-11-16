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

using Interpolations	# for Linear, Gridded, interpolate

using PhyloBits.PNtypes	# for e.g. HybridNetwork
using PhyloBits.TrUtils # for e.g. flat2
using PhyloBits.TreeTable	# for e.g. get_nonrootnodes_trdf
#using Distributed			# for e.g. Distributed.@spawn

print("...done.\n")


export area_of_areas_df_to_vectors, get_area_of_range, get_area_of_range_using_interpolator, update_Qij_e_vals!, update_Qij_d_vals!, get_elist_at_time_t!, update_Qij_e_vals_t!, update_Qij_d_vals_t!, get_dmat_at_time_t!, update_min_vdist_at_time_t_withp!, update_min_vdist_at_time_t!, get_mindist_between_pair_of_ranges, get_jmat_at_time_t!, update_Cijk_j_rates_t!, update_Cijk_j_rates!, update_Cijk_rates!, update_Cijk_v_rates!, update_Cijk_rates_sub_i_t!, update_mus_time_t!, update_QC_mats_time_t!, construct_QC_interpolators


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
	@inbounds @simd for i in 1:num_areas
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
#	@inbounds @simd for i in 1:length(p.setup.elist_t)
#		p.setup.elist_t[i] = p.setup.elist_base[i] * p.setup.area_of_areas[i]^p.bmo.est[p.setup.bmo_rows.u_e]
#	end
	get_elist_at_time_t!(p)
	
	# Update the Qmat, using elist_t
	#prtQp(p)
	#Rnames(p.p_indices)
	# PRE-ALLOCATE THIS FOR SPEED
	#e_rows = (1:length(p.p_indices.Qarray_event_types))[p.p_indices.Qarray_event_types .== "e"]
#	@inbounds @simd for i in 1:length(p.setup.e_rows)
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
	update_Qij_e_vals_t!(p)
	
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


function get_elist_at_time_t!(p)
	@inbounds @simd for i in 1:length(p.setup.elist_t)
		p.setup.elist_t[i] = p.setup.elist_base[i] * p.setup.area_of_areas[i]^p.bmo.est[p.setup.bmo_rows.u_e]
	end
end

function update_Qij_e_vals_t!(p)
	@inbounds @simd for i in 1:length(p.setup.e_rows)
		p.params.Qij_vals_t[p.setup.e_rows[i]] = p.setup.elist_t[p.setup.losses[p.setup.e_rows[i][1]] ][1]
	end
end

function update_dmat_at_time_t_t!(p)
	@inbounds @simd for i in 1:length(p.setup.dmat_t)
		p.setup.dmat_t[i] = p.setup.dmat_base[i] * p.setup.dispersal_multipliers_mat[i]^p.bmo.est[p.setup.bmo_rows.w] * p.setup.distmat[i]^p.bmo.est[p.setup.bmo_rows.x] * p.setup.envdistmat[i]^p.bmo.est[p.setup.bmo_rows.n] * p.setup.distmat2[i]^p.bmo.est[p.setup.bmo_rows.x2] * p.setup.distmat3[i]^p.bmo.est[p.setup.bmo_rows.x3]
	end
end

# Update the vicariance distances (e.g. minimum distance between two ranges) at time t
# (NOTE: dmat_t must be updated FIRST)
# 
# Here, we assume the minimum distance is relevant. Another option would be e.g. average distances
function update_min_vdist_at_time_t_withp!(p; mindist=1.0e9)
	@inbounds @simd for i in 1:length(p.setup.v_rows)
		# Get the ranges (pretty small)
		#jrange = p.setup.states_list[p.p_indices.Carray_jvals[setup.v_rows[i]]]
		#krange = p.setup.states_list[p.p_indices.Carray_kvals[setup.v_rows[i]]]

		# Save the minimum distances
		p.setup.vicdist_t[i] = get_mindist_between_pair_of_ranges(p, p.setup.states_list[p.p_indices.Carray_jvals[setup.v_rows[i]]], p.setup.states_list[p.p_indices.Carray_kvals[setup.v_rows[i]]]; mindist=mindist)
	end
end

function update_min_vdist_at_time_t!(vicdist_t, v_rows, distmat, states_list, Carray_jvals, Carray_kvals; mindist=1.0e9)
	@inbounds @simd for i in 1:length(v_rows)
		# Get the ranges (pretty small)
		#jrange = p.setup.states_list[p.p_indices.Carray_jvals[setup.v_rows[i]]]
		#krange = p.setup.states_list[p.p_indices.Carray_kvals[setup.v_rows[i]]]

		# Save the minimum distances
		vicdist_t[i] = get_mindist_between_pair_of_ranges(distmat, states_list[Carray_jvals[v_rows[i]]], states_list[Carray_kvals[v_rows[i]]]; mindist=mindist)
	end
end


function get_mindist_between_pair_of_ranges(distmat, jrange, krange; mindist=1.0e9)
	# Create a tight inner loop; check if each dist is lower
	@inbounds for j in 1:length(jrange)
		@inbounds for k in 1:length(krange)
			if distmat[jrange[j], krange[k]] < mindist
				mindist = distmat[jrange[j], krange[k]]
			end # END if
		end # end k
	end # end j
	
	return(mindist)
end

function update_Qij_d_vals_t!(p)
	@inbounds @simd for i in 1:length(p.setup.d_drows)
		p.params.Qij_vals_t[p.setup.d_drows[i]] += p.setup.dmat_t[p.setup.d_froms[i], p.setup.d_tos[i]]
	end
end


# 2022-09-26:
# 
# Update jmat_t, modifying jmat_base (which will typically be 1s)
# NOTE CAREFULLY: IN v10, jmat_t will modify the "val" of Cmat, i.e. the rate of
# jump dispersal as set AFTER the j_wt has been converted to a rate.
#
# This will match BioGeoBEARS in the default case, but not necessarily with
# changing distances etc.
# 
# (There is no other way to do it, any other method would require re-doing the j_wt-to-rates
#  calculation in each timestep t)
#
# (This will mean the meaning of j_wt is complex, but then it always was!)
function get_jmat_at_time_t!(p)
	@inbounds @simd for i in 1:length(p.setup.jmat_t)
		p.setup.jmat_t[i] = p.setup.jmat_base[i] * p.setup.dispersal_multipliers_mat[i]^p.bmo.est[p.setup.bmo_rows.w] * p.setup.distmat[i]^p.bmo.est[p.setup.bmo_rows.x] * p.setup.envdistmat[i]^p.bmo.est[p.setup.bmo_rows.n] * p.setup.distmat2[i]^p.bmo.est[p.setup.bmo_rows.x2] * p.setup.distmat3[i]^p.bmo.est[p.setup.bmo_rows.x3]
	end
end

function update_Cijk_j_rates_t!(p)
	p.params.Cijk_rates_t .= p.params.Cijk_rates # initialize at default rates 
	p.params.Cijk_rates_t[p.setup.j_rows] .= 0.0 # initialize j values at 0.0
	# Change the (relative) j rates, nothing else
	@inbounds @simd for i in 1:length(p.setup.j_jrows)
		# Add up the jmat_t modifiers for this dispersal event
		# Here, we are doing it fractionally by starting p.params.Cijk_rates[p.setup.j_jrows[i]]
		# The rate for a j event will be Cijk_rates * jmat_t
		# We take the Cijk_rates (not t-dependent), and give back a sum of the distance-weighted rates
		p.params.Cijk_rates_t[p.setup.j_jrows[i]] += p.params.Cijk_rates[p.setup.j_jrows[i]] * p.setup.jmat_t[p.setup.j_froms[i], p.setup.j_tos[i]] / p.setup.j_numdispersals[i]
	end	
end


function update_Cijk_j_rates!(p)
	get_jmat_at_time_t!(p)
	p.params.Cijk_rates_t[p.setup.j_rows] .= 0.0 # zero out Cijk_rates_t
	update_Cijk_j_rates_t!(p) # modify the (post-weights) jump dispersal rates
	update_Cijk_rates_sub_i_t!(p) # updates the pre-allocation of Cijk_rates_t, in Cijk_rates_sub_i_t, to simplify core SSE simd
end


function update_Cijk_rates!(p)
	get_jmat_at_time_t!(p)
	p.params.Cijk_rates_t[p.setup.j_rows] .= 0.0 # zero out Cijk_rates_t
	update_Cijk_j_rates_t!(p) # modify the (post-weights) jump dispersal rates
	update_Cijk_v_rates!(p)		# modify the (post-weights) vicariance rates
	update_Cijk_rates_sub_i_t!(p) # updates the pre-allocation of Cijk_rates_t, in Cijk_rates_sub_i_t, to simplify core SSE simd
end



# Vicariance
# Modify the default Cijk rates for vicariance
function update_Cijk_v_rates!(p)
	# Zero out rates_t before adding
	#p.params.Cijk_rates_t[p.setup.v_rows] .= 0.0
	@inbounds @simd for i in 1:length(p.setup.v_rows)
		# Modify the base vicariance rate by multiplying by mindistance ^ xv (where xv is positive)
		p.params.Cijk_rates_t[p.setup.v_rows[i]] = p.params.Cijk_rates[p.setup.v_rows[i]] * p.setup.vicdist_t[i]^p.bmo.est[p.setup.bmo_rows.xv]
	end
end

# Transfer Cijk rates to Cijk_rates_sub_i
function update_Cijk_rates_sub_i_t!(p)
	# Update the Cijk_rates_sub_i_t (where anc==i)
	@inbounds @simd for i in 1:length(p.setup.states_list)
		p.p_TFs.Cijk_rates_sub_i_t[i] .= p.params.Cijk_rates_t[p.p_TFs.Ci_eq_i[i]]
	end
end









"""
# 1. Update parameter 'd' based on a distance matrix (NOTE: MUST BE ALREADY DETERMINED for time 't')
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
	update_dmat_at_time_t_t!(p)
	
	#p.setup.amat_t .= p.setup.amat_base .* p.setup.dispersal_multipliers_mat.^p.bmo.est[p.setup.bmo_rows.w] .* p.setup.distmat.^p.bmo.est[p.setup.bmo_rows.x] .* p.setup.envdistmat.^p.bmo.est[p.setup.bmo_rows.n] .* p.setup.distmat2.^p.bmo.est[p.setup.bmo_rows.x2] .* p.setup.distmat3.^p.bmo.est[p.setup.bmo_rows.x3]
	
	# Update the Qmat, using elist_t
	#prtQp(p)
	#Rnames(p.p_indices)
	# PRE-ALLOCATE THIS FOR SPEED
	#e_rows = (1:length(p.p_indices.Qarray_event_types))[p.p_indices.Qarray_event_types .== "e"]
	
	#p.params.Qij_vals[p.setup.d_rows] .= 0.0
	p.params.Qij_vals_t[p.setup.d_rows] .= 0.0
	#@inbounds @simd for i in 1:length(p.setup.d_drows)
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
	update_Qij_d_vals_t!(p)
	
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


function update_mus_time_t!(p, t)
  #mu = p.params.mu_vals # base extinction rate of each range
  #mu_t = p.params.mu_t_vals # mu_t = mu at time t
  
  # Populate changing mus with time
  @inbounds @simd for i in 1:n
  	# total_area = get_area_of_range(tval, state_as_areas_list, area_of_areas_interpolator)
  	p.params.mu_t_vals[i] = p.params.mu_vals[i] * get_area_of_range(t, p.states_as_areas_lists[i], p.interpolators.area_of_areas_interpolator(t))^p.bmo.est[p.setup.bmo_rows.u_mu]
  end
  # Correct "Inf" max_extinction_rates
  p.params.mu_t_vals[p.params.mu_t_vals .> max_extinction_rate] .= max_extinction_rate

end



"""
Update the Qarray and Carray vectors at time t
"""
function update_QC_mats_time_t!(p, t)
 	# The row of bmo that refers to "u", the effect of area on extinction rate
 	# u_row = (1:Rnrow(bmo))[bmo.rownames .== "u"][]
 	max_extinction_rate = p.setup.max_extinction_rate

 	# Get the area of areas at time t
	p.setup.area_of_areas .= p.interpolators.area_of_areas_interpolator(t)
 	
	update_mus_time_t!(p, t)  
  
  # Get the e_vals for the Qij matrix, at time t
  update_Qij_e_vals!(p)
  # (updates p.params.Qij_vals)

 # Get the d_vals for the Qij matrix, at time t
  # 1. Update the distance matrices etc.
  p.setup.distmat .= p.interpolators.distances_interpolator(t)
  
  
  # Update the vicariance minimum distance 
  p.setup.vicdist_t .= p.interpolators.vicariance_mindists_interpolator(t)
    
  # Using the current t's distmat, etc. update the dmat_t, then 
  # propagate through the Q matrix
  update_Qij_d_vals!(p)
  
  # Using the current t's distmat, etc. update the jmat_t, then
  # propagate through the C matrix
  # Replaces: update_Cijk_j_rates!(p)
	# updates j events
  # updates vicariance also
  update_Cijk_rates!(p)

end # END function update_QC_mats_time_t!(p, t)





# Also interpolate mus
function construct_QC_interpolators(p, tvals)
	# Construct interpolators for Q_vals_t and C_rates_t
	mu_vals_by_t = [Vector{Float64}(undef, length(p.params.mu_t_vals)) for _ = 1:length(tvals)]
	Q_vals_by_t = [Vector{Float64}(undef, length(p.params.Qij_vals_t)) for _ = 1:length(tvals)]
	C_rates_by_t = [Vector{Float64}(undef, length(p.params.Cijk_rates_t)) for _ = 1:length(tvals)]

	for i in 1:length(tvals)
		# Zero out the rates
		Q_vals_by_t[i] .= 0.0
		C_rates_by_t[i] .= 0.0
	
		# Update the rates
		update_QC_mats_time_t!(p, tvals[i])
	
		# Save these rates
		mu_vals_by_t[i] = p.params.mu_t_vals
		Q_vals_by_t[i] .= p.params.Qij_vals_t
		C_rates_by_t[i] .= p.params.Cijk_rates_t
	end
	
	mu_vals_by_t
	Q_vals_by_t
	C_rates_by_t

	mu_vals_interpolator = interpolate((tvals,), mu_vals_by_t, Gridded(Linear()));
	Q_vals_interpolator = interpolate((tvals,), Q_vals_by_t, Gridded(Linear()));
	C_rates_interpolator = interpolate((tvals,), C_rates_by_t, Gridded(Linear()));
	
	interpolators = (times_for_dists_interpolator=p.interpolators.times_for_dists_interpolator, times_for_SSE_interpolators=p.interpolators.times_for_SSE_interpolators, distances_interpolator=p.interpolators.distances_interpolator, area_of_areas_interpolator=p.interpolators.area_of_areas_interpolator, vicariance_mindists_interpolator=p.interpolators.vicariance_mindists_interpolator, mu_vals_interpolator=mu_vals_interpolator, Q_vals_interpolator=Q_vals_interpolator, C_rates_interpolator=C_rates_interpolator)
	
	# Copy the tuple p
	p2 = (n=p.n, params=p.params, p_indices=p.p_indices, p_TFs=p.p_TFs, uE=p.uE, terms=p.terms, setup=p.setup, states_as_areas_lists=p.setup.states_list, interpolators=interpolators, use_distances=p.use_distances, bmo=p.bmo)
	
	return(p2)
end # END function construct_QC_interpolators(p, tvals)




end # END TimeDep