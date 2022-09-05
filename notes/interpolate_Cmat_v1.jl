
#######################################################
# Interpolate a Carray matrix
#
# Notes: Saving interpolated function into a separate file in Julia
# https://stackoverflow.com/questions/58499330/saving-interpolated-function-into-a-separate-file-in-julia
#######################################################

using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using PhyBEARS
using DataFrames
using CSV

trfn = "/GitHub/PhyBEARS.jl/data/Psychotria_tree.newick"
tr = readTopology(trfn)

lgdata_fn = "/GitHub/PhyBEARS.jl/data/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
include_null_range = false
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "birthRate"] .= 0.1
bmo.est[bmo.rownames .== "deathRate"] .= 0.01
bmo.est[bmo.rownames .== "d"] .= 0.02
bmo.est[bmo.rownames .== "e"] .= 0.02
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.11
bmo.est[:] = bmo_updater_v1(bmo);

# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Ds_v5, Es_tspan) = inputs;

solver_options.solver = CVODE_BDF(linear_solver=:GMRES);
solver_options.save_everystep = true;
solver_options.abstol = 1e-12;
solver_options.reltol = 1e-12;

p_Ds_v5 = inputs.p_Ds_v5;
p_Ds_v5_updater_v1!(p_Ds_v5, inputs);

prtCp(p_Ds_v5)


#######################################################
# OK let's make a function that interpolates different rates
# based on changing distances / areas
#######################################################
times = [0.0, 4.0, 20.0, 22.0, 30.0, 80.0, 600.0]
# current NZ area = 268,021 km2
areas = [1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 10.0]

area_interpolator = interpolate((times,), areas, Gridded(Linear()))
tvals = seq(0.0, 40, 1.0)
area_interpolator(tvals)


#######################################################
# Area interpolator for a series of states
#######################################################
areas_list = inputs.setup.areas_list
states_list = inputs.setup.states_list
area_of_areas_file = "/GitHub/PhyBEARS.jl/notes/area_of_areas_NZ_Oz_v1.txt"
area_of_areas_table = CSV.read(area_of_areas_file, DataFrame; delim="\t")
area_of_areas_table = readtable(area_of_areas_file)

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




# Now make an extinction_rate_interpolator
function get_extinction_rate_multiplier(tval, uval, state_as_areas_list, area_of_areas_interpolator)
	extinction_rate_multiplier = 1.0 * (get_area_of_range(tval, state_as_areas_list, area_of_areas_interpolator) ^ uval)
	# Error check
	extinction_rate_multiplier = minimum([extinction_rate_multiplier, 10000.0])
	return(extinction_rate_multiplier)
end

get_extinction_rate_multipliers = x -> get_extinction_rate_multiplier(x, uval, state_as_areas_list, area_of_areas_interpolator)

# https://docs.juliahub.com/JLD2/O1EyT/0.1.13/
cd("/GitHub/PhyBEARS.jl/src")
@save "list_of_interpolators.jld2" list_of_interpolators

# Reloading
using JLD2
using Interpolations # for Interpolations.scale, Interpolations.interpolate
@load "list_of_interpolators.jld2" list_of_interpolators # <- It loads to "list_of_interpolators"...has to be this

# Find location of a module
pathof(PhyBEARS)
interp_jld2_path = join([pathof(PhyBEARS), ])

itp = interpolate((x,), y, Gridded(Linear()))

