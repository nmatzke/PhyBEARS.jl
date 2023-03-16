#######################################################
# Tutorial: On an example simulated dataset
#######################################################

using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using DataFrames
using CSV

using PhyBEARS
#using PhyBEARS.Parsers


# Change the working directory as needed
wd = "/GitHub/PhyBEARS.jl/sims/sunkNZ_v1/"
cd(wd)

# This simulation has 148 living species
trfn = "living_tree_noNodeLabels.newick"
tr = readTopology(trfn)
trdf = prt(tr);
oldest_possible_age = 125.0

lgdata_fn = "geog_living.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
include_null_range = false
numareas = Rncol(geog_df)-1
max_range_size = numareas
n = numstates_from_numareas(numareas, max_range_size, include_null_range)

# DEC-type SSE model on Hawaiian Psychotria
# We are setting "j" to 0.0, for now -- so, no jump dispersal
bmo = construct_BioGeoBEARS_model_object();
bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "deathRate"] .= 0.9*ML_yule_birthRate(tr);
bmo.est[bmo.rownames .== "d"] .= 0.034;
bmo.est[bmo.rownames .== "e"] .= 0.028;
bmo.est[bmo.rownames .== "a"] .= 0.0;
bmo.est[bmo.rownames .== "j"] .= 0.1;
bmo.est[bmo.rownames .== "u"] .= -1.0;
bmo.min[bmo.rownames .== "u"] .= -2.5;
bmo.max[bmo.rownames .== "u"] .= 0.0;

bmo.type[bmo.rownames .== "j"] .= "free";
bmo.type[bmo.rownames .== "u"] .= "free";
bmo.type[bmo.rownames .== "birthRate"] .= "free";
bmo.type[bmo.rownames .== "deathRate"] .= "birthRate";



# Set up the model
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=include_null_range, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

bmo.est[:] = bmo_updater_v2(bmo, inputs.setup.bmo_rows);

inputs.setup.txt_states_list
