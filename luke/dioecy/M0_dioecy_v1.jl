using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames						# for DataFrame()
using DelimitedFiles				# for readdlm()
using NLopt									# seems to be the best gradient-free, box-constrained								

# List each PhyBEARS code file prefix here
using PhyloBits.TrUtils			# for e.g. numstxt_to_df()
using PhyloBits.TreeTable
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.ModelLikes # e.g. setup_DEC_SSE2
using PhyBEARS.Uppass

"""
# Run with:
cd("/GitHub/PhyBEARS.jl/luke/dioecy/")
include("M0_dioecy_v1.jl")
"""

cd("/GitHub/PhyBEARS.jl/luke/dioecy/")

trfn = "tree.newick"
tr = readTopology(trfn)

lgdata_fn = "geog.txt"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
rowSums_df(geog_df)
maximum(rowSums(geog_df[:,2]))
minimum(rowSums(geog_df))

# Add trait to bmo
geotrait_fn = "geog+trait.txt"
geotrait_df = Parsers.getranges_from_LagrangePHYLIP(geotrait_fn)

bmo = construct_BioGeoBEARS_model_object();

# Set up the model
numareas = 2;
inputs = ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;
(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs

prtCp(p_Ds_v5)
prtQp(p_Ds_v5)

# Add trait to bmo
num_trait_states = 2

rownames = ["t12", "t21", "m1", "m2", "m3", "m4"]
type_vec = ["fixed", "fixed", "fixed", "fixed", "fixed", "fixed"]
init_vec = [0.01, 0.01, 0.0, 0.0, 0.0, 0.0]
min_vec = [1.0e-12, 1.0e-12, 0.0, 0.0, 0.0, 0.0]
max_vec = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
est_vec = [0.01, 0.01, 0.0, 0.0, 0.0, 0.0]
note_vec = ["devel", "devel", "devel", "devel", "devel", "devel"]
desc_vec = ["base transition rate from H (hermaphroditic/monoecious) to D (dioecious)", "base transition rate from D (dioecious) to H (hermaphroditic/monoecious)", "multiplier on t12 when in range L", "multiplier on t21 when in range L", "multiplier on t12 when in range M", "multiplier on t21 when in range M"]

bmo_add = DataFrames.DataFrame(rownames=rownames, type=type_vec, init=init_vec, min=min_vec, max=max_vec, est=est_vec, note=note_vec, desc=desc_vec)


# Add to the state space

# Better: include the trait as another state, 


# In this case, we are multiplying the geography state space by 2
# Start by duplicating everything
 mu_vals
 :mu_vals_t
 :psi_vals
 :psi_vals_t
 :Qij_vals
 :Qij_vals_t
 :Cijk_weights
 :Cijk_probs
 :Cijk_rates
 :Cijk_vals
 :Cijk_rates_t
 :row_weightvals





bmo = Rrbind(bmo, bmo_add)

rownames = [];

bmo = DataFrames.DataFrame(rownames=rownames, type=type_vec, init=init_vec, min=min_vec, max=max_vec, est=est_vec, note=note_vec, desc=desc_vec)


for i in 1:num_trait_states
	rownames = 
end