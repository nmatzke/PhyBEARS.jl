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

lgdata_fn = "geog.txt"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
rowSums_df(geog_df)
maximum(rowSums(geog_df[:,2]))
minimum(rowSums(geog_df))

