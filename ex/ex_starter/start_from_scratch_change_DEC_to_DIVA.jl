# Example: go from a plain Julia startup to tables for
# Q (anagenetic) and C (cladogenetic) event rates 
julia --startup-file=no

# Installation
using Pkg

# Install these once (and any needed dependencies)
"""
Pkg.add(PackageSpec(url="https://github.com/nmatzke/PhyloBits.jl")) # for readTopology from PhyloNetworks, prt() for tree tables, etc.

Pkg.add(PackageSpec(url="https://github.com/nmatzke/PhyBEARS.jl"))  # formerly BioGeoJulia
"""

using Test
using DataFrames
using Dates									# for e.g. Dates.now(), DateTime
using Distributed						# for e.g. @spawnat
using Combinatorics					# for e.g. combinations()
using StatsBase
using NLopt									# seems to be the best gradient-free, box-constrained								
using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
using DataFrames  # for DataFrame
using DifferentialEquations

# Load functions
# List each PhyBEARS code file prefix here
using PhyBEARS.BGExample
using PhyloBits.PNreadwrite  # for readTopology etc.
using PhyloBits.TrUtils
using PhyBEARS.StateSpace
using PhyloBits.TreeTable		# for prt, nodetimes
using PhyBEARS.TreePass
using PhyBEARS.Parsers
using PhyBEARS.SSEs
using PhyBEARS.ModelLikes
using PhyBEARS.Optimizers

# Checking number of threads, cores and workers
using Distributed
Distributed.nprocs()
using Hwloc
Hwloc.num_physical_cores()
Hwloc.num_virtual_cores()
Distributed.nprocs()
numthreads = Base.Threads.nthreads()

# Load files
lgdata_fn = "/GitHub/PhyBEARS.jl/data/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Psychotria tree
tr = readTopology("((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);")

# DEC model on Hawaiian Psychotria
bmo = construct_BioGeoBEARS_model_object()
bmo.est[bmo.rownames .== "birthRate"] .= 0.32881638319078066
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.03505038
bmo.est[bmo.rownames .== "e"] .= 0.02832370
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
birthRate = 0.32881638319078066
numareas = Rncol(geog_df)-1

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
#inputs = ModelLikes.setup_DEC_SSE(numareas, tr; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, in_params=in_params)
root_age_mult=1.5; max_range_size=NaN; include_null_range=false; max_range_size=NaN
inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, solver_options, p_Es_v5, Es_tspan) = inputs;

prtQp(p_Es_v5)
prtCp(p_Es_v5)

Rnames(p_Es_v5)



