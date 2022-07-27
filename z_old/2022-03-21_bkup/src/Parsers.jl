#######################################################
# Parsers
#######################################################
module Parsers
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/


using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA           # for lsoda()
using Sundials        # for CVODE_BDF(linear_solver=:GMRES)
using DifferentialEquations
using Distributed
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloNetworks
#using Plots						# for plot
using DataFrames          # for DataFrame()
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
using BioGeoJulia.StateSpace
using BioGeoJulia.TreeTable
using BioGeoJulia.SSEs

print("\nBioGeoJulia: loading Parsers.jl")


# (1) List all function names here:
export getranges_from_LagrangePHYLIP, tipranges_to_tiplikes

#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("Parsers.jl")
"""
#######################################################



"""
# Goal: Read in various text files,
# starting with a C++Lagrange / BioGeoBEARS-format
# geography file, e.g.:
#
# =================
# 19	4 (K O M H)
# P_mariniana_Kokee2	1000
# P_mariniana_Oahu	0100
# P_mariniana_MauiNui	0010
# ...
# =================
#
# This is basically a "PHYLIP"-formatted file.


include("/GitHub/BioGeoJulia.jl/notes/Parsers.jl")
import .Parsers
lgdata_fn = "/GitHub/BioGeoJulia.jl/Rsrc/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
"""


"""
Function to read in Lagrange/BioGeoBEARS-type PHYLIP-formatted
 geography input files.

lgdata_fn = lagrange-style geography input file
block_allQs = give error if a row has all-question-marks
"""

function getranges_from_LagrangePHYLIP(lgdata_fn; block_allQs=true)
	#lgdata_fn = "/GitHub/BioGeoJulia.jl/Rsrc/Psychotria_geog.data"

	# Obtain a file handle
	fhandle = open(lgdata_fn)

	# Read the lines to an array
	lines = readlines(fhandle)
	close(fhandle)
	
	# Count the species read in
	spnum = 0
	
	# We have to assume that the first line is a header
	for i in 1:length(lines)
		# If the line is blank, skip it
		line = lines[i]
		if (length(line) == 0)
			continue
		end
		#println(i, ": ", line)

		# Parse the header line
		if (i == 1)
			# Split on "("
			parts = split(line, "(")
			# remove whitespace, convert to Integer using parse.()
			parts1 = strip(parts[1])
			global (numtaxa, numareas) = parse.(Int, split(strip(parts1)))

			# Part 2: remove trailing ")", parse the rest to areas list
			parts2 = strip(strip(strip(parts[2]), [')']))
			# Split on whitespace to get the list of areas
			global areas_list = split(parts2)

			# Check that the given number matches the actual number of areas
			if (length(areas_list) != numareas)
				txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). In your input file,\n'", lgdata_fn, "',\nthe given number of areas in line 1 (numareas=", numareas, ") does not equal the length of the areas_list (", parts2, ").\nPlease correct your input geography file and retry. Have a nice day."])
				error(txt)
			end # end error check
	 
			# Set up the output matrix
			global sp_names = collect(repeat([""], numtaxa))
			tmp_areas = collect(repeat(["0"], numareas))
			global sp_areas = collect(repeat([tmp_areas], numtaxa))

			continue # Go to next round of for-loop without finishing this one
		end # END if (i == 1)

		# All of the other lines, fill in sp_names and sp_areas
		words = strip.(split(line))
		spnum = spnum + 1
	
		# Error check
		if (spnum > numtaxa)
			txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). While reading in species #", spnum, ", you exceeded the limit declared in the first line of your input geography file, where numtaxa=", numtaxa, " species. Please correct your input geography file and retry. Have a nice day."])
			error(txt)
		end

		
		sp_names[spnum] = words[1]
		# Split areanums to "1", "0", etc., parse to Integers
		#areanums_for_this_species = parse.(Int, split(words[2], ""))
		areanums_for_this_species = split(words[2], "")
		#print(areanums_for_this_species)
		# [:] avoids creating a linked reference
		sp_areas[spnum] = areanums_for_this_species[:]

	end

	# Error check
	if (spnum != numtaxa)
		txt = paste0(["STOP ERROR in getranges_from_LagrangePHYLIP(). While reading in species, only ", spnum, " taxa were read in, however the first line of your input geography file declared there would be numtaxa=", numtaxa, " species. Please correct your input geography file and retry. Have a nice day."])
		error(txt)
	end

	# DataFrame of the area presence/absences
	# Flatten to vector, reshape to array, transpose to get
	# back original inputs in matrix form
	sp_areas2 = permutedims(reshape(flat2(sp_areas), numareas, numtaxa))
	sp_areas2_Matrix = convert(Matrix, sp_areas2)
	sp_areas3_Matrix = Rcbind(sp_names, sp_areas2_Matrix)
	#geog_df = convert(DataFrame, sp_areas3_Matrix)
	geog_df = DataFrame(sp_areas3_Matrix, :auto)
	new_names = Rrbind(["tipnames"], areas_list)
	geog_df = DataFrames.rename(geog_df, new_names)

	geog_df[!,:tipnames]

	return geog_df
end # END getranges_from_LagrangePHYLIP


#######################################################
# Put the tip ranges into the likelihoods
#######################################################
function tipranges_to_tiplikes(inputs, geog_df)
	dfnames = names(geog_df)
	area_column_nums = 2:length(dfnames)
	areas_txt_list = dfnames[area_column_nums]
	numareas = length(inputs.setup.areas_list)
	
	# Check if the number of areas in the geography file matches the number in the geog_df
	if (inputs.setup.numareas != numareas)
		txt = paste0(["STOP ERROR in tipranges_to_tiplikes(): inputs.setup.numareas=", numareas, ", but the number of areas in geog_df is ", numareas, ". Please fix and re-run."])
		error(txt)
	end
	
# 	statenums = collect(1:length(states_list))
# 	observed_statenums = collect(repeat([0], nrow(geog_df)))
	trdf_nodenums = collect(1:nrow(inputs.trdf))

	# Convert the observed geography into a list of states comparable to states_list
	geog_observed_list = []  # empty array

	for i in 1:nrow(geog_df)
		tmprow = geog_df[i,area_column_nums]
		# The function "parse" assumes a *string* input, not an Int or Float
		tmpnums = parse.(Int, flat2(tmprow)) # convert area 1/0s to Integer
		range_as_areanums = inputs.setup.areas_list[tmpnums .== 1]
		# Compare observed range_as_areanums to full states_list
		TF = [range_as_areanums] .== inputs.setup.states_list
		if (sum(TF) != 1)
			txt = paste0(["STOP ERROR: An observed range in your geography file, from tipname '", geog_df[i,:tipnames], "', is  found either 0 or >1 times in the list of states in inputs.setup.states_list. Printing range_as_areanums, then inputs.setup.states_list"])
			print(txt)
			print("\nrange_as_areanums (this is the observed range that was not found):\n")
			print(range_as_areanums)
			print("\nstates_list:\n")
			print(inputs.setup.states_list)
			error(txt)
		end
	
		# Yay, a single match to the states_list was found!
		# Convert to a state index number
		inputs.setup.observed_statenums[i] = inputs.setup.statenums[TF][1]
	end

	# Go through the geography file tips, match each one to the correct node of trdf,
	# then update the tip likelihoods at that node.

	for i in 1:nrow(geog_df)
		spname = geog_df[i,:tipnames]
		TF = spname .== inputs.trdf[!,:nodeName]
		nodeNum = trdf_nodenums[TF][1]
	
		# Input likelihoods of 1 for the observed state, 0 otherwise
		inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0         # zero out
		inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0     # zero out
		inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
		inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
		inputs.res.sumLikes_at_node_at_branchTop[nodeNum] = 1.0
	end

	return inputs
end # END function tipranges_to_tiplikes ()


end # ENDING Parsers