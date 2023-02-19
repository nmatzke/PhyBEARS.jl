#######################################################
# Parsers
#######################################################
module Parsers
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/


print("PhyBEARS: loading Parsers.jl dependencies...")
#using BenchmarkTools # for @time
using DelimitedFiles	# for writedlm()
using CSV							# for CSV.write()
using Interpolations	# for Linear, Gridded, interpolate
using InvertedIndices # for Not
using LSODA           # for lsoda()
using Sundials        # for CVODE_BDF(linear_solver=:GMRES)
using DifferentialEquations
using Base.Threads  # <-- this is good for @spawn, Distributed.@spawn is BAD, it produces Futures that have to be scheduled etc]
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloBits
#using Plots						# for plot
using DataFrames          # for DataFrame()
using PhyloBits.TrUtils # for flat2() (similar to unlist) and offdiag()
using PhyBEARS.StateSpace
using PhyloBits.TreeTable
using PhyBEARS.SSEs
using PhyBEARS.TimeDep # for update_min_vdist_at_time_t!

print("...done.\n")


# (1) List all function names here:
export extract_first_integer_from_string, states_list_to_R_cmd, getranges_from_LagrangePHYLIP, tipranges_to_tiplikes, modify_tiplikes_sampling_fossils_v7!, check_tr_geog_tip_labels, parse_distances_fn, parse_areas_fn, parse_numbers_list_fn, parse_times_fn, files_to_interpolators, model_to_text_v12

#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/PhyBEARS.jl/notes/")
include("Parsers.jl")
"""
#######################################################


"""
states_list = [[1], [2], [1,2]]
area_names = ["A", "B"]
txt_states_list = states_list_to_txt(states_list, area_names)
"""
# See function in library PhyBEARS.StateSpace
# function states_list_to_txt(states_list, area_names)
# 	txt_states_list = collect(repeat([""], length(states_list)))
# 	for i in 1:length(states_list)
# 		txt_states_list[i] = paste0(area_names[states_list[i]])
# 	end
# 	return(txt_states_list)
# end # END function states_list_to_txt(states_list, area_names)





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


include("/GitHub/PhyBEARS.jl/notes/Parsers.jl")
import .Parsers
lgdata_fn = "/GitHub/PhyBEARS.jl/data/Psychotria/Psychotria_geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)
"""
function extract_first_integer_from_string(str)
	tmp = [filter(isdigit, collect(s)) for s in str]
	tmp2 = tmp[length.(tmp) .> 0]
	tmp3 = parse.(Int, string(tmp2[1][]))
	return(tmp3)
end

"""
states_list = Vector{Any}[[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
states_list_to_R_cmd(states_list; outfn="states_list.R")
"""

function states_list_to_R_cmd(states_list; outfn="")
	tmpstr = string(states_list)
	tmpstr2 = replace(tmpstr, "Vector{Any}["=>"states_list = list(") 
	tmpstr3 = replace(tmpstr2, "]]"=>"])") 
	tmpstr4 = replace(tmpstr3, "["=>"c(") 
	tmpstr5 = replace(tmpstr4, "]"=>")") 
	
	if outfn != ""
		open(outfn, "w") do io
			write(io, tmpstr5)
			write(io, "\n")
		end # END open,do
	end # END if
	
	return(tmpstr5)
end




"""
Function to read in Lagrange/BioGeoBEARS-type PHYLIP-formatted
 geography input files.

lgdata_fn = lagrange-style geography input file
block_allQs = give error if a row has all-question-marks


# Read in a DataFrame from script
geog_df = DataFrame(AbstractVector[["sp4", "sp5", "sp6", "sp7"], [1, 1, 0, 0], [0, 0, 1, 1]], DataFrames.Index(Dict(:A => 2, :tipnames => 1, :B => 3), [:tipnames, :A, :B]))

# Output to screen:
julian_dput(geog_df)


# Read from file
wd = "/GitHub/PhyBEARS.jl/test/apes_SSE/"
cd(wd)

lgdata_fn = "geog.data"
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);

"""

function getranges_from_LagrangePHYLIP(lgdata_fn; block_allQs=true)
	#lgdata_fn = "/GitHub/PhyBEARS.jl/data/Psychotria/Psychotria_geog.data"

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
function tipranges_to_tiplikes(inputs, geog_df; fossils_older_than=1.0e-5)
	# Error check
	taxa = inputs.trdf.taxa
	tipnames = sort(taxa[inputs.trdf.nodeType .== "tip"])
	check_tr_geog_tip_labels(tipnames, geog_df)
	
	
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
		tmpnums = parse.(Int, string.(flat2(tmprow))) # convert area 1/0s to Integer
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





"""
#######################################################
Modify tiplikes by:
#######################################################
psi, sampling rate for serially-sampled fossils
rho (sampling_f)
rho (tipsamp_f)


res.sampling_f		# rho sampling probabilities for each state
res.tipsamp_f    	# sampling probabilities for each tip (a priori) -- the problem with this
									# is that it modifies Ei at the tip, for all states. Different tip branches
									# thus can't share the same Es calculation. Better to use states with different
									# rho to control for this (e.g. just a 'poorly sampled' and 'well-sampled' 
									# e.g. forests vs. not). But, leaving it in for experimental purposes.
inputs.setup.fossil_TF # is each node a fossil? If so, D multiplied by psi
inputs.setup.direct_TF # is each node a direct ancestor? 
												# NOTE: fossil direct ancestors should be represented by hooks, NOT
												#       direct ancestors, which are for recording/hypothesizing ancestral
												#       states. (Hooks are ultrashort branches that are not treated
												#       as speciation events.)

inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0         # zero out
inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum] .= 0.0     # zero out
inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum][inputs.setup.observed_statenums[i]] = 1.0
inputs.res.sumLikes_at_node_at_branchTop
"""
# Time-constant version of psi
function modify_tiplikes_sampling_fossils_v7!(inputs, p_Ds_v5)
	# Error checks
	taxa = inputs.trdf.taxa
	tipnames = sort(taxa[inputs.trdf.nodeType .== "tip"])
	check_tr_geog_tip_labels(tipnames, geog_df)

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

	# Go through the geography file tips, match each one to the correct node of trdf,
	# then update the tip likelihoods at that node.

	for i in 1:nrow(geog_df)
		spname = geog_df[i,:tipnames]
		TF = spname .== inputs.trdf[!,:nodeName]
		nodeNum = trdf_nodenums[TF][1]
		node_age = trdf.node_age[nodeNum]
		# Check if it's a fossil
		if inputs.setup.fossil_TF[nodeNum] == true
			if inputs.trdf.hook[nodeNum] == false
				# An m-type fossil (not a hook)
				# Di = observed state/range for fossil * psi for each state * Ei(node_age) for each state
				inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .* p_Ds_v5.params.psi_vals .* p_Ds_v5.sol_Es_v5(node_age)
				# Ei for node is just normal Ei(t)
			else
				# A k-type fossil (a hook, ie a direct ancestor, but making it a hook tip for 
				#   ease of storing/representation)
				# Di = observed state/range for fossil * psi for each state * Ei(node_age) for each state
				# Unlike Beaulieu & O'Meara, we assume there is trait information for the fossil
				inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .* p_Ds_v5.params.psi_vals
				# Ei for node is just normal Ei(t)
			end # END if inputs.trdf.hook[nodeNum] == false
		else
			# For living tips:
			inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .= inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] .* inputs.res.sampling_f[nodeNum]
			# (In the E's calculation, which should have been done already, 
			#  the tip Es will be, instead of 0.0... : u0 .= 1 .- inputs.res.sampling_f[nodeNum]
			
		end # END if inputs.setup.fossil_TF[nodeNum] == true

	inputs.res.sumLikes_at_node_at_branchTop[nodeNum] = sum(inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum])
	inputs.res.normlikes_at_each_nodeIndex_branchTop[nodeNum] .= inputs.res.likes_at_each_nodeIndex_branchTop[nodeNum] ./ inputs.res.sumLikes_at_node_at_branchTop[nodeNum]
		
	end
	
end




# Check that the tree and the geog file have matching tip labels

function check_tr_geog_tip_labels(tr::PhyloBits.PNtypes.HybridNetwork, geog_df::DataFrame)
	tipnames = sort(tr.names)
	geognames = sort(geog_df.tipnames)
	
	if (length(tipnames) != length(geognames))
		txt1 = "ERROR in check_tr_geog_tip_labels(tr, geog_df): the tree and the geography file must have the same number of species/OTUs."
		txt2 = paste0(["Instead, the tree has ", length(tipnames), " tips, and the geography file has ", length(geog_df.tipnames), ". Fix these input files and re-run."])
		println(txt1)
		println(txt2)
		error(paste0([txt1, " ", txt2]))
	end
	
	# Otherwise, check that they match
	TF = all(tipnames .== geognames)
	if TF == false
		txt = "ERROR in check_tr_geog_tip_labels(tr, geog_df): the tree and the geography species/OTU names do not match. Printing them to screen."
		println(txt)
		display(Rcbind(tipnames, geognames))
		error(txt)
	end
	return()
end


# Alternative version
function check_tr_geog_tip_labels(tipnames::Vector{String}, geog_df::DataFrame)
	#tipnames = sort(tr.names)
	tipnames = sort(tipnames)
	geognames = sort(geog_df.tipnames)
	
	if (length(tipnames) != length(geognames))
		txt1 = "ERROR in check_tr_geog_tip_labels(tr, geog_df): the tree and the geography file must have the same number of species/OTUs."
		txt2 = paste0(["Instead, the tree has ", length(tipnames), " tips, and the geography file has ", length(geog_df.tipnames), ". Fix these input files and re-run."])
		println(txt1)
		println(txt2)
		error(paste0([txt1, " ", txt2]))
	end
	
	# Otherwise, check that they match
	TF = all(tipnames .== geognames)
	if TF == false
		txt = "ERROR in check_tr_geog_tip_labels(tr, geog_df): the tree and the geography species/OTU names do not match. Printing them to screen."
		println(txt)
		display(Rcbind(tipnames, geognames))
		error(txt)
	end
	return()
end




"""
Read a text file containing whitespace-delimited distance matrices

(with a blank line between them)

fn = "/GitHub/PhyBEARS.jl/files/distances_changing_v1_wEND.txt";
distmats = parse_distances_fn(fn)

fn = "/GitHub/PhyBEARS.jl/files/distances_changing_v1.txt";
distmats = parse_distances_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/area_of_areas_changing_v1_wEND.txt";
area_vectors = parse_areas_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/area_of_areas_changing_v1.txt";
area_vectors = parse_areas_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1.txt";
times = parse_times_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1.txt";
times = parse_numbers_list_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1_no0time.txt";
times = parse_times_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1_no0time.txt";
times = parse_numbers_list_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_wEND_v1.txt";
times = parse_times_fn(fn)
"""
function parse_distances_fn(fn)
	lines = readlines(fn);
	
	# Figure out how many blocks there are
	words = Any[]
	numblocks = 0
	list_of_numlines_per_block = Any[]
	numlines = 0
	count_one = 0
	#lastword = ""
	for i in 1:length(lines)
		# Try to parse, if not, save the string
		try
			words = parse.(Float64, split(lines[i]))
		catch
			words = lines[i]
		end
		
		if (length(words) == 0)
			numblocks = numblocks + 1
			count_one = count_one + 1
			push!(list_of_numlines_per_block, numlines)
			numlines = 0
			if count_one > 1
				txt = "ERROR in parse_distances_fn(): file has 2 or more consecutive blank lines. Fix and re-run."
				error(txt)
			end
		end

		# Catch "END"
		if (uppercase.(lines[i]) == "END")
			lastword = "END"
		 	break
		end

		# If the last line has numbers
		if (length(words) > 0) && (i == length(lines))
			numblocks = numblocks + 1
			numlines = numlines + 1
			push!(list_of_numlines_per_block, numlines)
		end


		if length(words) > 0
			# Not a blank line, reset to 0
			count_one = 0
			numlines = numlines + 1
		end
	end
	
	# If it didn't end with an END, add 1 to numblocks
	#if lastword == ""
	#	numblocks = numblocks + 1
	#end
	
	# OK, now list_of_numlines_per_block has the number
	# of lines per block
	# Check if they match
	
	if all(unique(list_of_numlines_per_block) .== list_of_numlines_per_block) == false
		txt = "ERROR in parse_distances_fn(): all of the distance matrices must have the same number of lines. Instead, the blocks have the number of lines printed below. Fix and re-run."
		print("\n\n")
		print(txt)
		print("\nlist_of_numlines_per_block = \n")
		print(list_of_numlines_per_block)
		error(txt)
	end
	
	# If all error checks passed, fill in the matrices
	numrows = unique(list_of_numlines_per_block)[1]
	numcols = unique(list_of_numlines_per_block)[1]
	num_matrices = length(list_of_numlines_per_block)
	distmats = [Matrix{Float64}(undef, numrows, numcols) for _ = 1:num_matrices]
	
	blocknum = 1
	j = 1
	for i in 1:length(lines)
		# Try to parse, if not, save the string
		try
			words = parse.(Float64, split(lines[i]))
		catch
			words = lines[i]
		end

		# Catch "END"
		if (uppercase.(lines[i]) == "END")
		 break
		end

		if length(words) == 0
			j = 1
			blocknum = blocknum + 1
		end

		if (length(words) == 0) && (i == length(lines))
			break
		end
	
	
	
		if length(words) > 0
		# Not a blank line, fill in data
			if length(words) == numcols
				distmats[blocknum][j,:] .= words
				j = j+1
			else
				txt = paste0(["ERROR in parse_distances_fn(): On line ", string(i), " of '", fn, "', there were ", string(length(words)), " numbers. There need to be ", string(numcols), "."])
				error(txt)
			end
		end
	end
	
	return(distmats)	
end # END function parse_distances_fn(fn)


"""
Read a text file containing whitespace-delimited vectors

(withOUT a blank line between them)

fn = "/GitHub/PhyBEARS.jl/files/distances_changing_v1_wEND.txt";
distmats = parse_distances_fn(fn)

fn = "/GitHub/PhyBEARS.jl/files/distances_changing_v1.txt";
distmats = parse_distances_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/area_of_areas_changing_v1_wEND.txt";
area_vectors = parse_areas_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/area_of_areas_changing_v1.txt";
area_vectors = parse_areas_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1.txt";
times = parse_times_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1.txt";
times = parse_numbers_list_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1_no0time.txt";
times = parse_times_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1_no0time.txt";
times = parse_numbers_list_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_wEND_v1.txt";
times = parse_times_fn(fn)
"""
function parse_areas_fn(fn)
	lines = readlines(fn);
	# Replace tabs with space
	#lines = readlines(IOBuffer(replace(read(fn), UInt8('\t') => UInt8(' '))));
	
	# Figure out how many blocks there are
	words = Any[]
	list_of_numcols = Any[]
	numlines = 0
	for i in 1:length(lines)
		# Try to parse, if not, save the string
		try
			words = parse.(Float64, split(lines[i]))
		catch
			words = lines[i]
		end
		
		if length(words) == 0
			break
		end

		# Catch "END"
		if (uppercase.(lines[i]) == "END")
		 break
		end

		if length(words) > 0
			push!(list_of_numcols, length(words))
			numlines = numlines + 1
		end
	end
	
	# OK, now list_of_numlines_per_block has the number
	# of lines per block
	# Check if they match
	
	if all(unique(list_of_numcols) .== list_of_numcols) == false
		txt = "ERROR in parse_areas_fn(): all of the numeric rows must have the same number of entries. Instead, the rows have the number of entries printed below. Fix and re-run."
		print("\n\n")
		print(txt)
		print("\nlist_of_numcols = \n")
		print(list_of_numcols)
		error(txt)
	end
	
	# If all error checks passed, fill in the matrices
	numrows = length(list_of_numcols)
	#println(list_of_numcols)
	numcols = unique(list_of_numcols)[1]
	area_vectors = [Vector{Float64}(undef, numcols) for _ = 1:numrows]
	
	for i in 1:length(lines)
		# Try to parse, if not, save the string
		try
			words = parse.(Float64, split(lines[i]))
		catch
			words = lines[i]
		end

		if length(words) == 0
			break
		end

		# Catch "END"
		if (uppercase.(lines[i]) == "END")
			break
		end
		
		#print(words)
		
		if length(words) > 0
			#print(words)
			area_vectors[i] .= words
		end
	end
	
	return(area_vectors)	
end # END function parse_areas_fn(fn)


"""
fn = "/GitHub/PhyBEARS.jl/files/distances_changing_v1_wEND.txt";
distmats = parse_distances_fn(fn)

fn = "/GitHub/PhyBEARS.jl/files/distances_changing_v1.txt";
distmats = parse_distances_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/area_of_areas_changing_v1_wEND.txt";
area_vectors = parse_areas_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/area_of_areas_changing_v1.txt";
area_vectors = parse_areas_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1.txt";
times = parse_times_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1.txt";
times = parse_numbers_list_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1_no0time.txt";
times = parse_times_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1_no0time.txt";
times = parse_numbers_list_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_wEND_v1.txt";
times = parse_times_fn(fn)
"""
function parse_numbers_list_fn(fn)
	lines = readlines(fn);
	
	# Remove any lines with "END"
	TF1 = uppercase.(lines) .== "END"
	TF2 = lines .== ""
	TF = (TF1 .+ TF2) .> 0
	nums_to_delete = (1:length(TF))[TF]
	deleteat!(lines, nums_to_delete)
	
	nums = parse.(Float64, lines)
	return(nums)	
end # END function parse_times_fn(fn)

"""
fn = "/GitHub/PhyBEARS.jl/files/distances_changing_v1_wEND.txt";
distmats = parse_distances_fn(fn)

fn = "/GitHub/PhyBEARS.jl/files/distances_changing_v1.txt";
distmats = parse_distances_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/area_of_areas_changing_v1_wEND.txt";
area_vectors = parse_areas_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/area_of_areas_changing_v1.txt";
area_vectors = parse_areas_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1.txt";
times = parse_times_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1.txt";
times = parse_numbers_list_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1_no0time.txt";
times = parse_times_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_v1_no0time.txt";
times = parse_numbers_list_fn(fn)

fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/times_wEND_v1.txt";
times = parse_times_fn(fn)
"""
function parse_times_fn(fn)
	times = parse_numbers_list_fn(fn)
	if times[1] != 0.0
		pushfirst!(times, 0.0)
	end
	return(times)	
end # END function parse_times_fn(fn)



#######################################################
# Read distances etc. files to interpolators
#######################################################
function files_to_interpolators(files, numareas, states_list, v_rows, Carray_jvals, Carray_kvals, trdf; oldest_possible_age=1000.0, fraction_of_minimum_nonzero_distance=0.5)
	
	if files.distances_fn != ""
		times = parse_times_fn(files.times_fn)
	else
		times = [0.0, oldest_possible_age]
	end

	times_for_SSE_interpolators = sort(unique(reduce(vcat, [times, sort(unique(trdf.node_age))])));

	
	#######################################################
	# Distances interpolator for a series of distances
	#######################################################
	if files.distances_fn != ""
		distmats = parse_distances_fn(files.distances_fn)
		# files.distances2_fn
		# files.distances3_fn
		# files.envdistances_fn
		# files.manual_dispersal_multipliers_fn
		TF = length(times) == length(distmats)
		if TF == false
			txt = "STOP ERROR in files_to_interpolators(). distances_fn must have the same number of entries as the times_fn file."
			error(txt)
		end
		
		# Check if the distances matrix is the correct size
		TF = dim(distmats[1])[1] == numareas
		if TF == false
			txt = paste0(["STOP ERROR in files_to_interpolators(). numareas=", numareas, ", but the distance matrices in files.distances_fn='", files.distances_fn, "' have dim(distmats[1])[1]=", dim(distmats[1])[1], ". These must match. Edit your inputs and re-try."])
			error(txt)
		end
		
		#######################################################
		# NOTES ON HAVING DISTANCES MODIFY RATES
		#
		# 1. Units don't matter. All that matters is the 
		#    proportions that one assigns to different distances.
		#
		# 2. Therefore, we could just divide by the minimum or 
		#    maximum distance. However, this doesn't work
		#		 when the minimum distance is 0.0
		#
		# 3. A FULL solution would have a separate parameter for
		#    dispersal-between-connected or dispersal-between-not
		#    connected, i.e. a multiplier. This could be done with 
		#    the dispersal multipliers file if desired.
		#
		# 4. However, just for the distances, let's do it like
		#    this. The goal is that the smallest distance gets
		#    a distance-based dispersal multiplier of 1.0 -- thus
		#    the base-rate of dispersal does NOT change, regardless
		#    of the value of "x". 
		#
		#    4a. *If* the smallest non-diagonal distance is positive,
		#        i.e. it's not 0.0, then divide all distances by this
		#        distance. Now all distances are just ratios of
		#        this distance.
		#
		#    4b. *If* the smallest non-diagonal distance is 0.0,
		# 			 divide all distances by:
		#        min_rel_dist = (minimum_nonzero_distance * fraction_of_minimum_nonzero_distance)
		#
		#        The default: fraction_of_minimum_nonzero_distance = 0.5.
		#
		#        E.g., if x=-1, then by default, dispersal over the minimum_nonzero_distance
		#        would have 50% of the rate of dispersal over the 0.0 physical distance
		#        
		#        This fraction_of_minimum_nonzero_distance is arbitrary.  If users
		#        decide to use parameter u to penalize dispersal between non-touching 
		#        areas, fraction_of_minimum_nonzero_distance should be set to 1.0.
		#######################################################
		
		# Let's divide the distances by the maximum distance
		# dists_vector = collect(Iterators.flatten(vec.(distmats)))
		# dists_vector2 = dists_vector[dists_vector .> 0.0]
		# maxval = maximum(dists_vector2)

		# Find the minimum off-diagonal in the distance matrices
		dists_vector = collect(Iterators.flatten(vec.(offdiag.(distmats))))
		minval = minimum(dists_vector)
		minval_nonzero = minimum(dists_vector[dists_vector .> 0.0])
		
		# Matrix of true/false for offdiagonal
		offdiag_TF = make_offdiag_TF(dim(distmats[1])[1])
		
		if minval == 0.0
			# Replace any off-diagonal zeros with
			# inval_nonzero * fraction_of_minimum_nonzero_distance
			for i in 1:length(distmats)
				offdiag_and_zeroTF = (offdiag_TF + (distmats[i] .== 0.0)) .== 2
				
				distmats[i][offdiag_and_zeroTF] .= minval_nonzero * fraction_of_minimum_nonzero_distance
			end
			# This is the new minval
			minval = minval_nonzero * fraction_of_minimum_nonzero_distance
		end
		
		for i in 1:length(distmats)
			distmats[i] .= distmats[i] ./ minval
		end
	
		distances_interpolator = interpolate((times,), distmats, Gridded(Linear()));
	else
		distmats = [Array{Float64}(undef, numareas, numareas) for _ = 1:length(times)]
		distances_interpolator = interpolate((times,), distmats, Gridded(Linear()));
	end
	
	#######################################################
	# Area of areas
	#######################################################
	if files.area_of_areas_fn != ""
		area_of_areas = parse_areas_fn(files.area_of_areas_fn)
		TF = length(times) == length(area_of_areas)
		if TF == false
			txt = "STOP ERROR in files_to_interpolators(). area_of_areas_fn must have the same number of entries as the times_fn file."
			error(txt)
		end

		# Check if the distances matrix is the correct size
		TF = length(area_of_areas[1]) == numareas
		if TF == false
			txt = paste0(["STOP ERROR in files_to_interpolators(). numareas=", numareas, ", but the area_of_areas matrices in files.area_of_areas_fn='", files.area_of_areas_fn, "' have length(area_of_areas[1])=", length(area_of_areas[1]), ". These must match. Edit your inputs and re-try."])
			error(txt)
		end

	
		# Let's divide the distances by the maximum area
		areas_vector = collect(Iterators.flatten(vec.(area_of_areas)))
		maxval = maximum(areas_vector)

		for i in 1:length(area_of_areas)
			area_of_areas[i] .= area_of_areas[i] ./ maxval
		end
	
		area_of_areas_interpolator = interpolate((times,), area_of_areas, Gridded(Linear()))
	else
		area_of_areas = [Vector{Float64}(undef, numareas) for _ = 1:length(times)]
		for i in 1:length(area_of_areas)
			area_of_areas[i] .= 1.0
		end
		area_of_areas_interpolator = interpolate((times,), area_of_areas, Gridded(Linear()))
	end

	
	#######################################################
	# Vicariance - minimum distances interpolator
	#######################################################
	#######################################################
	# Construct vicariance minimum-distances interpolator
	# (for the v_rows of the matrix)
	#######################################################
	# A Vector of # times Vectors, each of length v_rows
	# THIS AVOIDS THE LINKED VECTORS ISSUE
	changing_mindists = [Vector{Float64}(undef, length(v_rows)) for _ = 1:length(times)]

	# Populate the minimum distances
	for i in 1:length(times)
		distmat = distances_interpolator(times[i])
		update_min_vdist_at_time_t!(changing_mindists[i], v_rows, distmat, states_list, Carray_jvals, Carray_kvals; mindist=1.0e9)
	end

	vicariance_mindists_interpolator = interpolate((times,), changing_mindists, Gridded(Linear()));

	interpolators = (times_for_dists_interpolator=times, times_for_SSE_interpolators=times_for_SSE_interpolators, distances_interpolator=distances_interpolator, area_of_areas_interpolator=area_of_areas_interpolator, vicariance_mindists_interpolator=vicariance_mindists_interpolator)
	
	return(interpolators)
end # END function files_to_interpolators()

#######################################################
# Output a time-varying model to text files
#######################################################

"""
outfns = ["setup_df.txt",
"timepoints.txt", 
"mu_vals_by_t.txt", 
"Qvals_by_t.txt",
"Crates_by_t.txt",
"Qarray.txt",
"Carray.txt",
"area_names.txt",
"states_list.R"
]
"""

function model_to_text_v12(p_Ds_v12, timepoints; prefix="")
	area_names = p_Ds_v12.setup.area_names
	# Initialize list of output filenames
	num_outfns = 9	# 8 filenames for now
	outfns = Vector{String}(undef, num_outfns) 
	
	TF = (prefix == "") || (endswith(prefix,"_"))
	if TF == false
		prefix2 = paste0([prefix, "_"])
	else
		prefix2 = paste0([prefix, ""])
	end
	
	# Run settings
	numareas = length(p_Ds_v12.setup.areas_list)
	numstates = length(p_Ds_v12.setup.states_list)
	max_range_size = p_Ds_v12.setup.max_range_size
	include_null_range = p_Ds_v12.setup.include_null_range
	setup_df = DataFrame(numareas=numareas, numstates=numstates, max_range_size=max_range_size, include_null_range=include_null_range)
	outfns[1] = paste0([prefix2, "setup_df.txt"]);
	CSV.write(outfns[1], setup_df; delim="\t")

	t = 0.0
	num_mu_vals = length(p_Ds_v12.interpolators.mu_vals_interpolator(t));
	numQvals = length(p_Ds_v12.interpolators.Q_vals_interpolator(t));
	numCrates = length(p_Ds_v12.interpolators.C_rates_interpolator(t));
	#Qvals_by_t = [Vector{Float64}(undef, numQvals) for _ = 1:length(timepoints)];
	#Crates_by_t = [Vector{Float64}(undef, numCrates) for _ = 1:length(timepoints)];
	mu_vals_by_t = Matrix{Float64}(undef, num_mu_vals, length(timepoints));
	Qvals_by_t = Matrix{Float64}(undef, numQvals, length(timepoints));
	Crates_by_t = Matrix{Float64}(undef, numCrates, length(timepoints));

	for i in 1:length(timepoints)
		mu_vals_by_t[:,i] .= p_Ds_v12.interpolators.mu_vals_interpolator(timepoints[i])
		Qvals_by_t[:,i] .= p_Ds_v12.interpolators.Q_vals_interpolator(timepoints[i])
		Crates_by_t[:,i] .= p_Ds_v12.interpolators.C_rates_interpolator(timepoints[i])
	end

	
	outfns[2] = paste0([prefix2, "timepoints.txt"]);
	open(outfns[2], "w") do io
		writedlm(io, timepoints, '\n')
	end;
	outfns[3] = paste0([prefix2, "mu_vals_by_t.txt"]);
	open(outfns[3], "w") do io
		writedlm(io, mu_vals_by_t, '\t')
	end;
	outfns[4] = paste0([prefix2, "Qvals_by_t.txt"]);
	open(outfns[4], "w") do io
		writedlm(io, Qvals_by_t, '\t')
	end;
	outfns[5] = paste0([prefix2, "Crates_by_t.txt"]);
	open(outfns[5], "w") do io
		writedlm(io, Crates_by_t, '\t')
	end;
	# Write the time-changing mus to a Matrix

	# Write the rest of the Qarray (i,j)
	outfns[6] = paste0([prefix2, "Qarray.txt"]);
	CSV.write(outfns[6], prtQp(p_Ds_v12); delim="\t")
	# Write the rest of the Carray (i,j,k)
	outfns[7] = paste0([prefix2, "Carray.txt"]);
	CSV.write(outfns[7], prtCp(p_Ds_v12); delim="\t")
	#moref(outfn)
	
	outfns[8] = paste0([prefix2, "area_names.txt"]);
	open(outfns[8], "w") do io
		for i=1:length(area_names)
			write(io, area_names[i])
			write(io, "\n")
		end
	end # END open,do
	
	outfns[9] = paste0([prefix2, "states_list.R"]);
	states_list_to_R_cmd(p_Ds_v12.setup.states_list; outfn=outfns[9])
	
	return(outfns)
end # END function model_to_text_v12(p_Ds_v12, timepoints; prefix="")

end # ENDING Parsers
