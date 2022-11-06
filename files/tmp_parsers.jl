"""
Construct interpolators for Qmat, Carray

Inputs: a "" that includes distmats
"""
function construct_interpolators(inputs)
	
end # END function construct_interpolators


"""
Read a text file containing whitespace-delimited distance matrices

(with a blank line between them)

fn = "/GitHub/PhyBEARS.jl/files/distances_changing_v1_wEND.txt";
"""
function parse_distances_fn(fn)
	lines = readlines(fn);
	
	# Figure out how many blocks there are
	numblocks = 0
	list_of_numlines_per_block = Any[]
	numlines = 0
	count_one = 0
	for i in 1:length(lines)
		# Try to parse, if not, save the string
		try
			words = parse.(Float64, split(lines[i]))
		catch
			words = lines[i]
		end
		
		if length(words) == 0
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
		 break
		end

		if length(words) > 0
			# Not a blank line, reset to 0
			count_one = 0
			numlines = numlines + 1
		end
	end
	
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
	distmats = [Matrix{Float64}(undef, numcols, numrows) for _ = 1:num_matrices]
	
	blocknum = 1
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

include("/Users/nickm/GitHub/PhyBEARS.jl/files/tmp_parsers.jl")
fn = "/Users/nickm/GitHub/PhyBEARS.jl/files/area_of_areas_changing_v1_wEND.txt";
area_vectors = parse_areas_fn4(fn)
"""
function parse_areas_fn4(fn)
	lines = readlines(fn);
	# Replace tabs with space
	#lines = readlines(IOBuffer(replace(read(fn), UInt8('\t') => UInt8(' '))));
	
	# Figure out how many blocks there are
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
	println(list_of_numcols)
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
			print(words)
			area_vectors[i] .= words
		end
	end
	
	return(area_vectors)	
end # END function parse_areas_fn(fn)




