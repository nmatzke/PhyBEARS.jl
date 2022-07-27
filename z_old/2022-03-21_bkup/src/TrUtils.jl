module TrUtils
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/

using DataFrames
using Plots  # for savefig
#using RCall

print("\nBioGeoJulia: loading TrUtils.jl")

export getwd, Rgetwd, setwd, recursive_find, include_jls, source, dim, Rdim, seq, Rchoose, Rcbind, Rrbind, paste, paste0, type, class, Rclass, slashslash, addslash, df_to_Rdata, Reval, Rdput, Rnames, Rtypes, ont, saveopen, Rnrow, Rncol, Rsize, Rorder, headLR, flat2, single_element_array_to_scalar, headf, moref, scr2str, lagrange_to_tip

# Basic checks during Julia startup
#function hello_world_TrUtils()
#	display("TrUtils.hello_world_TrUtils() says hello on 2020-07-06_abc")
#end


# R-like utilities, and other short functions

# Handy aliases

# Reload BioGeoJulia
function re()
	# Remove and re-install
	Pkg.rm("BioGeoJulia")
	Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
	@eval using BioGeoJulia
	@eval using BioGeoJulia.MaxentInterp
	@eval using BioGeoJulia.BGExample
	@eval using BioGeoJulia.TrUtils
	@eval using BioGeoJulia.TreeTable
	@eval using BioGeoJulia.StateSpace
	@eval using BioGeoJulia.Parsers
	@eval using BioGeoJulia.SSEs
	@eval using BioGeoJulia.TreePass

	# Refresh Revise's look
	atreplinit() do repl
	try
	@eval using Revise
	@async Revise.wait_steal_repl_backend()
	catch
	end
	end
end
 

# getwd
function getwd()
	pwd()
end

function Rgetwd()
	pwd()
end

# setwd
function setwd(path=expanduser("~"))
	cd(path)
end


# Find all code *.jl files in a package

"""
package_path = "/GitHub/BioGeoJulia.jl"
recursive_find(package_path)
include_jls(package_path)
"""

function recursive_find(package_path)
	# Look in the source directory, "src"
	srcpath = joinpath(package_path, "src")
	srcfiles = readdir(srcpath)
	jl_files = []  # empty array
	for fn in srcfiles
		TF = endswith(fn, ".jl")
		if (TF == true)
			tmpfn = joinpath(srcpath, fn)
			push!(jl_files, tmpfn)
		end
	end
	
	# Look in the tests directory, "test"
# 	srcpath = joinpath(package_path, "test")
# 	srcfiles = readdir(srcpath)
# 	for fn in srcfiles
# 		TF = endswith(fn, ".jl")
# 		if (TF == true)
# 			tmpfn = joinpath(srcpath, fn)
# 			push!(jl_files, tmpfn)
# 		end
# 	end
	
	return jl_files
end

"""
package_path = "/GitHub/BioGeoJulia.jl"
recursive_find(package_path)
include_jls(package_path)
"""
function include_jls(package_path)
	srcfiles = recursive_find(package_path)
	for fn in srcfiles
		include(fn)
	end
end

# source
function source(str)
	include(str)
end


# dimensions
"""
A = ones(3,3)
dim(A)
Rdim(A)
dim(A)[1]
dim(A)[2]
"""
function dim(A)
	size(A)
end

# dimensions
function Rdim(A)
	size(A)
end


# seq
function seq(from, to, by=1)
	return(collect(from:by:to))
end


# choose
# n choose k
function Rchoose(n,k)
	return(binomial(n,k))
end


# cbind()
function Rcbind(A...)
	hcat(A...)
end

# rbind
# c()
# concatenate
function Rrbind(A...)
	vcat(A...)
end

# paste
function paste(array_of_strings; delim)
	newtxt = join(array_of_strings, delim)
	return(newtxt)
end

# paste0
function paste0(array_of_strings; delim="")
	newtxt = join(array_of_strings, delim)
	return(newtxt)
end

# type
function type(obj)
	typeof(obj)
end

# class
# Returns a plain-test version of the type/class
function class(obj)
	string(typeof(obj))
end

# Rclass
# Returns a plain-test version of the type/class
function Rclass(obj)
	string(typeof(obj))
end


# Convert any multiple slashes to single slashes
function slashslash(txt)
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	txt = replace(txt, "//" => "/")
	return(txt)
end

# Add a slash to the end of a string, if it is not there
# (Unless the txt string is "", then return "")
function addslash(txt)
	if (txt == "")
		return("")
	end
	if (endswith(txt, "/") == false)
		txt = join([txt, "/"])
	end
	return(txt)
end


# Save a julia DataFrame to an R data.frame
# as an Rdata file that can be easily loaded.
# Source: https://stackoverflow.com/questions/28084403/saving-julia-dataframe-to-read-in-r-using-hdf5/57903489#57903489

#######################################################
# NOTE: YOU MAY HAVE TO PASTE THE FUNCTION DEFINITION
# OF df_to_Rdata(), ALONG WITH THE "using" STATEMENTS,
# INTO YOUR MAIN ENVIRONMENT FOR IT TO WORK!!!
#######################################################

#using DataFrames
#using RCall
function df_to_Rdata(df; fn="dfjulia.RData", path=expanduser("~"))
	# Create the output fn
	pathfn = slashslash(join([addslash(path), fn], ""))
		
	# R environment in a session started from Julia
	g = globalEnv  # This requires "using RCall" previously, or you get:
								 # ERROR: UndefVarError: globalEnv not defined
	reval(rparse("dfls <- NULL"))

	# add columns one at a time converting Julia vectors to R-types via   RCall.sexp
	#  https://github.com/JuliaStats/RCall.jl/blob/master/src/sexp.jl
	for cnm in DataFrames._names(df)
		g[:colcnm] = sexp(convert(Array, df[!,cnm]))
		reval(rparse("dfls\$$cnm <- colcnm"))
	end
	reval(rparse("df <- data.frame(dfls)"))
	
	# Make and run the command to save the .Rdata file
	# (the .Rdata file will load to object df in R)
	txt = join(["save(file='", pathfn, "', df)"], "")
	reval(rparse(txt))
	return(pathfn)
end


# Input an object from a string representation from repr()
# dget, load, eval()
# eval
"""
tmpmatrix = [3 1; 3 2; 5 3; 5 4; 7 5; 7 6]
tmpstr = repr(tmpmatrix)
tmpstr2 = eval(Meta.parse(tmpstr))
tmpstr2
"""
function Reval(tmpstr)
	eval(Meta.parse(tmpstr))
end

# Output an object to a string representation with repr()
# dput, dump, str()
"""
tmpstr = "[[1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]"
tmpstr2 = tmpstr
states_list = Reval(tmpstr)
@test Rdput(states_list) == tmpstr2
"""
function Rdput(item)
	tmpstr = repr(item)
end


# fields / "names" of an object
# https://stackoverflow.com/questions/41687418/how-to-get-fields-of-a-julia-object
"""
obj = construct_Res()
Rnames(obj)
Rtypes(obj)
Rcbind(Rnames(obj), Rtypes(obj))
ont(obj)
"""
function Rnames(obj)
	flat2(fieldnames(typeof(obj)))
end


"""
obj = construct_Res()
Rnames(obj)
Rtypes(obj)
Rcbind(Rnames(obj), Rtypes(obj))
ont(obj)
"""
function Rtypes(obj)
	tmpnames = Rnames(obj)
	types = collect(repeat([DataType],length(tmpnames))) # empty array
	for i in 1:length(tmpnames)
		types[i] = typeof(getfield(obj, tmpnames[i]))		
	end
	return types
end




# ont = object names and types
"""
obj = construct_Res()
Rnames(obj)
Rtypes(obj)
Rcbind(Rnames(obj), Rtypes(obj))
ont(obj)
"""
function ont(obj)
	Rcbind(Rnames(obj), Rtypes(obj))
end



# Send the just-done plot to PDF, and open
function saveopen(fn)
	print("Saving ", """'""", fn, """'""")
	savefig(fn)
	cmdtxt = join(["`", "open ", fn, "`"], "")
	print("""Running """, cmdtxt, "\n")
	run(`open $fn`)
end


function Rnrow(obj)
	return size(obj)[1]
end

function Rncol(obj)
	return size(obj)[2]
end

function Rsize(obj)
	return size(obj)
end

function Rorder(obj)
	return sortperm(obj)
end

# Print the left and rightmost columns of a table
function headLR(df, num_startcols=4, num_endcols=4)
	ncols = Rncol(df)
	startcols = collect(1:num_startcols)
	endcols = collect((ncols-(num_endcols-1)):ncols)
	colnums = flat2(collect([startcols, endcols]))
	colnums
	print(df[:,colnums])
end


# Flattens an array of arrays into a vector
# Similar to R's unlist()
function flat2(arr)
    rst = Any[]
    grep(v) = for x in v
        if isa(x, Array) grep(x) else push!(rst, x) end
    end
    grep(arr)
    rst
end


# Convert a single-element array to scalar
# Julia often produces single-element arrays. 
# To convert to scalar, just take item [1]
# https://stackoverflow.com/questions/39079428/1-element-array-to-scalar-in-julia
function single_element_array_to_scalar(tmparray)
	if length(tmparray) != 1
		txt = ["STOP ERROR in single_element_array_to_scalar().\nThe input 'tmparray' has to be of length 1, but length(tmparray)=", string(length(tmparray)), ".\nPrinting input tmparray...\n"]
		errortxt = join(txt, "")
		println(errortxt)
		print(tmparray)
		error(errortxt)
	end
	
	# If check passed, go ahead.
	tmpscalar = tmparray[1]
	return tmpscalar
end



# Print the file to screen, with line numbers
function headf(fn; numlines=5)
	open(fn, "r") do f
		for (i,ln) in enumerate(eachline(f))
			if i > numlines
				break
			end
			println("$i $ln")
		end
	end
end

# Print the file to screen, with line numbers
function moref(fn)
	open(fn, "r") do f
	 for (i,ln) in enumerate(eachline(f))
		 println("$i $ln")
	 end
	end
end


function scr2str(obj)
	io = IOBuffer()
	show(io, "text/plain", obj)
	str = String(take!(io))
	return str
end

function lagrange_to_tip(inputs, geog_df)
	dfnames = names(geog_df)
	area_column_nums = 2:length(dfnames)
	areas_txt_list = dfnames[area_column_nums]
	numareas = length(inputs.setup.areas_list)
	areas_list = collect(1:numareas)
	maxareas = numareas
	include_null_range = true
	states_list = areas_list_to_states_list(areas_list, maxareas, include_null_range)
	statenums = collect(1:length(states_list))
	observed_statenums = collect(repeat([0], nrow(geog_df)))
	trdf_nodenums = collect(1:nrow(trdf))

	for i in 1:nrow(geog_df)
		tmprow = geog_df[i,area_column_nums]
		tmpnums = parse.(Int, flat2(tmprow))
		range_as_areanums = inputs.setup.areas_list[tmpnums .== 1]
		# Compare observed range_as_areanums to full states_list
		TF = [range_as_areanums] .== inputs.setup.states_list
		if (sum(TF) != 1)
			txt = paste0(["STOP ERROR: An observed range in your geography file, from tipname '", geog_df[i,:tipnames], "', is not found in the list of states in inputs.setup.states_list. Printing range_as_areanums, then inputs.setup.states_list"])
			print(txt)
			print("\nrange_as_areanums (this is the observed range that was not found):\n")
			print(range_as_areanums)
			print("\nstates_list:\n")
			print(inputs.setup.states_list)
			error(txt)
		end
	
		# Yay, a single match to the states_list was found!
		# Convert to a state index number
		observed_statenums[i] = statenums[TF][1]
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
end # END function lagrange_to_tip(inputs, geog_df)




end # end of module