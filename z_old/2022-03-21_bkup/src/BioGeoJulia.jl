module BioGeoJulia
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/
export hello_BioGeoJulia, add_one_BioGeoJulia

print("\nBioGeoJulia: loading BioGeoJulia.jl")

# List each BioGeoJulia code file here
# NOTE: LOAD THE DEPENDENCY .jl FILES *FIRST*, or you get "not recognized" errors

include("BGExample.jl")			# default examples
include("TrUtils.jl")			# basic utility functions 
include("MaxentInterp.jl")	# preconstructed interpolator for weighting rangesize of smaller daughter
include("TreeTable.jl")			# for prt() tree tables (DFs), bd_liks(), etc.
include("StateSpace.jl")	# set up lists of areas and states (geographic ranges)
include("SSEs.jl")				# SSE calculations with various amounts of speed optimization
include("Parsers.jl")			# Parsers to read e.g. geography file
include("TreePass.jl")		# downpass and uppass through the phylogeny; prt() etc.


#######################################################
# Put functions here
#######################################################
"""
    hello(who::String)

Return "BioGeoJulia says, hi `who`".
"""
hello_BioGeoJulia(who::String) = "hello_BioGeoJulia() says '$who'"

"""
    add_one_BioGeoJulia(x::Number)

Return `x + 1`.
"""
add_one_BioGeoJulia(x::Number) = x + 1

end # Ends the module command
