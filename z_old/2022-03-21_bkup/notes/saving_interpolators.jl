
#######################################################
# Saving interpolated function into a separate file in Julia
# https://stackoverflow.com/questions/58499330/saving-interpolated-function-into-a-separate-file-in-julia
#######################################################
# https://docs.juliahub.com/JLD2/O1EyT/0.1.13/
cd("/GitHub/BioGeoJulia.jl/src")
@save "list_of_interpolators.jld2" list_of_interpolators

# Reloading
using JLD2
using Interpolations # for Interpolations.scale, Interpolations.interpolate
@load "list_of_interpolators.jld2" list_of_interpolators # <- It loads to "list_of_interpolators"...has to be this

# Find location of a module
pathof(BioGeoJulia)
interp_jld2_path = join([pathof(BioGeoJulia), ])

itp = interpolate((x,), y, Gridded(Linear()))

