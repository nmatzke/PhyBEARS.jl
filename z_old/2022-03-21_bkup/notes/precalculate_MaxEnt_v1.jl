#######################################################
# Pre-calculate MaxEnt solutions for any reasonable n
# (number of areas), so that we don't have to 
# have a numerical solver step, nor use the
# annoying SCS package that keeps changing
# syntax
#######################################################

#######################################################
# Example use of maximum entropy on discrete case
# https://github.com/JuliaOpt/Convex.jl/issues/64
#######################################################
using Distributions  # for quantile
using Convex				 # for Convex.entropy(), maximize()
using SCS						 # for SCSSolve, solve (maximize(entropy()))
using Interpolations # for Interpolations.scale, Interpolations.interpolate


maxent_constraint_01 = 0.0

points_to_interpolate_at = seq(0.0, 1.0, 0.01)


numareas_in_ancestor = collect(1:20)
maxent_interp = []

print("\n\nConstructing interpolator #: ")
for j in 1:length(numareas_in_ancestor)
	print(paste0([j, " "]))
	
	list_of_vectors = []
	for i in 1:length(points_to_interpolate_at)

		#total_numareas = 6
		#n = total_numareas
		n = numareas_in_ancestor[j]
		x = 1:n
	
		maxent_constraint_01 = points_to_interpolate_at[i]
	
		#discrete_values_padded = cat(0, collect(x)[:], n+1; dims=1)
		discrete_values_padded = collect(x)
		maxent_constraint = quantile(discrete_values_padded, maxent_constraint_01)
		#maxent_constraint = 1

		# Vector of probabilities that must sum to 1
		p = Variable(length(x))
		probability_contraints = [0 <= p, p <= 1, sum(p) == 1];
		feature_constraints = sum(p'*x) == maxent_constraint

		# base or prior (e.g. Uniform)
		#h = pdf.(Uniform(1, n), x)

		# This solution updates p (look in p.values)
		problem = Convex.maximize(Convex.entropy(p), 0 <= p, p <= 1, sum(p) == 1, feature_constraints)

		# worked 2020-09-29 at work, failed 2020-09-30 at home
		#sol = Convex.solve!(problem, SCS.SCSSolver(verbose=0))

		# worked e.g. July 2020 (version issue I guess)
		# worked at home 2020-09-30
		# worked at home & laptop, 2022-03-11
		sol = Convex.solve!(problem, SCS.Optimizer(verbose=0))
		# ...assumes Pkg.add(name="SCS", version="0.9.0")
		# ...assumes Pkg.add(name="Convex", version="0.14.18")

		maxent_result = abs.(round.(p.value; digits=4))
		push!(list_of_vectors, maxent_result)
	end # for i in 1:length(points_to_interpolate_at)
	
	# Make the interpolator from one vector to the next
	itp = interpolate((points_to_interpolate_at,), list_of_vectors, Gridded(Linear()))

	push!(maxent_interp, itp)
end # for j in 1:length(numareas_in_ancestor)

maxent_interp


#######################################################
# Saving interpolated function into a separate file in Julia
# https://stackoverflow.com/questions/58499330/saving-interpolated-function-into-a-separate-file-in-julia
#######################################################
# https://docs.juliahub.com/JLD2/O1EyT/0.1.13/
cd("/GitHub/BioGeoJulia.jl/src")
@save "maxent_interp.jld2" maxent_interp

# Reloading
using JLD2
using Interpolations # for Interpolations.scale, Interpolations.interpolate
@load "maxent_interp.jld2" maxent_interp # <- It loads to "maxent_interp"...has to be this

# Find location of a module
pathof(BioGeoJulia)
path0 = replace(pathof(BioGeoJulia), "BioGeoJulia.jl/"=>"tmpabc.jl/")
path1 = replace(path0, "BioGeoJulia.jl"=>"")
path2 = replace(path1, "tmpabc.jl/"=>"BioGeoJulia.jl/")

interp_jld2_path = join([path2, "maxent_interp.jld2"])
@load interp_jld2_path maxent_interp # 
length(maxent_interp)



itp = interpolate((x,), y, Gridded(Linear()))



using Interpolations, Plots

julia> y = [0.0, 0.1, 0.85, 0.9, 1.0];

julia> x = range(0, stop = 1, length = 5);

julia> itp = interpolate((x,), y, Gridded(Linear()))

using Interpolations, Plots
t = 0:.2:1
x, y = 2sin.(π*t), cos.(π*t)
itp = Interpolations.scale(interpolate([x y], (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)
tfine = 0:.01:1
xs, ys = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine]
x0, y0 = itp(0.5,1), itp(0.5,2) # interpolate point at t=0.5
plot(xs, ys, aspect_ratio=1, label="BSpline", title="Interpolations.jl parametric BSpline")
scatter!(x, y, label="input points")
scatter!([x0], [y0], ms=5, mc=:red, label="Interpolated") 