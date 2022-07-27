#######################################################
# work precision diagrams
#######################################################

module WorkPrecision

print("\n\nStarting module 'WorkPrecision'...loading dependencies...\n")
using BenchmarkTools # for @time
using InvertedIndices # for Not
using LSODA
using DifferentialEquations
using Distributed
using Random					# for MersenneTwister()
using Dates						# for e.g. DateTime, Dates.now()
using PhyloNetworks
#using Plots						# for plot
using DataFrames          # for DataFrame()

using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
#Pkg.add(PackageSpec(url="https://github.com/JuliaDiffEq/deSolveDiffEq.jl"))
using deSolveDiffEq 
# https://docs.juliadiffeq.org/stable/solvers/ode_solve/index.html


using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs

# (1) List all function names here:
export say_hello3, workprecision

#######################################################
# Temporary file to store functions under development
#
# Start with:
# 
# Setup:

"""
cd("/GitHub/BioGeoJulia.jl/notes/")
include("tst3.jl")

cd("/GitHub/BioGeoJulia.jl/notes/")
include("WorkPrecision.jl")

"""
#######################################################


#######################################################
# (2) write the functions here
#######################################################

say_hello3() = println("Hello dude2!")


"""
using Pkg
using StatsPlots
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
#Pkg.add(PackageSpec(url="https://github.com/JuliaDiffEq/deSolveDiffEq.jl"))
using deSolveDiffEq 
# https://docs.juliadiffeq.org/stable/solvers/ode_solve/index.html

cd("/GitHub/BioGeoJulia.jl/notes/")
include("ModelLikes.jl")
import .ModelLikes
#using .Tmp

using Profile     # for @profile
using DataFrames  # for DataFrame
using PhyloNetworks
using BioGeoJulia.TrUtils # for flat2() (similar to unlist)
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs

using RCall       # for df_to_Rdata, reval, g = globalEnv

using Random
Random.seed!(123)

# 2D Linear ODE
function f(du,u,p,t)
  @inbounds for i in eachindex(u)
    du[i] = 1.01*u[i]
  end
end
function f_analytic(u₀,p,t)
  u₀*exp(1.01*t)
end
tspan = (0.0,10.0)
prob = ODEProblem(ODEFunction(f,analytic=f_analytic),rand(100,100),tspan)

abstols = 1.0 ./ 10.0 .^ (3:13)
reltols = 1.0 ./ 10.0 .^ (0:10);



setups = [Dict(:alg=>DP5())
          Dict(:alg=>ode45())
          Dict(:alg=>ARKODE(Sundials.Explicit(),etable=Sundials.DORMAND_PRINCE_7_4_5))
          Dict(:alg=>Tsit5())]
solnames = ["OrdinaryDiffEq";"ODE";"Sundials ARKODE";"OrdinaryDiffEq Tsit5"]

df = ModelLikes.workprecision(prob, setups, abstols, reltols, solnames=solnames, save_everystep=false,numruns=1)




"""

function workprecision(prob, setups, abstols=1.0 ./ 10.0 .^ (3:13), reltols=1.0 ./ 10.0 .^ (0:10); solnames=string.(collect(1:length(setups))), save_everystep=false, numruns=1, dt=0.01)
	
	# First column is the tolerances
	results_matrix = abstols
	
	ival = 0
	# Go through the dictionary of setups
	for setup in setups
		#tmpdict = setups[i]
		solver = setup[:alg]
		ival += 1
		
		calctimes = collect(repeat([0.0], length(abstols)))
		
		print(paste0(["\n\nRunning ", solnames[ival], "...\n"]))
		
		for j in 1:length(abstols)
			if (j == 1)
				# Preliminary run to do compiling
				if (startswith(solnames[ival], "E"))
					sol = solve(prob, solver, save_everystep=save_everystep, abstol=abstols[j], reltol=reltols[j], dt=dt);
				else
					sol = solve(prob, solver, save_everystep=save_everystep, abstol=abstols[j], reltol=reltols[j]);
					#print(sol)
					#print("\n\n")
				end
			end
			
			diagnostics = collect(repeat([Dates.now()], 3))
			diagnostics[1] = Dates.now()		

			if (startswith(solnames[ival], "E"))
				sol = solve(prob, solver, save_everystep=save_everystep, abstol=abstols[j], reltol=reltols[j], dt=dt);
			else
				sol = solve(prob, solver, save_everystep=save_everystep, abstol=abstols[j], reltol=reltols[j]);
			end
			
			# Final run diagnostics
			diagnostics[2] = Dates.now()
			diagnostics[3] = diagnostics[2]-diagnostics[1]

			if ((sol.retcode == :Success) || (sol.retcode == :Default))
				total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
			else
				total_calctime_in_sec = NaN
			end
			calctimes[j] = total_calctime_in_sec
			
		end
		results_matrix = hcat(results_matrix, calctimes)
	end
	
	df = DataFrame(results_matrix)
	solnames2 = deepcopy(solnames)
	prepend!(solnames2, ["abstol"])
	rename!(df, solnames2)
	
	return df
end # End function workprecision



end # End Module
