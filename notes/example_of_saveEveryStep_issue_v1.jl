#######################################################
# 2021-09-02, Example 1
# Replicable example showing how solve() can give different
# results with saveeverystep=true or saveeverystep==false
# 
# This is apparently intended behavior, see:
# https://github.com/SciML/DifferentialEquations.jl/issues/796
# 
# ...just always use save_everystep=true to avoid the problems!
# ...or saveat=
#######################################################


#######################################################
# I stumbled on this problem over a year ago so it's not 
# just a recent version.  But my current setup:
#
# Mac OSX 10.13.6
#
# julia --version
# julia version 1.6.2
#
# Versions:
using Pkg
pkgs = Pkg.installed();

pkgs["DifferentialEquations"]
# v"6.17.2"
pkgs["OrdinaryDiffEq"]
# v"5.60.0"
pkgs["Sundials"]
# v"4.5.3"
pkgs["DiffEqDevTools"]
# v"2.27.2"
pkgs["ODEInterfaceDiffEq"]
# v"3.10.1"
pkgs["ODE"]
# v"2.13.0"
pkgs["LSODA"]
# v"0.7.0"
#######################################################


#######################################################
# Dependencies
#######################################################

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DataFrames  # for DataFrame
using DifferentialEquations
using OrdinaryDiffEq, Sundials, DiffEqDevTools, ODEInterfaceDiffEq, ODE, LSODA


#######################################################
# Parameters
#######################################################
# Number of states
n = 2

# Parameter values
mu_vals = [0.0, 0.0]
Qij_vals = [0.0, 0.0]
Cijk_weights = [1.0, 1.0]
Cijk_vals = [0.2222222, 0.2222222]
params = (mu_vals=mu_vals, Qij_vals=Qij_vals, Cijk_weights=Cijk_weights, Cijk_vals=Cijk_vals)

# Indices
Qarray_ivals = [1, 2]
Qarray_jvals = [2, 1]
Qarray_event_types = ["a", "a"]
Carray_ivals = [1, 2]
Carray_jvals = [1, 2]
Carray_kvals = [1, 2]
Carray_event_types = ["y", "y"]

p_indices = (Qarray_ivals=Qarray_ivals, Qarray_jvals=Qarray_jvals, Qarray_event_types=Qarray_event_types, Carray_ivals=Carray_ivals, Carray_jvals=Carray_jvals, Carray_kvals=Carray_kvals, Carray_event_types=Carray_event_types)

Qi_eq_i = Any[Bool[1, 0], Bool[0, 1]]
Ci_eq_i = Any[Bool[1, 0], Bool[0, 1]]
Qi_sub_i = Any[[1], [2]]
Qj_sub_i = Any[[2], [1]]
Ci_sub_i = Any[[1], [2]]
Cj_sub_i = Any[[1], [2]]
Ck_sub_i = Any[[1], [2]]

p_TFs = (Qi_eq_i=Qi_eq_i, Ci_eq_i=Ci_eq_i, Qi_sub_i=Qi_sub_i, Qj_sub_i=Qj_sub_i, Ci_sub_i=Ci_sub_i, Cj_sub_i=Cj_sub_i, Ck_sub_i=Ck_sub_i)

# Starting values
uE = [0.0, 0.0]

p = (n=n, params=params, p_indices=p_indices, p_TFs=p_TFs, uE=uE)



#######################################################
# Differential equation
#######################################################
diffeq_Es = (du,u,p,t) -> begin

  # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Qij_vals = p.params.Qij_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Qarray_ivals = p.p_indices.Qarray_ivals
	Qarray_jvals = p.p_indices.Qarray_jvals
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
	# Pre-calculated solution of the Es
	uE = p.uE
	
	two = 1.0
  @inbounds for i in 1:n
		# Calculation of "D" (likelihood of tip data)
		du[i] = -(sum(Cijk_vals[Carray_ivals .== i]) + sum(Qij_vals[Qarray_ivals .== i]) + mu[i])*u[i] +
			(sum(Qij_vals[Qarray_ivals .== i] .* u[(Qarray_jvals[Qarray_ivals .== i])]))	
  end
end


# Timespan to solve over
Es_tspan = (0.0, 1.4)

# starting state
uE = [0.0, 1.0]

# Define the problem
prob_Es_v5 = DifferentialEquations.ODEProblem(diffeq_Es, uE, Es_tspan, p)




##################################################################################
# EXAMPLES OF DIFFERING RESULTS BETWEEN save_everystep=false AND save_everystep=true
##################################################################################

# The right answer is:
exp(-0.2222222)
# 0.8007374207109728
#
#   are provided by 
# Tsit5, save_everystep=true
# CVODE, save_everystep=true
#
# ...but not save_everystep=false
# 
# LSODA, save_everystep=false is closer but not really good, although typically on other problems
# LSODA save_everystep=true is closer but again not really good.
# 


# Tsit5
sol_Es_v5_Tsit5_savefalse = solve(prob_Es_v5, Tsit5(), save_everystep=false, abstol=1e-9, reltol=1e-9);
sol_Es_v5_Tsit5_savefalse(1.0)
# 0.0
# 0.8090232075241014

sol_Es_v5_Tsit5_savetrue = solve(prob_Es_v5, Tsit5(), save_everystep=true, abstol=1e-9, reltol=1e-9);
sol_Es_v5_Tsit5_savetrue(1.0)
# 0.0
# 0.8007374207059241


# CVODE
sol_Es_v5_cvode_savefalse = solve(prob_Es_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=false, abstol=1e-9, reltol=1e-9);
sol_Es_v5_cvode_savefalse(1.0)
# 0.0
# 0.8090232023651034

sol_Es_v5_cvode_savetrue = solve(prob_Es_v5, CVODE_BDF(linear_solver=:GMRES), save_everystep=true, abstol=1e-9, reltol=1e-9);
sol_Es_v5_cvode_savetrue(1.0)
# 0.0
# 0.8007374145719871


# LSODA
sol_Es_v5_lsoda_savefalse = solve(prob_Es_v5, lsoda(), save_everystep=false, abstol=1e-9, reltol=1e-9);
sol_Es_v5_lsoda_savefalse(1.0)
# 0.0
# 0.8090232070918716

sol_Es_v5_lsoda_savetrue = solve(prob_Es_v5, lsoda(), save_everystep=true, abstol=1e-9, reltol=1e-9);
sol_Es_v5_lsoda_savetrue(1.0)
# 0.0
# 0.8012403954152403

