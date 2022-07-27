


# GPU-Accelerated ODE Solving in R with Julia, the Language of Libraries
# August 24 2020 in Differential Equations, Julia, Programming, R, Uncategorized | Tags: diffeqr, differentialequations, gpu, high-performance, jit, r | Author: Christopher Rackauckas
# 
# https://www.stochasticlifestyle.com/gpu-accelerated-ode-solving-in-r-with-julia-the-language-of-libraries/?fbclid=IwAR0I2sH0f-civuWs1nZ2MRD9nOGtWIv4gdICWKAqP00fdZj8Nw1LphzhZ8k
# 


using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("DiffEqGPU")

using DifferentialEquations
function lorenz(du,u,p,t)
  du[1] = p[1]*(u[2]-u[1])
  du[2] = u[1]*(p[2]-u[3]) - u[2]
  du[3] = u[1]*u[2] - p[3]*u[3]
end
u0 = [1.0,1.0,1.0]
tspan = (0.0,100.0)
p = [10.0,28.0,8/3]
prob = ODEProblem(lorenz,u0,tspan,p)
sol = solve(prob,saveat=1.0)


@time for i in 1:100 solve(prob,saveat=1.0) end


using DifferentialEquations, DiffEqGPU
function lorenz(du,u,p,t)
  du[1] = p[1]*(u[2]-u[1])
  du[2] = u[1]*(p[2]-u[3]) - u[2]
  du[3] = u[1]*u[2] - p[3]*u[3]
end
u0 = [1.0,1.0,1.0]
tspan = (0.0,100.0)
p = [10.0,28.0,8/3]
prob = ODEProblem(lorenz,u0,tspan,p)
prob_func = (prob,i,repeat) -> remake(prob,u0=rand(3).*u0,p=rand(3).*p)
monteprob = EnsembleProblem(prob, prob_func = prob_func, safetycopy=false)
sol = solve(monteprob,Tsit5(),EnsembleSerial(),trajectories=100_000,saveat=1.0f0)


@time sol = solve(monteprob,Tsit5(),EnsembleSerial(),trajectories=100_000,saveat=1.0f0)

# 23.037315 seconds (16.10 M allocations: 2.022 GiB, 12.14% gc time)
# EnsembleSolution Solution of length 100000 with uType:
# ODESolution{Float64,2,Array{Array{Float64,1},1},Nothing,Nothing,Array{Float64,1},Array{Array{Array{Float64,1},1},1},ODEProblem{Array{Float64,1},Tuple{Float64,Float64},true,Array{Float64,1},ODEFunction{true,typeof(lorenz),LinearAlgebra.UniformScaling{Bool},Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing},Base.Iterators.Pairs{Union{},Union{},Tuple{},NamedTuple{(),Tuple{}}},DiffEqBase.StandardODEProblem},Tsit5,OrdinaryDiffEq.InterpolationData{ODEFunction{true,typeof(lorenz),LinearAlgebra.UniformScaling{Bool},Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing},Array{Array{Float64,1},1},Array{Float64,1},Array{Array{Array{Float64,1},1},1},OrdinaryDiffEq.Tsit5Cache{Array{Float64,1},Array{Float64,1},Array{Float64,1},OrdinaryDiffEq.Tsit5ConstantCache{Float64,Float64}}},DiffEqBase.DEStats}



# (Crash)
# @time sol = solve(monteprob,Tsit5(),EnsembleGPUArray(),trajectories=100_000,saveat=1.0f0)





# 
# R doesn't reach quite the level of Julia here, and if you profile you'll see it's because the `prob_func`, i.e. the function
# that tells you which problems to solve, is still a function written in R and this becomes the bottleneck as the computation
# becomes faster and faster. Thus you will get closer and closer to the Julia speed with longer and harder ODEs, but it still
# means there's work to be done. Another detail is that the Julia code is able to be further accelerated by using 32-bit
# numbers. Let's see that in action:


using DifferentialEquations, DiffEqGPU
function lorenz(du,u,p,t)
  du[1] = p[1]*(u[2]-u[1])
  du[2] = u[1]*(p[2]-u[3]) - u[2]
  du[3] = u[1]*u[2] - p[3]*u[3]
end
u0 = Float32[1.0,1.0,1.0]
tspan = (0.0f0,100.0f0)
p = Float32[10.0,28.0,8/3]
prob = ODEProblem(lorenz,u0,tspan,p)
prob_func = (prob,i,repeat) -> remake(prob,u0=rand(Float32,3).*u0,p=rand(Float32,3).*p)
monteprob = EnsembleProblem(prob, prob_func = prob_func, safetycopy=false)
sol = solve(monteprob,Tsit5(),EnsembleSerial(),trajectories=100_000,saveat=1.0f0)
@time sol = solve(monteprob,Tsit5(),EnsembleSerial(),trajectories=100_000,saveat=1.0f0)


# 23.542183 seconds (16.00 M allocations: 1.622 GiB, 15.66% gc time)
# EnsembleSolution Solution of length 100000 with uType:
# ODESolution{Float32,2,Array{Array{Float32,1},1},Nothing,Nothing,Array{Float32,1},Array{Array{Array{Float32,1},1},1},ODEProblem{Array{Float32,1},Tuple{Float32,Float32},true,Array{Float32,1},ODEFunction{true,typeof(lorenz),LinearAlgebra.UniformScaling{Bool},Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing},Base.Iterators.Pairs{Union{},Union{},Tuple{},NamedTuple{(),Tuple{}}},DiffEqBase.StandardODEProblem},Tsit5,OrdinaryDiffEq.InterpolationData{ODEFunction{true,typeof(lorenz),LinearAlgebra.UniformScaling{Bool},Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing,Nothing},Array{Array{Float32,1},1},Array{Float32,1},Array{Array{Array{Float32,1},1},1},OrdinaryDiffEq.Tsit5Cache{Array{Float32,1},Array{Float32,1},Array{Float32,1},OrdinaryDiffEq.Tsit5ConstantCache{Float32,Float32}}},DiffEqBase.DEStats}




# 
# I tried this code and code further than I thought I would, but the key issue is that Mac have AMD GPUs (Graphical Processing Units), whereas DiffEqGPU is programmed for NVIDIA GPUs, I think a Windows/Linux thing.  So, no GPUs for us in the short-term.
# You sent
# a few seconds ago
# 
# But, if we ever really need to crunch, we might be able to try this:
# 
# https://discourse.julialang.org/t/opencl-status-in-julia/34442/4
# ==========
# I can certainly appreciate the desire to not be locked into NVIDIA. I am the same boat with Mac. I require DifferentialEquations.jl’s ensemble GPU support which, as far as I’m aware, only works on NVIDIA GPUs. My solution is to setup a cloud VM and connect to a remote process in the VM with JUNO. It works very well.
# =========