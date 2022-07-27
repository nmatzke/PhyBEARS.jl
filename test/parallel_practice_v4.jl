"""
Note:
JULIA_NUM_THREADS must be defined before starting julia; defining it in startup.jl is too late in the startup process.

https://docs.julialang.org/en/v1.0/manual/environment-variables/#JULIA_NUM_THREADS-1
"""


# Good page on Base.Threads.@spawn (INSTEAD OF Distributed.@spawn)
# Notes on Multithreading with Julia
# Eric Aubanel
# June 2020
# (revised August 2020)
# http://www.cs.unb.ca/~aubanel/JuliaMultithreadingNotes.html


"""

> JULIA_NUM_THREADS=20 julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.3.0-rc1.0 (2019-08-18)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

Threads.nthreads()
8
"""

import Base.Threads.@spawn  # <-- this is good, Distributed.@spawn is BAD, it produces Futures that have to be scheduled etc]
Threads.nthreads()

using Dates			 # for DateTime()
using Random     # for MersenneTwister()
using Statistics     # for mean()
using BenchmarkTools # for @benchmark
using DataFrames
using PhyloNetworks
using Base.Threads			# for e.g. Base.Threads.@spawn  (NOTE: Distributed.@spawn is different, deprecated, and BAD - returns a Future)


function countloop2(num_iterations, current_nodeIndex)
	txt = join(["\nThread ID: ", string(Threads.threadid())], "")
	print(txt)
	x = 0.0
	random_number_generator = MersenneTwister(current_nodeIndex);

	for i in 1:num_iterations
	   x = x + (randn(random_number_generator, 1)[1] / num_iterations)
	end
	txt = join(["\nThread ID: ", string(Threads.threadid()), "...done."], "")	
	return(x)
end


countloop_num_iterations = 10000000
y = countloop2(countloop_num_iterations, 1)
y1 = Base.Threads.@spawn countloop2(countloop_num_iterations, 1)
y2 = Base.Threads.@spawn countloop2(countloop_num_iterations, 1)
y3 = Base.Threads.@spawn countloop2(countloop_num_iterations, 1)
y4 = Base.Threads.@spawn countloop2(countloop_num_iterations, 1)
y5 = Base.Threads.@spawn countloop2(countloop_num_iterations, 1)

