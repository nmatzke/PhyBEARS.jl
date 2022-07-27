module ExampleModule
using DataFrames
using Plots  # for savefig
#using RCall

export example1, example2

# Example functions to provide a model workflow for
# developing new functions

function example1()
	print("\n\nRunning example1()\n")
end

function example2()
	print("\n\nRunning example2 again()\n")	
end


end # end of module