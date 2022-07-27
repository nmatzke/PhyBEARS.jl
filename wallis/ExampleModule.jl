module ExampleModule
using DataFrames
using Plots  # for savefig
#using RCall

export example1, example2, example3, PUTTHISHERE

# Example functions to provide a model workflow for
# developing new functions

function example1()
	print("\n\nRunning example1()\n")
end

function example2()
	print("\n\nRunning example2 again()\n")	
end

function example3()
	print("\n\nRunning example3 again again()\n")	
end

function PUTTHISHERE()
	for a in 5:15
		if a > 10
			println("a is greater than 10")
		elseif a < 10
			println("a is less than 10")
		else
			println("a is equal to 10")
		end
	end
end

end # end of module
