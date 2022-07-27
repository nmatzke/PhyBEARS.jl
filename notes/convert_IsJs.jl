
# Just convert is, js to a single-number index for a matrix
"""
include("/GitHub/PhyBEARS.jl/notes/convert_IsJs.jl")

A=[1 2 3; 4 5 6; 7 8 9]
# See how the single-number indices count down the columns
A[1:9]
n = dim(A)[1]
is= [[1, 1, 1, 1, 1, 1, 1, 1],
 [2, 2, 2, 2, 2, 2, 2, 2],
 [3, 3, 3, 3, 3, 3, 3, 3]]

js= [[1, 1, 1, 1, 1, 1, 1, 1],
 [2, 2, 2, 2, 2, 2, 2, 2],
 [3, 3, 3, 3, 3, 3, 3, 3]]

indices = convert_is_js_to_single_index2(is, js, n)
"""

function convert_is_js_to_single_index2(is, js, n)
	if length(is) != length(js)
		txt = "STOP ERROR in convert_is_js_to_single_index(is, js, n). The lengths of 'is' and 'js' must be identical. Check these inputs and re-run."
		throw(txt)
	end
	
	
	tmpvec = Any[]
	indices = collect(repeat([tmpvec], length(is)))
	for tmpi in 1:length(indices)
		if maximum(is[tmpi])[1] > n
			txt = "STOP ERROR in convert_is_js_to_single_index(is, js, n). The maximum of 'is' cannot be larger than 'n'.  Check these inputs and re-run."
			throw(txt)
		end
	
		if maximum(js[tmpi])[1] > n
			txt = "STOP ERROR in convert_is_js_to_single_index(is, js, n). The maximum of 'js' cannot be larger than 'n'.  Check these inputs and re-run."
			throw(txt)
		end

		#indices[i] .= ((js[i].-1) * n) + is[i]
		single_num_indices = is[tmpi] .+ (js[tmpi] .-1) .* n
		indices[tmpi] = single_num_indices
		# Adding the [1]'s mean that this error is avoided:
		# ERROR: MethodError: no method matching -(::Vector{Int64}, ::Int64)
		# ...which comes from mistmatched types
	end # END for i in 1:length(indices)
	
	return indices
end # END function convert_is_js_to_single_index(is, js, n)
