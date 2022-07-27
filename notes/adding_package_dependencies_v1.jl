
#
# This worked on 2022-03-12 with old versions:
	# ...assumes Pkg.add(name="SCS", version="0.9.0")
	# ...assumes Pkg.add(name="Convex", version="0.14.18")

# Check versions with
Pkg.status("SCS")




# Adding LSODA to package environment:

cd("/GitHub/PhyBEARS.jl")
]

pkg> comes up...

# ACTIVATE, **THEN** add packages
activate .
add LoopVectorization
add LSODA
add Distributions
#add Convex
#add SCS
add PhyloNetworks		# for trees etc.
add LinearAlgebra		# for factorize()
resolve

# Changes to .toml files appear on GitHub
# Commit & Push


backspace to exit pkg>


;;j to:

exit()




pwd

julia



