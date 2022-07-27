# Pkg.rm("BioGeoJulia")
# Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))

using DataFrames
using PhyloNetworks
using BioGeoJulia

using BioGeoJulia.TrUtils
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs





#great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
great_ape_newick_string = "((chimp:1,human:1):1,gorilla:2);"
tro = readTopology(great_ape_newick_string)
tro
# HybridNetwork, Rooted Network
# 4 edges
# 5 nodes: 3 tips, 0 hybrid nodes, 2 internal tree nodes.
# tip labels: chimp, human, gorilla
# ((chimp:1.0,human:1.0):1.0,gorilla:2.0);



rootnodenum = tro.root
trdf = prt(tro, rootnodenum)
trdf






length(tro.node)
# 5

tro.edge[1]


newnode, newedge = PhyloNetworks.breakedge!(tro.edge[1], tro);

newnode
# PhyloNetworks.Node:
#  number:4
#  attached to 2 edges, numbered: 1 5

newedge # new edge 5 goes from node -3 and 4 (new)
# PhyloNetworks.Edge:
#  number:5
#  length:0.5
#  attached to 2 node(s) (parent first): -3 4


length(tr.node) # one more than before
# 6

writeTopology(tro) # note extra pair of parentheses around chimp
"((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"




tro = readTopology(great_ape_newick_string)
tr = readTopology("((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);")


tro
tr

# Node numbers
tro.node
tr.node

# Edge numbers
tro.edge
tr.edge

prt(tro, tr.root)
# Works

prt(tr, tr.root)
# Error in prt which means to me that:

# 1. We need to add the edge numbers to the prt table
# 2. Need to change prt to allow 2-degree nodes

# ERROR: BoundsError: attempt to access 6-element Array{Int64,1} at index [7]
# Stacktrace:
#  [1] setindex! at ./array.jl:847 [inlined]
#  [2] get_nodenumbers_above_node(::HybridNetwork, ::Int64, ::Array{Int64,1}, ::Int64; indexNum_table::Array{Int64,2}) at /Users/nickm/.julia/packages/BioGeoJulia/uDruG/src/TreePass.jl:54
#  [3] get_nodenumbers_above_node(::HybridNetwork, ::Int64, ::Array{Int64,1}, ::Int64; indexNum_table::Array{Int64,2}) at /Users/nickm/.julia/packages/BioGeoJulia/uDruG/src/TreePass.jl:57 (repeats 4 times)
#  [4] #get_nodenumbers_above_node#1 at /Users/nickm/.julia/packages/BioGeoJulia/uDruG/src/TreePass.jl:58 [inlined]
#  [5] get_LR_uppass_nodeIndexes(::HybridNetwork) at /Users/nickm/.julia/packages/BioGeoJulia/uDruG/src/TreePass.jl:312
#  [6] get_node_heights(::HybridNetwork) at /Users/nickm/.julia/packages/BioGeoJulia/uDruG/src/TreePass.jl:873
#  [7] get_node_ages(::HybridNetwork) at /Users/nickm/.julia/packages/BioGeoJulia/uDruG/src/TreePass.jl:910
#  [8] prt(::HybridNetwork, ::Int64, ::Bool) at /Users/nickm/.julia/packages/BioGeoJulia/uDruG/src/TreePass.jl:521
#  [9] prt(::HybridNetwork, ::Int64) at /Users/nickm/.julia/packages/BioGeoJulia/uDruG/src/TreePass.jl:499
#  [10] top-level scope at REPL[139]:3

