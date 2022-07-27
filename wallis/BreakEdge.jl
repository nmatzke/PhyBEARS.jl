using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using PhyloNetworks					# most maintained, emphasize; for HybridNetwork
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames

using LinearAlgebra  # for "I" in: Matrix{Float64}(I, 2, 2)
										 # https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using DataFrames  # for DataFrame
using DifferentialEquations
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA


# List each PhyBEARS code file prefix here
using PhyBEARS.Example
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.TrUtils
using PhyBEARS.SSEs

include("/GitHub/PhyBEARS.jl/src/TreePass.jl")
import .TreePass

# Repeat calculation in Julia
include("/GitHub/PhyBEARS.jl/notes/ModelLikes.jl")
import .ModelLikes

########################



test = readTopology("(((chimp:1,human:1):1,gorilla:2):1,orang:3);")
length(test.node)
test.edge[2]

newnode, newedge = PhyloNetworks.breakedge!(test.edge[2], test)
length(test.node)
writeTopology(test)

in_params = (birthRate=0.2, deathRate=0.1, d_val=0.0, e_val=0.0, a_val=0.0, j_val=0.0)
numareas = 2
n = 4

numareas=2
tr=readTopology("((chimp:1,human:1):1,gorilla:2);")
inputs = ModelLikes.setup_DEC_SSE(numareas, test; root_age_mult=1.5, max_range_size=NaN, include_null_range=false, in_params=in_params)
(setup, res, trdf, solver_options, p_Ds_v5, Es_tspan) = inputs

# getting and issue here. 
# Error reads: ERROR: BoundsError: attempt to access Int64
#  at index [2]

# Also cannot convert tr.edge[2] to an float. is it trying to look for an Int but only finding a floater? hmm

"""
PhyloNetwork's practice code reads as: 

But even in their practice, it breaks at newnode, newedge

net = readTopology("(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));")
length(net.node)
net.edge[4]
newnode, newedge = PhyloNetworks.breakedge!(net.edge[4], net)
length(net.node)
newedge # new edge 21 goes from node -8 and 11 (new)
net.edge[4] # original edge 4 now goes from node 11 (new) to 3
writeTopology(net)

"""

"""
Lets have a look at PhyloNetwork's Code for Phylonetwork.breakedge!
"""

function breakedge!(edge::Edge, net::HybridNetwork)
    pn = getParent(edge) # parent node
    # new child edge = old edge, same hybrid attribute
    removeEdge!(pn,edge)
    removeNode!(pn,edge)
    max_edge = maximum(e.number for e in net.edge)
    max_node = maximum(n.number for n in net.node)
    newedge = Edge(max_edge+1) # create new parent (tree) edge
    newnode = Node(max_node+1,false,false,[edge,newedge]) # tree node
    setNode!(edge,newnode) # newnode comes 2nd, and parent node along 'edge'
    edge.isChild1 = true
    setNode!(newedge,newnode) # newnode comes 1st in newedge, but child node
    newedge.isChild1 = true
    setEdge!(pn,newedge)
    setNode!(newedge,pn) # pn comes 2nd in newedge
    if edge.length == -1.0
        newedge.length = -1.0
    else
        edge.length /= 2
        newedge.length = edge.length
    end
    newedge.containRoot = edge.containRoot
    pushEdge!(net,newedge)
    pushNode!(net,newnode)
    return newnode, newedge
end



"""
Find function definitions for:
    setEdge!   --> Found in PhyloNetworks.auxillary.jl
    setNode!   --> Found in PhyloNetworks.auxillary.jl
    getParent  --> getParents (with an s) found in PhyloNetworks.manipulateNet.jl
                   singular getParents found in PhyloNetworks.auxillary.jl
    removeEdge! --> Found in PhyloNetworks.auxillary.jl
    removeNode! --> Found in PhyloNetworks.auxillary.jl
    pushEdge --> Found in PhyloNetworks.auxillary.jl 
    pushNode --> Found in PhyloNetworks.auxillary.jl

    etc etc etc. all are undefined when run alone?
    WHYYYYY

See if Nick knows on Friday?

We know what they MEAN, but I cant figure out where it's trying to pull an Int64 from?
"""

# Walking through breakedge function

edge = net.edge[4]
net

using PhyloNetworks

@inline function getParents(node::Node)
    parents = Node[]
    for e in node.edge
            if node == getChild(e)
                push!(parents, getParent(e))
            end
    end
    return parents
end

# I end up with "Node not defined?"
# It will not let me initialize a blank array? Node[], x[], nada works?
# Julia...what's going on girl?

pn = getParents(edge) # parent node
# new child edge = old edge, same hybrid attribute

removeEdge!(pn,edge)
removeNode!(pn,edge)
max_edge = maximum(e.number for e in net.edge)
max_node = maximum(n.number for n in net.node)
newedge = Edge(max_edge+1) # create new parent (tree) edge
newnode = Node(max_node+1,false,false,[edge,newedge]) # tree node
setNode!(edge,newnode) # newnode comes 2nd, and parent node along 'edge'
edge.isChild1 = true
setNode!(newedge,newnode) # newnode comes 1st in newedge, but child node
newedge.isChild1 = true
setEdge!(pn,newedge)
setNode!(newedge,pn) # pn comes 2nd in newedge
if edge.length == -1.0
    newedge.length = -1.0
else
    edge.length /= 2
    newedge.length = edge.length
end
newedge.containRoot = edge.containRoot
pushEdge!(net,newedge)
pushNode!(net,newnode)
return newnode, newedge

