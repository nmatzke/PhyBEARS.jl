module TreePass
__precompile__(false)  # will cause using / import to load it directly into the 
                       # current process and skip the precompile and caching. 
                       # This also thereby prevents the module from being 
                       # imported by any other precompiled module.
                       # https://docs.julialang.org/en/v1/manual/modules/
using BioGeoJulia.TrUtils # for e.g. flat2
using BioGeoJulia.SSEs 
#using BioGeoJulia.Flow 
using DataFrames			# for e.g. DataFrame()
using PhyloNetworks		# for e.g. readTopology()
using Dates						# for e.g. DateTime, Dates.now()
using Distributed			# for e.g. @spawn
using Random					# for MersenneTwister()
using DifferentialEquations # for ODEProblem
using LSODA						# for lsoda()
export get_nodenumbers_above_node, get_postorder_nodenumbers_above_node, initialize_edgematrix, get_pruningwise_postorder_edgematrix, get_LR_uppass_edgematrix, get_LR_downpass_edgematrix, get_LR_uppass_nodeIndexes, get_LR_downpass_nodeIndexes, get_Rnodenums, get_nodeIndex_PNnumber, get_nodeIndex_from_PNnumber, prt, get_taxa_descending_from_each_node, isTip_TF, get_NodeIndexes_from_edge, get_NodeIndex_df_by_tree_edges, get_node_heights, get_node_ages, SolverOpt, construct_SolverOpt, Res, construct_Res, count_nodes_finished, nodeOp_average_likes, nodeOp, nodeOp_Cmat, nodeOp_ClaSSE_v5!, branchOp_example, branchOp_ClaSSE_Ds_v5, branchOp, setup_inputs_branchOp_ClaSSE_Ds_v5, countloop, iterative_downpass!, iterative_downpass_nonparallel_ClaSSE_v5!, iterative_downpass_nonparallel!






# Get the 2 nodeIndexes descending from a node; iterates for uppass
# (assumes a binary, PhyloNetwork, rooted tree)
# Reverse for downpass

"""
using DataFrames
using PhyloNetworks

# For Nick's editing (ignore)
include("/GitHub/BioGeoJulia.jl/notes/TreePassO.jl")

#######################################################
# Typical bifurcating (binary) tree
#######################################################
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreePassO.get_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table )


#######################################################
# Tree with a 2-degree node inside a branch
#######################################################
include("/GitHub/BioGeoJulia.jl/notes/TreePassO.jl")

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreePassO.get_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table )
"""

function get_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum; indexNum_table)
  if (tr.node[rootnodenum].leaf != true)
  	# * A typical internal node will be attached to 3 edges
  	#   (left descendant, right descendant, ancestor edge)
  	# * A root node will be attached to 2 edges (left, right descendant edges
  	# * A degree-2 (mid-branch) node will have 1 descendant edge

  	if (length(tr.node[rootnodenum].edge) > 3)
  		txt = join(["STOP ERROR in get_nodenumbers_above_node(): tr.node[rootnodenum=", string(rootnodenum), "] has more than 3 edges. Probably this is a multifurcating node, which is not allowed."])
  		error(txt)
  	end

  	
  	# Is the current node a typical internal node?
  	typicalTF = length(tr.node[rootnodenum].edge) == 3
  	
  	# Is the current node the root?
  	root_PNnumber = tr.node[tr.root].number  # PhyloNetworks node number of root
  	# rootnodenum = current node being examined
  	current_PNnumber = tr.node[rootnodenum].number
  	rootTF = root_PNnumber == current_PNnumber
		
		# If typical or root, proceed
		typical_or_root_TF = typicalTF || rootTF
  	if (typical_or_root_TF == true)
			# Left descendant edge
			one_edge = tr.node[rootnodenum].edge[1]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			left_dec_PNnumber = anc_decPNnumbers[2]
			left_dec_nodeIndex = get_nodeIndex_from_PNnumber(left_dec_PNnumber, indexNum_table=indexNum_table)
		
			one_edge = tr.node[rootnodenum].edge[2]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			right_dec_PNnumber = anc_decPNnumbers[2]
			right_dec_nodeIndex = get_nodeIndex_from_PNnumber(right_dec_PNnumber, indexNum_table=indexNum_table)
		
			# Then, iterate through left and right clades
			#println(rootnodenum)
			nodeIndex_array[iterNum] = rootnodenum
			#print(nodeIndex_array)
			iterNum = iterNum + 1
			(nodeIndex_array, iterNum) = get_nodenumbers_above_node(tr, right_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			(nodeIndex_array, iterNum) = get_nodenumbers_above_node(tr, left_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			return (nodeIndex_array, iterNum)
		else # It must be a nontypical node (i.e. a 2-degree, inside branch)
			# Check number of edges
			num_edges_attached_to_this_node = length(tr.node[rootnodenum].edge)
			if (num_edges_attached_to_this_node != 2)
				txt = ["ERROR in get_nodenumbers_above_node(): node ", string(rootnodenum), " has ", string(num_edges_attached_to_this_node), " edges attached. It should have exactly 2. Exiting with error."]
				error(join(txt, ""))
			end # END if (num_edges_attached_to_this_node != 2)
			
			one_edge = tr.node[rootnodenum].edge[1]

			# Ancestor node PNnumber for the descendant edge
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			dec_PNnumber = anc_decPNnumbers[2]
			dec_nodeIndex = get_nodeIndex_from_PNnumber(dec_PNnumber, indexNum_table=indexNum_table)
			
			# Ancestor edge not needed
			
			# Then, iterate through single descendant clade
			nodeIndex_array[iterNum] = rootnodenum
			iterNum = iterNum + 1
			(nodeIndex_array, iterNum) = get_nodenumbers_above_node(tr, dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			return (nodeIndex_array, iterNum)
		end # END if (typical_or_root_TF == true)
  else
  	# Leaf node:
  	#println(rootnodenum)
  	nodeIndex_array[iterNum] = rootnodenum
  	#print(nodeIndex_array)
  	iterNum = iterNum + 1
  	return (nodeIndex_array, iterNum)
  end

	# Error check
	txt = ["ERROR in get_nodenumbers_above_node(): You shouldn't be able to get to the last line of this function."]
	error(join(txt, ""))

end # END get_nodenumbers_above_node


"""
using DataFrames
using PhyloNetworks

# For Nick's editing (ignore)
include("/GitHub/BioGeoJulia.jl/notes/TreePassO.jl")

#######################################################
# Typical bifurcating (binary) tree
#######################################################
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreePassO.get_postorder_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table)


#######################################################
# Tree with a 2-degree node inside a branch
#######################################################
include("/GitHub/BioGeoJulia.jl/notes/TreePassO.jl")

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreePassO.get_postorder_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table )
"""

function get_postorder_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum; indexNum_table)
  if (tr.node[rootnodenum].leaf != true)
  	# * A typical internal node will be attached to 3 edges
  	#   (left descendant, right descendant, ancestor edge)
  	# * A root node will be attached to 2 edges (left, right descendant edges
  	# * A degree-2 (mid-branch) node will have 1 descendant edge

  	if (length(tr.node[rootnodenum].edge) > 3)
  		txt = join(["STOP ERROR in get_postorder_nodenumbers_above_node(): tr.node[rootnodenum=", string(rootnodenum), "] has more than 3 edges. Probably this is a multifurcating node, which is not allowed."])
  		error(txt)
  	end
  	
  	# Is the current node a typical internal node?
  	typicalTF = length(tr.node[rootnodenum].edge) == 3
  	
  	# Is the current node the root?
  	root_PNnumber = tr.node[tr.root].number  # PhyloNetworks node number of root
  	# rootnodenum = current node being examined
  	current_PNnumber = tr.node[rootnodenum].number
  	rootTF = root_PNnumber == current_PNnumber
		
		# If typical or root, proceed
		typical_or_root_TF = typicalTF || rootTF
  	if (typical_or_root_TF == true)
			# Left descendant edge
			one_edge = tr.node[rootnodenum].edge[1]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			left_dec_PNnumber = anc_decPNnumbers[2]
			left_dec_nodeIndex = get_nodeIndex_from_PNnumber(left_dec_PNnumber, indexNum_table=indexNum_table)
		
			one_edge = tr.node[rootnodenum].edge[2]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			right_dec_PNnumber = anc_decPNnumbers[2]
			right_dec_nodeIndex = get_nodeIndex_from_PNnumber(right_dec_PNnumber, indexNum_table=indexNum_table)
		
			# Then, iterate through left and right clades
			#println(rootnodenum)
			nodeIndex_array[iterNum] = rootnodenum
			iterNum = iterNum + 1
			(nodeIndex_array, iterNum) = get_postorder_nodenumbers_above_node(tr, right_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			(nodeIndex_array, iterNum) = get_postorder_nodenumbers_above_node(tr, left_dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			#print(nodeIndex_array)

			return (nodeIndex_array, iterNum)
		else # It must be a nontypical node (i.e. a 2-degree, inside branch)
			# Check number of edges
			num_edges_attached_to_this_node = length(tr.node[rootnodenum].edge)
			if (num_edges_attached_to_this_node != 2)
				txt = ["ERROR in get_postorder_nodenumbers_above_node(): node ", string(rootnodenum), " has ", string(num_edges_attached_to_this_node), " edges attached. It should have exactly 2. Exiting with error."]
				error(join(txt, ""))
			end # END if (num_edges_attached_to_this_node != 2)
			
			one_edge = tr.node[rootnodenum].edge[1]

			# Ancestor node PNnumber for the descendant edge
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			dec_PNnumber = anc_decPNnumbers[2]
			dec_nodeIndex = get_nodeIndex_from_PNnumber(dec_PNnumber, indexNum_table=indexNum_table)
			
			# Ancestor edge not needed
			
			# Then, iterate through single descendant clade
			nodeIndex_array[iterNum] = rootnodenum
			iterNum = iterNum + 1
			(nodeIndex_array, iterNum) = get_postorder_nodenumbers_above_node(tr, dec_nodeIndex, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
			return (nodeIndex_array, iterNum)
		end # END if (typical_or_root_TF == true)
  else
  	# Leaf node:
  	#println(rootnodenum)
  	nodeIndex_array[iterNum] = rootnodenum
  	iterNum = iterNum + 1
  	#print(nodeIndex_array)
  	return (nodeIndex_array, iterNum)
  end

	# Error check
	txt = ["ERROR in get_postorder_nodenumbers_above_node(): You shouldn't be able to get to the last line of this function."]
	error(join(txt, ""))

end # END get_postorder_nodenumbers_above_node


function initialize_edgematrix(tr)
	#ancNodeIndex = collect(repeat([0], 2*(tr.numNodes-tr.numTaxa)))
  #decNodeIndex = collect(repeat([0], 2*(tr.numNodes-tr.numTaxa)))
  ancNodeIndex = collect(repeat([0], length(tr.edge)))
	decNodeIndex = collect(repeat([0], length(tr.edge)))
  edgematrix = hcat(ancNodeIndex,decNodeIndex)
  return(edgematrix)
end


# Returns an edgematrix, with
# 1st column gives ancestral nodeIndex for the edge
# 2nd column gives descendant nodeIndex for the edge
# The edges are in order:
#   pair of edges descending from the root
#   pairs following up the right branch above root
#   pairs following up the left branch above root
#
# Iterate backwards for a postorder traversal in "pruningwise" order 
# (as in APE phylo.reorder "pruningwise")
#

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,orangutan:12);"
great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,(orang1:1,orang2:1):11);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root





rootnodenum = tr.root
#tipNodeIndex_array = collect(repeat([0], length(tr.numTaxa)))
NodeIndex_array = collect(repeat([0], 2*(tr.numNodes-tr.numTaxa)));
iterNum = 0;
indexNum_table = get_nodeIndex_PNnumber(tr);
uppass_edgematrix = initialize_edgematrix(tr)
res = get_pruningwise_postorder_edgematrix(tr, rootnodenum);
uppass_edgematrix = res[1]

get_LR_uppass_edgematrix(tr)
get_LR_downpass_edgematrix(tr)

"""

function get_pruningwise_postorder_edgematrix(tr, rootnodenum, iterNum=1; edgematrix=initialize_edgematrix(tr), indexNum_table=get_nodeIndex_PNnumber(tr))
  if (tr.node[rootnodenum].leaf != true)
  	# * A typical internal node will be attached to 3 edges
  	#   (left descendant, right descendant, ancestor edge)
  	# * A root node will be attached to 2 edges (left, right descendant edges
  	# * A degree-2 (mid-branch) node will have 1 descendant edge
  	
  	if (length(tr.node[rootnodenum].edge) > 3)
  		txt = join(["STOP ERROR in get_pruningwise_postorder_edgematrix(): tr.node[rootnodenum=", string(rootnodenum), "] has more than 3 edges. Probably this is a multifurcating node, which is not allowed."])
  		error(txt)
  	end
  	
  	
  	# Is the current node a typical internal node?
  	typicalTF = length(tr.node[rootnodenum].edge) == 3
  	
  	# Is the current node the root?
  	root_PNnumber = tr.node[tr.root].number  # PhyloNetworks node number of root
  	# rootnodenum = current node being examined
  	current_PNnumber = tr.node[rootnodenum].number
  	rootTF = root_PNnumber == current_PNnumber
		
		# If typical or root, proceed
		typical_or_root_TF = typicalTF || rootTF
  	if (typical_or_root_TF == true)
			# Left descendant edge
			one_edge = tr.node[rootnodenum].edge[1]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			left_anc_PNnumber = anc_decPNnumbers[1]
			left_dec_PNnumber = anc_decPNnumbers[2]
			left_anc_nodeIndex = get_nodeIndex_from_PNnumber(left_anc_PNnumber, indexNum_table=indexNum_table)
			left_dec_nodeIndex = get_nodeIndex_from_PNnumber(left_dec_PNnumber, indexNum_table=indexNum_table)
		
			# Right descendant edge
			one_edge = tr.node[rootnodenum].edge[2]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			right_anc_PNnumber = anc_decPNnumbers[1]
			right_dec_PNnumber = anc_decPNnumbers[2]
			right_anc_nodeIndex = get_nodeIndex_from_PNnumber(right_anc_PNnumber, indexNum_table=indexNum_table)
			right_dec_nodeIndex = get_nodeIndex_from_PNnumber(right_dec_PNnumber, indexNum_table=indexNum_table)
		
			# Then, iterate through left and right clades
			#println(rootnodenum)
			edgematrix[iterNum,1] = right_anc_nodeIndex
			edgematrix[iterNum,2] = right_dec_nodeIndex
			iterNum = iterNum + 1
			edgematrix[iterNum,1] = left_anc_nodeIndex
			edgematrix[iterNum,2] = left_dec_nodeIndex
			iterNum = iterNum + 1
		
			(edgematrix, iterNum) = get_pruningwise_postorder_edgematrix(tr, right_dec_nodeIndex, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
			(edgematrix, iterNum) = get_pruningwise_postorder_edgematrix(tr, left_dec_nodeIndex, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
			#print(nodeIndex_array)
			return (edgematrix, iterNum)
		else
			# Single descendant edge
			one_edge = tr.node[rootnodenum].edge[1]
			anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
			anc_PNnumber = anc_decPNnumbers[1]
			dec_PNnumber = anc_decPNnumbers[2]
			anc_nodeIndex = get_nodeIndex_from_PNnumber(anc_PNnumber, indexNum_table=indexNum_table)
			dec_nodeIndex = get_nodeIndex_from_PNnumber(dec_PNnumber, indexNum_table=indexNum_table)

			# Then, iterate through the single descendant clade
			#println(rootnodenum)
			edgematrix[iterNum,1] = anc_nodeIndex
			edgematrix[iterNum,2] = dec_nodeIndex
			iterNum = iterNum + 1
			
			(edgematrix, iterNum) = get_pruningwise_postorder_edgematrix(tr, dec_nodeIndex, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
			return (edgematrix, iterNum)
		end # END if (typical_or_root_TF == true)
  else
  	#println(rootnodenum)
  	#iterNum = iterNum + 1
  	#nodeIndex_array[iterNum] = rootnodenum
  	#print(nodeIndex_array)
  	return (edgematrix, iterNum)
  end
	
	# Shouldn't get here
	error("get_pruningwise_postorder_edgematrix(): Shouldn't get here!")
end




# Get the node indexes for an uppass from the root to the tips
# Reverse for downpass (postorder)
# Uses get_nodenumbers_above_node() recursively to iterate up tree

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

uppass_edgematrix = get_LR_uppass_edgematrix(tr)
uppass_edgematrix
"""

function get_LR_uppass_edgematrix(tr)
	rootnodenum = tr.root
	iterNum = 0
	edgematrix = initialize_edgematrix(tr)
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_pruningwise_postorder_edgematrix(tr, rootnodenum, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
	uppass_edgematrix = res[1]
	return(uppass_edgematrix)
end


# Get the node indexes for an downpass from the root to the tips
# Reverse of get_LR_uppass_nodeIndexes()
# 

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

uppass_edgematrix = get_LR_uppass_edgematrix(tr)
uppass_edgematrix

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
downpass_edgematrix
"""
function get_LR_downpass_edgematrix(tr)
	rootnodenum = tr.root
	iterNum = 1
	edgematrix = initialize_edgematrix(tr)
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_pruningwise_postorder_edgematrix(tr, rootnodenum, iterNum, edgematrix=edgematrix, indexNum_table=indexNum_table)
	uppass_edgematrix = res[1]
	
	# Reverse
	numrows = size(uppass_edgematrix)[1]
	reverse_rownums = seq(numrows,1,-1)
	downpass_edgematrix = uppass_edgematrix[reverse_rownums,:]
	return(downpass_edgematrix)
end



# Preorder traversal (I think)
"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,(orang1:1,orang2:1):11);"
tr = readTopology(great_ape_newick_string)

# Preorder traversal
uppass_nodeIndexes = get_LR_uppass_nodeIndexes(tr)

# Reverse preorder traversal (not the same as postorder!)
downpass_nodeIndexes = get_LR_downpass_nodeIndexes(tr)
"""
function get_LR_uppass_nodeIndexes(tr)
	rootnodenum = tr.root
	iterNum = 1
	nodeIndex_array = collect(repeat([0], tr.numNodes))
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
	uppass_nodeIndexes = res[1]
	return(uppass_nodeIndexes)
end

# Reverse-Preorder traversal (I think)
"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,orangutan:12);"
#great_ape_newick_string = "(((human:6,(chimp1:0.5,bonobo:0.5):5.5):1,gorilla:7):5,(orang1:1,orang2:1):11);"
tr = readTopology(great_ape_newick_string)

# Preorder traversal
uppass_nodeIndexes = get_LR_uppass_nodeIndexes(tr)

# Reverse preorder traversal (not the same as postorder!)
downpass_nodeIndexes = get_LR_downpass_nodeIndexes(tr)
"""
function get_LR_downpass_nodeIndexes(tr)
	rootnodenum = tr.root
	iterNum = 1
	nodeIndex_array = collect(repeat([0], tr.numNodes))
	indexNum_table = get_nodeIndex_PNnumber(tr)

	res = get_nodenumbers_above_node(tr, rootnodenum, nodeIndex_array, iterNum, indexNum_table=indexNum_table)
	downpass_nodeIndexes = reverse(res[1])
	return(downpass_nodeIndexes)
end


# Get Rnodenums
# They seem to be just 
# sort(tip_PNnumbers) + sort(abs(internal_PNnumbers))

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
df

# Get the edges, get the nodeIndexes corresponding to these
edge = tr.edge
edge

one_edge = edge[1]

edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
edge_df

uppass_edgematrix = get_LR_uppass_edgematrix(tr)
uppass_edgematrix

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
downpass_edgematrix


# Get the R node numbers, append to indexNum_table
Rnodenums_in_indexNum_table_order = get_Rnodenums(tr, indexNum_table)
indexNum_table2 = hcat(indexNum_table, Rnodenums_in_indexNum_table_order)
indexNum_table2

indexNum_table3 = indexNum_table2[sortperm(indexNum_table2[:,3]),:]
indexNum_table3

"""


function get_Rnodenums(tr, indexNum_table)
	numnodes = length(tr.node)
	Rnodenums = collect(1:numnodes)
	tipsTF = indexNum_table[:,2] .> 0
	internalTF = tipsTF .== false
	
	numtips = sum(tipsTF)
	numinternal = sum(internalTF)
	
	
	PNnumbers_in_Rnodenums_order = collect(repeat([0], numnodes))
	PNnumbers_in_Rnodenums_order[1:numtips] = sort(indexNum_table[:,2][tipsTF])
	
	PNnumbers_in_Rnodenums_order[(numtips+1):(numtips+numinternal)] = reverse(sort(indexNum_table[:,2][internalTF]))
	
	
	tmpmat = hcat(Rnodenums, PNnumbers_in_Rnodenums_order)
	
	# Match to indexNum_table
	indices = collect(1:numnodes)
	Rnodenums_in_indexNum_table_order = collect(repeat([0], numnodes))
	for i in 1:numnodes
		current_PNnumber = indexNum_table[i,2] 
		TF = tmpmat[:,2] .== current_PNnumber
		rownum_in_tmpmat = indices[TF][1]
		Rnodenums_in_indexNum_table_order[i] = tmpmat[rownum_in_tmpmat,1]
	end
	return(Rnodenums_in_indexNum_table_order)
end





"""
using DataFrames
using PhyloNetworks
using PhyloPlots
include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

"""
function get_nodeIndex_PNnumber(tr)
	# Get a numNode x 2 table
	# Index, then ".number"
	numnodes = length(tr.node)
	indexNum_table = Array{Int}(undef, numnodes, 2)
	for i in 1:numnodes
		indexNum_table[i,1] = i
		indexNum_table[i,2] = tr.node[i].number
	end
	return(indexNum_table)
end



# Go from a PhyloNetwork Node Number (PNnumber) to the node index 
# (i.e., the index of that node in the list of nodes)
function get_nodeIndex_from_PNnumber(PNnumber; indexNum_table)
	TF01 = indexNum_table[:,2] .== PNnumber
	# Adding the [1] returns a scalar
	nodeIndex = indexNum_table[:,1][TF01][1]
	return nodeIndex
end



"""
using DataFrames
using PhyloNetworks

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
trdf


using DataFrames
using PhyloNetworks

# For Nick's editing (ignore)
include("/GitHub/BioGeoJulia.jl/notes/TreePassO.jl")

#######################################################
# Typical bifurcating (binary) tree
#######################################################
great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreePassO.get_postorder_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table)


#######################################################
# Tree with a 2-degree node inside a branch
#######################################################
include("/GitHub/BioGeoJulia.jl/notes/TreePassO.jl")

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr

nodeIndex_array = collect(repeat([0], tr.numNodes))

iterNum = 1

indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

TreePassO.get_postorder_nodenumbers_above_node(tr, tr.root, nodeIndex_array, iterNum, indexNum_table=indexNum_table )

"""
# Return a DataFrame with the edge numbers
function prt(tr, rootnodenum=tr.root, get_taxa_by_node=true)
	#using DataFrames
	numnodes = length(tr.node)
	# number of digits for internal node numbers
	numdigits = length(digits(numnodes))
	if (numdigits < 2)
		numdigits = 2
	end
	
	# Initialize the dataframe
	trdf = DataFrames.DataFrame(nodeIndex=collect(1:numnodes), PNnumber=collect(repeat([0],numnodes)))
	
	# Fill in the PNnumber node numbers
	indexNum_table = get_nodeIndex_PNnumber(tr)
	trdf[!, :PNnumber] = indexNum_table[:,2]
	trdf
	
	# Add the R node numbers
	Rnodenums = get_Rnodenums(tr, indexNum_table)
	trdf[!,:Rnodenums] = Rnodenums
	trdf
	
	
	# Add the node ages
	node_age = get_node_ages(tr)
	trdf[!,:node_age] = node_age
	
	# Add branch lengths
	brlen = collect(repeat([0.0], numnodes))
	ancNodeIndex = collect(repeat([0], numnodes))
	leftNodeIndex = collect(repeat([0], numnodes))
	rightNodeIndex = collect(repeat([0], numnodes))
	nodeName = collect(repeat([""], numnodes))
	nodeType = collect(repeat([""], numnodes))
	
	edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
	edge_df_rownums = collect(1:Rnrow(edge_df))
	for i in 1:Rnrow(trdf)
		nodeIndex = trdf[i,:nodeIndex]
		TF = edge_df[:,:edge_decNodeIndex] .== nodeIndex
		if (sum(TF) > 0)
			edge_df_rownum = edge_df_rownums[TF][1]
			brlen[i] = edge_df[:,:edge_length][edge_df_rownum]
			ancNodeIndex[i] = edge_df[:,:edge_ancNodeIndex][edge_df_rownum]
		else
			brlen[i] = -999.0
			ancNodeIndex[i] = -999
		end
		
		# Left and right descendant nodeIndexes
		TF_edgerows_descend_from_decNodeIndex = edge_df[:,:edge_ancNodeIndex] .== nodeIndex
		
		# Direct-ancestor node
		if (sum(TF_edgerows_descend_from_decNodeIndex) == 1)
			leftNodeIndex[i] = edge_df[TF_edgerows_descend_from_decNodeIndex,:edge_decNodeIndex][1]
			rightNodeIndex[i] = -99999 # No right descendant, as it's a direct ancestor
			nodeType[i] = "directAnc"  # da = direct ancestor node

			# Node names
			tmp_nodeName = tr.node[nodeIndex].name
			if (tmp_nodeName == "")
				internal_nodeIndex_as_string = lpad(string(nodeIndex), numdigits, "0")
				nodeName[i] = join(["da", internal_nodeIndex_as_string], "")
			else
				nodeName[i] = tmp_nodeName
			end
		end # END if (sum(TF_edgerows_descend_from_decNodeIndex) == 1)
		
		# Bifurcating node
		if (sum(TF_edgerows_descend_from_decNodeIndex) == 2)
			# Internal node
			leftNodeIndex[i] = edge_df[TF_edgerows_descend_from_decNodeIndex,:edge_decNodeIndex][1]
			rightNodeIndex[i] = edge_df[TF_edgerows_descend_from_decNodeIndex,:edge_decNodeIndex][2]
			nodeType[i] = "internal"
			if (nodeIndex == tr.root)
				nodeType[i] = "root"
			end # END if (nodeIndex == tr.root)
			
			# Node names
			tmp_nodeName = tr.node[nodeIndex].name
			if (tmp_nodeName == "")
				internal_nodeIndex_as_string = lpad(string(nodeIndex), numdigits, "0")
				nodeName[i] = join(["in", internal_nodeIndex_as_string], "")
			else
				nodeName[i] = tmp_nodeName
			end
		else
			# Tip node
			leftNodeIndex[i] = -999
			rightNodeIndex[i] = -999
			nodeType[i] = "tip"
			nodeName[i] = tr.node[nodeIndex].name
		end # END if (sum(TF_edgerows_descend_from_decNodeIndex) == 2)
	end # END for i in 1:Rnrow(trdf)
	
	# Add the fields
	trdf[!,:brlen] = brlen
	trdf[!,:ancNodeIndex] = ancNodeIndex
	trdf[!,:leftNodeIndex] = leftNodeIndex
	trdf[!,:rightNodeIndex] = rightNodeIndex
	trdf[!,:nodeName] = nodeName
	trdf[!,:nodeType] = nodeType
	
	if get_taxa_by_node == true
		downpass_edgematrix = get_LR_downpass_edgematrix(tr)
		taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=downpass_edgematrix)
		trdf[!,:taxa] = taxa
	end
	
	return(trdf)
end



# Get tipnames descending from each node
# The list is comma-delimited and alphabetical
"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=get_LR_downpass_edgematrix(tr))


#######################################################
# Tree with a 2-degree node inside a branch
#######################################################
include("/GitHub/BioGeoJulia.jl/notes/TreePassO.jl")

great_ape_newick_string = "((human:1.0,(chimp:0.5):0.5):1.0,gorilla:2.0);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root

# Get trdf WITHOUT "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, false)

# Not all columns print now that trdf is big, let's print just the left and right
# Get trdf WITH "taxa" field (species descending from each node)
trdf = prt(tr, rootnodenum, true);
headLR(trdf)

downpass_edgematrix = get_LR_downpass_edgematrix(tr)
taxa = get_taxa_descending_from_each_node(tr, trdf, downpass_edgematrix=get_LR_downpass_edgematrix(tr))
"""

function get_taxa_descending_from_each_node(tr, trdf; downpass_edgematrix=get_LR_downpass_edgematrix(tr))
	# Get sizes
	numnodes = length(tr.node)
	node_indices_in_downpass_edgematrix_tmp = flat2(downpass_edgematrix)
	node_indices_in_downpass_edgematrix = node_indices_in_downpass_edgematrix_tmp[node_indices_in_downpass_edgematrix_tmp .!== 0]
	numIndices = length(sort(unique(node_indices_in_downpass_edgematrix)))
	
	# Remove 0s from edge matrix, as they constitute "fake" right branches
	
	
	
	# Error check
	if (numnodes != numIndices)
		txt = ["ERROR in get_taxa_descending_from_each_node(): the number of nodes in tr (", string(numnodes), ") does not match the number of rows in downpass_edgematrix (", string(numIndices), ")."]
		error(join(txt, ""))
	end
	
	taxa = collect(repeat([""], numnodes))
	
	# Assumes a binary tree, and a nodeIndexes df in downpass order, from
	# downpass_edgematrix = get_LR_downpass_edgematrix(tr)
	
	# Step through the downpass_edgematrix, in pairs
	edgematrix_rows_to_visit = collect(1:2:Rnrow(downpass_edgematrix))
	
	numrows_in_downpass_edgematrix = Rnrow(downpass_edgematrix)
	i = 1
	while i <= numrows_in_downpass_edgematrix
		j = i+1
		nodeIndex_ancL = downpass_edgematrix[i,1]
		nodeIndex_ancR = downpass_edgematrix[j,1]
		
		# If direct ancestor node
		if (nodeIndex_ancL != nodeIndex_ancR)
			nodeIndex_left = downpass_edgematrix[i,2]
			nodeIndex_right = -99999
			nodeIndex_anc = nodeIndex_ancL
			
			# Get the taxa
			tmp_taxa = [""]
			if (tr.node[nodeIndex_left].leaf == true)
				taxa[nodeIndex_left] = tr.node[nodeIndex_left].name
				tmp_taxa[1] = tr.node[nodeIndex_left].name
			else
				tmp_taxa[1] = taxa[nodeIndex_left]
			end

			# Put the taxa in the ancestral node
			ancNodeIndex_of_Left = nodeIndex_ancL
			ancNodeIndex = ancNodeIndex_of_Left
			taxa_unordered = flat2(split.(tmp_taxa, ","))
			taxa_ordered = sort(taxa_unordered)
			taxa_at_ancNode = join(taxa_ordered, ",")
			taxa[ancNodeIndex] = taxa_at_ancNode
		
			# Advance i by just 1 (because this was a node with a single descendant)
			i = i+1
		end # END if (nodeIndex_ancL != nodeIndex_ancR) # (end direct ancestor node)

		# Bifurcating node, 2 adjacent edges
		if (nodeIndex_ancL == nodeIndex_ancR)
			nodeIndex_anc = nodeIndex_ancL
			nodeIndex_left = downpass_edgematrix[i,2]
			nodeIndex_right = downpass_edgematrix[j,2]

			tmp_taxa = collect(repeat([""], 2))
			if (tr.node[nodeIndex_left].leaf == true)
				taxa[nodeIndex_left] = tr.node[nodeIndex_left].name
				tmp_taxa[1] = tr.node[nodeIndex_left].name
			else
				tmp_taxa[1] = taxa[nodeIndex_left]
			end
			if (tr.node[nodeIndex_right].leaf == true)
				taxa[nodeIndex_right] = tr.node[nodeIndex_right].name
				tmp_taxa[2] = tr.node[nodeIndex_right].name
			else
				tmp_taxa[2] = taxa[nodeIndex_right]
			end
		
			# Join the two lists and put in the ancestral node
			ancNodeIndex_of_Left = nodeIndex_ancL
			ancNodeIndex_of_Right = nodeIndex_ancR
			if (ancNodeIndex_of_Left != ancNodeIndex_of_Right)
				txt = ["ERROR in get_taxa_descending_from_each_node(): ancNodeIndex_of_Left must match ancNodeIndex_of_Right in trdf, but doesn't.\nAncestor of nodeIndex_left ", string(nodeIndex_left), " is ", string(ancNodeIndex_of_Left), ", but\nancestor of nodeIndex_right ", string(nodeIndex_right), " is ", string(ancNodeIndex_of_Right), "\n"]
				msg = join(txt, "")
				println(msg)
				error(msg)
			end
			# No error, so ancNodeIndex was found:
			ancNodeIndex = ancNodeIndex_of_Left
			taxa_unordered = flat2(split.(tmp_taxa, ","))
			taxa_ordered = sort(taxa_unordered)
			taxa_at_ancNode = join(taxa_ordered, ",")
			taxa[ancNodeIndex] = taxa_at_ancNode
			
			# Advance i by 2 (because this was a node with 2 descendants)
			i = i+2
		end # END if (nodeIndex_ancL == nodeIndex_ancR) # (bifurcating node)
	end # END while
	
	
	return(taxa)
end









# Is the node a tip?
function isTip_TF(tr, nodeIndex)
	# Returns true if it's a tip (a "leaf"), false if not
	return(tr.node[nodeIndex].leaf)
end

# Get the nodeIndex ancestral to 2 nodes



# Elaborate function to MAKE SURE we are getting the ancestor & descendant PNnumbers correct

"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

# Get the number (index I think!) of the root node, make a DataFrame of the 
# nodeIndex and PNnumber
rootnodenum = tr.root
trdf = prt(tr, rootnodenum)
df

# Get the edges, get the nodeIndexes corresponding to these
edge = tr.edge
edge

one_edge = edge[1]

# Every edge has only 2 nodes. 
for i in 1:length(edge)
	one_edge = edge[i]
	anc_decPNnumbers = get_NodeIndexes_from_edge(one_edge)
	print(anc_decPNnumbers)
end
"""

function get_NodeIndexes_from_edge(one_edge)
	# Make sure there are only 2 nodes for each edge.
	numnodec_on_this_edge = length(one_edge.node)
	if ( numnodec_on_this_edge != 2 )
		array_of_strings = ["\nError in get_NodeIndexes_from_edge(): this function assumes each edge has only 2 nodes, ancestor & descendant.\nInstead, this edge gives length(edge.node)=", string(numnodec_on_this_edge), ".\nPrinting edge to screen...\n"]
		txt = join(array_of_strings, "")
		println(txt)
		print(edge)
		error(txt)
	end
	
	# Initialize
	tmp_nodeIndices = collect(repeat([0], numnodec_on_this_edge))
	tmp_nodePNnumbers = collect(repeat([0], numnodec_on_this_edge))
	ancPNnumber = 0
	decPNnumber = 0
	
	# Record the PN (PhyloNetworks) node numbers
	tmp_nodePNnumbers[1] = one_edge.node[1].number
	tmp_nodePNnumbers[2] = one_edge.node[2].number
	tmp_nodeIndices_on_edge = collect(1:2)
	
	# Positive node numbers are tips; negative node numbers are internal

	# Declare error if both node numbers are the same
	PNnumber_eq_0_TF = tmp_nodePNnumbers .== 0
	if (sum(PNnumber_eq_0_TF) > 0)
		error("Error in get_NodeIndexes_from_edge(): both node numbers attached to this edge are the same. This should be impossible.")
	end

	
	# Declare error if both node numbers are positive (both are tips)
	positiveTF = tmp_nodePNnumbers .> 0
	if (sum(positiveTF) == 2)
		error("Error in get_NodeIndexes_from_edge(): both node numbers attached to this edge are tips (PhyloNetworks node number > 0). This should be impossible.")
	end
	
	# If one PNnumber is positive and one negative, then you have a descendant tip, and ancestor internal node
	if (sum(positiveTF) == 1)
		anc_nodeIndex_on_edge = tmp_nodeIndices_on_edge[positiveTF .== false]
		ancPNnumber = tmp_nodePNnumbers[anc_nodeIndex_on_edge] # internal nodenum
		dec_nodeIndex_on_edge = tmp_nodeIndices_on_edge[positiveTF .== true]
		decPNnumber = tmp_nodePNnumbers[dec_nodeIndex_on_edge]
	end
	
	# If both are negative, then you have 2 internal nodes. The one closest to 0 is the root
	if (sum(positiveTF) == 0)
		min_absval = -1 * minimum(abs.(tmp_nodePNnumbers))
		matches_min_TF = tmp_nodePNnumbers .== min_absval
		anc_nodeIndex_on_edge = tmp_nodeIndices_on_edge[matches_min_TF .== true]
		dec_nodeIndex_on_edge = tmp_nodeIndices_on_edge[matches_min_TF .== false]
		ancPNnumber = tmp_nodePNnumbers[anc_nodeIndex_on_edge] # ancestral internal nodenum
		decPNnumber = tmp_nodePNnumbers[dec_nodeIndex_on_edge] # ancestral internal nodenum
	end
	
	anc_decPNnumbers = [ancPNnumber[1], decPNnumber[1]]
	return(anc_decPNnumbers)
end


"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr
tr.root

# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
indexNum_table

edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
"""
function get_NodeIndex_df_by_tree_edges(tr; indexNum_table)
	# Get the ancestral and descendant node numbers for each edge
	num_edges = length(tr.edge)
	edgeIndex = collect(1:num_edges)
	edge_ancNodeIndex = collect(repeat([0], num_edges))
	edge_decNodeIndex = collect(repeat([0], num_edges))
	edge_ancPNnumber = collect(repeat([0], num_edges))
	edge_decPNnumber = collect(repeat([0], num_edges))
	edge_length = collect(repeat([0.0], num_edges))

	# Initialize the dataframe
	edge_df = DataFrames.DataFrame(edgeIndex=edgeIndex, edge_ancNodeIndex=edge_ancNodeIndex, edge_decNodeIndex=edge_decNodeIndex, edge_ancPNnumber=edge_ancPNnumber, edge_decPNnumber=edge_decPNnumber, edge_length=edge_length)
	edge_df

	# Every edge has only 2 nodes. 
	for i in 1:length(tr.edge)
		one_edge = tr.edge[i]
		anc_dec_PNnumbers = get_NodeIndexes_from_edge(one_edge)
	
		# Use a named keyword argument, indexNum_table, to provide the translation from PNnumber to node index
		anc_dec_nodeIndices = get_nodeIndex_from_PNnumber.(anc_dec_PNnumbers, indexNum_table=indexNum_table)
		edge_df[i, :edge_ancNodeIndex] = anc_dec_nodeIndices[1]
		edge_df[i, :edge_decNodeIndex] = anc_dec_nodeIndices[2]
		edge_df[i, :edge_ancPNnumber] = anc_dec_PNnumbers[1]
		edge_df[i, :edge_decPNnumber] = anc_dec_PNnumbers[2]
		edge_df[i, :edge_length] = one_edge.length
	end

	return(edge_df)
end


"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)
"""
function get_node_heights(tr)
	indexNum_table = get_nodeIndex_PNnumber(tr)
	uppass_nodeIndexes = get_LR_uppass_nodeIndexes(tr)
	edge_df = get_NodeIndex_df_by_tree_edges(tr, indexNum_table=indexNum_table)
	cumulative_height_at_each_node = collect(repeat([0.0], length(tr.node)))
	
	print("\nuppass_nodeIndexes:\n")
	print(uppass_nodeIndexes)
	print("\n")
	print("\nedge_df:\n")
	print(edge_df)
	
	# Iterate up through the nodes from the root
	for i in 1:length(uppass_nodeIndexes)
		print("\n")
		print("\n")
		print(i)
		
		nodeIndex = uppass_nodeIndexes[i]
		if (nodeIndex == tr.root)
			cumulative_height_at_each_node[nodeIndex] = 0.0
		else
			print("\nedge_df[:,:edge_decNodeIndex]:\n")
			print(edge_df[:,:edge_decNodeIndex])
			print("\nnodeIndex:\n")
			print(nodeIndex)
			ancTF = edge_df[:,:edge_decNodeIndex] .== nodeIndex
			print("\nancTF:\n")
			print(ancTF)
			anc_nodeIndex = edge_df[:,:edge_ancNodeIndex][ancTF][1]
			# Get previous age
			previous_height_above_root = cumulative_height_at_each_node[anc_nodeIndex]
			length_of_this_branch = edge_df[:,:edge_length][ancTF][1]
			current_height = length_of_this_branch + previous_height_above_root
			cumulative_height_at_each_node[nodeIndex] = current_height
		end
	end
	return(cumulative_height_at_each_node)
end




"""
using DataFrames
using PhyloNetworks
using PhyloPlots

include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)
"""
function get_node_ages(tr)
	cumulative_height_at_each_node = get_node_heights(tr)
	tree_height = maximum(cumulative_height_at_each_node)
	node_age = tree_height .- cumulative_height_at_each_node
	return node_age
end
















#######################################################
# Parallel operations on binary trees
#######################################################
# Threaded downpass that spawns new processes when the 2 nodes above are done.

# Solver options structure
# "mutable" means you can change the values referred to by the keys
mutable struct SolverOpt
	solver::Any
	save_everystep::Bool
	abstol::Float64
	reltol::Float64
end

function construct_SolverOpt()
	solver = Tsit5()
	save_everystep = false
	abstol = 1e-9
	reltol = 1e-9
	solver_options = SolverOpt(solver, save_everystep, abstol, reltol)
	return solver_options
end


# Results structure
struct Res
	# The states can be 
	# "not_ready" (value at branch top NOT available)
	# "ready_for_nodeOp" (values at branches above are ready)
	# "ready_for_branchOp" (value at branch top available)
	# "calculating" (value at branch bottom being calculated)
	# "done" (value at branch bottom available)
	node_state::Array{String}
	node_Lparent_state::Array{String}
	node_Rparent_state::Array{String}
	
	# Tree structure
	root_nodeIndex::Int64
	numNodes::Int64
	uppass_edgematrix::Array{Int64}
	sumLikes_at_node_at_branchTop::Array{Float64}
	lnL_at_node_at_branchTop::Array{Float64}
	lq_at_branchBot::Array{Float64}
	like_at_branchBot::Array{Float64}

	thread_for_each_nodeOp::Array{Int64}
	thread_for_each_branchOp::Array{Int64}
	
	# Calculation stats
	calc_spawn_start::Array{DateTime}
	calc_start_time::Array{DateTime}
	calc_end_time::Array{DateTime}
	calc_duration::Array{Float64}
	calctime_iterations::Array{Float64}
	Es_at_each_nodeIndex_branchTop::Array{Array{Float64,1},1}
	Es_at_each_nodeIndex_branchBot::Array{Array{Float64,1},1}
	fakeX0s_at_each_nodeIndex_branchBot::Array{Array{Float64,1},1}
	likes_at_each_nodeIndex_branchTop::Array{Array{Float64,1},1}
	normlikes_at_each_nodeIndex_branchTop::Array{Array{Float64,1},1}
	likes_at_each_nodeIndex_branchBot::Array{Array{Float64,1},1}
	normlikes_at_each_nodeIndex_branchBot::Array{Array{Float64,1},1}
end

# Construct a default, simple results structure
# (likes_at_each_nodeIndex_branchTop is an array)
function construct_Res_old()
	n = 1 # number of states
	node_state = ["ready_for_branchOp", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready"]
	node_Lparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	node_Rparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	root_nodeIndex = 7
	numNodes = 7
	uppass_edgematrix = [7 6; 7 5; 5 4; 5 3; 3 2; 3 1]
	likes_at_each_nodeIndex_branchTop = collect(repeat([1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0], n))
	likes_at_each_nodeIndex_branchBot = collect(repeat([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], n))
	sumLikes_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lnL_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lq_at_branchBot = collect(repeat([0.0], numNodes))
	like_at_branchBot = collect(repeat([0.0], numNodes))

	thread_for_each_nodeOp = collect(repeat([0], numNodes))
	thread_for_each_branchOp = collect(repeat([0], numNodes))

	calc_spawn_start = collect(repeat([Dates.now()], numNodes))
	calc_start_time = collect(repeat([Dates.now()], numNodes))
	calc_end_time = collect(repeat([Dates.now()], numNodes))
	calc_duration = collect(repeat([0.0], numNodes))

	calctime_iterations = [0.0, 0.0]

	res = Res(node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, likes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, sumLikes_at_node_at_branchTop, lnL_at_node_at_branchTop, lq_at_branchBot, like_at_branchBot, thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations)
	return res
end


# Construct a default, simple results structure
# (likes_at_each_nodeIndex_branchTop is an array of arrays)
function construct_Res()
	n = 1 # number of states
	numNodes = 7  # number of nodes
	node_state = ["ready_for_branchOp", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready", "ready_for_branchOp", "not_ready"]
	node_Lparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	node_Rparent_state = ["NA", "NA", "not_ready", "NA", "not_ready", "NA", "not_ready"]
	root_nodeIndex = 7
	uppass_edgematrix = [7 6; 7 5; 5 4; 5 3; 3 2; 3 1]
	likes_OneNode = collect(repeat([0.0], n))
# 	likes_at_each_nodeIndex_branchTop = repeat([likes_OneNode], numNodes)
# 	likes_at_each_nodeIndex_branchBot = repeat([likes_OneNode], numNodes)
# 	normlikes_at_each_nodeIndex_branchTop = repeat([likes_OneNode], numNodes)
# 	normlikes_at_each_nodeIndex_branchBot = repeat([likes_OneNode], numNodes)
	Es_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	Es_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	fakeX0s_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	likes_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	likes_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	normlikes_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	normlikes_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	
	default_likes_at_each_nodeIndex_branchTop = [1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0]
	for i in 1:length(likes_at_each_nodeIndex_branchTop)
		likes_at_each_nodeIndex_branchTop[i][1] = default_likes_at_each_nodeIndex_branchTop[i]
	end
	#typeof(likes_at_each_nodeIndex_branchTop)
	sumLikes_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lnL_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lq_at_branchBot = collect(repeat([0.0], numNodes))
	like_at_branchBot = collect(repeat([0.0], numNodes))


	thread_for_each_nodeOp = collect(repeat([0], 7))
	thread_for_each_branchOp = collect(repeat([0], 7))

	calc_spawn_start = collect(repeat([Dates.now()], 7))
	calc_start_time = collect(repeat([Dates.now()], 7))
	calc_end_time = collect(repeat([Dates.now()], 7))
	calc_duration = collect(repeat([0.0], 7))

	calctime_iterations = [0.0, 0.0]

	res = Res(node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, sumLikes_at_node_at_branchTop, lnL_at_node_at_branchTop, lq_at_branchBot, like_at_branchBot, thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations, Es_at_each_nodeIndex_branchTop, Es_at_each_nodeIndex_branchBot, fakeX0s_at_each_nodeIndex_branchBot, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, normlikes_at_each_nodeIndex_branchBot)
	return res
end





"""
# Load a simple tree, see a simple list of nodeIndexes
using DataFrames
using PhyloNetworks
include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr.root
# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)

"""
function construct_Res(tr::HybridNetwork)
	root_nodeIndex = tr.root
	numNodes = tr.numNodes
	uppass_edgematrix = get_LR_uppass_edgematrix(tr)
	
	# Give tip nodeIndexes their nodeNodex as the "likelihood"
	indexNum_table = get_nodeIndex_PNnumber(tr)
	tipsTF = indexNum_table[:,2] .> 0
	tipLikes = indexNum_table[tipsTF,2] * 1.0

	# Set up an array of length nstates (n), to hold the likelihoods for each node
	n = 1
	blank_states = collect(repeat([0.0], n))
# 	likes_at_each_nodeIndex_branchTop = collect(repeat([blank_states], numNodes))
# 	likes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
# 	normlikes_at_each_nodeIndex_branchTop = collect(repeat([blank_states], numNodes))
# 	normlikes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
	Es_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	Es_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	fakeX0s_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)

	likes_at_each_nodeIndex_branchTop = collect(repeat([collect(repeat([0.0], n))], numNodes))
	likes_at_each_nodeIndex_branchBot = collect(repeat([collect(repeat([0.0], n))], numNodes))
	normlikes_at_each_nodeIndex_branchTop = collect(repeat([collect(repeat([0.0], n))], numNodes))
	normlikes_at_each_nodeIndex_branchBot = collect(repeat([collect(repeat([0.0], n))], numNodes))

	# Put in the tip node numbers as the fake likelihoods
	function f(numNodes, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, fakeX0s_at_each_nodeIndex_branchBot, tipsTF)
		j = 0
		for i in 1:numNodes
			if (tipsTF[i] == true)
				j = j+1
				# Transfer from 1D array to 1D array
				likes_at_each_nodeIndex_branchTop[i] = [tipLikes[j]]
				normlikes_at_each_nodeIndex_branchTop[i] = [tipLikes[j] / sum(tipLikes[j])]
				fakeX0s_at_each_nodeIndex_branchBot[i] = [tipLikes[j] / sum(tipLikes[j])]
			end
		end
	end
	# Run function f()
	f(numNodes, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, fakeX0s_at_each_nodeIndex_branchBot, tipsTF)
	
	likes_at_each_nodeIndex_branchTop
	sumLikes_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lnL_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lq_at_branchBot = collect(repeat([0.0], numNodes))
	like_at_branchBot = collect(repeat([0.0], numNodes))
	
	# Fill in the node_states
	node_state = collect(repeat(["not_ready"], numNodes))
	node_state[tipsTF] .= "ready_for_branchOp"
	node_Lparent_state = collect(repeat(["not_ready"], numNodes))
	node_Rparent_state = collect(repeat(["not_ready"], numNodes))
	node_Lparent_state[tipsTF] .= "NA"
	node_Rparent_state[tipsTF] .= "NA"
	
	# Initialize with zeros for the other items
	#likes_at_each_nodeIndex_branchBot = collect(repeat([0.0], numNodes))
	likes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
	thread_for_each_nodeOp = collect(repeat([0], numNodes))
	thread_for_each_branchOp = collect(repeat([0], numNodes))
	
	calc_spawn_start = collect(repeat([Dates.now()], numNodes))
	calc_start_time = collect(repeat([Dates.now()], numNodes))
	calc_end_time = collect(repeat([Dates.now()], numNodes))
	calc_duration = collect(repeat([0.0], numNodes))

	calctime_iterations = [0.0, 0.0]
	number_of_whileLoop_iterations = [0]	

	# Initialize res object
	res = Res(node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, sumLikes_at_node_at_branchTop, lnL_at_node_at_branchTop, lq_at_branchBot, like_at_branchBot, thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations, Es_at_each_nodeIndex_branchTop, Es_at_each_nodeIndex_branchBot, fakeX0s_at_each_nodeIndex_branchBot, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, normlikes_at_each_nodeIndex_branchBot)
	return res
end





"""
# Load a simple tree, see a simple list of nodeIndexes
using DataFrames
using PhyloNetworks
include("/drives/Dropbox/_njm/__julia/julia4Rppl_v1.jl")

great_ape_newick_string = "(((human:6,chimpanzee:6):1,gorilla:7):5,orangutan:12);"
tr = readTopology(great_ape_newick_string)
tr.root
# Get a table with the index numbers of the nodes
indexNum_table = get_nodeIndex_PNnumber(tr)
n=10 # number of states in the state space

"""
function construct_Res(tr::HybridNetwork, n)
	root_nodeIndex = tr.root
	numNodes = tr.numNodes
	uppass_edgematrix = get_LR_uppass_edgematrix(tr)
	
	
	# Set up an array of length nstates (n), to hold the likelihoods for each node
	blank_states = collect(repeat([0.0], n))
# 	likes_at_each_nodeIndex_branchTop = collect(repeat([blank_states], numNodes))
# 	likes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
# 	normlikes_at_each_nodeIndex_branchTop = collect(repeat([blank_states], numNodes))
# 	normlikes_at_each_nodeIndex_branchBot = collect(repeat([blank_states], numNodes))
	Es_at_each_nodeIndex_branchTop = repeat([collect(repeat([0.0], n))], numNodes)
	Es_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	fakeX0s_at_each_nodeIndex_branchBot = repeat([collect(repeat([0.0], n))], numNodes)
	likes_at_each_nodeIndex_branchTop = collect(repeat([collect(repeat([0.0], n))], numNodes))
	likes_at_each_nodeIndex_branchBot = collect(repeat([collect(repeat([0.0], n))], numNodes))
	normlikes_at_each_nodeIndex_branchTop = collect(repeat([collect(repeat([0.0], n))], numNodes))
	normlikes_at_each_nodeIndex_branchBot = collect(repeat([collect(repeat([0.0], n))], numNodes))

	# Give tip nodeIndexes a likelihood of 1 at all states
	indexNum_table = get_nodeIndex_PNnumber(tr)
	tipsTF = indexNum_table[:,2] .> 0
	tipnums = seq(1, length(tipsTF), 1)[tipsTF]
	
# 	put_in_fake_tipLikes = false
# 	if put_in_fake_tipLikes == true
# 		for i in 1:length(tipnums)
# 			tipLikes = collect(repeat([1.0], n))
# 			likes_at_each_nodeIndex_branchTop[tipnums[i]] = tipLikes
# 			normlikes_at_each_nodeIndex_branchTop[i] = [tipLikes[j] / sum(tipLikes[j])]
# 		end
# 	end
	sumLikes_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lnL_at_node_at_branchTop = collect(repeat([0.0], numNodes))
	lq_at_branchBot = collect(repeat([0.0], numNodes))
	like_at_branchBot = collect(repeat([0.0], numNodes))
	
	# Fill in the node_states
	node_state = collect(repeat(["not_ready"], numNodes))
	node_state[tipsTF] .= "ready_for_branchOp"
	node_Lparent_state = collect(repeat(["not_ready"], numNodes))
	node_Rparent_state = collect(repeat(["not_ready"], numNodes))
	node_Lparent_state[tipsTF] .= "NA"
	node_Rparent_state[tipsTF] .= "NA"
	
	# Initialize with zeros for the other items
	thread_for_each_nodeOp = collect(repeat([0], numNodes))
	thread_for_each_branchOp = collect(repeat([0], numNodes))
	
	calc_spawn_start = collect(repeat([Dates.now()], numNodes))
	calc_start_time = collect(repeat([Dates.now()], numNodes))
	calc_end_time = collect(repeat([Dates.now()], numNodes))
	calc_duration = collect(repeat([0.0], numNodes))

	calctime_iterations = [0.0, 0.0]
	number_of_whileLoop_iterations = [0]	

	# Initialize res object
	res = Res(node_state, node_Lparent_state, node_Rparent_state, root_nodeIndex, numNodes, uppass_edgematrix, sumLikes_at_node_at_branchTop, lnL_at_node_at_branchTop, lq_at_branchBot, like_at_branchBot, thread_for_each_nodeOp, thread_for_each_branchOp, calc_spawn_start, calc_start_time, calc_end_time, calc_duration, calctime_iterations, Es_at_each_nodeIndex_branchTop, Es_at_each_nodeIndex_branchBot, fakeX0s_at_each_nodeIndex_branchBot, likes_at_each_nodeIndex_branchTop, normlikes_at_each_nodeIndex_branchTop, likes_at_each_nodeIndex_branchBot, normlikes_at_each_nodeIndex_branchBot)
	return res
end




# Convert this res object to a DataFrame
#function res_to_df(res)
#	
#end

function count_nodes_finished(node_state)
	sum(node_state .== "done")
end

# Average the downpass likelihoods at a node
# (default)
function nodeOp_average_likes(tmp1, tmp2)
	nodeData_at_top = (tmp1 + tmp2)/2
	return(nodeData_at_top)
end

# Use the Cmat to combine the likelihoods
nodeOp_Cmat = (tmpDs; tmp1, tmp2, p_Ds_v5) -> begin
	p = p_Ds_v5
#	hcat(p.p_indices.Carray_ivals, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals, p.params.Cijk_vals)
	
	# Go through each Ci (ancestral state index)

 # Possibly varying parameters
  n = p.n
  mu = p.params.mu_vals
  Cijk_vals = p.params.Cijk_vals
	
	# Indices for the parameters (events in a sparse anagenetic or cladogenetic matrix)
	Carray_ivals = p.p_indices.Carray_ivals
	Carray_jvals = p.p_indices.Carray_jvals
	Carray_kvals = p.p_indices.Carray_kvals
	
# 	print("\n")
# 	print("\n")
# 	print("Running nodeOp_Cmat:\n\n")
	
	# Calculate likelihoods of states just before speciation
  @inbounds for i in 1:n
  	# These are the i's, j's, and k's FOR AN ANCESTOR I
		Ci_sub_i = p.p_TFs.Ci_sub_i[i]
		Cj_sub_i = p.p_TFs.Cj_sub_i[i]
		Ck_sub_i = p.p_TFs.Ck_sub_i[i]
		
		# This is the TFs for an ancestor i - NEEDED FOR FETCHING Cij_vals!!
		Ci_eq_i = p.p_TFs.Ci_eq_i[i]

# 		Calculation of "D" (likelihood of tip data)
# 		for state i=1, multiply
# 		
# 		All the speciation rates where the ancestor = i; length is length(Ci_sub_i)
# 		Cijk_vals[Ci_sub_i]
# 		
# 		All the ancestors (i's) where the ancestor = i; length is length(Ci_sub_i)
# 		Carray_ivals[Ci_sub_i]
# 
# 		All the left descendant states (j's) where the ancestor just before speciation = i; length is sum(Ci_sub_i)
# 		Carray_jvals[Ci_sub_i]
# 		
# 		All the right descendant states (k's) where the ancestor just before speciation = i; length is sum(Ci_sub_i)
# 		Carray_kvals[Ci_sub_i]
# 		
# 		Coming down from left branch, contributing to likelihood of ancestor state i;
# 		resulting length is sum(Ci_sub_i)
# 		tmp1[Carray_jvals[Ci_sub_i]]
# 
# 		Coming down from left branch, contributing to likelihood of ancestor state i;
# 		resulting length is sum(Ci_sub_i)
# 		tmp2[Carray_kvals[Ci_sub_i]]
		
		# Equivalent:
		# inputs.p_Ds_v5.params.Cijk_vals[inputs.p_Ds_v5.p_TFs.Ci_eq_i[i]]
		# 6-element Array{Float64,1}:
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		# 
		# julia> inputs.p_Ds_v5.p_TFs.Ci_sub_i[i]
		# 6-element Array{Int64,1}:
		#  3
		#  3
		#  3
		#  3
		#  3
		#  3
		# 
		# julia> inputs.p_Ds_v5.params.Cijk_vals[ inputs.p_Ds_v5.p_TFs.Ci_sub_i[i]]
		# 6-element Array{Float64,1}:
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333
		#  0.03333333333333333		
		
		# Parameter values for these events with nonzero rates
		tmpDs[i] = sum(Cijk_vals[Ci_eq_i] .* tmp1[Cj_sub_i] .* tmp2[Ck_sub_i])
# 		print(tmpDs[i])
# 		print("\n")
# 		print(Cijk_vals)
# 		print("\n")
# 		print(Ci_sub_i)
# 		print("\n")
# 		print(Cj_sub_i)
# 		print("\n")
# 		print(Ck_sub_i)
# 		print("\n")
# 		print(tmp1)
# 		print("\n")
# 		print(tmp2)
# 		print("\n")
  end
  return(tmpDs)
end


	
	




# Combine likelihoods from above
function nodeOp(current_nodeIndex, res; nodeOp_function=nodeOp_average_likes)
	res.node_state[current_nodeIndex] = "calculating_nodeOp"
	uppass_edgematrix = res.uppass_edgematrix
	
	# Record the thread, for kicks
	tmp_threadID = Threads.threadid()
	res.thread_for_each_nodeOp[current_nodeIndex] = tmp_threadID
	TF = uppass_edgematrix[:,1] .== current_nodeIndex
	if (sum(TF) == 2)
		# Get likelihoods from above (iterates up to tips)
		parent_nodeIndexes = uppass_edgematrix[TF,2]

		# Crucial
# 		tmp1 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
# 		tmp2 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		tmp1 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		tmp2 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]

		# Check that data are actually available
		if (sum(tmp1) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp1) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		if (sum(tmp2) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp2) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		#nodeData_at_top = tmp1 + tmp2
		#nodeData_at_top = (tmp1 + tmp2)/2
		nodeData_at_top = nodeOp_function(tmp1, tmp2)
		
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] = nodeData_at_top
		
		# Check if it's the root node
		if (current_nodeIndex == res.root_nodeIndex)
			res.node_state[current_nodeIndex] = "done"
		else
			res.node_state[current_nodeIndex] = "ready_for_branchOp"
		end
		return()
	elseif (sum(TF) == 0)
	  # If a tip
	  txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): shouldn't be run on a tip node."], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	else
	  txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(TF) should be 0 or 2"], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	end
	txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): shouldn't get here."], "")
	print("\n")
	print(txt)
	print("\n")
	return(error(txt))
end # END function nodeOp(current_nodeIndex, res; nodeOp_function=nodeOp_average_likes)


# Node operation, combining probabilities from above (assumed to be fast)
function nodeOp_ClaSSE_v5!(current_nodeIndex, res; p_Ds_v5)
	res.node_state[current_nodeIndex] = "calculating_nodeOp"
	uppass_edgematrix = res.uppass_edgematrix
	
	# Record the thread, for kicks
	tmp_threadID = Threads.threadid()
	res.thread_for_each_nodeOp[current_nodeIndex] = tmp_threadID
	TF = uppass_edgematrix[:,1] .== current_nodeIndex
	if (sum(TF) == 2)
		# Get likelihoods from above (iterates up to tips)
		parent_nodeIndexes = uppass_edgematrix[TF,2]
		#tmp1 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		#tmp2 = res.likes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		# Crucial
		# The likelihoods were taken out to produce normlikes, stored in "res.lq_at_branchBot"
		tmp1 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[1]]
		tmp2 = res.normlikes_at_each_nodeIndex_branchBot[parent_nodeIndexes[2]]
		combined_branch_lnLs = res.lq_at_branchBot[parent_nodeIndexes[1]] + res.lq_at_branchBot[parent_nodeIndexes[2]]
		
		
		# Check that data are actually available
		if (sum(tmp1) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp1) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		if (sum(tmp2) == 0.0)
			txt = join(["Error in nodeOp(current_nodeIndex=", string(current_nodeIndex), "): sum(tmp2) == 0.0, indicating data at parent nodes actually not available."], "")
			res.node_state[current_nodeIndex] = txt
			print("\n")
			print(txt)
			print("\n")
			return(error(txt))
		end

		#nodeData_at_top = tmp1 + tmp2
		#nodeData_at_top = (tmp1 + tmp2)/2
		#nodeData_at_top = nodeOp_function(tmp1, tmp2)
# 		print("\n\n")
# 		
# 		print("\ncurrent_nodeIndex:\n")
# 		print(current_nodeIndex)

		tmpDs = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] # Placeholder
# 		print("\ntmpDs:\n")
# 		print(tmpDs)
		nodeData_at_top = nodeOp_Cmat(tmpDs, tmp1=tmp1, tmp2=tmp2, p_Ds_v5=p_Ds_v5)
		# Somehow adding .+ 0.0 individualizes the assignment!
		
		sum_likes_at_node = sum(nodeData_at_top)
		#sum_likes_at_node = 1.0
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex] = (nodeData_at_top .+ 0.0)
		res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
# 		print("\n\ncurrent_nodeIndex:")
# 		print(current_nodeIndex)
# 		print("\n")
		res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex] = ((nodeData_at_top .+ 0.0) ./ sum_likes_at_node)
		
		res.sumLikes_at_node_at_branchTop[current_nodeIndex] = sum_likes_at_node + 0.0
		res.lnL_at_node_at_branchTop[current_nodeIndex] = log(sum_likes_at_node) + 0.0
# 		print("\nnodeData_at_top:\n")
# 		print(nodeData_at_top)
# 
# 		print("\nres.likes_at_each_nodeIndex_branchTop[current_nodeIndex]:\n")
# 		print(res.likes_at_each_nodeIndex_branchTop[current_nodeIndex])
# 
# 		print("\nres.likes_at_each_nodeIndex_branchTop:\n")
# 		print(res.likes_at_each_nodeIndex_branchTop)

		
		# Check if it's the root node
		if (current_nodeIndex == res.root_nodeIndex)
			res.node_state[current_nodeIndex] = "done"
		else
			res.node_state[current_nodeIndex] = "ready_for_branchOp"
		end
		return(res)  # 2020-08-04_NJM
	elseif (sum(TF) == 0)
	  # If a tip
	  txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): shouldn't be run on a tip node."], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	else
	  txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): sum(TF) should be 0 or 2"], "")
	  print("\n")
	  print(txt)
	  print("\n")
		return(error(txt))
	end # End if (sum(TF) == 2)
	txt = join(["Error in nodeOp_ClaSSE_v5(current_nodeIndex=", string(current_nodeIndex), "): shouldn't get here."], "")
	print("\n")
	print(txt)
	print("\n")
	return(error(txt))
end # END function nodeOp_ClaSSE_v5!(current_nodeIndex, res; p_Ds_v5)






# Calculate down a branch
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_example(current_nodeIndex, res; num_iterations=10000000)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Example slow operation
	y = countloop(num_iterations, current_nodeIndex)

	nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	nodeData_at_bottom = nodeData_at_top / 2.0

	return(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time)
end




# Calculate Ds down a branch
#
# Modifies branchOp to do Ds calculation down a branch
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp_ClaSSE_Ds_v5(current_nodeIndex, res; u0, tspan, p_Ds_v5, solver_options=solver_options)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	# Example slow operation
	#y = countloop(num_iterations, current_nodeIndex)
# 	print("\n")
# 	print("branchOp_ClaSSE_Ds_v5: d: ")
# 	print(p_Ds_v5.params.Qij_vals[1])
# 	print("branchOp_ClaSSE_Ds_v5: e: ")
# 	print(p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)])
# 	print("\n")
	
	prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
	sol_Ds = solve(prob_Ds_v5, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol)
	
# 	print(sol_Ds[length(sol_Ds)])
# 	print("\n")
	
	#nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
	#nodeData_at_bottom = nodeData_at_top / 2.0
	#nodeData_at_bottom = sol_Ds.u
	
	return(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time)
end








"""
	Set up inputs
"""

function setup_inputs_branchOp_ClaSSE_Ds_v5(u0, tspan, p_Ds_v5; solver="Tsit5()", 
				 save_everystep="false", abstol="1e-9", reltol="1e-9")
	
	prob_str = "prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)"
	solve_str = join(["sol_Ds = solve(prob_Ds_v5, ", solver, ", save_everystep=", save_everystep, ", abstol=", abstol, ", reltol=", reltol, ")"])
	store_str = "nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)]"
	
	# Assemble the NamedTuple of inputs
	inputs = (u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver=solver, save_everystep=save_everystep, abstol=abstol, reltol=reltol, prob_str=prob_str, solve_str=solve_str, store_str=store_str)
	return inputs
end


# Calculate Ds down a branch
#
# Modifies branchOp to do Ds calculation down a branch
#
# This function can read from res, but writing to res is VERY BAD as 
# it created conflicts apparently when there were more @spawns than cores
# Do all the writing to res in the while() loop
function branchOp(current_nodeIndex, res, inputs)
	calc_start_time = Dates.now()
	spawned_nodeIndex = current_nodeIndex
	tmp_threadID = Threads.threadid()
	
	
	# The old practice input was an Int64
	if (typeof(inputs) != Int64)
		u0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
		tspan = inputs.tspan
		p_Ds_v5 = inputs.p_Ds_v5
		#solver = inputs.solver
		#save_everystep = inputs.save_everystep
		#abstol = inputs.abstol
		#reltol = inputs.reltol


		# Example slow operation
		#y = countloop(num_iterations, current_nodeIndex)
		#	prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
		#	sol_Ds = solve(prob_Ds_v5, Tsit5(), save_everystep=true, abstol = 1e-9, reltol = 1e-9)
		#print(u0)
		
		# It DOESN'T really work to set the functions by text;
		# global variable problems, etc.
		# prob_Ds_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Ds_v5, u0, tspan, p_Ds_v5)
		# sol_Ds = solve(prob_Ds_v5, Tsit5(), save_everystep=false, abstol=1e-9, reltol=1e-9)
		# nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)]
		eval(Meta.parse(inputs.prob_str))
		eval(Meta.parse(inputs.solve_str))
		eval(Meta.parse(inputs.store_str))
		
	else
		nodeData_at_top = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
		nodeData_at_bottom = nodeData_at_top / 2.0
		#nodeData_at_bottom = sol_Ds
	end

	
	return(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time)
end




function countloop(num_iterations, current_nodeIndex)
	x = 0.0
	random_number_generator = MersenneTwister(current_nodeIndex);

	for i in 1:num_iterations
	   x = x + (randn(random_number_generator, 1)[1] / num_iterations)
	end
	return(x)
end







"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf, p_Ds_v5, solver_options=construct_SolverOpt(), max_iterations=10^10, return_lnLs=false)
	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()
	
	# Re-set the node states for new downpass
	TF = trdf[:,:nodeType] .== "tip"
	res.node_state[TF] .= "ready_for_branchOp"
	res.node_state[TF .== false] .= "not_ready"
	res.node_Lparent_state[TF] .= "NA"
	res.node_Lparent_state[TF .== false] .= "not_ready"
	res.node_Rparent_state[TF] .= "NA"
	res.node_Rparent_state[TF .== false] .= "not_ready"
		
	# Setup
	current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "done"
				return()
			end
	
			# Retrieve the inputs for the calculation down the branch
			
			# Use the RAW likelihoods (don't normalize)
			u0 = res.likes_at_each_nodeIndex_branchTop[current_nodeIndex]
			#u0 = u0 ./ (sum(u0))
			
			# Use the NORMALIZED (rescaled to sum to 1) likelihoods
			# Doesn't work -- claims an interpolation error for going beyond range
			# branchOp on current_nodeIndex=4ERROR: LoadError: Solution interpolation 
			# cannot extrapolate past the final timepoint. Either solve on a longer 
			# timespan or use the local extrapolation from the integrator interface.
			#u0 = res.normlikes_at_each_nodeIndex_branchTop[current_nodeIndex]
			
			brlen = trdf[current_nodeIndex, :brlen]
			age_branchtop = trdf[current_nodeIndex, :node_age]
			age_branchbot = age_branchtop + brlen
			tspan = (age_branchtop, age_branchbot)
			#p_Ds_v5 = inputs.p_Ds_v5

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			#print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
#			if (parallel_TF == true)
#				push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
#			else
			tmp_results = branchOp_ClaSSE_Ds_v5(current_nodeIndex, res, u0=u0, tspan=tspan, p_Ds_v5=p_Ds_v5, solver_options=solver_options)
			#tmp_results = branchOp(current_nodeIndex, res, num_iterations)
			push!(tasks, tmp_results)			 # Add results to "tasks"
#			end
			push!(tasks_fetched_TF, false) # Add a "false" to tasks_fetched_TF
		end # END for current_nodeIndex in indexes_ready
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			if (tasks_fetched_TF[i] == false)
				#if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
					(tmp_threadID, sol_Ds, spawned_nodeIndex, calc_start_time) = tasks[i]
					nodeData_at_bottom = sol_Ds.u[length(sol_Ds.u)] .+ 0.0
# 					end
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
# 					print("\n\n12345\n\n")
# 					print("res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex]:\n")
# 					print(res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex])
# 					print("\n\nnodeData_at_bottom:\n")
# 					print(nodeData_at_bottom)
# 					print("\n\n12345\n\n")

					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
					res.lq_at_branchBot[spawned_nodeIndex] = log(sum_nodeData_at_bottom)
					res.like_at_branchBot[spawned_nodeIndex] = sum_nodeData_at_bottom

					# Get the ancestor nodeIndex
					uppass_edgematrix = res.uppass_edgematrix
					TF = uppass_edgematrix[:,2] .== spawned_nodeIndex
					parent_nodeIndex = uppass_edgematrix[TF,1][1]

					# Get the left daughter nodeIndex (1st in the uppass_edgematrix)
					edge_rows_TF = uppass_edgematrix[:,1] .== parent_nodeIndex
					left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
					right_nodeIndex = uppass_edgematrix[edge_rows_TF,2][2]

					# Update the state of the parent_node's daughters
					if (spawned_nodeIndex == left_nodeIndex)
						res.node_Lparent_state[parent_nodeIndex] = "ready"
					end
					if (spawned_nodeIndex == right_nodeIndex)
						res.node_Rparent_state[parent_nodeIndex] = "ready"
					end

					# Update the state of the current node
					res.node_state[spawned_nodeIndex] = "done"
				#end
			end
		end
	
		# Update which nodes have had both parents complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = res.node_Rparent_state .== "ready"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"
	
		# Update nodes when the branches above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			# Combine the downpass branch likelihoods
			#nodeOp(current_nodeIndex, res, nodeOp_function=nodeOp_average_likes)
			res = nodeOp_ClaSSE_v5!(current_nodeIndex, res, p_Ds_v5=p_Ds_v5)
			# (updates res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass_nonparallel(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
			error(txt)
		end
		
		# Test for concluding the while loop
		are_we_done && break
	end # END while(are_we_done == false)
	
	# This breaks it for some reason:
	# ERROR: setfield! immutable struct of type Res cannot be changed
	#global res.number_of_whileLoop_iterations = iteration_number

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	Julia_sum_lq = sum(res.lq_at_branchBot[1:(length(res.lq_at_branchBot)-1)])

	# Add the root probabilities

	# Assuming diversitree options:
	# root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
	# i.e., the root state probs are just the root_Ds/sum(root_Ds)
	d_root_orig = res.likes_at_each_nodeIndex_branchTop[length(res.likes_at_each_nodeIndex_branchTop)]
	root_stateprobs = d_root_orig/sum(d_root_orig)
	rootstates_lnL = log(sum(root_stateprobs .* d_root_orig))
	Julia_total_lnLs1 = Julia_sum_lq + rootstates_lnL
	
	if return_lnLs == true
		#txt = paste0(["d=", p_Ds_v5.params.Qij_vals[1], ",	e=", p_Ds_v5.params.Qij_vals[length(p_Ds_v5.params.Qij_vals)], ",	Julia_sum_lq=", round(Julia_sum_lq; digits=3), ", rootstates_lnLB=", round(rootstates_lnL; digits=3), ",	Julia_total_lnLs1B=", Julia_total_lnLs1])
		#print(txt) 
		#print("\n")
		return(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1)
	else
		return(total_calctime_in_sec, iteration_number)
	end
	
	# shouldn't get here
	return NaN
end # END iterative_downpass_nonparallel_ClaSSE_v5!



















"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
"""
function iterative_downpass!(res; max_iterations=10^10, num_iterations=10000000)

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	if (parallel_TF == false)
		txt = "Error in iterative_downpass!(): This function probably requires multiple threads operating to consistently compile. Try starting julia with e.g. 'JULIA_NUM_THREADS=8 julia'."
		error(txt)
	end

	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()
	
	# Setup
	current_nodeIndex = res.root_nodeIndex
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "done"
				return()
			end

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
			push!(tasks_fetched_TF, false)
		end
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			if (tasks_fetched_TF[i] == false)
				if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
					(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
					
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
					res.lq_at_branchBot[spawned_nodeIndex] = log(sum_nodeData_at_bottom)
					res.like_at_branchBot[spawned_nodeIndex] = sum_nodeData_at_bottom

					# Get the ancestor nodeIndex
					uppass_edgematrix = res.uppass_edgematrix
					TF = uppass_edgematrix[:,2] .== spawned_nodeIndex
					parent_nodeIndex = uppass_edgematrix[TF,1][1]

					# Get the left daughter nodeIndex (1st in the uppass_edgematrix)
					edge_rows_TF = uppass_edgematrix[:,1] .== parent_nodeIndex
					left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
					right_nodeIndex = uppass_edgematrix[edge_rows_TF,2][2]

					# Update the state of the parent_node's daughters
					if (spawned_nodeIndex == left_nodeIndex)
						res.node_Lparent_state[parent_nodeIndex] = "ready"
					end
					if (spawned_nodeIndex == right_nodeIndex)
						res.node_Rparent_state[parent_nodeIndex] = "ready"
					end

					# Update the state of the current node
					res.node_state[spawned_nodeIndex] = "done"
				end
			end
		end
	
		# Update which nodes have had both parents complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = res.node_Rparent_state .== "ready"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"
	
		# Update nodes when the branches above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			nodeOp(current_nodeIndex, res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
			error(txt)
		end
		
		# Test for concluding the while loop
		are_we_done && break
	end
	
	# This breaks it for some reason:
	# ERROR: setfield! immutable struct of type Res cannot be changed
	#global res.number_of_whileLoop_iterations = iteration_number

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	return(total_calctime_in_sec, iteration_number)
end # END iterative_downpass!




"""
Iterate through the "res" object many times to complete the downpass, spawning jobs along the way
Non-parallel version (no istaskdone, etc.)
"""
function iterative_downpass_nonparallel!(res; max_iterations=10^10, num_iterations=10000000)
	diagnostics = collect(repeat([Dates.now()], 3))
	diagnostics[1] = Dates.now()
	
	# Setup
	current_nodeIndex = res.root_nodeIndex

	# Check number of threads
	numthreads = Threads.nthreads()
	parallel_TF = numthreads > 1
	tasks = Any[]
	tasks_fetched_TF = Bool[]
	are_we_done = false

	iteration_number = 0
	while(are_we_done == false)
		iteration_number = iteration_number+1
		# As long as all the nodes are not done,
		# check for "ready" nodes
		# When they finish, change to "done"
		indexes_ready = findall(res.node_state .== "ready_for_branchOp")
		for current_nodeIndex in indexes_ready
			# Before spawning, do some checks
			res.node_state[current_nodeIndex] = "calculating_branchOp"
			# Check for root; no calculation on root branch for now
			if current_nodeIndex == res.root_nodeIndex
				res.node_state[current_nodeIndex] = "done"
				return()
			end

			# Spawn a branch operation, and a true-false of whether they are fetched
			res.calc_spawn_start[current_nodeIndex] = Dates.now()
			print(join(["\nbranchOp on current_nodeIndex=", string(current_nodeIndex)], ""))
# 			if (parallel_TF == true)
# 				push!(tasks, @spawn branchOp(current_nodeIndex, res, num_iterations=num_iterations))
# 			else
				tmp_results = branchOp(current_nodeIndex, res, num_iterations)
				push!(tasks, tmp_results)
# 			end
			push!(tasks_fetched_TF, false)
		end
	
		# Check which jobs are done, fetch them, and update status of that node
		num_tasks = length(tasks)
		for i in 1:num_tasks
			if (tasks_fetched_TF[i] == false)
				#if (istaskdone(tasks[i]) == true)
					# Get the results
					calc_end_time = Dates.now()
# 					if (parallel_TF == true)
# 						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = fetch(tasks[i])
# 					else
						(tmp_threadID, nodeData_at_bottom, spawned_nodeIndex, calc_start_time) = tasks[i]
# 					end
					# Store run information
					res.calc_start_time[spawned_nodeIndex] = calc_start_time
					res.calc_end_time[spawned_nodeIndex] = calc_end_time
					res.calc_duration[spawned_nodeIndex] = (calc_end_time - calc_start_time).value / 1000.0
					tasks_fetched_TF[i] = true
					
					# Record information
					res.thread_for_each_branchOp[spawned_nodeIndex] = tmp_threadID
					sum_nodeData_at_bottom = sum(nodeData_at_bottom)
					res.likes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = nodeData_at_bottom .+ 0.0
					res.normlikes_at_each_nodeIndex_branchBot[spawned_nodeIndex] = (nodeData_at_bottom .+ 0.0) ./ sum_nodeData_at_bottom
					res.lq_at_branchBot[spawned_nodeIndex] = log(sum_nodeData_at_bottom)
					res.like_at_branchBot[spawned_nodeIndex] = sum_nodeData_at_bottom

					
					# Get the ancestor nodeIndex
					uppass_edgematrix = res.uppass_edgematrix
					TF = uppass_edgematrix[:,2] .== spawned_nodeIndex
					parent_nodeIndex = uppass_edgematrix[TF,1][1]

					# Get the left daughter nodeIndex (1st in the uppass_edgematrix)
					edge_rows_TF = uppass_edgematrix[:,1] .== parent_nodeIndex
					left_nodeIndex = uppass_edgematrix[edge_rows_TF,2][1]
					right_nodeIndex = uppass_edgematrix[edge_rows_TF,2][2]

					# Update the state of the parent_node's daughters
					if (spawned_nodeIndex == left_nodeIndex)
						res.node_Lparent_state[parent_nodeIndex] = "ready"
					end
					if (spawned_nodeIndex == right_nodeIndex)
						res.node_Rparent_state[parent_nodeIndex] = "ready"
					end

					# Update the state of the current node
					res.node_state[spawned_nodeIndex] = "done"
				#end
			end
		end
	
		# Update which nodes have had both parents complete
		TF1 = res.node_state .== "not_ready"
		TF2 = res.node_Lparent_state .== "ready"
		TF3 = res.node_Rparent_state .== "ready"
		TF = (TF1 + TF2 + TF3) .== 3
		res.node_state[TF] .= "ready_for_nodeOp"
	
		# Update nodes when the branches above finish
		indexes_ready = findall(res.node_state .== "ready_for_nodeOp")
		for current_nodeIndex in indexes_ready
			# Spawn a node operation
			#push!(tasks, @spawn nodeOp(current_nodeIndex, res))
			nodeOp(current_nodeIndex, res)
		end
	
		# Check if we are done?
		are_we_done = count_nodes_finished(res.node_state) >= res.numNodes
		
		# Error trap
		if (iteration_number >= max_iterations)
			txt = join(["Error in iterative_downpass_nonparallel(): iteration_number ", string(iteration_number), " exceeded max_iterations. Probably your loop is not concluding, or you have a massively huge tree or slow calculation, and need to set max_iterations=Inf."], "")
			error(txt)
		end
		
		# Test for concluding the while loop
		are_we_done && break
	end
	
	# This breaks it for some reason:
	# ERROR: setfield! immutable struct of type Res cannot be changed
	#global res.number_of_whileLoop_iterations = iteration_number

	print_num_iterations = false
	if print_num_iterations
		txt = join(["\nFinished at iteration_number ", string(iteration_number), "."], "")
		print(txt)
		print("\n")
	end
	
	# Final run diagnostics
	diagnostics[2] = Dates.now()
	diagnostics[3] = diagnostics[2]-diagnostics[1]
	total_calctime_in_sec = (diagnostics[2]-diagnostics[1]).value / 1000
	
	res.calctime_iterations[1] = total_calctime_in_sec
	res.calctime_iterations[2] = iteration_number / 1.0
	
	return(total_calctime_in_sec, iteration_number)
end # END iterative_downpass_nonparallel!





















end # end of module