using DataFrames
using PhyloNetworks
using PhyloPlots

#include("/GitHub/PhyBEARS.jl/notes/TreePassO.jl")


#######################################################
# Standard bifurcating tree
#######################################################
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

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)

#######################################################
# Tree with a 2-degree node inside a branch
#######################################################

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

# Make a trdf, a DataFrame holding all of the key tree information
cumulative_height_at_each_node = get_node_heights(tr)
node_age = get_node_ages(tr)
