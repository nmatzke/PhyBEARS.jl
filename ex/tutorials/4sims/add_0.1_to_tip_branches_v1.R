#######################################################
# Adding 0.1 to each tip branch
# (this is kind of a hack; it would be better for
#  the simulator to go to e.g. 51 species, then cut
#  the 51st species)
#######################################################

wd = "/GitHub/PhyBEARS.jl/ex/tutorials/4sims/ss8_sim_001/"
setwd(wd)

treefilename = "tree.newick"

tr = read.tree(treefilename)

tip_nodes = 1:(length(tr$tip.label))

trtable = prt(tr, printflag=FALSE)
trtable

branch_numbers_to_change = trtable$parent_br[tip_nodes]

tr2 = tr

tr2$edge.length[branch_numbers_to_change] = tr2$edge.length[branch_numbers_to_change] + 0.1

trtable2 = prt(tr2, printflag=FALSE)