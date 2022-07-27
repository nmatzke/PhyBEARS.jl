#######################################################
# Reorder a castor simulation from e.g. simulate_tdsse2()
# to the default APE node order
#
# Procedure: 
# 1. prints simulation$tree to a newick file,
#    then reads it back in
#
# 2. Reorders the simulation$tip_states and simulation$node_states to match 
#
# 3. Adds simulation$states, to put all the states in one place
#######################################################
reorder_castor_sim_to_default_ape_node_order <- function(simulation)
	{
	# Make tree table from simulation tree
	simtree = simulation$tree
	simtable = prt(simtree, printflag=FALSE, get_tipnames=TRUE)

	# Write simtree to Newick; read back in to get standard APE newick format
	newtree = read.tree(file="", text=write.tree(simulation$tree, file=""))
	# Write tree table for new tree
	newtable = prt(newtree, printflag=FALSE, get_tipnames=TRUE)

	match_simtree_in_new2 = match(x=newtable$tipnames, table=simtable$tipnames)
	match_simtree_in_new2

	match_simtree_in_new2sub_tips = match_simtree_in_new2[match_simtree_in_new2 <= length(simtree$tip.label)]
	match_simtree_in_new2sub_nodes = match_simtree_in_new2[match_simtree_in_new2 > length(simtree$tip.label)]-length(simtree$tip.label)

	simstates = c(as.numeric(simulation$tip_states), simulation$node_states)

	simstates
	simstates[match_simtree_in_new2]

	simulation$node_states
	simulation$node_states[match_simtree_in_new2sub_nodes]
	
	simulation$tip_states
	simulation$tip_states[match_simtree_in_new2sub_tips]
	
	# Re-save the re-ordered simulation
	simulation2 = simulation
	simulation2$tree = newtree
	simulation2$node_states = simulation$node_states[match_simtree_in_new2sub_nodes]
	simulation2$tip_states = simulation$tip_states[match_simtree_in_new2sub_tips]
	names(simulation2$tip_states) = 1:length(simulation2$tip_states)
	names(simulation2$tip_states)
	simulation2$states = simstates[match_simtree_in_new2]
	return(simulation2)
	}

