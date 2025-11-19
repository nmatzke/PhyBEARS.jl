library(cladoRcpp)
library(BioGeoBEARS)
library(ape)
library(castor)

# Set the random number seed, to make it repeatable
set.seed(54321)

time_grid = seq(0,10,0.1) # About halfway through
A = get_random_mk_transition_matrix(Nstates=4, rate_model="ER", max_rate=0.1)

# Make it more like a DEC model (anagenetic)
d_rate = 0.06 # range expansion
e_rate = 0.01 # range contraction

# 4 states are:
# null, A, B, AB
A[,] = 0
A[2,4] = d_rate # A->AB
A[3,4] = d_rate # B->AB
A[2,1] = e_rate # A->null
A[3,1] = e_rate # B->null
A[4,2] = e_rate # AB->A
A[4,3] = e_rate # AB->B
A

# In Q transition matrices, the diagonals = -sum(off-diagonal for that row)
diag(A) = 0.0
A
diag(A) = -rowSums(A)
A


# Cladogenetic part of the DEC model
# At speciation, we have:
# Specify probabilities of different events, given that speciation has occurred



#######################################################
# 1st regime: vicariance, subset sympatry not allowed - BAYAREALIKE model
#######################################################

# Sympatry
# null->null,null
# A -> A,A  # 100
# B -> B,B	# 100
# Vicariance
# AB -> A,B	# 0
# AB -> B,A	# 0
# Subset sympatry (speciation within widespread ancestor)
# AB -> AB,A	# 0
# AB -> A,AB	# 0
# AB -> AB,B	# 0
# AB -> B,AB	# 0
# AB -> AB, AB # 1.0

transition_matrix_C = get_random_mk_transition_matrix(Nstates=4, rate_model="ER", max_rate=0.1)
transition_matrix_C[,] = 0.0
transition_matrix_C
transition_matrix_C[1,1] = 1.0 # null->null
transition_matrix_C[2,2] = 1.0 # A->A
transition_matrix_C[3,3] = 1.0 # B->B
transition_matrix_C[4,] = 0.0
transition_matrix_C[4,1] = 0.0
transition_matrix_C[4,2] = 0.0 # AB->A
transition_matrix_C[4,3] = 0.0 # AB->B
transition_matrix_C[4,4] = 6/6 # AB->AB
transition_matrix_C1 = transition_matrix_C
transition_matrix_C

# Rows of transition_table_C
transition_table_C = NULL
tmprow = c(1,1,1,1.0)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(2,2,2,1.0)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(3,3,3,1.0)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(4,2,3,0)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(4,2,4,0)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(4,3,4,0)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(4,4,4,1.0)
transition_table_C = rbind(transition_table_C, tmprow)
transition_table_C

# convert from 1-based to 0-based indices (turns out this is NOT needed for the C++ code)
transition_table_indices_C = transition_table_C[,1:3] - 1
transition_table_probs_C1 = transition_table_C[,4]


#######################################################
# 2nd regime: vicariance, subset sympatry are allowed - DEC model
#######################################################

# Sympatry
# null->null,null
# A -> A,A  # 100
# B -> B,B	# 100
# Vicariance
# AB -> A,B	# 1/6
# AB -> B,A	# 1/6
# Subset sympatry (speciation within widespread ancestor)
# AB -> AB,A	# 1/6
# AB -> A,AB	# 1/6
# AB -> AB,B	# 1/6
# AB -> B,AB	# 1/6

transition_matrix_C = get_random_mk_transition_matrix(Nstates=4, rate_model="ER", max_rate=0.1)
transition_matrix_C[,] = 0.0
transition_matrix_C
transition_matrix_C[1,1] = 1.0 # null->null
transition_matrix_C[2,2] = 1.0 # A->A
transition_matrix_C[3,3] = 1.0 # B->B
transition_matrix_C[4,] = 0.0
transition_matrix_C[4,1] = 0.0
transition_matrix_C[4,2] = 2/6 # AB->A
transition_matrix_C[4,3] = 2/6 # AB->B
transition_matrix_C[4,4] = 2/6 # AB->AB
transition_matrix_C2 = transition_matrix_C

# Rows of transition_table_C
transition_table_C = NULL
tmprow = c(1,1,1,1.0)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(2,2,2,1.0)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(3,3,3,1.0)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(4,2,3,1/3)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(4,2,4,1/3)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(4,3,4,1/3)
transition_table_C = rbind(transition_table_C, tmprow)
tmprow = c(4,4,4,0.0)
transition_table_C = rbind(transition_table_C, tmprow)
transition_table_C

# convert from 1-based to 0-based indices (turns out this is NOT needed for the C++ code)
transition_table_indices_C = transition_table_C[,1:3] - 1
transition_table_probs_C2 = transition_table_C[,4]






# Make the arrays for the parameters
transition_matrix_C_array = array(data=0.0, dim=c(nrow(transition_matrix_C1), ncol(transition_matrix_C), length(time_grid)))

# Probs encoded in a num_clado_events * 2 * num_times array
# the "2" just duplicates, ensures we have a 3D array
transition_table_probs_C_matrix = array(data=0.0, dim=c(length(transition_table_probs_C1), 2, length(time_grid)))

# Fill in the transition_matrix_C and transition_table_probs
for (i in 1:length(time_grid))
	{
	if (time_grid[i] < 5.0)
		{
		transition_matrix_C_array[,,i] = transition_matrix_C1
		transition_table_probs_C_matrix[,1,i] = transition_table_probs_C1
		transition_table_probs_C_matrix[,2,i] = transition_table_probs_C1
		} else {
		transition_matrix_C_array[,,i] = transition_matrix_C2		
		transition_table_probs_C_matrix[,1,i] = transition_table_probs_C2
		transition_table_probs_C_matrix[,2,i] = transition_table_probs_C2
		}
	}


parameters = list(birth_rates         = 0.3,
                  death_rates         = 0.0,
                  transition_matrix_A = A,
                  transition_matrix_C = transition_matrix_C_array,
                  transition_table_indices_C = transition_table_indices_C,
                  transition_table_probs_C = transition_table_probs_C_matrix )
 

simulation = simulate_tdsse2( Nstates = 4, 
                            parameters = parameters, 
														splines_degree      = 1,
                            start_state = 2,
                            max_tips = 50,
                            time_grid = time_grid,
                            include_birth_times=TRUE,
                            include_death_times=TRUE,
                            coalescent=FALSE)
plot(simulation$tree); axisPhylo(); mtext(text="Millions of years ago (Ma)", side=1, line=2)
simulation$Ntransitions_A
simulation$Ntransitions_C

simtree = simulation$tree
simtable = prt(simtree, printflag=FALSE, get_tipnames=TRUE)
newtree = read.tree(file="", text=write.tree(simulation$tree, file=""))
newtable = prt(newtree, printflag=FALSE, get_tipnames=TRUE)

match_simtree_in_new1 = match(x=simtable$tipnames, table=newtable$tipnames)
match_simtree_in_new1

match_simtree_in_new2 = match(x=newtable$tipnames, table=simtable$tipnames)
match_simtree_in_new2


match_simtree_in_new1sub = match_simtree_in_new1[match_simtree_in_new1 > length(simtree$tip.label)]-length(simtree$tip.label)
match_simtree_in_new1sub
match_simtree_in_new2sub = match_simtree_in_new2[match_simtree_in_new2 > length(simtree$tip.label)]-length(simtree$tip.label)
match_simtree_in_new2sub

simstates = c(as.numeric(simulation$tip_states), simulation$node_states)

simstates
simstates[match_simtree_in_new1]
simstates[match_simtree_in_new2]

simulation$node_states
simulation$node_states[match_simtree_in_new1sub]
simulation$node_states[match_simtree_in_new2sub]


#######################################################
# Make colors
#######################################################
library(cladoRcpp)
areas = c("A", "B")
max_range_size = 2
include_null_range = TRUE
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
states_list_0based

possible_ranges_list_txt = areas_list_to_states_list_new(areas=areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
possible_ranges_list_txt

# Make the list of ranges by node
ranges_list = NULL



for (i in 1:length(states_list_0based))
    {    
    if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
        {
        tmprange = "_"
        } else {
        tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
        }
    ranges_list = c(ranges_list, tmprange)
    }

# Look at the ranges list
ranges_list




colors_matrix = get_colors_for_numareas(length(areas))
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index=states_list_0based, plot_null_range=include_null_range)
colors_list_for_states[length(colors_list_for_states)] = "cyan"   # usually "all areas" is white, 
																																		# but with just AB, I like orange

# Get the colors by node
statetxt_by_node = rep("_", times=length(simstates))
for (i in 1:length(statetxt_by_node))
	{
	statetxt_by_node[i] = ranges_list[simstates[i]]
	}
	
statetxt_by_node
cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates=statetxt_by_node)
cols_byNode


pdffn = "compare_castor_simtree_newtree_v1.pdf"
pdf(file=pdffn, width=6, height=12)


# Plot the original simtree
tipnums = 1:length(simtree$tip.label)
nodenums = (length(simtree$tip.label)+1):(length(simtree$tip.label)+simtree$Nnode)
plot(simtree, show.tip.label=FALSE)
axisPhylo()
mtext(text="Millions of years ago (Mega-annum, Ma)", side=1, line=2.5)
title("Plotting simtree & states from simulate_tdsse2()")
tiplabels(text=statetxt_by_node[tipnums], tip=tipnums, bg=cols_byNode[tipnums])
nodelabels(text=statetxt_by_node[nodenums], node=nodenums, bg=cols_byNode[nodenums])

# Try the re-orderings #2
statetxt_by_node2 = statetxt_by_node[match_simtree_in_new2]
cols_byNode2 = cols_byNode[match_simtree_in_new2]

plot(newtree, show.tip.label=FALSE)
axisPhylo()
mtext(text="Millions of years ago (Mega-annum, Ma)", side=1, line=2.5)
title("simulate_tdsse2() tree after reordering attempt #2\nto default APE format -- CORRECT!")
tiplabels(text=statetxt_by_node2[tipnums], tip=tipnums, bg=cols_byNode2[tipnums])
nodelabels(text=statetxt_by_node2[nodenums], node=nodenums, bg=cols_byNode2[nodenums])

# Try the re-orderings #1
statetxt_by_node1 = statetxt_by_node[match_simtree_in_new1]
cols_byNode1 = cols_byNode[match_simtree_in_new1]

plot(newtree, show.tip.label=FALSE)
axisPhylo()
mtext(text="Millions of years ago (Mega-annum, Ma)", side=1, line=2.5)
title("simulate_tdsse2() tree after reordering attempt #1\nto default APE format")
tiplabels(text=statetxt_by_node1[tipnums], tip=tipnums, bg=cols_byNode1[tipnums])
nodelabels(text=statetxt_by_node1[nodenums], node=nodenums, bg=cols_byNode1[nodenums])

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


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
	} # END reorder_castor_sim_to_default_ape_node_order


# This assumes the list of states can be generated by
# rcpp_areas_list_to_states_list(); a manual states list 
# would need another function
plot_castor_simulation <- function(simulation, areas=2, max_range_size=2, include_null_range=TRUE)
	{
	defaults='
	areas = c("A", "B")
	max_range_size = 2
	include_null_range = TRUE
	'
	
	states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
	states_list_0based

	possible_ranges_list_txt = areas_list_to_states_list_new(areas=areas,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
	possible_ranges_list_txt

	# Make the list of ranges by node
	ranges_list = NULL

	for (i in 1:length(states_list_0based))
			{    
			if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
					{
					tmprange = "_"
					} else {
					tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
					}
			ranges_list = c(ranges_list, tmprange)
			}

	# Look at the ranges list
	ranges_list

	colors_matrix = get_colors_for_numareas(length(areas))
	colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index=states_list_0based, plot_null_range=include_null_range)
	colors_list_for_states[length(colors_list_for_states)] = "cyan"   # usually "all areas" is white, 
																																			# but with just AB, I like orange

	# Get the colors by node
	statetxt_by_node = rep("_", times=length(simstates))
	for (i in 1:length(statetxt_by_node))
		{
		statetxt_by_node[i] = ranges_list[simstates[i]]
		}
	statetxt_by_node
	
	cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates=statetxt_by_node)
	cols_byNode

	# Plot the original simtree
	tipnums = 1:length(simtree$tip.label)
	nodenums = (length(simtree$tip.label)+1):(length(simtree$tip.label)+simtree$Nnode)
	plot(simtree, show.tip.label=FALSE)
	axisPhylo()
	mtext(text="Millions of years ago (Mega-annum, Ma)", side=1, line=2.5)
	title("Plotting simtree & states from simulate_tdsse2()")
	tiplabels(text=statetxt_by_node[tipnums], tip=tipnums, bg=cols_byNode[tipnums])
	nodelabels(text=statetxt_by_node[nodenums], node=nodenums, bg=cols_byNode[nodenums])

	return(cols_byNode)
	} # END plot_castor_simulation

simulation2 = reorder_castor_sim_to_default_ape_node_order(simulation)

# Get the colors by node
statetxt_by_node3 = rep("_", times=length(simulation2$states))
for (i in 1:length(statetxt_by_node))
	{
	statetxt_by_node3[i] = ranges_list[simulation2$states[i]]
	}
	
statetxt_by_node3
cols_byNode3 = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates=statetxt_by_node3)
cols_byNode3


plot(simulation2$tree, show.tip.label=FALSE)
axisPhylo()
mtext(text="Millions of years ago (Mega-annum, Ma)", side=1, line=2.5)
title("simulate_tdsse2() tree after reordering with\ncastor_sim_to_default_ape_node_order()")
tiplabels(text=statetxt_by_node3[tipnums], tip=tipnums, bg=cols_byNode3[tipnums])
nodelabels(text=statetxt_by_node3[nodenums], node=nodenums, bg=cols_byNode3[nodenums])


#plot_castor_simulation(simulation2, areas=2, max_range_size=2, include_null_range=TRUE)

