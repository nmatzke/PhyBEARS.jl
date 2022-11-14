library(cladoRcpp)
library(BioGeoBEARS)
library(ape)
library(castor)

# Setup
numstates = 8
wd = "/GitHub/PhyBEARS.jl/ex/siminf_v12a/"
setwd(wd)
outfns = c(
"timepoints.txt", 
"mu_vals_by_t.txt",
"Qvals_by_t.txt",
"Crates_by_t.txt",
"Qarray.txt",
"Carray.txt")


# Set the random number seed, to make it repeatable
set.seed(54321)
# Set the max_simulation time -- whatever changing distances you have will be
# extended with this timepoint, i.e. the final geography/rates will continue
# to this timepoint
max_simulation_time = NULL
max_simulation_time = 100.0
max_tips = NULL

# Read timepoints
# Read in Q/A matrix, populate one and array
# Read in C matrix, populate one and array

# Load the files
#time_grid = seq(0,10,0.1) # About halfway through
time_grid = c(read.table(outfns[1], header=FALSE))[[1]]

# Change from time reading backwards to time reading forwards
colnums = rev(1:length(time_grid))
time_grid = -1*(max(time_grid) - time_grid) - min(-1*(max(time_grid) - time_grid))
time_grid
mu_vals_by_t = as.matrix(read.table(outfns[2], header=FALSE))[,colnums]
Qvals_by_t = as.matrix(read.table(outfns[3], header=FALSE))[,colnums]
Crates_by_t = as.matrix(read.table(outfns[4], header=FALSE))[,colnums]
Qarray = read.table(outfns[5], header=TRUE)
Carray = read.table(outfns[6], header=TRUE)

# Add the final time at end of forward simulation, if needed
if (max(time_grid) < max_simulation_time)
	{
	time_grid = c(time_grid, max_simulation_time)
	mu_vals_by_t = cbind(mu_vals_by_t, mu_vals_by_t[,ncol(mu_vals_by_t)])
	Qvals_by_t = cbind(Qvals_by_t, Qvals_by_t[,ncol(Qvals_by_t)])
	Crates_by_t = cbind(Crates_by_t, Crates_by_t[,ncol(Crates_by_t)])
	}

# Set a maximum rate for extreme cases
max_rate = 10.0

# Enforce maximum rate on Crates_by_t
Crates_by_t[Crates_by_t > max_rate] = max_rate
mu_vals_by_t[mu_vals_by_t > max_rate] = max_rate

# Produce the A transition matrix / array
A = array(data=0.0, dim=c(numstates,numstates,length(time_grid)))
for (i in 1:length(time_grid))
	{
	for (j in 1:nrow(Qarray))
		{
		A[Qarray$i[j],Qarray$j[j],i] = Qvals_by_t[j,i]
		}
	
	# Enforce maximum rate
	A[A > max_rate] = max_rate
	
	# Set the diagonal
	diag(A[,,i]) = 0.0
	diag(A[,,i]) = -rowSums(A[,,i])
	}

# All rows in A sum to 0.0!
round(apply(X=A, MARGIN=3, rowSums), digits=10)


ignore = '
A = get_random_mk_transition_matrix(Nstates=numstates, rate_model="ER", max_rate=0.1)

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
' # END Ignore

# Cladogenetic part of the DEC model
# At speciation, we have:
# Specify probabilities of different events, given that speciation has occurred
# (ADD THESE UP to provide the speciation rates / lambdas)

# transition_table_C: columns are i, j, k, prob
# transrates_table_C: columns are i, j, k, rate # (later, we will add up these rates to get the total lambda by state)
transition_table_C = NULL
transrates_table_C = NULL
transition_table_probs_C_matrix = array(data=0.0, dim=c(nrow(Crates_by_t), 4, length(time_grid)))
transition_table_rates_C_matrix = array(data=0.0, dim=c(nrow(Crates_by_t), 4, length(time_grid)))
dim(transition_table_probs_C_matrix)
dim(transition_table_rates_C_matrix)

rates_sums_by_t = NULL
for (i in 1:length(time_grid))
	{
	tmprates = as.data.frame(cbind(Carray$i, Carray$j, Carray$k, Crates_by_t[,i]), stringsAsFactors=FALSE)
	names(tmprates) = c("i", "j", "k", "rates_t")
	# Convert rates to probabilities
	rates_sums = aggregate(tmprates$rates_t, by=list(tmprates$i), FUN=sum)
	# Column names
	names(rates_sums) = c("i", "rates_sum")
	
	
	# Make sure rates_sums are ordered by 'i' (ancestral state index)
	tmporder = order(rates_sums$i)
	rates_sums = rates_sums[tmporder,]
	
	# Check for null range; add row for "i" if needed
	if (rates_sums$i[1] == 2)
		{
		tmprow = as.data.frame(matrix(data=c(1,0.0), nrow=1), stringsAsFactors=FALSE)
		names(tmprow) = c("i", "rates_sum")
		rates_sums = rbind(tmprow, rates_sums)
		}
	
	rates_sums_by_t	= cbind(rates_sums_by_t, rates_sums$rates_sum)

	# Convert rates to probabilities
	tmpprobs = tmprates
	names(tmpprobs) = c("i", "j", "k", "probs_t")
	for (q in 1:nrow(rates_sums))
		{
		tmpi = rates_sums$i[q]
		TF = tmprates$i == tmpi

		tmpprobs$probs_t[TF] = tmprates$rates_t[TF] / rates_sums$rates_sum[q]
		}
	
	transition_table_rates_C_matrix[,,i] = as.matrix(tmprates)
	transition_table_probs_C_matrix[,,i] = as.matrix(tmpprobs)
	}


# Square transition matrix (not used, but sanity check)
transition_matrix_C = array(data=0.0, dim=c(numstates,numstates,1))
transition_matrix_C_array = array(data=0.0, dim=c(numstates,numstates,length(time_grid)))
for (i in 1:length(time_grid))
	{
	for (j in 1:nrow(Carray))
		{
		if (Carray$i[j] == Carray$j[j])
			{
			transition_matrix_C_array[Carray$i[j],Carray$j[j],i] = transition_matrix_C_array[Carray$i[j],Carray$j[j],i] + Crates_by_t[j,i] / rates_sums_by_t[Carray$i[j],i]
			} else {
			transition_matrix_C_array[Carray$i[j],Carray$j[j],i] = transition_matrix_C_array[Carray$i[j],Carray$j[j],i] + Crates_by_t[j,i] / rates_sums_by_t[Carray$i[j],i] / 2
			transition_matrix_C_array[Carray$i[j],Carray$k[j],i] = transition_matrix_C_array[Carray$i[j],Carray$k[j],i] + Crates_by_t[j,i] / rates_sums_by_t[Carray$i[j],i] / 2
			}
		# Ignore, rates are already doubled if needed
		# if (Carray$i[j] != Carray$j[j])
		}
	# DO NOT set the diagonal on the cladogenetic transition rates
	#diag(transition_matrix_C_array[,,i]) = 0.0
	#diag(transition_matrix_C_array[,,i]) = -rowSums(transition_matrix_C_array[,,i])

	# Null range correction
	if (rowSums(transition_matrix_C_array[,,i])[1] == 0.0)
		{
		transition_matrix_C_array[1,1,i] = 1.0
		}
	}
transition_matrix_C = transition_matrix_C_array[,,dim(transition_matrix_C_array)[3]]
transition_matrix_C
rowSums(transition_matrix_C)


# Check probs
for (i in 1:length(time_grid))
	{
	tmpprobs = transition_table_probs_C_matrix[,,i]
	probs_sum = aggregate(tmpprobs[,4], by=list(tmpprobs[,1]), FUN=sum)
	tmprates = transition_table_rates_C_matrix[,,i]
	rates_sum = aggregate(tmprates[,4], by=list(tmprates[,1]), FUN=sum)
	sums_by_i = cbind(rates_sum, probs_sum[,2])
	names(sums_by_i) = c("i", "rates_sum", "probs_sum")
	print(sums_by_i)
	}


# convert from 1-based to 0-based indices (turns out this is NOT needed for the C++ code)
# transition_table_indices_C = transition_table_C[,1:3] - 1
transition_table_indices_C = transition_table_rates_C_matrix[,1:3,1] - 1
transition_table_indices_C_matrix = array(data=transition_table_indices_C, dim=c(dim(transition_table_indices_C),length(time_grid)))
# transition_table_probs_C1 = transition_table_C[,4]
transition_table_probs_C1 = transition_table_probs_C_matrix[,4,1]






# Make the arrays for the parameters

# Probs encoded in a num_clado_events * 2 * num_times array
# the "2" just duplicates, ensures we have a 3D array
transition_table_probs_C = array(data=0.0, dim=c(nrow(transition_table_indices_C), 2, length(time_grid)))

# Fill in the transition_matrix_C and transition_table_probs
for (i in 1:length(time_grid))
	{
	transition_table_probs_C[,1,i] = transition_table_probs_C_matrix[,4,i]
	transition_table_probs_C[,2,i] = transition_table_probs_C_matrix[,4,i]
	}


parameters = list(birth_rates         = rates_sums_by_t,
                  death_rates         = mu_vals_by_t,
                  transition_matrix_A = A,
                  transition_matrix_C = transition_matrix_C_array,
                  transition_table_indices_C = transition_table_indices_C,
                  transition_table_probs_C = transition_table_probs_C )
 

simulation = simulate_tdsse2( Nstates = numstates, 
                            parameters = parameters, 
														splines_degree      = 1,
                            start_state = 2,
                            max_tips = max_tips,
                            max_time = max_time,
                            time_grid = time_grid,
                            include_birth_times=TRUE,
                            include_death_times=TRUE,
                            coalescent=FALSE
)
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


plot_castor_simulation(simulation2, areas=2, max_range_size=2, include_null_range=TRUE)

