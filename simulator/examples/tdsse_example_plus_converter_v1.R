library(BioGeoBEARS)
library(ape)
library(castor)

time_grid = seq(0,10,0.1) # About halfway through
A = get_random_mk_transition_matrix(Nstates=4, rate_model="ER", max_rate=0.1)

# Make it more like a DEC model (anagenetic)
d_rate = 0.034 # range expansion
e_rate = 0.028 # range contraction

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
                            time_grid = time_grid)
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

