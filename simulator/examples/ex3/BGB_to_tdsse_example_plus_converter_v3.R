#######################################################
# Go from a BioGeoBEARS run object, to a castor TDSSE simulation, to converter
#######################################################

wd = "/GitHub/PhyBEARS.jl/simulator/examples/ex3/"
setwd(wd)

# Dependencies
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(ape)
library(castor)

# Set the random number seed, to make it repeatable
set.seed(54321)

#######################################################
# Load an example BioGeoBEARS run object
#######################################################
trfn = "Psychotria_5.2.newick"
geogfn = "Psychotria_geog.data"
tr = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
numareas = ncol(tipranges@df)
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = numareas
include_null_range = TRUE
numstates = numstates_from_numareas(numareas=numareas, maxareas=numareas, include_null_range=include_null_range)

# Make it like Psychotria DEC+J model (anagenetic)
d_rate = 0.000001 # range expansion
e_rate = 0.000001 # range contraction


BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = d_rate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = e_rate

returned_mats = get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object=BioGeoBEARS_run_object, numstates=numstates, include_null_range=include_null_range, timeperiod_i=1)
names(returned_mats)

Carray_df = get_Cevent_probs_df_from_mats(mats=returned_mats, include_null_range=include_null_range)


time_grid = seq(0,20,0.1) # About halfway through
#A = get_random_mk_transition_matrix(Nstates=4, rate_model="ER", max_rate=0.1)
A = returned_mats$Qmat


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
                            max_tips = 51,
                            time_grid = time_grid,
                            include_birth_times=TRUE,
                            include_death_times=TRUE,
                            coalescent=FALSE)
plot(simulation$tree); axisPhylo(); mtext(text="Millions of years ago (Ma)", side=1, line=2)
simulation$Ntransitions_A
simulation$Ntransitions_C


# Helper functions for castor
source("/GitHub/PhyBEARS.jl/Rsrc/castor_helpers.R")

wd = "/GitHub/PhyBEARS.jl/simulator/examples/ex1/"
setwd(wd)

write_out_original_castor_simfiles(simulation, wd, fn="rawsim")

# Remove the last tip from the simulation
simulation2 = remove_last_tip_from_simulation(simulation)
plot(simulation2$tree)

# Write out files with tip 51 removed
write_out_original_castor_simfiles(simulation2, wd, fn="rawsim_wo_tip51")

# Reorder to be more readable
# Write out A WHOLE BUNCH OF FILES, e.g. geography etc.
simulation3 = reorder_castor_sim_to_default_ape_node_order(simulation2)

# External inputs for making geography files etc.
area_names = c("A","B","C","D")
states_list = rcpp_areas_list_to_states_list(areas=area_names, maxareas=length(area_names), include_null_range=TRUE)

# Reorder, and also write out the reordered files plus geography etc.
write_out_reordered_castor_simfiles(simulation3, wd, area_names=c("A","B","C","D"), states_list=states_list)

# Look at the files in Finder
system(paste0("open ", wd))

