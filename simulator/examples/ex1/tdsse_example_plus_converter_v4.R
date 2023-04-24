
# Dependencies
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(ape)
library(castor)

# Helper functions for castor
source("/GitHub/PhyBEARS.jl/Rsrc/castor_helpers.R")

wd = "/GitHub/PhyBEARS.jl/simulator/examples/ex1/"
setwd(wd)

# Set the random number seed, to make it repeatable
set.seed(543210)
set.seed(54321)

##########################################################
# External inputs for making geography files etc.
##########################################################
# This example simulation is 2 AREAS, 4 STATES
# Areas: A, B
# States: null, A, B, AB

area_names = c("A","B")
max_range_size = length(area_names)
include_null_range = TRUE

states_list = rcpp_areas_list_to_states_list(areas=area_names, maxareas=max_range_size, include_null_range=include_null_range)

numareas = length(area_names)
numstates = length(states_list)

time_grid = seq(0,10,0.1) # Change about halfway through
##########################################################

# Start with a random rate matrix between the 4 states
# (we will edit this later)
A = get_random_mk_transition_matrix(Nstates=4, rate_model="ER", max_rate=0.1)

# Make it more like a DEC model (anagenetic)
d_rate = 0.06 # range expansion
e_rate = 0.01 # range contraction
j_rate = 0.0 # Not used -- if you want to use it, use in transition_rate_C etc.

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
                            max_tips = 51,
                            time_grid = time_grid,
                            include_birth_times=TRUE,
                            include_death_times=TRUE,
                            coalescent=FALSE)
plot(simulation$tree); axisPhylo(); mtext(text="Millions of years ago (Ma)", side=1, line=2)
simulation$Ntransitions_A
simulation$Ntransitions_C

orig_sim = simulation

simulation2 = reorder_castor_sim_to_default_ape_node_order(simulation)

# This proves we know node states match to R
# plot_castor_simulation - assumes simulation object has been reordered
plot_castor_simulation(simulation2, areas=area_names, max_range_size=max_range_size, include_null_range=include_null_range)

write_out_original_castor_simfiles(simulation, wd, fn="rawsim")
plot(simulation$tree)
simulation$node_states
simulation$tip_states

dorder = postorder_nodes_phylo4_return_table2(simulation$tree)
preorder = dorder[rev(order(dorder$nodenums)),]
preorder

# Remove the last tip from the simulation
source('~/HD/GitHub/PhyBEARS.jl/Rsrc/castor_helpers.R', chdir = TRUE)
simulation3 = remove_last_tip_from_simulation(simulation2)
plot(simulation3$tree)

plot_castor_simulation(simulation3, areas=area_names, max_range_size=max_range_size, include_null_range=include_null_range)




# Write out files with tip 51 removed
write_out_original_castor_simfiles(simulation3, wd, fn="rawsim_wo_tip51")

# Reorder to be more readable
# Write out A WHOLE BUNCH OF FILES, e.g. geography etc.
simulation4 = reorder_castor_sim_to_default_ape_node_order(simulation3)


# Reorder, and also write out the reordered files plus geography etc.
write_out_reordered_castor_simfiles(simulation4, wd, area_names=area_names, states_list=states_list)

plot_castor_simulation(simulation4, areas=area_names, max_range_size=max_range_size, include_null_range=include_null_range)


moref("geog_living.data")

# Look at the files in Finder
system(paste0("open ", wd))



# Plot the damn thing to double-check!
dev.off()
dev.off()


#######################################################
# PDF plots via BioGeoBEARS -- MUCH SIMPLER to just
# use plot_castor_simulation()
#######################################################
pdffn = "simbiogeog_on_tree_wFossils_noNodeLabels.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - full tree with fossils
#######################################################
trfn = "tree_wFossils_noNodeLabels.newick"
geogfn = "geog_wFossils.data"
simstates_fn = "simstates_all.txt"

tr = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(geogfn)

analysis_titletxt = paste0("Castor biogeog simulation plotted on ", trfn, "\n(full tree, with fossils)")

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = include_null_range
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

BioGeoBEARS_run_object$allow_null_tips = TRUE # allows null-range tips

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = d_rate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = e_rate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = j_rate

# Calculate (fake) ancestral states for a single parameter guess
res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=TRUE, skip_optim_option="return_all")

res$optim_result 

# Replace fake ancestral states with "true" simulated ancestral states
simstates = read.table(simstates_fn)

for (i in 1:nrow(simstates))
	{
	res$ML_marginal_prob_each_state_at_branch_top_AT_node[i,] = 0.0
	res$ML_marginal_prob_each_state_at_branch_top_AT_node[i,simstates[i,]] = 1.0
	}
res$ML_marginal_prob_each_state_at_branch_top_AT_node

# Setup for plot
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=include_null_range, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=include_null_range, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - living-only tree, *WITHOUT* fossils
#######################################################
trfn = "tree_wFossils_noNodeLabels.newick"
geogfn = "geog_living.data"
simstates_fn = "simstates_living.txt"

tr = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(geogfn)

analysis_titletxt = paste0("Castor biogeog simulation plotted on ", trfn, "\n(living tree, *without* fossils)")

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = include_null_range
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

BioGeoBEARS_run_object$allow_null_tips = TRUE # allows null-range tips

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = d_rate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = e_rate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = j_rate

# Calculate (fake) ancestral states for a single parameter guess
res = bears_optim_run(BioGeoBEARS_run_object, skip_optim=TRUE, skip_optim_option="return_all")

# Replace fake ancestral states with "true" simulated ancestral states
simstates = read.table(simstates_fn)

for (i in 1:nrow(simstates))
	{
	res$ML_marginal_prob_each_state_at_branch_top_AT_node[i,] = 0.0
	res$ML_marginal_prob_each_state_at_branch_top_AT_node[i,simstates[i,]] = 1.0
	}
res$ML_marginal_prob_each_state_at_branch_top_AT_node

# Setup for plot
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=include_null_range, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=include_null_range, tr=tr, tipranges=tipranges)



dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it





