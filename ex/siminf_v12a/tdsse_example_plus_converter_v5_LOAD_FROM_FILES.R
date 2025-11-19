library(cladoRcpp)
library(BioGeoBEARS)
library(ape)

# modified version of castor

# Install modified "castor" package in R
install='
install.packages(c('naturalsort', 'RSpectra'))
R CMD install ~/GitHub/PhyBEARS.jl/simulator/castor_1.7.2.000004.tar.gz
'
# install.packages(pkgs="~/GitHub/PhyBEARS.jl/simulator/castor_1.7.2.000004.tar.gz", lib="/Library/Frameworks/R.framework/Resources/library/", repos=NULL, type="source")

# install.packages("~/GitHub/PhyBEARS.jl/simulator/castor_1.7.2.000004.tar.gz
library(castor)

# for: reorder_castor_sim_to_default_ape_node_order(simulation)
source("~/GitHub/PhyBEARS.jl/Rsrc/castor_helpers.R")

wd = "~/GitHub/PhyBEARS.jl/ex/siminf_v12a/"
setwd(wd)
simfns = c("setup_df.txt",
"timepoints.txt", 
"mu_vals_by_t.txt", 
"Qvals_by_t.txt",
"Crates_by_t.txt",
"Qarray.txt",
"Carray.txt",
"area_names.txt",
"states_list.R")


simulation2 = simulate_tdsse2_for_timeperiod(wd, start_state=2, max_simulation_time=12, min_tips=50, max_tips=50, simfns=default_simfns(), seedval=54321, max_rate=10.0, numtries=250)


# Setup
setup='
numareas = 3
numstates = 8
max_range_size = 3
include_null_range = TRUE
'

# Read in the setup
setup_df = read.table("setup_df.txt", header=TRUE)
numareas = setup_df$numareas
numstates = setup_df$numstates
max_range_size = setup_df$max_range_size
if (tolower(setup_df$include_null_range) == "true")
	{
	include_null_range = TRUE
	} else {
	include_null_range = FALSE
	}

# User sets these here!
# Set the random number seed, to make it repeatable
seedval = 54321
set.seed(seedval)
# Set the max_simulation time -- whatever changing distances you have will be
# extended with this timepoint, i.e. the final geography/rates will continue
# to this timepoint
start_state = 8 # number of the starting state
max_simulation_time = 15.0 # Set to 0 if the user doesn't want to set a max simulation time
max_tips = NULL

# Read timepoints
# Read in Q/A matrix, populate one and array
# Read in C matrix, populate one and array

# Load the files
#time_grid = seq(0,10,0.1) # About halfway through
time_grid = c(read.table(simfns[2], header=FALSE))[[1]]

# Change from time reading backwards to time reading forwards
colnums = rev(1:length(time_grid))
time_grid = -1*(max(time_grid) - time_grid) - min(-1*(max(time_grid) - time_grid))
time_grid
mu_vals_by_t = as.matrix(read.table(simfns[3], header=FALSE))[,colnums]
Qvals_by_t = as.matrix(read.table(simfns[4], header=FALSE))[,colnums]
Crates_by_t = as.matrix(read.table(simfns[5], header=FALSE))[,colnums]
Qarray = read.table(simfns[6], header=TRUE)
Carray = read.table(simfns[7], header=TRUE)

# Add the final time at end of forward simulation, if needed
if (max(time_grid) < max_simulation_time)
	{
	time_grid = c(time_grid, max_simulation_time)
	mu_vals_by_t = cbind(mu_vals_by_t, mu_vals_by_t[,ncol(mu_vals_by_t)])
	Qvals_by_t = cbind(Qvals_by_t, Qvals_by_t[,ncol(Qvals_by_t)])
	Crates_by_t = cbind(Crates_by_t, Crates_by_t[,ncol(Crates_by_t)])
	}

# Set a maximum rate for extreme cases
max_rate = 30.0

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


# Cladogenetic part of the DEC model
# At speciation, we have:
# Specify probabilities of different events, given that speciation has occurred
# (ADD THESE UP to provide the speciation rates / lambdas)

# transition_table_C: columns are i, j, k, prob
# transition_table_rates_C_matrix: columns are i, j, k, rate 
# (later, we will add up these rates to get the total lambda by state)
transition_table_C = NULL
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
		if ((Carray$i[j] == Carray$j[j]) && (Carray$i[j] == Carray$k[j]))
			{
			# state i = state j = state k
			transition_matrix_C_array[Carray$i[j],Carray$j[j],i] = transition_matrix_C_array[Carray$i[j],Carray$j[j],i] + Crates_by_t[j,i] / rates_sums_by_t[Carray$i[j],i]
			} else {
			# Otherwise, you have two different "anagenetic" events, with the probability split between them
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
 

simulation
simulation = simulate_tdsse2( Nstates = numstates, 
                            parameters = parameters, 
														splines_degree      = 1,
                            start_state = 8,
                            max_tips = max_tips,
                            max_time = max_simulation_time,
                            time_grid = time_grid,
                            include_birth_times=TRUE,
                            include_death_times=TRUE,
                            coalescent=FALSE)
simulation


source("~/GitHub/PhyBEARS.jl/Rsrc/castor_helpers.R")

# Check for successful simulation
if (simulation$success == TRUE)
	{
	simulation2 = reorder_castor_sim_to_default_ape_node_order(simulation)
	print(unname(simulation2$tip_states))
	plot_castor_simulation(simulation2, areas=areas, max_range_size=length(areas), include_null_range=TRUE)
	} else {
	simulation2 = NULL
	}

