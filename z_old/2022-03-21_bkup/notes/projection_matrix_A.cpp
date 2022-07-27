Calculation of A(t), the linear dynamics at time t, given precalculated Es and model parameters

# castor:

# Returns a matrix A, for time t

# Elsewhere, castor says:

	// model parameters
	std::vector<double> transition_rates; // 2D array of size Nstates x Nstates, in row-major format, listing Markov transition rates between states. transition_rates[r,c] is the transition rate r-->c. Non-diagonal entries must be positive, the sum of each row must be zero.
	
	...so the diagonal is negative, like a Q matrix

// get_LinearDynamics_A
	// provide matrix encoding linear rates of change of the current state X, i.e. return A(t), where:
	//   dX/dt = A(t)*X(t) (if inverse==false)
	// or:
	//   dX/dt = X(t)*A(t) (if inverse==true)
	// note that, in principle, A may also depend on the current state X, i.e. A=A(t,X(t))
	// The returned A must be in row-major format
	void getLinearDynamics (double age, std::vector<double> &A) const{
		const MuSSEstateE current_E = E(age);
		// The mapping A is essentially the transition_matrix, plus some additional terms on the diagonal
		A = transition_rates;
		for(long r=0; r<Nstates; ++r){
			A[r*Nstates+r] += - (speciation_rates[r]+extinction_rates[r]) + 2*speciation_rates[r]*current_E[r]; // add more terms to diagonal
		}
		if(inverse) A *= -1;
	}


# Diversitree:

# Combine the lambdas and Qs across into a new Q matrix,
# for the off-diagnoals
# the instantaneous rate of change of the likelihoods
# On the diagonal, subtract colsums(A) (giving instantaneous 
# rate of nothing happening
# And add sum(lambdas for that state) - mu(lambdas for that state)

# diversitree:::projection.matrix.classe
projection.matrix.classe <- function(pars, k) 
	{
	A <- matrix(0, nrow = k, ncol = k)
	nsum <- k * (k + 1)/2  # nsum: number of descendant state-pairs for a single state
	
	# number of lambda parameters = nsum * k = k * k * (k+1)/2 
	# For each ancestral state, there are a pair of descendant states
	# BUT, we only want combinations, not permutations
	# for 4 states, there are 6 combinations of different pairs:
	# gtools::combinations(n=4,r=2)
	# ...and 4 identical pairs; 6+4 = q0
	# For 4 ancestral states, there are 4*10 = 40 lambdas 
	
	kseq <- seq_len(k)
	pars.lam <- pars[seq(1, nsum * k)]  # params 1-40
	pars.mu <- pars[seq(nsum * k + 1, (nsum + 1) * k)]  # pars 41-44 (40+1, 11*4)
	pars.q <- pars[seq((nsum + 1) * k + 1, length(pars))] # pars 45-56
	
	# Build a matrix of ancestor / left / right triplets
	# 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3
  #       3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4
	# 1 1 1 1 2 2 2 3 3 4 1 1 1 1 2 2 2 3 3 4 1 1 1
  #       1 2 2 2 3 3 4 1 1 1 1 2 2 2 3 3 4
  # 1 2 3 4 2 3 4 3 4 4
	idx.lam <- cbind(rep(kseq, each = nsum), rep(rep(kseq, times = seq(k, 1, -1)), k), unlist(lapply(kseq, function(i) i:k)))
	
	# Build the Q matrix
	# 2 3 4 1 3 4 1 2 4 1 2 3
  # 1 1 1 2 2 2 3 3 3 4 4 4
	idx.q <- cbind(unlist(lapply(kseq, function(i) (kseq)[-i])), rep(kseq, each = k - 1))
	
	# Go through all the lambdas, sum the transition rates for each Q entry
	for (n in seq_len(nsum * k))
		{
		r <- idx.lam[n, ]
		A[r[2], r[1]] <- A[r[2], r[1]] + pars.lam[n]
		A[r[3], r[1]] <- A[r[3], r[1]] + pars.lam[n]
		}
	A[idx.q] <- A[idx.q] + pars.q # Add the Q rates
	diag(A) <- 0
	# Add to diagonals
	diag(A) <- -colSums(A) + unlist(lapply(kseq, function(i) sum(pars.lam[seq((i-1) * nsum + 1, i * nsum)]) - pars.mu[i]))
	A
	}

add_to_diags <- function(i)
	{
	start_lambda_index = (i-1) * nsum + 1
	stop_lambda_index = i * nsum
	lambda_indices = seq(start_lambda_index, stop_lambda_index)
	sum(pars.lam[lambda_indices] - pars.mu[i]
	}
