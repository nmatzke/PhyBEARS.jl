
#######################################################
# Likelihood equation in the birthdeath function
# (derived by pulling apart the birthdeath() function from ape)
# This version stores all of the pieces of the calculation, for comparison
#######################################################
bd_liks <- function(tr, birthRate=1.0, deathRate=0.0)
	{
	ex='
	trtxt = "((((((((P_hawaiiensis_WaikamoiL1:0.9665748366,P_mauiensis_Eke:0.9665748366):0.7086257935,(P_fauriei2:1.231108298,P_hathewayi_1:1.231108298):0.4440923324):0.1767115552,(P_kaduana_PuuKukuiAS:1.851022399,P_mauiensis_PepeAS:1.851022399):0.0008897862802):0.3347375986,P_kaduana_HawaiiLoa:2.186649784):0.302349378,(P_greenwelliae07:1.132253042,P_greenwelliae907:1.132253042):1.35674612):1.689170274,((((P_mariniana_MauiNui:1.99490084,P_hawaiiensis_Makaopuhi:1.99490084):0.7328279804,P_mariniana_Oahu:2.72772882):0.2574151709,P_mariniana_Kokee2:2.985143991):0.4601084855,P_wawraeDL7428:3.445252477):0.732916959):0.7345185743,(P_grandiflora_Kal2:2.480190277,P_hobdyi_Kuia:2.480190277):2.432497733):0.2873119899,((P_hexandra_K1:2.364873976,P_hexandra_M:2.364873976):0.4630447802,P_hexandra_Oahu:2.827918756):2.372081244);"
	tr = read.tree(file="", text=trtxt)
	branching.times(tr)
	birthRate=0.3288164   # ML birthRate for Psychotria tree
	deathRate=0.0					# ML deathRate for Psychotria tree
	bd = bd_liks(tr, birthRate, deathRate)
	bd
	
	trtxt = "((chimp:1,human:1):1,gorilla:2);"
	tr = read.tree(file="", text=trtxt)
	branching.times(tr)
	birthRate=1.0
	deathRate=0.0
	bd = bd_liks(tr, birthRate, deathRate)
	bd

	trtxt = "((chimp:1,human:1):1,gorilla:2);"
	tr = read.tree(file="", text=trtxt)
	branching.times(tr)
	birthRate=1.0
	deathRate=0.999
	bd = bd_liks(tr, birthRate, deathRate)
	bd



	# Getting the birthRate and deathRate from
	# a = deathRate / birthRate	# relative death rate
	# r = birthRate - deathRate	# net diversification rate
	BD =  birthdeath(tr)
	BD
	names(BD)

	# Calculate the birthRate and deathRate from the outputs
	x1 = unname(BD$para["d/b"])
	x2 = unname(BD$para["b-d"])
	deathRate = (x2*x1) / (1-x1)
	birthRate = deathRate+x2
	c(birthRate, deathRate)
	'
	
	
	a = deathRate / birthRate	# relative death rate
	r = birthRate - deathRate	# net diversification rate

	N <- length(tr$tip.label)
	nb_node = tr$Nnode - 1
	sptimes <- c(NA, branching.times(tr)) # NA so the number of times equals number of tips?
	x = sptimes
	# a = "d/b"
	# r = "b-d"

	# dev <- function(a=0.1, r=0.2, N, x, return_deviance=FALSE)
	# 	{
	if (r < 0 || a > 1)
		{
		return(-1e+100)
		}
	
	lnl_topology = lfactorial(tr$Nnode)
	lnl_numBirths = nb_node * log(r)
	lnl_Births_above_root = r * sum(sptimes[3:N])
	
	lnl_numtips_wOneMinusDeathRate = N * log(1 - a)
	# Interpretation: more tips are less likely, if relativeDeathRate is >0
	# If relativeDeathRate = 1, a=0, and lnl=-Inf... 
	#    CLEARLY WRONG EXCEPT IN A MAXIMUM LIKELIHOOD CONTEXT!!!
	# If relativeDeathRate = 0, a=0, and lnl=0, i.e. any number of tips is equiprobable
	
	lnl_branching_times = -2 * sum(log(exp(r * sptimes[2:N]) - a))
	# For each observed branching,
	# netDiversificationRate * timeOfEvent <- take exponent of that ; this means recorded events are less likely in the past
	# <- subtract "a", a constant (relativeDeathRate)
	#
	# This would be a straight likelihood as:
	# 1/
	# (exp(r*branching_time)-a)^2
	#
	# Sum the logs of these
	#
	# If deathRate = 0
	# lnl_branching_times = -2 * sum(log(exp(birthRate*sptimes[2:N]) - 0))
	# lnl_branching_times = -2 * sum(log( exp(birthRate*sptimes[2:N]) )
	# lnl_branching_times = -2 * sum( birthRate*sptimes[2:N] )
	#
	# Note: sum(X) = 9 = total branchlength of tr
	# In BD:
	# -2*sum(sptimes[2:N]) = -12
	# sum(sptimes[3:N]) = 3
	# So, lnl_branching_times + lnl_Births_above_root = yule's -lambda * X
	lnL = lnl_topology + lnl_numBirths + lnl_Births_above_root + lnl_numtips_wOneMinusDeathRate + lnl_branching_times
	dev =  -2 * lnL
	
	bd = NULL
	bd$tr = tr
	bd$birthRate = birthRate
	bd$deathRate = deathRate
	bd$relative_deathRate = a
	bd$net_diversification_rate = r
	bd$dev = dev
	bd$lnl_topology = lnl_topology
	bd$lnl_numBirths = lnl_numBirths
	bd$lnl_Births_above_root = lnl_Births_above_root
	bd$lnl_numtips_wOneMinusDeathRate = lnl_numtips_wOneMinusDeathRate
	bd$lnl_branching_times = lnl_branching_times
	bd$lnL = lnL
	
	return(bd)
	}






# Get total LnL and branch LnL from ClaSSE output
get_classe_LnLs <- function(classe_res)
	{
	
	branch_LnL = sum(attr(classe_res, "intermediates")$lq)
	ttl_LnL = classe_res
	attributes(ttl_LnL) = NULL
	
	return(c(ttl_LnL, branch_LnL))
	}



BGBres_into_classe_params <- function(res, classe_params, birthRate=0.2)
	{
	mats = get_Qmat_COOmat_from_res(res, numstates=ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node), include_null_range=TRUE, max_range_size=res$inputs$max_range_size, timeperiod_i=1)
	numstates = length(mats$states_list)
	include_null_range = res$inputs$include_null_range
	
	# BioGeoBEARS Cevent weights into DF
	Carray_df = get_Cevent_probs_df_from_mats(mats, include_null_range=include_null_range)
	
	# Diversitree params into a df
	lambda_ijk_df = classe_lambdas_to_df(classe_params, k=numstates)
	
	# Put the Cevent_weights * birthRate into diversitree table
	lambda_ijk_df = get_Cevent_lambdas_from_BGB_Carray(lambda_ijk_df, Carray_df, birthRate=birthRate)
	lambda_ijk_df[lambda_ijk_df$lambda != 0,]

	# Convert the diversitree table back to classe_params
	lambdas_to_put_in_params = rownames(lambda_ijk_df[lambda_ijk_df$lambda != 0,])
	indices_in_classe_params = match(x=lambdas_to_put_in_params, table=names(classe_params))
	classe_params[indices_in_classe_params] = lambda_ijk_df[lambdas_to_put_in_params,"lambda"]
	classe_params[153]
	classe_params["lambda020202"]
	classe_params["lambda060203"]
	
	
	# Now do the Qmat
	Qij_df = classe_Qs_to_df(classe_params, k=numstates)
	Qmat = mats$Qmat
	Qij_df = get_Qijs_from_BGB_Qarray(Qij_df, Qmat)
	Qs_to_put_in_params = rownames(Qij_df[Qij_df$q != 0,])
	indices_in_classe_params = match(x=Qs_to_put_in_params, table=names(classe_params))
	classe_params[indices_in_classe_params] = Qij_df[Qs_to_put_in_params,"q"]
	
	
	# And the mus	
	# Because mu (the lineage extinction rate) is always 0.0 in BioGeoBEARS models
	# (i.e., a pure-birth Yule process is being assumed)
	classe_params[grepl(pattern="mu", x=names(classe_params))] = 0.0 
	
	return(classe_params)
	} # BGBres_into_classe_params <- function(res, classe_params, birthRate=0.2)

# k= number of states
# Note that diversitree classe includes only e.g. q123, not q 132
classe_Qs_to_df <- function(classe_params, k=3)
	{
	ex='
	k = 3
	classe_params = c(lambda111 = 0.2, lambda112 = 0, lambda113 = 0, lambda122 = 0, 
lambda123 = 0, lambda133 = 0, lambda211 = 0, lambda212 = 0, lambda213 = 0, 
lambda222 = 0.2, lambda223 = 0, lambda233 = 0, lambda311 = 0, 
lambda312 = 0.0666666666666667, lambda313 = 0.0666666666666667, 
lambda322 = 0, lambda323 = 0.0666666666666667, lambda333 = 0, 
mu1 = 0.1, mu2 = 0.1, mu3 = 0.1, q12 = 0, q13 = 0, q21 = 0, q23 = 0, 
q31 = 0, q32 = 0)
	
	Qij_df = classe_qs_to_df(classe_params=classe_params, k=3)
	Qij_df
	' # end example
	
	# Number of qs per state (e.g., for 3 states, this is 6 qs per state
	#nsum <- k * (k + 1)/2

	# Convert the classe_params qs to a table
	Q_vals = classe_params[grepl(pattern="q", x=names(classe_params))]
	Q_names = names(Q_vals)
	
	# Get i, j, k indices
	Qij_txt = gsub(pattern="q", replacement="", x=Q_names)
	
	if (k <= 9)
		{
		ijs_vector = unlist(sapply(X=Qij_txt, FUN=strsplit, split="", USE.NAMES=FALSE))
		Qij_mat = matrix(as.numeric(ijs_vector), ncol=k, byrow=TRUE)
		Qij_mat
		} else {
		# Split by 2
		is_vec = unlist(sapply(X=Qij_txt, FUN=substr, start=1, stop=2, USE.NAMES=FALSE))
		js_vec = unlist(sapply(X=Qij_txt, FUN=substr, start=3, stop=4, USE.NAMES=FALSE))
		Qij_mat = cbind(is_vec, js_vec)
		}
	
	Qij_df = as.data.frame(cbind(Qij_mat, Q_vals), stringsAsFactors=FALSE)
	names(Qij_df) = c("i", "j", "q")
	# Convert qs to numeric
	Qij_df$q = as.numeric(as.character(Qij_df$q))
	return(Qij_df)
	} # END classe_Qs_to_df <- function(classe_params, k=3)


get_Qijs_from_BGB_Qarray <- function(Qij_df, Qmat)
	{
	for (r in 1:nrow(Qij_df))
		{
		# r = 153
		# Qij_df[153,]
		#               i  j  k lambda
		# lambda020202 02 02 02      0
		i = as.numeric(Qij_df$i[r])
		j = as.numeric(Qij_df$j[r])
	
		# Don't include the diagonals from Q
		if (i == j)
			{
			next()
			}

		Qij_df$q[r] = Qij_df$q[r] + Qmat[i,j]
		} # END for (r in nrow(Qij_df))

	Qij_df[Qij_df$q != 0,]
	return(Qij_df)
	} # END get_Qijs_from_BGB_Qarray <- function(Qij_df, Qmat, birthRate=1.0)


# k= number of states
# Note that diversitree classe includes only e.g. lambda123, not lambda 132
classe_lambdas_to_df <- function(classe_params, k=3)
	{
	ex='
	k = 3
	classe_params = c(lambda111 = 0.2, lambda112 = 0, lambda113 = 0, lambda122 = 0, 
lambda123 = 0, lambda133 = 0, lambda211 = 0, lambda212 = 0, lambda213 = 0, 
lambda222 = 0.2, lambda223 = 0, lambda233 = 0, lambda311 = 0, 
lambda312 = 0.0666666666666667, lambda313 = 0.0666666666666667, 
lambda322 = 0, lambda323 = 0.0666666666666667, lambda333 = 0, 
mu1 = 0.1, mu2 = 0.1, mu3 = 0.1, q12 = 0, q13 = 0, q21 = 0, q23 = 0, 
q31 = 0, q32 = 0)
	
	lambda_ijk_df = classe_lambdas_to_df(classe_params=classe_params, k=3)
	lambda_ijk_df
	' # end example
	
	# Number of lambdas per state (e.g., for 3 states, this is 6 lambdas per state
	nsum <- k * (k + 1)/2

	# Convert the classe_params lambdas to a table
	lambda_vals = classe_params[grepl(pattern="lambda", x=names(classe_params))]
	lambda_names = names(lambda_vals)
	
	# Get i, j, k indices
	lambda_ijk_txt = gsub(pattern="lambda", replacement="", x=lambda_names)
	
	if (k <= 9)
		{
		ijks_vector = unlist(sapply(X=lambda_ijk_txt, FUN=strsplit, split="", USE.NAMES=FALSE))
		lambda_ijk_mat = matrix(as.numeric(ijks_vector), ncol=k, byrow=TRUE)
		lambda_ijk_mat
		} else {
		# Split by 2
		is_vec = unlist(sapply(X=lambda_ijk_txt, FUN=substr, start=1, stop=2, USE.NAMES=FALSE))
		js_vec = unlist(sapply(X=lambda_ijk_txt, FUN=substr, start=3, stop=4, USE.NAMES=FALSE))
		ks_vec = unlist(sapply(X=lambda_ijk_txt, FUN=substr, start=5, stop=6, USE.NAMES=FALSE))
		lambda_ijk_mat = cbind(is_vec, js_vec, ks_vec)
		}
	
	lambda_ijk_df = as.data.frame(cbind(lambda_ijk_mat, lambda_vals), stringsAsFactors=FALSE)
	names(lambda_ijk_df) = c("i", "j", "k", "lambda")
	# Convert lambdas to numeric
	lambda_ijk_df$lambda = as.numeric(as.character(lambda_ijk_df$lambda))
	return(lambda_ijk_df)
	} # END classe_lambdas_to_df <- function(classe_params, k=3)


get_Cevent_lambdas_from_BGB_Carray <- function(lambda_ijk_df, Carray_df, birthRate=1.0)
	{
	for (r in 1:nrow(lambda_ijk_df))
		{
		# r = 153
		# lambda_ijk_df[153,]
		#               i  j  k lambda
		# lambda020202 02 02 02      0
		i = as.numeric(lambda_ijk_df$i[r])
		j = as.numeric(lambda_ijk_df$j[r])
		k = as.numeric(lambda_ijk_df$k[r])
	
		iTF = Carray_df$i == i
		jTF = Carray_df$j == j
		kTF = Carray_df$k == k
		TF = (iTF + jTF + kTF) == 3
		if (sum(TF) == 1)
			{
			lambda_ijk_df$lambda[r] = lambda_ijk_df$lambda[r] + (Carray_df$prob[TF] * birthRate)
			}
	
		# Don't repeat the search when j==k (left and right states are the same)
		if (j == k)
			{
			next()
			}
	
		# Repeat, switching j and k states
		jTF = Carray_df$j == k
		kTF = Carray_df$k == j
		TF = (iTF + jTF + kTF) == 3
		if (sum(TF) == 1)
			{
			lambda_ijk_df$lambda[r] = lambda_ijk_df$lambda[r] + (Carray_df$prob[TF] * birthRate)
			} # END if (sum(TF) == 1)
		} # END for (r in nrow(lambda_ijk_df))

	lambda_ijk_df[lambda_ijk_df$lambda != 0,]
	return(lambda_ijk_df)
	} # END get_Cevent_lambdas_from_BGB_Carray <- function(lambda_ijk_df, Carray_df)



# Calculate the sum of the log computed likelihoods at each node
# i.e., the likelihoods of the speciation events, assuming normalized
# likelihoods at the branch bottoms above each node
# "base" is "t(base)", actually
get_sum_log_computed_likes_at_each_node <- function(tr, base, lq, classe_params)
	{
	# Number of lambdas per state (e.g., for 3 states, this is 6 lambdas per state
	k = ncol(base) / 2
	nsum <- k * (k + 1)/2
	
	Ds_cols = (k+1):(2*k)
	base_likes = apply(X=base[,Ds_cols], MARGIN=2, FUN="*", exp(lq))
	base_normlikes = base_likes / rowSums(base_likes)
	
	# Get a data.frame tabulating the lambdas
	lambda_ijk_df = classe_lambdas_to_df(classe_params=classe_params, k=k)
	lambda_ijk_df$i = as.numeric(as.character(lambda_ijk_df$i))
	lambda_ijk_df$j = as.numeric(as.character(lambda_ijk_df$j))
	lambda_ijk_df$k = as.numeric(as.character(lambda_ijk_df$k))
	lambda_ijk_df$lambda = as.numeric(as.character(lambda_ijk_df$lambda))
	
	
	# Go through the ancestral states
	computed_likelihoods_at_each_node_just_before_speciation = matrix(0.0, nrow=nrow(base), ncol=k)
	
	# Reorder the edge matrix into pruningwise order
	# This is CRUCIAL!!
	tr2 <- reorder(tr, "pruningwise")
	num_internal_nodes = tr$Nnode
	
	# DEFINE DOWNPASS THROUGH THE BRANCHES	
	i = 1
	edges_to_visit = seq(from=1, by=2, length.out=num_internal_nodes)
	
	for (i in edges_to_visit)
		{
		# Get the node numbers at the tips of these two edges		
		j = i+1
		left_desc_nodenum <- tr2$edge[i, 2]
		right_desc_nodenum <- tr2$edge[j, 2]
		# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
		anc <- tr2$edge[i, 1]
		
		# For this anc node, go through the states and sum the likes
		tmp_likes_AT_node = rep(0.0, times=k)
		for (l in 1:k) # l = ancestral state number
			{
			i_TF = lambda_ijk_df$i == l
			j_ne_k_TF = lambda_ijk_df$j != lambda_ijk_df$k
			rows_use_lambda_div2_TF = (i_TF + j_ne_k_TF) == 2
			rows_use_lambda_div1_TF = (i_TF + rows_use_lambda_div2_TF) == 1
			lambda_ijk_df[rows_use_lambda_div1_TF,]
			lambda_ijk_df[rows_use_lambda_div2_TF,]
			
			# Skip e.g. null-range states
			if (sum(rows_use_lambda_div1_TF) > 0)
				{
				# Sum likes where the daughters the same
				ind = rows_use_lambda_div1_TF
				lcol = lambda_ijk_df$j[ind]
				rcol = lambda_ijk_df$k[ind]
				tmp_likes_AT_node[l] = sum(lambda_ijk_df$lambda[ind] * base_normlikes[left_desc_nodenum,lcol] * base_normlikes[right_desc_nodenum,rcol])
				} # END if (sum(rows_use_lambda_div1_TF) > 0)

			if (sum(rows_use_lambda_div2_TF) > 0)
				{
				# Sum likes where the daughters are NOT the same
				# (divide lambda by 2, but use twice)
				ind = rows_use_lambda_div2_TF
				lcol = lambda_ijk_df$j[ind]
				rcol = lambda_ijk_df$k[ind]
				# Left, then right
				tmp_likes_AT_node[l] = tmp_likes_AT_node[l] + sum(lambda_ijk_df$lambda[ind]/2 * base_normlikes[left_desc_nodenum,lcol] * base_normlikes[right_desc_nodenum,rcol])
				tmp_likes_AT_node[l] = tmp_likes_AT_node[l] + sum(lambda_ijk_df$lambda[ind]/2 * base_normlikes[left_desc_nodenum,rcol] * base_normlikes[right_desc_nodenum,lcol])
				} # END if (sum(rows_use_lambda_div1_TF) > 0)
			} # END for (l in 1:k)
		
		computed_likelihoods_at_each_node_just_before_speciation[anc,] = tmp_likes_AT_node
		} # END for (i in edges_to_visit) 
	
	return(computed_likelihoods_at_each_node_just_before_speciation)
	} # END get_sum_log_computed_likes_at_each_node <- function(tr, base, lq, classe_params)



# Convert a "res" object from ClaSSE to a 
# BioGeoBEARS-like set of matrices
claSSE_res_to_prt_OLD <- function(res)
	{
	# Branch top values ("initial")
	init = attr(res,"intermediates")$init
	
	# Branch bottom values ("base")
	base = attr(res,"intermediates")$base
	
	# lqs = log-likelihoods at each branch bottom
	lq = attr(res,"intermediates")$lq
	q = exp(attr(res,"intermediates")$lq)
	
	vals = t(attr(res2, "intermediates")$vals)	# Es and Ds at the root	
	E_indices = 1:nstates
	d_root_orig = vals[-E_indices]							# Original D likelihoods at root

	# Assumes bifurcating tree
	numstates = nrow(init) / 2
	numnodes = ncol(init) # internal plus tip nodes
	numTips = (ncol(init) + 1) / 2
	numInternal = numTips - 1
	
	Es_atNode_branchTop = matrix(data=0, ncol=numstates, nrow=numnodes)
	Es_atNode_branchBot = matrix(data=0, ncol=numstates, nrow=numnodes) 
	likes_at_each_nodeIndex_branchTop = matrix(data=0, ncol=numstates, nrow=numnodes)
	likes_at_each_nodeIndex_branchBot = matrix(data=0, ncol=numstates, nrow=numnodes) 
	normlikes_at_each_nodeIndex_branchTop = matrix(data=0, ncol=numstates, nrow=numnodes)
	normlikes_at_each_nodeIndex_branchBot = matrix(data=0, ncol=numstates, nrow=numnodes)
	
	Ecols = 1:numstates
	Dcols = (numstates+1):(2*numstates)
	Es_atNode_branchTop = (t(init))[,Ecols]
	Es_atNode_branchBot = (t(base))[,Ecols]
	likes_at_each_nodeIndex_branchTop = (t(init))[,Dcols]
	normlikes_at_each_nodeIndex_branchBot = (t(base))[,Dcols]
	normlikes_at_each_nodeIndex_branchTop = likes_at_each_nodeIndex_branchTop / rowSums(likes_at_each_nodeIndex_branchTop)
	likes_at_each_nodeIndex_branchBot = normlikes_at_each_nodeIndex_branchBot * q
	
	# Assign each branch bottom 
	trtable = prt(tr, printflag=FALSE)
	
	Es_atNode_branchTop
	Es_atNode_branchBot
	likes_at_each_nodeIndex_branchTop
	likes_at_each_nodeIndex_branchBot
	normlikes_at_each_nodeIndex_branchTop
	normlikes_at_each_nodeIndex_branchBot
	}





# Convert a "res" object from diversitree 
# ClaSSE to a 
# BioGeoBEARS-like set of matrices
claSSE_res_to_prt <- function(res, tr, classe_params)
	{
	ex='
	res1t = structure(-6.22837508651605, intermediates = list(init = structure(c(0, 
0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0.153452947030521, 
0.153452947030521, 0.153452947030521, 0, 0, 0.0326226539925451, 
0.0868935659418745, 0.0868935659418745, 0.0868935659418745, 0, 
0, 0.0333321365204794), .Dim = 6:5), base = structure(c(0.0868935659418745, 
0.0868935659418745, 0.0868935659418745, 0, 0.994007973162888, 
0.00599202683711178, 0.0868935659418746, 0.0868935659418746, 
0.0868935659418745, 0.994007973162888, 0, 0.00599202683711177, 
0.153452946974297, 0.153452946974297, 0.153452946974297, 0, 0.978679619776353, 
0.0213203802236467, NA, NA, NA, NA, NA, NA, 0.153452947030521, 
0.153452947030521, 0.153452947030521, 0, 0, 1), .Dim = 6:5), 
    lq = c(-0.275795607421526, -0.275795607421526, -0.511628046092027, 
    0, -3.68502440109251), vals = c(0.153452947030521, 0.153452947030521, 
    0.153452947030521, 0, 0, 0.0326226539925451), branchLnL = -4.74824366202759, 
    root.p = c(0, 0, 1)), vals = c(0.153452947030521, 0.153452947030521, 
0.153452947030521, 0, 0, 0.0326226539925451))
	res = res1t
	
	trstr = "((chimp:1,human:1):1,gorilla:2);"
	tr = read.tree(file="", text=trstr)
	
	k = 3
	classe_params = c(lambda111 = 0.2, lambda112 = 0, lambda113 = 0, lambda122 = 0, 
lambda123 = 0, lambda133 = 0, lambda211 = 0, lambda212 = 0, lambda213 = 0, 
lambda222 = 0.2, lambda223 = 0, lambda233 = 0, lambda311 = 0, 
lambda312 = 0.0666666666666667, lambda313 = 0.0666666666666667, 
lambda322 = 0, lambda323 = 0.0666666666666667, lambda333 = 0, 
mu1 = 0.1, mu2 = 0.1, mu3 = 0.1, q12 = 0, q13 = 0, q21 = 0, q23 = 0, 
q31 = 0, q32 = 0)
	
	lambda_ijk_df = classe_lambdas_to_df(classe_params=classe_params, k=3)
	lambda_ijk_df
	
	classe_res_dfs = claSSE_res_to_prt(res, tr, classe_params)
	classe_res_dfs
	names(classe_res_dfs)
	' # END ex
	
	# Branch top values ("initial")
	init = t(attr(res,"intermediates")$init)
	
	# Branch bottom values ("base")
	base = t(attr(res,"intermediates")$base)
	
	# lqs = log-likelihoods at each branch bottom
	lq = attr(res,"intermediates")$lq
	q = exp(attr(res,"intermediates")$lq)
	
	vals = t(attr(res2, "intermediates")$vals)	# Es and Ds at the root	
	E_indices = 1:k
	d_root_orig = vals[-E_indices]							# Original D likelihoods at root

	# Assumes bifurcating tree
	numstates = ncol(init) / 2
	rootnode = length(tr$tip.label) + 1
	numnodes = nrow(init) # internal plus tip nodes
#	numTips = (nrow(init) + 1) / 2
	numTips = length(tr$tip.label)
	numInternal = numTips - 1
	
	# Intitializing
	Es_atNode_branchTop = matrix(data=0, ncol=numstates, nrow=numnodes)
	Es_atNode_branchBot = matrix(data=0, ncol=numstates, nrow=numnodes) 
	likes_at_each_nodeIndex_branchTop = matrix(data=0, ncol=numstates, nrow=numnodes)
	likes_at_each_nodeIndex_branchBot = matrix(data=0, ncol=numstates, nrow=numnodes) 
	normlikes_at_each_nodeIndex_branchTop = matrix(data=0, ncol=numstates, nrow=numnodes)
	normlikes_at_each_nodeIndex_branchBot = matrix(data=0, ncol=numstates, nrow=numnodes)
	
	# Filling in
	Ecols = 1:numstates
	Dcols = (numstates+1):(2*numstates)
	Es_atNode_branchTop = init[,Ecols]
	Es_atNode_branchBot = base[,Ecols]
	likes_at_each_nodeIndex_branchTop = init[,Dcols]
	normlikes_at_each_nodeIndex_branchBot = base[,Dcols]
	normlikes_at_each_nodeIndex_branchTop = likes_at_each_nodeIndex_branchTop / rowSums(likes_at_each_nodeIndex_branchTop)
	likes_at_each_nodeIndex_branchBot = normlikes_at_each_nodeIndex_branchBot * q
	
	# Likelihoods just after nodeOp, before normalization
	computed_likelihoods_at_each_node_just_before_speciation = get_sum_log_computed_likes_at_each_node(tr, base, lq, classe_params)
	log(rowSums(computed_likelihoods_at_each_node_just_before_speciation))
	TF = is.finite(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation)))
	sum_log_nodelikes = sum(lq) + sum(log(rowSums(computed_likelihoods_at_each_node_just_before_speciation))[TF])
sum_log_nodelikes

	# LnLs1
	# If root=ROOT.OBS, root.p=NULL, condition.surv=FALSE
	d_root_orig = likes_at_each_nodeIndex_branchTop[rootnode,]						# Original D likelihoods at root
	root.p = d_root_orig/sum(d_root_orig)
	loglik_s1 = log(sum(root.p * d_root_orig)) + sum(lq)
	loglik_s1

	# LnLs1t
	# If root=ROOT.OBS, root.p=NULL, condition.surv=TRUE
	root.p = d_root_orig/sum(d_root_orig)
	# MuSSE/ClaSSE
	pars = classe_params
	nsum <- k * (k + 1)/2
	lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
	i <- seq_len(k)
	e.root <- Es_atNode_branchTop[rootnode,]
	d.root <- d_root_orig/sum(root.p * lambda * (1 - e.root)^2)
	loglik_s1t = log(sum(root.p * d.root)) + sum(lq)
	loglik_s1t

	R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = sum_log_nodelikes


# R_result_branch_lnL = -4.748244
# R_result_total_LnLs1 = -8.170992
# R_result_total_LnLs1t = -6.228375
# R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = -11.57223

	R_result_branch_lnL = sum(lq)
	R_result_total_LnLs1 = loglik_s1
	R_result_total_LnLs1t = loglik_s1t
	R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = sum_log_nodelikes
	
	
	classe_res_dfs = NULL
	classe_res_dfs$init = t(attr(res,"intermediates")$init)
	classe_res_dfs$base = t(attr(res,"intermediates")$init)
	classe_res_dfs$lq = lq
	classe_res_dfs$Es_atNode_branchTop = Es_atNode_branchTop
	classe_res_dfs$Es_atNode_branchBot = Es_atNode_branchBot
	classe_res_dfs$likes_at_each_nodeIndex_branchTop = likes_at_each_nodeIndex_branchTop
	classe_res_dfs$likes_at_each_nodeIndex_branchBot = likes_at_each_nodeIndex_branchBot
	classe_res_dfs$normlikes_at_each_nodeIndex_branchTop = normlikes_at_each_nodeIndex_branchTop
	classe_res_dfs$normlikes_at_each_nodeIndex_branchBot = normlikes_at_each_nodeIndex_branchBot
	classe_res_dfs$computed_likelihoods_at_each_node_just_before_speciation = computed_likelihoods_at_each_node_just_before_speciation
	classe_res_dfs$R_result_branch_lnL = R_result_branch_lnL
	classe_res_dfs$R_result_total_LnLs1 = R_result_total_LnLs1
	classe_res_dfs$R_result_total_LnLs1t = R_result_total_LnLs1t
	classe_res_dfs$R_result_sum_log_computed_likelihoods_at_each_node_x_lambda = R_result_sum_log_computed_likelihoods_at_each_node_x_lambda
	
	names(classe_res_dfs)
	return(classe_res_dfs)
	} # END claSSE_res_to_prt <- function(res, tr, classe_params)





#diversitree:::all.branches.matrix
all.branches.matrix <- function(pars, cache, initial.conditions, branches, preset = NULL) 
{
    len <- cache$len
    depth <- cache$depth
    children <- cache$children
    #order <- cache$order[-length(cache$order)]
    order <- cache$order
    root <- cache$root
    n <- length(len)
    lq <- rep(0, n)
    n.tip <- cache$n.tip
    y <- cache$y
    branch.init <- branch.base <- matrix(NA, cache$info$ny, n)
    if (!is.null(preset)) {
        lq[preset$target] <- preset$lq
        branch.base[, preset$target] <- preset$base
    }
    if (is.null(names(y))) {
        for (x in y) {
            if (!is.null(x)) {
                idx <- x$target
                branch.init[, idx] <- x$y
                ans <- branches(x$y, x$t, pars, 0, idx)
                lq[idx] <- ans[[1]]
                branch.base[, idx] <- ans[[2]]
            }
        }
    }
    else {
        tip.t <- y$t
        tip.target <- y$target
        tip.y <- branch.init[, tip.target] <- y$y
        for (i in seq_along(tip.t)) {
            idx <- tip.target[i]
            ans <- branches(tip.y[, i], tip.t[i], pars, 0, idx)
            lq[idx] <- ans[[1]]
            branch.base[, idx] <- ans[[2]]
        }
    }
    
    # NJM ADD additional things to log
    ans_lists = NULL
    branches_function_list = NULL
    countval = 0
    for (i in order) {
    	# Count so you can treat the root specially
    	countval = countval + 1
        y.in <- initial.conditions(branch.base[, children[i, ]], pars, depth[i], i)
    	
    	# Normalize to sum to 1, if you are at the root state
    	# THIS CONVERTS LIKELIHOODS TO ROOT.OBS
    	#if (countval != length(cache$order))
    	#	{
    	#	y.in = y.in / sum(y.in)
	    #   }
   		#y.in = y.in / sum(y.in)
        if (!all(is.finite(y.in))) 
            stop("Bad initial conditions: calculation failure along branches?")
        branch.init[, i] <- y.in
        ans <- branches(y.in, len[i], pars, depth[i], i)
        lq[i] <- ans[[1]]
        branch.base[, i] <- ans[[2]]

	    # NJM ADD additional things to log
	    ans_lists[[i]] = ans
        branches_function_list[[i]] = body(branches)
    }
    y.in <- initial.conditions(branch.base[, children[root, ]], 
        pars, depth[root], root)
    branch.init[, root] <- y.in
    
    # This sets:
    # 
    # $init
    # $lq
    # $vals
    # $ans_lists
    # 
    list(init = branch.init, base = branch.base, lq = lq, vals = y.in, ans_lists=ans_lists, branches_function_list=branches_function_list)
}







#######################################################
# Convert the log-likelihoods from DEC, DEC+J etc.
# (from the Lagrange or BioGeoBEARS programs) into 
# the equivalent lnL from Diversitree's ClaSSE.
# 
# These log-likelihoods will be obtained from the 
# ClaSSE settings below, assuming you have the same
# phylogeny and tip data between BioGeoBEARS and ClaSSE.
#
# (Note that tip data, i.e. geographic ranges, have to be 
#  input differently between diversitree and BioGeoBEARS)
# 
# Under these settings (basically, equal root-state frequencies)
# the *differences* in log-likelihood will be identical between:
#
# 1. BioGeoBEARS DEC and DEC+J
# 2. Diversitree ClaSSE-DEC and ClaSSE-DEC+J
# 
# Exact matches can be obtained with the 
# classe settings below.
# 
# (see "compare_BGB_diversitree_DEC_v1.R" for a demonstration)
#
#
# Specifically, these setups in the script:
# 
##########################
# res2: root=ROOT.FLAT, root.p=NULL, condition.surv=FALSE
##########################
tmptxt='
res2 = classe_Kstates(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
' # END tmptxt

##########################
# res3: root=ROOT.GIVEN, root.p=root_probs, condition.surv=FALSE
##########################
# All probabilities equal, except null range has prob=0
tmptxt='
root_probs_equal = rep(1, times=numstates)
root_probs_equal[sum(include_null_range)] = 0
root_probs_equal = root_probs_equal / sum(root_probs_equal)
root_probs = root_probs_equal

res3 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
' # END tmptxt

##########################
# res5: root=ROOT.GIVEN, root.p=root_probs, condition.surv=FALSE
##########################
# All states, except null range, get "probability" 1
# (i.e., ignore root state frequencies, like DEC-type models)
tmptxt='
root_probs_single = rep(1, times=numstates)
root_probs_single[sum(include_null_range)] = 0

res5 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
' # END tmptxt

######################################################
# Exact matches can also be obtained if condition.surv=TRUE,
# on the above settings
######################################################
# 
# These are the meanings of the different "res1", "res1t" etc. abbreviations:
# 
# res1 = classe_Kstates(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
# res2 = classe_Kstates(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
# root_probs = root_probs_equal
# res3 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
# root_probs = root_probs_biased
# res4 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
# root_probs = root_probs_single
# res5 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
# res6 = classe_Kstates(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
# 
# res1t = classe_Kstates(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# res2t = classe_Kstates(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# root_probs = root_probs_equal
# res3t = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
# root_probs = root_probs_biased
# res4t = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
# root_probs = root_probs_single
# res5t = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
# res6t = classe_Kstates(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# 

convert_BGB_lnL_to_ClaSSE <- function(res, tr=NULL, root_probs_biased=NULL)
	{
	# Load tree from stored trfn (tree filename), if needed
	if (is.null(tr) == TRUE)
		{
		tr = try(read.tree(res$inputs$trfn))
		
		# Error trap
		if (class(tr) == "try-error")
			{
			txt = paste0("STOP ERROR in convert_BGB_lnL_to_ClaSSE(): the tree file '", es$inputs$trfn, "' could not be read into an R 'phylo' object by ape's read.tree(). Try providing the read-in tree directly through the 'tr=' argument.")
			
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if (class(tr) == "try-error")
		} # END if (is.null(tr) == TRUE)
	
	
	# Error trap: tree must be ultrametric
	if (is.ultrametric(tr) == FALSE)
		{
		txt = "STOP ERROR in convert_BGB_lnL_to_ClaSSE(): is.ultrametric(tr) returned FALSE. convert_BGB_lnL_to_ClaSSE()'s calculations only make sense if you have an ultrametric tree, i.e. all of the tips survive to the present."
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (is.ultrametric(tr) == FALSE)
	
	
	#######################################################
	# Matching diversitree and BioGeoBEARS
	#######################################################
	# BioGeoBEARS stores the likelihoods at each node, 
	# including the root node.
	#
	# But diversitree stores the likelihoods at branch bottoms,
	# and treats the root node differently, depending on 
	# various user options.

	# Get basic info from BioGeoBEARS res object:
	include_null_range = res$inputs$include_null_range
	numstates = ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)

	
	
	# A. First, we will need the tree likelihood under a pure-birth (yule) process,
	# as well as the Maximum Likelihood estimate of the birthRate
	
	# birthRate = yule(tr)$lambda  # The ML lambda from Yule. 
	                             # Equals (#speciations-1)/tree_length
	birthRate = (tr$Nnode-1)/sum(tr$edge.length)
	deathRate = 0
	
	# Get the pieces of the birthdeath log-likelihood
	bd_ape = bd_liks(tr, birthRate=birthRate, deathRate=deathRate)
	
	
	# B. Let's start by summing the BioGeoBEARS likelihoods
	# but exclude the root node.
	root_nodenum = length(tr$tip.label) + 1
	sumBGBlike_not_root = sum(log(res$computed_likelihoods_at_each_node[-root_nodenum]))
	sumBGBlike_not_root
	
	# C. Let's take the sum of the branch-bottom likelihoods from the birth-death
	# process.
	# We could do this with the "bd_liks" function
	# (which takes the lnL pieces from ape::birthdeath()
	tmptxt='
	bd_ape = bd_liks(tr, birthRate=birthRate, deathRate=deathRate)
	bd_ape$lnl_numBirths + bd_ape$lnl_Births_above_root + bd_ape$lnl_branching_times
  bd_ape$lnl_numBirths + -(tr$Nnode-1)
  '
  # But, we can do it from scratch:
  sum_branchBot_likes = sum(-birthRate * trtable$edge.length, na.rm=TRUE) + (tr$Nnode-1)*log(birthRate) 
  
  # Add the lnL of root speciation event, -1 for extra node
  1-log(1/birthRate)
  -(log(1/birthRate) - 1)
  
  # D. Put them all together to get the "branch_lnL" -- matches the
  # "branch_lnLs", i.e. sum(lq) from diversitree's claSSE
	branch_lnL = sumBGBlike_not_root + sum_branchBot_likes + (1-log(1/birthRate))
  
  # This matches the "branch lnLs" from diversitree claSSE:
  # lq is the Ds log-likelihood summed at each branch bottom, and extracted
  # 
  tmptxt = '
  res1 = classe_Kstates(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
  init = t(attr(res1, "intermediates")$init)
	base = t(attr(res1, "intermediates")$base)
  Ds_cols = (numstates+1):(2*numstates)
	lq = t(attr(res1, "intermediates")$lq)
  log(rowSums(base[,Ds_cols] * exp(c(lq))))
	' # END tmptxt
	
	# res0: Adding BioGeoBEARS likelihood to a Yule process
	#       +SFs (state frequencies/base frequencies
	bd_ape = bd_liks(tr, birthRate=birthRate, deathRate=deathRate)
	bgb_plus_Yule = bgb1 + bd_ape$lnL
	bgb_plus_Yule_minus_topo = bgb1 + bd_ape$lnL - bd_ape$lnl_topology
	bgb_plus_Yule_minus_topo2 = bgb1 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate)
	# this last matches claSSE branch_lnLs
	
	
  # To get to the res1_lnL, we can also add the root state likelihoods:
	BGBlnL_at_root = log(res$computed_likelihoods_at_each_node[root_nodenum]) - 1
	d_root_orig_BGB = res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[root_nodenum,] * exp(BGBlnL_at_root)
	d_root_orig_BGB
	sum(d_root_orig_BGB)

  tmptxt = '
	vals = t(attr(res1, "intermediates")$vals)	# Es and Ds at the root
	E_indices = 1:numstates
	d_root_orig_diversitree = vals[-E_indices]
	d_root_orig_diversitree
	sum(d_root_orig_diversitree)
	' # END tmptxt

	#####################################
  # 1. res1 for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	root.p = d_root_orig_BGB/sum(d_root_orig_BGB)
	rootlikes = log(sum(root.p * d_root_orig_BGB))
	res1_rootlikes = rootlikes

	claSSE_res1_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
  

	#####################################
  # 2. res2 for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# Match BioGeoBEARS to diverstree res2 (ROOT.FLAT state frequencies)
	root.p = rep(1/numstates, times=numstates)
	rootlikes = log(sum(root.p * d_root_orig_BGB))
	res2_rootlikes = rootlikes


	claSSE_res2_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes


	#####################################
  # 3. res3 for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# Match BioGeoBEARS to diverstree res3 (all equal, except null range -- ROOT.GIVEN)
	# Set up various assumptions about the root state probabilities
	# All probabilities equal, except null range has prob=0
	root_probs_equal = rep(1, times=numstates)
	root_probs_equal[sum(include_null_range)] = 0
	root_probs_equal = root_probs_equal / sum(root_probs_equal)
	root.p = root_probs_equal
	rootlikes = log(sum(root.p * d_root_orig_BGB))
	res3_rootlikes = rootlikes

	claSSE_res3_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes

	#####################################
  # 4. res4 for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# SKIP? Not really relevant (requires unequal root state probabilities)

	# Root probs highly biased towards the last state
	if (is.null(root_probs_biased) == TRUE)
		{
		root_probs_biased = rep(0.01, times=numstates)
		root_probs_biased[sum(include_null_range)] = 0
		root_probs_biased[length(root_probs_biased)] = 0.01 * (numstates-include_null_range)
		root_probs_biased = root_probs_biased / sum(root_probs_biased)
		}

	root_probs = root_probs_biased
	root.p = root_probs
	rootlikes = log(sum(root.p * d_root_orig_BGB))
	res4_rootlikes = rootlikes

	claSSE_res4_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes


	#####################################
  # 5. res5 for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# Match BioGeoBEARS to diverstree res5 (all 1s, except null range -- ROOT.GIVEN)
	# All states, except null range, get "probability" 1
	# (i.e., ignore root state frequencies, like DEC-type models)
	root_probs_single = rep(1, times=numstates)
	root_probs_single[sum(include_null_range)] = 0
	root.p = root_probs_single
	rootlikes = log(sum(root.p * d_root_orig_BGB))
	res5_rootlikes = rootlikes

	claSSE_res5_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes

	#####################################
  # 6. res6 for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# Match BioGeoBEARS to diverstree res6 (equilibrium frequencies for root states)
	# SKIP: Requires classe inputs and functions
	# res6: Equilibrium root frequencies
	# If root=ROOT.EQUI, condition.surv=FALSE

	# We have to extract the classe_params
	classe_params = BGBres_construct_classe_states_plus_params(res, tr=tr)
	
	# Project the ClaSSE model onto an instantaneous rate matrix, A
	k = numstates
	A = projection.matrix.classe(pars=classe_params, k) 

	# Calculate equilibrium frequencies by eigenvectors
	evA <- eigen(A)
	i <- which(evA$values == max(evA$values))
	equilibrium_root_freqs = evA$vectors[, i]/sum(evA$vectors[, i])
	
	rootlikes = log(sum(equilibrium_root_freqs * d_root_orig_BGB))
	res6_rootlikes = rootlikes
	
	claSSE_res6_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
	

	#######################################################
	# For the condition on survival=TRUE condition, we 
	# multiply the d_roots * root.p * lambdas * (1-e.root)^2
	# to get new d_roots.
	#
	# However, for BioGeoBEARS models, all lambdas=birthRate, and
	# all e.roots are 0.  So we are just multiplying by
	# birthRate
	#######################################################

	#####################################
  # 1. res1t for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	root.p = d_root_orig_BGB/sum(d_root_orig_BGB)
	d_root_orig_BGB_condsurv = d_root_orig_BGB / sum(root.p * birthRate)
	rootlikes = log(sum(root.p * d_root_orig_BGB_condsurv))
	res1t_rootlikes = rootlikes

	claSSE_res1t_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
  claSSE_res1t_lnL

	#####################################
  # 2. res2t for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# Match BioGeoBEARS to diverstree res2 (ROOT.FLAT state frequencies)
	root.p = rep(1/numstates, times=numstates)
	root.p[sum(include_null_range)] = 0
	root.p = root.p / sum(root.p)
	d_root_orig_BGB_condsurv = d_root_orig_BGB / sum(root.p * birthRate)
	rootlikes = log(sum(root.p * d_root_orig_BGB_condsurv))
	res2t_rootlikes = rootlikes

	claSSE_res2t_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
	claSSE_res2t_lnL

	#####################################
  # 3. res3t for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# Match BioGeoBEARS to diverstree res3t (all equal, except null range -- ROOT.GIVEN)
	# Set up various assumptions about the root state probabilities
	# All probabilities equal, except null range has prob=0
	root_probs_equal = rep(1, times=numstates)
	root_probs_equal[sum(include_null_range)] = 0
	root_probs_equal = root_probs_equal / sum(root_probs_equal)
	root.p = root_probs_equal
	d_root_orig_BGB_condsurv = d_root_orig_BGB / sum(root.p * birthRate)
	rootlikes = log(sum(root.p * d_root_orig_BGB_condsurv))
	res3t_rootlikes = rootlikes

	claSSE_res3t_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
	claSSE_res3t_lnL

	#####################################
  # 4. res4t for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# SKIP? Not really relevant (requires unequal root state probabilities)

	# Root probs highly biased towards the last state
	if (is.null(root_probs_biased) == TRUE)
		{
		root_probs_biased = rep(0.01, times=numstates)
		root_probs_biased[sum(include_null_range)] = 0
		root_probs_biased[length(root_probs_biased)] = 0.01 * (numstates-include_null_range)
		root_probs_biased = root_probs_biased / sum(root_probs_biased)
		}

	root_probs = root_probs_biased
	root.p = root_probs
	d_root_orig_BGB_condsurv = d_root_orig_BGB / sum(root.p * birthRate)
	rootlikes = log(sum(root.p * d_root_orig_BGB_condsurv))
	res4t_rootlikes = rootlikes

	claSSE_res4t_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
	claSSE_res4t_lnL

	#####################################
  # 5. res5t for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# Match BioGeoBEARS to diverstree res5t (all 1s, except null range -- ROOT.GIVEN)
	# All states, except null range, get "probability" 1
	# (i.e., ignore root state frequencies, like DEC-type models)
	root_probs_single = rep(1, times=numstates)
	root_probs_single[sum(include_null_range)] = 0
	root.p = root_probs_single
	d_root_orig_BGB_condsurv = d_root_orig_BGB / sum(root.p * birthRate)
	rootlikes = log(sum(root.p * d_root_orig_BGB_condsurv))
	res5t_rootlikes = rootlikes

	claSSE_res5t_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes


	#####################################
  # 6. res6t for diversitree claSSE, from BioGeoBEARS starting point
	#####################################
	# Match BioGeoBEARS to diverstree res6t (equilibrium frequencies for root states)
	# SKIP: Requires classe inputs and functions
	# res6: Equilibrium root frequencies
	# If root=ROOT.EQUI, condition.surv=FALSE

	# We have to extract the classe_params
	classe_params = BGBres_construct_classe_states_plus_params(res, tr=tr)
	
	# Project the ClaSSE model onto an instantaneous rate matrix, A
	k = numstates
	A = projection.matrix.classe(pars=classe_params, k) 

	# Calculate equilibrium frequencies by eigenvectors
	evA <- eigen(A)
	i <- which(evA$values == max(evA$values))
	equilibrium_root_freqs = evA$vectors[, i]/sum(evA$vectors[, i])
	equilibrium_root_freqs_noNull = equilibrium_root_freqs
	equilibrium_root_freqs_noNull[sum(include_null_range)] = 0
	
	d_root_orig_BGB_condsurv = d_root_orig_BGB / sum(equilibrium_root_freqs_noNull * birthRate)
	rootlikes = log(sum(equilibrium_root_freqs_noNull * d_root_orig_BGB_condsurv))
	res6t_rootlikes = rootlikes
	
	claSSE_res6t_lnL = sumBGBlike_not_root + sum_branchBot_likes - (log(1/birthRate) - 1) + rootlikes
	claSSE_res6t_lnL
	
	
	# Assemble results
	row_names = c("res1", "res2", "res3", "res4", "res5", "res6", "res1t", "res2t", "res3t", "res4t", "res5t", "res6t")
	
	BGB_lnLs = rep(res$total_loglikelihood, times=length(row_names))
	branch_lnLs = rep(branch_lnL, times=length(row_names))
	sumBGBlnL_notRoots = rep(sumBGBlike_not_root, times=length(row_names))
	sumBranchBot_lnL = rep(sum_branchBot_likes, times=length(row_names))
	one_minus_lnBirthRate = rep((1-log(1/birthRate)), times=length(row_names))
	ln_rootlikes = c(res1_rootlikes, res2_rootlikes, res3_rootlikes, res4_rootlikes, res5_rootlikes, res6_rootlikes, res1t_rootlikes, res2t_rootlikes, res3t_rootlikes, res4t_rootlikes, res5t_rootlikes, res6t_rootlikes)
	BGB_to_classe_lnLs = c(claSSE_res1_lnL, claSSE_res2_lnL, claSSE_res3_lnL, claSSE_res4_lnL, claSSE_res5_lnL, claSSE_res6_lnL, claSSE_res1t_lnL, claSSE_res2t_lnL, claSSE_res3t_lnL, claSSE_res4t_lnL, claSSE_res5t_lnL, claSSE_res6t_lnL)

	resmat = cbind(BGB_lnLs, branch_lnLs, sumBGBlnL_notRoots, sumBranchBot_lnL, one_minus_lnBirthRate, ln_rootlikes, BGB_to_classe_lnLs)
	resdf = adf(resmat)
	row.names(resdf) = row_names
	resdf
	} # END convert_BGB_lnL_to_ClaSSE <- function(res, tr=NULL, root_probs_biased=NULL)


















#######################################################
# Take a BioGeoBEARS results (res) object and return
# the classe_params
#######################################################
BGBres_construct_classe_states_plus_params <- function(res, tr=NULL)
	{

runtxt = '
# Load simple example tree (newick format, must be ultrametric, i.e. 
# all the tips come to the present)
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
tr = read.tree(trfn)

# Load geography data
geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))

library(BioGeoBEARS)
max_range_size = 4
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isnt much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
BioGeoBEARS_run_object$printlevel = 0					 # skip printing to screen
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
include_null_range = BioGeoBEARS_run_object$include_null_range

# Run the Maximum Likelihood optimization
res = bears_optim_run(BioGeoBEARS_run_object)
classe_params = BGBres_construct_classe_states_plus_params(res, tr=tr)
' # END runtxt



	# Set up a diversitree ClaSSE model from BioGeoBEARS

	# Load tree from stored trfn (tree filename), if needed
	if (is.null(tr) == TRUE)
		{
		tr = try(read.tree(res$inputs$trfn))
		
		# Error trap
		if (class(tr) == "try-error")
			{
			txt = paste0("STOP ERROR in convert_BGB_lnL_to_ClaSSE(): the tree file '", es$inputs$trfn, "' could not be read into an R 'phylo' object by ape's read.tree(). Try providing the read-in tree directly through the 'tr=' argument.")
			
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if (class(tr) == "try-error")
		} # END if (is.null(tr) == TRUE)
	
	
	# Error trap: tree must be ultrametric
	if (is.ultrametric(tr) == FALSE)
		{
		txt = "STOP ERROR in convert_BGB_lnL_to_ClaSSE(): is.ultrametric(tr) returned FALSE. convert_BGB_lnL_to_ClaSSE()'s calculations only make sense if you have an ultrametric tree, i.e. all of the tips survive to the present."
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (is.ultrametric(tr) == FALSE)

	
	# BioGeoBEARS assumes a Yule birthRate
	# birthRate = yule(tr)$lambda  # The ML lambda from Yule. 
	                             # Equals (#speciations-1)/tree_length
	birthRate = (tr$Nnode-1)/sum(tr$edge.length)
	deathRate = 0


	# Get the tip statenums
	numtips = length(tr$tip.label)
	numstates = ncol(res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
	tip_statenums = rep(0, times=numtips)
	for (i in 1:numtips)
		{ # Find the "1" (the observed state for each tip)
		TF = res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[i,] == 1
		tip_statenums[i] = (1:numstates)[TF]
		}
	tip_statenums
	names(tip_statenums) = tr$tip.label
	states = tip_statenums

	# Set the sampling rate to 1 for each state
	sampling.f = rep(1, times=numstates)
	k = numstates

	# Create the ClaSSE likelihood function for k states
	# (strict=FALSE means that some states in the state space can be absent from the tips)
	classe_Kstates = make.classe(tree=tr, states=states, k=k, sampling.f=sampling.f, strict=FALSE)

	# The names of the ClaSSE parameters:
	# Note that diverstree creates ALL the possible parameters, which gets
	# ridiculous quickly, e.g. 
	# 4 areas = 16 geographic range states = 2432 parameters in param_names
	param_names = argnames(classe_Kstates)
	length(param_names)
	length(param_names[grepl(pattern="lambda", x=param_names)]) # 2176 speciation rates
	length(param_names[grepl(pattern="mu", x=param_names)])     #   16 extinction rates
	length(param_names[grepl(pattern="q", x=param_names)])      #  240 Q transition rates

	# Most parameters will be zero
	classe_params = rep(0, times=length(param_names))
	names(classe_params) = param_names
	head(classe_params)
	tail(classe_params)

	# Make a data.frame to match up with the BioGeoBEARS Carray_df
	lambda_ijk_df = classe_lambdas_to_df(classe_params, k=numstates)
	head(lambda_ijk_df)

	# Fill in the params from the BioGeoBEARS "res" results
	classe_params = BGBres_into_classe_params(res, classe_params, birthRate=birthRate)

	return(classe_params)
	} # END BGBres_construct_classe_states_plus_params <- function(res, tr=NULL)








#######################################################
# convert_BGB_to_BGB_Yule_SFs
# 
# This is a simpler version, that just presents the 
# short and simple version of these conversions from
# BioGeoBEARS likelihoods to claSSE likelihoods.
# 
# Convert the log-likelihoods from DEC, DEC+J etc.
# (from the Lagrange or BioGeoBEARS programs) into 
# the equivalent lnL from Diversitree's ClaSSE.
# 
# These log-likelihoods will be obtained from the 
# ClaSSE settings below, assuming you have the same
# phylogeny and tip data between BioGeoBEARS and ClaSSE.
#
# (Note that tip data, i.e. geographic ranges, have to be 
#  input differently between diversitree and BioGeoBEARS)
# 
# Under these settings (basically, equal root-state frequencies)
# the *differences* in log-likelihood will be identical between:
#
# 1. BioGeoBEARS DEC and DEC+J
# 2. Diversitree ClaSSE-DEC and ClaSSE-DEC+J
# 
# Exact matches can be obtained with the 
# classe settings below.
# 
# (see "compare_BGB_diversitree_DEC_v1.R" for a demonstration)
#
#
# Specifically, these setups in the script:
# 
##########################
# res2: root=ROOT.FLAT, root.p=NULL, condition.surv=FALSE
##########################
tmptxt='
res2 = classe_Kstates(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
' # END tmptxt

##########################
# res3: root=ROOT.GIVEN, root.p=root_probs, condition.surv=FALSE
##########################
# All probabilities equal, except null range has prob=0
tmptxt='
root_probs_equal = rep(1, times=numstates)
root_probs_equal[sum(include_null_range)] = 0
root_probs_equal = root_probs_equal / sum(root_probs_equal)
root_probs = root_probs_equal

res3 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
' # END tmptxt

##########################
# res5: root=ROOT.GIVEN, root.p=root_probs, condition.surv=FALSE
##########################
# All states, except null range, get "probability" 1
# (i.e., ignore root state frequencies, like DEC-type models)
tmptxt='
root_probs_single = rep(1, times=numstates)
root_probs_single[sum(include_null_range)] = 0

res5 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
' # END tmptxt

######################################################
# Exact matches can also be obtained if condition.surv=TRUE,
# on the above settings
######################################################
# 
# These are the meanings of the different "res1", "res1t" etc. abbreviations:
# 
# res1 = classe_Kstates(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
# res2 = classe_Kstates(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
# root_probs = root_probs_equal
# res3 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
# root_probs = root_probs_biased
# res4 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
# root_probs = root_probs_single
# res5 = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=FALSE)
# res6 = classe_Kstates(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
# 
# res1t = classe_Kstates(pars=classe_params, root=ROOT.OBS, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# res2t = classe_Kstates(pars=classe_params, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# root_probs = root_probs_equal
# res3t = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
# root_probs = root_probs_biased
# res4t = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
# root_probs = root_probs_single
# res5t = classe_Kstates(pars=classe_params, root=ROOT.GIVEN, root.p=root_probs, intermediates=TRUE, condition.surv=TRUE)
# res6t = classe_Kstates(pars=classe_params, root=ROOT.EQUI, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# 

convert_BGB_to_BGB_Yule_SFs <- function(res, tr=NULL)
	{
	runtxt='
# Load simple example tree (newick format, must be ultrametric, i.e. 
# all the tips come to the present)
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
tr = read.tree(trfn)

# Load geography data
geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))

library(BioGeoBEARS)
max_range_size = 4
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isnt much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
BioGeoBEARS_run_object$printlevel = 0					 # skip printing to screen
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
include_null_range = BioGeoBEARS_run_object$include_null_range

# Run the Maximum Likelihood optimization
res = bears_optim_run(BioGeoBEARS_run_object)
DEC_converted_lnLs = convert_BGB_to_BGB_Yule_SFs(res, tr=tr)
DEC_converted_lnLs


BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isnt much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
BioGeoBEARS_run_object$printlevel = 0					 # skip printing to screen
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
include_null_range = BioGeoBEARS_run_object$include_null_range

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.001

# Run the Maximum Likelihood optimization
res = bears_optim_run(BioGeoBEARS_run_object)
DECj_converted_lnLs = convert_BGB_to_BGB_Yule_SFs(res, tr=tr)
DECj_converted_lnLs

DECj_converted_lnLs - DEC_converted_lnLs
	' # END runtxt
	
	
	
	# Load tree from stored trfn (tree filename), if needed
	if (is.null(tr) == TRUE)
		{
		tr = try(read.tree(res$inputs$trfn))
		
		# Error trap
		if (class(tr) == "try-error")
			{
			txt = paste0("STOP ERROR in convert_BGB_lnL_to_ClaSSE(): the tree file '", es$inputs$trfn, "' could not be read into an R 'phylo' object by ape's read.tree(). Try providing the read-in tree directly through the 'tr=' argument.")
			
			cat("\n\n")
			cat(txt)
			cat("\n\n")
			stop(txt)
			} # END if (class(tr) == "try-error")
		} # END if (is.null(tr) == TRUE)
	
	
	# Error trap: tree must be ultrametric
	if (is.ultrametric(tr) == FALSE)
		{
		txt = "STOP ERROR in convert_BGB_lnL_to_ClaSSE(): is.ultrametric(tr) returned FALSE. convert_BGB_lnL_to_ClaSSE()'s calculations only make sense if you have an ultrametric tree, i.e. all of the tips survive to the present."
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		stop(txt)
		} # END if (is.ultrametric(tr) == FALSE)
	
	
	#######################################################
	# Matching diversitree and BioGeoBEARS
	#######################################################
	# BioGeoBEARS stores the likelihoods at each node, 
	# including the root node.
	#
	# But diversitree stores the likelihoods at branch bottoms,
	# and treats the root node differently, depending on 
	# various user options.

	# Get basic info from BioGeoBEARS res object:
	include_null_range = res$inputs$include_null_range
	numstates = ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)

	
	
	# A. First, we will need the tree likelihood under a pure-birth (yule) process,
	# as well as the Maximum Likelihood estimate of the birthRate
	
	# birthRate = yule(tr)$lambda  # The ML lambda from Yule. 
	                             # Equals (#speciations-1)/tree_length
	birthRate = (tr$Nnode-1)/sum(tr$edge.length)
	deathRate = 0
	
	# Get the pieces of the birthdeath log-likelihood
	bd_ape = bd_liks(tr, birthRate=birthRate, deathRate=deathRate)
	
	bd_ape$lnl_topology
	bd_ape$lnl_numBirths
	bd_ape$lnl_Births_above_root
	bd_ape$lnl_numtips_wOneMinusDeathRate
	bd_ape$lnl_branching_times
	bd_ape$lnL

	bgb1 = sum(log(res$computed_likelihoods_at_each_node[-root_nodenum]))
	bgb2 = sum(log(res$computed_likelihoods_at_each_node))
	BGB_lnL = bgb2
	bgb_root_lnL = sum(log(res$computed_likelihoods_at_each_node[root_nodenum]))
	equal_root_prob = log(1/numstates)
	equal_root_prob2 = log(1/(numstates-include_null_range)) 

	bggb_plus_Yule = bgb2 + bd_ape$lnL
	bggb_plus_Yule_minus_root = bgb1 + bd_ape$lnL
	bggb_plus_Yule_minus_root_topology = bgb1 + bd_ape$lnL - bd_ape$lnl_topology

	# Matches classe branch_lnL
	bggb_plus_Yule_1br_minus_root_topology = bgb1 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate)
	classe_branch_lnL = bggb_plus_Yule_1br_minus_root_topology

	# res1 match
	res1_bgb_classe = bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL +  log(sum(d_root_orig_BGB*d_root_orig_BGB/sum(d_root_orig_BGB)))

	# res2 match
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology - log(1/birthRate) + equal_root_prob
	res2_bgb_classe = bgb2 + bd_ape$lnL - bd_ape$lnl_topology - log(1/birthRate) + equal_root_prob
	
	
	# res3 match
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology - log(1/birthRate) + equal_root_prob2
	res3_bgb_classe = bgb2 + bd_ape$lnL - bd_ape$lnl_topology - log(1/birthRate) + equal_root_prob2

	# res5 match
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) 
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology - log(1/birthRate) 
	res5_bgb_classe = 

	# res1t match
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL +  log(sum(d_root_orig_BGB*d_root_orig_BGB/sum(d_root_orig_BGB))) + log(1/birthRate)
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1 - bgb_root_lnL +  log(sum(d_root_orig_BGB*d_root_orig_BGB/sum(d_root_orig_BGB)))
	res1t_bgb_classe = bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1 - bgb_root_lnL +  log(sum(d_root_orig_BGB*d_root_orig_BGB/sum(d_root_orig_BGB)))
	
	# res2t match
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/birthRate)
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2
	res2t_bgb_classe = bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2


	# res3t match
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/birthRate)
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2
	res3t_bgb_classe = bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2

	# res5t match
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + 1-log(1/birthRate) - bgb_root_lnL + log(sum(d_root_orig_BGB)) + equal_root_prob2 + log(1/(birthRate))
	bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 
	res5t_bgb_classe = bgb2 + bd_ape$lnL - bd_ape$lnl_topology + equal_root_prob2 

	
	# Assemble results
	resmat = rbind(BGB_lnL, bggb_plus_Yule, bggb_plus_Yule_minus_root, bggb_plus_Yule_minus_root_topology, bggb_plus_Yule_1br_minus_root_topology, classe_branch_lnL, res1_bgb_classe, res2_bgb_classe, res3_bgb_classe, res5_bgb_classe, res1t_bgb_classe, res2t_bgb_classe, res3t_bgb_classe, res5t_bgb_classe)
	converted_lnLs = adf(resmat)

	row_names = c("BGB_lnL", "bggb_plus_Yule", "bggb_plus_Yule_minus_root", "bggb_plus_Yule_minus_root_topology", "bggb_plus_Yule_1br_minus_root_topology", "classe_branch_lnL", "res1", "res2", "res3", "res5", "res1t", "res2t", "res3t", "res5t")
	colnames(converted_lnLs) = "lnL"
	row.names(converted_lnLs) = row_names
	converted_lnLs
	return(converted_lnLs)
	} # END convert_BGB_lnL_to_ClaSSE <- function(res, tr=NULL, root_probs_biased=NULL)



# Save, and return, the DEC and DEC+J lnLs for Lagrange/BioGeoBEARS,
# and ClaSSE with various options.
#
# See the script compare_BGB_diversitree_DEC_v1.R for their generation.
#
Psychotria_DEC_DECj_lnLs <- function()
	{
	runthis ='
Psychotria_BGB_lnLs = Psychotria_DEC_DECj_lnLs()

DEC_lnL = Psychotria_BGB_lnLs$DEC_lnL
DECj_lnL = Psychotria_BGB_lnLs$DECj_lnL
DEC_all_lnLs = Psychotria_BGB_lnLs$DEC_all_lnLs
DECj_all_lnLs = Psychotria_BGB_lnLs$DECj_all_lnLs

# BioGeoBEARS differences
# EXACT MATCH
DEC_lnL - DECj_lnL
-13.59554
DEC_R_result_total_LnLs5 - DECj_R_result_total_LnLs5 # EXACT MATCH
-13.59554

# BioGeoBEARS differences
DEC_DECj_lnL_difference = DEC_lnL - DECj_lnL
DEC_DECj_lnL_difference

# Diversitree differences
res5_diffs = DEC_R_result_total_LnLs5 - DECj_R_result_total_LnLs5
res5_diffs
abs(res5_diffs - DEC_DECj_lnL_difference) < 0.01

branch_lnLs_diffs = DEC_all_lnLs$branch_LnL - DECj_all_lnLs$branch_LnL
branch_lnLs_diffs
abs(branch_lnLs_diffs - DEC_DECj_lnL_difference) < 0.01

all_lnLs_diffs = DEC_all_lnLs$ttl_LnL - DECj_all_lnLs$ttl_LnL
all_lnLs_diffs
abs(all_lnLs_diffs - DEC_DECj_lnL_difference) < 0.01


' # END runthis
	
	DEC_lnL = -34.54313
	DECj_lnL = -20.94759

	DEC_all_lnLs = structure(list(ttl_LnL = c(-72.602116, -74.33632, -74.271782, 
	-74.605042, -71.563731, -74.263777, -71.48986, -73.159526, -73.159526, 
	-73.492786, -73.159526, -73.079568), branch_LnL = c(-67.629498, 
	-67.629498, -67.629498, -67.629498, -67.629498, -67.629498, -67.629498, 
	-67.629498, -67.629498, -67.629498, -67.629498, -67.629498), 
			ObsDiff = c(-4.972618, -6.706822, -6.642284, -6.975544, -3.934234, 
			-6.634279, -3.860362, -5.530028, -5.530028, -5.863288, -5.530028, 
			-5.45007), LnLdiff = c(-3.8604, -5.5946, -5.53, -5.8633, 
			-2.822, -5.522, -2.7481, -4.4178, -4.4178, -4.751, -4.4178, 
			-4.3378), exp_ObsDiff = c(0.00692499, 0.00122254, 0.00130405, 
			0.00093446, 0.0195607, 0.00131453, 0.0210604, 0.00396588, 
			0.00396588, 0.00284188, 0.00396588, 0.004296), exp_LnLdiff = c(0.0210604, 
			0.00371801, 0.00396588, 0.00284188, 0.0594882, 0.00399775, 
			0.064049, 0.0120611, 0.0120611, 0.00864277, 0.0120611, 0.0130651
			)), .Names = c("ttl_LnL", "branch_LnL", "ObsDiff", "LnLdiff", 
	"exp_ObsDiff", "exp_LnLdiff"), row.names = c("LnLs1", "LnLs2", 
	"LnLs3", "LnLs4", "LnLs5", "LnLs6", "LnLs1t", "LnLs2t", "LnLs3t", 
	"LnLs4t", "LnLs5t", "LnLs6t"), class = "data.frame")



	DECj_all_lnLs = 
	structure(list(ttl_LnL = c(-58.837585, -60.740782, -60.676244, 
	-61.330746, -57.968194, -59.96035, -57.725329, -59.563988, -59.563988, 
	-60.218491, -59.563988, -58.848094), branch_LnL = c(-55.37332, 
	-55.37332, -55.37332, -55.37332, -55.37332, -55.37332, -55.37332, 
	-55.37332, -55.37332, -55.37332, -55.37332, -55.37332), ObsDiff = c(-3.464265, 
	-5.367462, -5.302924, -5.957426, -2.594874, -4.58703, -2.352009, 
	-4.190668, -4.190668, -4.845171, -4.190668, -3.474774), LnLdiff = c(-2.352, 
	-4.2552, -4.1907, -4.8452, -1.4826, -3.4748, -1.2398, -3.0784, 
	-3.0784, -3.7329, -3.0784, -2.3625), exp_ObsDiff = c(0.031296, 
	0.00466596, 0.00497702, 0.00258656, 0.0746553, 0.0101831, 0.0951777, 
	0.0151362, 0.0151362, 0.00786628, 0.0151362, 0.0309688), exp_LnLdiff = c(0.0951777, 
	0.0141902, 0.0151362, 0.00786628, 0.227043, 0.0309688, 0.289456, 
	0.0460323, 0.0460323, 0.023923, 0.0460323, 0.0941827)), .Names = c("ttl_LnL", 
	"branch_LnL", "ObsDiff", "LnLdiff", "exp_ObsDiff", "exp_LnLdiff"
	), row.names = c("LnLs1", "LnLs2", "LnLs3", "LnLs4", "LnLs5", 
	"LnLs6", "LnLs1t", "LnLs2t", "LnLs3t", "LnLs4t", "LnLs5t", "LnLs6t"
	), class = "data.frame")


	# BioGeoBEARS differences
	DEC_DECj_lnL_difference = DEC_lnL - DECj_lnL
	DEC_DECj_lnL_difference

	# Diversitree differences
	res5_diffs = DEC_R_result_total_LnLs5 - DECj_R_result_total_LnLs5
	res5_diffs
	
	branch_lnLs_diffs = DEC_all_lnLs$branch_LnL - DECj_all_lnLs$branch_LnL
	branch_lnLs_diffs
	
	all_lnLs_diffs = DEC_all_lnLs$ttl_LnL - DECj_all_lnLs$ttl_LnL
	all_lnLs_diffs
	
	cat("\n\nBioGeoBEARS DEC lnL:\n")
	cat(DEC_lnL)

	cat("\n\nBioGeoBEARS DEC+J lnL:\n")
	cat(DECj_lnL)
	
	cat("\n\nDiversitree claSSE lnLs: res5:\n")
	cat("DEC_R_result_total_LnLs5:\n")
	print(DEC_R_result_total_LnLs5)
	cat("\n\nDECj_R_result_total_LnLs5:\n")
	print(DECj_R_result_total_LnLs5)
	cat("\n\nDEC vs. DEC+J diffs in claSSE res5_diffs:\n")
	print(res5_diffs)
	print(abs(res5_diffs-DEC_DECj_lnL_difference) < 0.01)
	
	cat("\n\nDiversitree claSSE lnLs: branch_lnLs:\n")
	cat("DEC_all_lnLs$branch_LnL:\n")
	print(DEC_all_lnLs$branch_LnL)
	cat("\n\nDECj_all_lnLs$branch_LnL:\n")
	print(DECj_all_lnLs$branch_LnL)
	cat("\n\nDEC vs. DEC+J diffs in claSSE branch_lnLs:\n")
	print(branch_lnLs_diffs)
	print(abs(branch_lnLs_diffs-DEC_DECj_lnL_difference) < 0.01)

	cat("\n\nDiversitree claSSE lnLs: all_lnLs:\n")
	cat("DEC_all_lnLs$ttl_LnL:\n")
	print(DEC_all_lnLs$ttl_LnL)
	cat("\n\nDECj_all_lnLs$ttl_LnL:\n")
	print(DECj_all_lnLs$branch_LnL)
	cat("\n\nDEC vs. DEC+J diffs in claSSE branch_lnLs:\n")
	print(all_lnLs_diffs)
	print(abs(all_lnLs_diffs-DEC_DECj_lnL_difference) < 0.01)
	
	Psychotria_BGB_lnLs = list()
	Psychotria_BGB_lnLs$DEC_lnL = DEC_lnL
	Psychotria_BGB_lnLs$DECj_lnL = DECj_lnL
	Psychotria_BGB_lnLs$DEC_all_lnLs = DEC_all_lnLs
	Psychotria_BGB_lnLs$DECj_all_lnLs = DECj_all_lnLs
	
	runtxt='
	DEC_lnL = Psychotria_BGB_lnLs$DEC_lnL
	DECj_lnL = Psychotria_BGB_lnLs$DECj_lnL
	DEC_all_lnLs = Psychotria_BGB_lnLs$DEC_all_lnLs
	DECj_all_lnLs = Psychotria_BGB_lnLs$DECj_all_lnLs
	'
	
	return(Psychotria_BGB_lnLs)
	}


