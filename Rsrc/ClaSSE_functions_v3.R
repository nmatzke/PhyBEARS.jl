#######################################################
# Located in:
# /drives/GDrive/z_help/ClaSSE_LnLs/grok_classe_setup/
#######################################################


# Set your working directory to wherever you unzipped this directory
# wd = "/drives/GDrive/z_help/ClaSSE_LnLs/"
# setwd(wd)
# 
# library(diversitree)
# 
# Made with make.classe
# classe_2areas <- function(pars, condition.surv = TRUE, root = ROOT.OBS, root.p = NULL, intermediates = FALSE)  
# 	{
# 	Returns pars2 as the lambdas & mus, and the Qmat
#     pars2 <- f.pars(pars)
#     
#     This calculates the likelihoods
#     res <- all.branches(pars2, intermediates)
#     
#     This calculates the root state likelihoods, and combines
#     res etc. into the outputs, if intermediates == TRUE
#     rootfunc(res, pars, condition.surv, root, root.p, intermediates)
# 	}
# 
# BGB_DEC_LnL + BD_LnL
# -8.17333
# 
# BGB_DEC_LnL_e01 + BD_LnL
# -8.304785
# 
# BGB_DECj_LnL + BD_LnL
# -7.842103
# 
# 
# ClaSSE-DEC
# -13.5624
# 
# ClaSSE-DEC_e01
# -13.63116
# 
# ClaSSE-DECj
# -13.05974
# 
# ClaSSE_LnLs = c(-13.5624, -13.63116, -13.05974)
# ClaSSE_LnLs - min(ClaSSE_LnLs)
# 
# BGB_LnLs = c(-8.17333, -8.304785, -7.842103)
# BGB_LnLs - min(BGB_LnLs)
# 
# 
# 
# ClaSSE_calcs_DEC = classe_2areas(pars=classe_params_DEC, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# lq <- attr(ClaSSE_calcs_DEC, "intermediates")$lq
# exp(lq)
# sum(lq)
# t(attr(ClaSSE_calcs_DEC, "intermediates")$init)[,5:8]
# 
# ClaSSE_calcs_DEC_e01 = classe_2areas(pars=classe_params_DEC_e01, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# lq <- attr(ClaSSE_calcs_DEC_e01, "intermediates")$lq
# exp(lq)
# sum(lq)
# t(attr(ClaSSE_calcs_DEC_e01, "intermediates")$init)[,5:8]
# 
# ClaSSE_calcs_DECj = classe_2areas(pars=classe_params_DECj, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=TRUE)
# lq <- attr(ClaSSE_calcs_DECj, "intermediates")$lq
# exp(lq)
# sum(lq)
# t(attr(ClaSSE_calcs_DECj, "intermediates")$init)[,5:8]
# 
# 
# cbind(classe_params_DEC, classe_params_DEC_e01, classe_params_DECj)
# 
# 
# t(attr(ClaSSE_calcs_DEC, "intermediates")$init)[,5:8]
# 
# 
# 
# 
# 
# 
# tree=tr
# states
# k
# sampling.f=NULL
# strict=FALSE
# control=list()
# 
# cache <- diversitree:::make.cache.classe(tree, states, k, sampling.f, strict)
# initial.conditions <- diversitree:::make.initial.conditions.classe(k)
# all.branches <- diversitree:::make.all.branches.dtlik(cache, control, initial.conditions)
# rootfunc <- diversitree:::rootfunc.classe
# 
# 
# f.pars <- diversitree:::make.pars.classe(k)
# 

# diversitree/R/util.R
# https://github.com/richfitz/diversitree/blob/b24aa1417f99525bd08aaa55691306b831db157c/R/util.R
quadratic.roots <- function(a, b, c)
	{
	(-b + c(-1, 1) * sqrt(b*b - 4*a*c))/(2 * a)
	}

# Made with make.bisse
# bisse_2areas = make.bisse(tree=tr, states=states, sampling.f=sampling.f, strict=FALSE)
# printed to screen with dput(bisse_2areas)


bisse_2areas_default <- function(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL, intermediates=FALSE) 
	{
	diversitree:::check.pars.bisse(pars)
	preset <- diversitree:::branches.unresolved.bisse(pars, unresolved)
	ans <- all.branches(pars, intermediates, preset)
	rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
	}




# Get total LnL and branch LnL from ClaSSE output
get_classe_LnLs <- function(classe_res)
	{
	
	branch_LnL = sum(attr(classe_res, "intermediates")$lq)
	ttl_LnL = classe_res
	attributes(ttl_LnL) = NULL
	
	return(c(ttl_LnL, branch_LnL))
	}



# Made with make.classe
classe_2areas <- function(pars, condition.surv = TRUE, root = ROOT.OBS, root.p = NULL, intermediates = FALSE)  
	{
	# Returns pars2 as the lambdas & mus, and the Qmat
    pars2 <- f.pars(pars)
    
    # This calculates the likelihoods
    res <- all.branches(pars2, intermediates)
    
    # This calculates the root state likelihoods, and combines
    # res etc. into the outputs, if intermediates == TRUE
    rootfunc(res, pars, condition.surv, root, root.p, intermediates)
	}



# f.pars <- diversitree:::make.pars.classe(k)
f.pars <- function(pars, k)
	{
    check.pars.classe(pars, k)
    qmat[idx.qmat] <- pars[idx.q]
    diag(qmat) <- -rowSums(qmat)
    c(pars[idx.lm], qmat)
	}

# diversitree:::check.pars.classe
check.pars.classe <- function(pars, k)
	{
	# This is an algebraic simplification of
	# k * ((k * (k+1)) / 2) (k times (1 + 2 + 3 + ... + k))  (number of speciation params for each starting k, times k starts)
	# (assumes symmetry)
	# k (number of extinction params)
	# k * (k-1) (number of Q parameters, assuming symmetry)
	check.pars.nonnegative(pars, (k + 3) * k * k/2)
	}

# diversitree:::check.pars.nonnegative 
check.pars.nonnegative <- function(pars, npar) 
	{
    if (length(pars) != npar) 
        stop(sprintf("Incorrect parameter length: expected %d, got %d", 
            npar, length(pars)))
    if (any(!is.finite(pars)) || any(pars < 0)) 
        stop("Parameters must be non-negative and finite")
    pars
	}


#all.branches <- diversitree:::make.all.branches.dtlik(cache, control, initial.conditions)
all.branches <- function(pars, intermediates, preset = NULL)
	{
	branches <- diversitree:::make.branches.dtlik(cache$info, control)
	all.branches.matrix(pars, cache, initial.conditions, branches, preset)
	}



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
                branch.init[, idx] <- x$y									# $init
                ans <- branches(x$y, x$t, pars, 0, idx)		# list of ans
                lq[idx] <- ans[[1]]												# $lq
                branch.base[, idx] <- ans[[2]]						# $base
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
    ans_lists = NULL
    branches_function_list = NULL
    countval = 0
    for (i in order) {
    	# Count so you can treat the root specially
    	countval = countval + 1
        y.in <- initial.conditions(branch.base[, children[i, ]], pars, depth[i], i)
    	
    	# Normalize to sum to 1, if you are at the root state
    	# THIS CONVERTS LIKELIHOODS TO ROOT.OBS
    	if (countval != length(cache$order))
    		{
    		y.in = y.in / sum(y.in)
	        }
   		#y.in = y.in / sum(y.in)
        if (!all(is.finite(y.in))) 
            stop("Bad initial conditions: calculation failure along branches?")
        branch.init[, i] <- y.in
        ans <- branches(y.in, len[i], pars, depth[i], i)
        lq[i] <- ans[[1]]
        branch.base[, i] <- ans[[2]]
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


make.bisse <- function(tree, states, unresolved=NULL, sampling.f=NULL, nt.extra=10, strict=TRUE, control=list()) 
	{
	cache <- diversitree:::make.cache.bisse(tree, states, unresolved, sampling.f, nt.extra, strict)
	unresolved <- cache$unresolved
	all.branches <- diversitree:::make.all.branches.dtlik(cache, control, diversitree:::initial.conditions.bisse)
	rootfunc <- rootfunc.musse
	ll <- function(pars, condition.surv = TRUE, root = ROOT.OBS, root.p = NULL, intermediates = FALSE)
		{
		diversitree:::check.pars.bisse(pars)
		preset <- diversitree:::branches.unresolved.bisse(pars, unresolved)
		ans <- all.branches(pars, intermediates, preset)
		rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
		}
	class(ll) <- c("bisse", "dtlik", "function")
	ll
	}


# diversitree:::rootfunc.musse
rootfunc.musse <- function (res, pars, condition.surv, root, root.p, intermediates) 
	{
	vals <- res$vals
	lq <- res$lq
	k <- length(vals)/2
	i <- seq_len(k)
	d.root <- vals[-i]
	if (k == 2)
		{
		root.equi <- stationary.freq.bisse
		} else {
		root.equi <- NULL
		}
	root.p <- root.p.calc(d.root, pars, root, root.p, root.equi)
	if (condition.surv)
		{
		lambda <- pars[i]
		e.root <- vals[i]
		d.root <- d.root/sum(root.p * lambda * (1 - e.root)^2)
		}
	if (root == ROOT.ALL)
		{
		loglik <- log(d.root) + sum(lq)
		} else {
		loglik <- log(sum(root.p * d.root)) + sum(lq)
		}
	if (intermediates == TRUE)
		{
		res$root.p <- root.p
		attr(loglik, "intermediates") <- res
		attr(loglik, "vals") <- vals
		}
	loglik
	}

#diversitree:::stationary.freq.bisse
stationary.freq.bisse <- function (pars) 
	{
	lambda0 <- pars[1]
	lambda1 <- pars[2]
	mu0 <- pars[3]
	mu1 <- pars[4]
	q01 <- pars[5]
	q10 <- pars[6]
	g <- (lambda0 - mu0) - (lambda1 - mu1)
	eps <- (lambda0 + mu0 + lambda1 + mu1) * 1e-14
	if (abs(g) < eps)
		{
		if (q01 + q10 == 0)
			{
			p <- 0.5
			} else {
			p <- q10/(q01 + q10)
			}
		} else {
		roots <- quadratic.roots(g, q10 + q01 - g, -q10)
		roots <- roots[roots >= 0 & roots <= 1]
		if (length(roots) > 1)
			{
			p <- NA
			}
		else p <- roots
		}
	c(p, 1 - p)
	}



# diversitree:::root.p.calc
root.p.calc <- function(vals, pars, root, root.p=NULL, root.equi=NULL) 
	{
	if (!is.null(root.p) && root != ROOT.GIVEN) 
			warning("Ignoring specified root state")
	k <- length(vals)
	if (root == ROOT.FLAT) {
			p <- 1/k
	}
	else if (root == ROOT.EQUI) {
			if (is.null(root.equi)) 
					stop("Equilibrium root probability not possible with this method")
			p <- root.equi(pars)
	}
	else if (root == ROOT.OBS) {
			p <- vals/sum(vals)
	}
	else if (root == ROOT.GIVEN) {
			if (length(root.p) != length(vals)) 
					stop("Invalid length for root.p")
			p <- root.p
	}
	else if (root == ROOT.ALL) {
			p <- rep(1, k)
	}
	else {
			stop("Invalid root mode")
	}
	p
	}



#diversitree:::make.cache.musse
make.cache.musse <- function(tree, states, k, sampling.f = NULL, strict = TRUE) 
	{
	tree <- check.tree(tree)
	states <- check.states(tree, states, strict = strict, strict.vals = 1:k)
	cache <- make.cache(tree)
	cache$info <- make.info.musse(k, tree)
	cache$states <- states
	cache$sampling.f <- check.sampling.f(sampling.f, k)
	cache$y <- initial.tip.xxsse(cache)
	cache
	}



#diversitree:::check.tree
check.tree <- function(tree, ultrametric = TRUE, bifurcating = TRUE, node.labels = FALSE) 
	{
	if (!inherits(tree, "phylo")) 
			stop("'tree' must be a valid phylo tree")
	if (ultrametric && !is.ultrametric(tree)) 
			stop("'tree' must be ultrametric")
	if (any(tree$edge.length < 0)) 
			stop("Negative branch lengths in tree")
	if (bifurcating && (!is.binary(tree) || any(tabulate(tree$edge[, 
			1]) == 1))) 
			stop("'tree must be bifurcating (no polytomies or unbranched nodes)'")
	if (any(duplicated(tree$tip.label))) 
			stop("Tree contains duplicated tip labels")
	if (node.labels)
		{
		if (is.null(tree$node.label))
			{
			tree$node.label <- sprintf("nd%d", seq_len(tree$Nnode))
			}
		else if (any(duplicated(tree$node.label)))
			{
			stop("Tree contains duplicated node labels")
			}
		}
	return(tree)
	} # END check.tree <- function(tree, ultrametric = TRUE, bifurcating = TRUE, node.labels = FALSE) 


#diversitree:::check.states
check.states <- function(tree, states, allow.unnamed=FALSE, strict=FALSE, strict.vals=NULL, as.integer=TRUE) 
	{
	if (is.matrix(states))
		{
			if (inherits(tree, "clade.tree")) 
					stop("Clade trees won't work with multistate tips yet")
			n <- rowSums(states > 0)
			if (any(n == 0)) 
					stop(sprintf("No state found for taxa: %s", paste(names(n)[n == 
							0], collapse = ", ")))
			i.mono <- which(n == 1)
			i.mult <- which(n > 1)
			tmp <- diversitree:::matrix.to.list(states)
			names(tmp) <- rownames(states)
			states.mult <- lapply(tmp[i.mult], as.numeric)
			states <- rep(NA, length(tmp))
			names(states) <- names(tmp)
			states[i.mono] <- sapply(tmp[i.mono], function(x) which(x != 
					0))
			attr(states, "multistate") <- list(i = i.mult, states = states.mult)
		} # END if (is.matrix(states))

	if (is.null(names(states)))
		{
		if (allow.unnamed)
			{
			if (length(states) == length(tree$tip.label))
				{
				names(states) <- tree$tip.label
				warning("Assuming states are in tree$tip.label order")
				} else {
				stop(sprintf("Invalid states length (expected %d)", length(tree$tip.label)))
				}
		} else {
		stop("The states vector must contain names")
		}
	} # END if (is.null(names(states)))

	if (!all(tree$tip.label %in% names(states))) 
		stop("Not all species have state information")
	
	
	if (!is.null(strict.vals))
		{
		if (isTRUE(all.equal(strict.vals, 0:1))) 
				if (is.logical(states)) 
						states[] <- as.integer(states)
		if (strict)
			{
			if (!isTRUE(all.equal(sort(strict.vals), sort(unique(na.omit(states)))))) 
					stop("Because strict state checking requested, all (and only) ", 
						sprintf("states in %s are allowed", paste(strict.vals, 
							collapse = ", ")))
			} else {
			extra <- setdiff(sort(unique(na.omit(states))), strict.vals)
			if (length(extra) > 0) 
					stop(sprintf("Unknown states %s not allowed in states vector", paste(extra, collapse = ", ")))
			}
		if (as.integer && any(!is.na(states))) 
				states <- check.integer(states)
		} # END if (!is.null(strict.vals))

	if (inherits(tree, "clade.tree"))
		{
		spp.clades <- unlist(tree$clades)
		if (!all(spp.clades %in% names(states))) 
				stop("Species in 'clades' do not have states information")
		states[union(tree$tip.label, spp.clades)]
		} else {
		ret <- states[tree$tip.label]
		attr(ret, "multistate") <- attr(states, "multistate")
		ret
		} # END if (inherits(tree, "clade.tree"))
	} # END check.states <- function(tree, states, allow.unnamed=FALSE, strict=FALSE, strict.vals=NULL, as.integer=TRUE) 



#diversitree:::check.integer
check.integer <- function(x) 
	{
	if (is.null(x)) 
		stop("NULL argument for ", deparse(substitute(x)))
	nna <- !is.na(x)
	if (length(x) > 0 && !any(nna)) 
		stop("No non-NA values for ", deparse(substitute(x)))
	if (length(x) && max(abs(x[nna] - round(x[nna]))) > 1e-08) 
		stop("Non-integer argument for ", deparse(substitute(x)))
	storage.mode(x) <- "integer"
	return(x)
	} # END check.integer <- function(x) 





#diversitree:::make.cache
make.cache <- function (tree) 
	{
	if (inherits(tree, "phylo")) 
			class(tree) <- "phylo"
	edge <- tree$edge
	edge.length <- tree$edge.length
	idx <- seq_len(max(edge))
	n.tip <- length(tree$tip.label)
	tips <- seq_len(n.tip)
	root <- n.tip + 1
	is.tip <- idx <= n.tip
	children <- get.children(edge, n.tip)
	parent <- edge[match(idx, edge[, 2]), 1]
	order <- get.ordering(children, is.tip, root)
	len <- edge.length[match(idx, edge[, 2])]
	height <- branching.heights(tree)
	depth <- max(height) - height
	depth2 <- branching.depth(len, children, order, tips)
	i <- abs(depth - depth2) < 1e-08
	depth[i] <- depth2[i]
	if (is.ultrametric(tree)) 
			depth[tips] <- 0
	anc <- vector("list", max(order))
	for (i in c(rev(order[-length(order)]), tips)) anc[[i]] <- c(parent[i], 
			anc[[parent[i]]])
	ans <- list(tip.label = tree$tip.label, node.label = tree$node.label, 
			len = len, children = children, parent = parent, order = order, 
			root = root, n.tip = n.tip, n.node = tree$Nnode, tips = tips, 
			height = height, depth = depth, ancestors = anc, edge = edge, 
			edge.length = edge.length)
	ans
	} # END make.cache <- function (tree) 


#diversitree:::get.children
get.children <- function(edge, n.tip) 
	{
	x <- as.integer(edge[, 1])
	levels <- as.integer((n.tip + 1):max(edge[, 1]))
	f <- match(x, levels)
	levels(f) <- as.character(levels)
	class(f) <- "factor"
	children <- split(edge[, 2], f)
	names(children) <- NULL
	if (!all(unlist(lapply(children, length)) == 2)) 
			stop("Multifircations/unbranched nodes in tree - best get rid of them")
	rbind(matrix(NA, n.tip, 2), t(matrix(unlist(children), 2)))
	} # END get.children <- function(edge, n.tip) 

# diversitree:::get.ordering
get.ordering <- function(children, is.tip, root) 
	{
	todo <- list(root)
	i <- root
	repeat
		{
		kids <- children[i, ]
		i <- kids[!is.tip[kids]]
		if (length(i) > 0) 
			todo <- c(todo, list(i))
		else break
		}
	as.vector(unlist(rev(todo)))
	} # END get.ordering <- function(children, is.tip, root) 

# diversitree:::branching.heights
branching.heights <- function(phy) 
	{
	if (!inherits(phy, "phylo")) 
			stop("object \"phy\" is not of class \"phylo\"")
	phy <- reorder(phy, "cladewise")
	edge <- phy$edge
	n.node <- phy$Nnode
	n.tip <- length(phy$tip.label)
	ht <- numeric(n.node + n.tip)
	for (i in seq_len(nrow(edge))) ht[edge[i, 2]] <- ht[edge[i,1]] + phy$edge.length[i]
	names.node <- phy$node.label
	if (is.null(names.node)) 
			names.node <- (n.tip + 1):(n.tip + n.node)
	names(ht) <- c(phy$tip.label, names.node)
	ht
	} # END branching.heights <- function(phy) 



# diversitree:::branching.depth
branching.depth <- function(len, children, order, tips) 
	{
	depth <- numeric(nrow(children))
	depth[tips] <- 0
	for (i in order) depth[i] <- depth[children[i, 1]] + len[children[i,1]]
	depth
	} # END branching.depth <- function(len, children, order, tips) 




# diversitree:::make.info.musse
make.info.musse <- function(k, phy) 
	{
	list(name = "musse", name.pretty = "MuSSE", np = as.integer(k * 
			(k + 2)), argnames = default.argnames.musse(k), ny = as.integer(2 * 
			k), k = as.integer(k), idx.e = as.integer(1:k), idx.d = as.integer((k + 
			1):(2 * k)), derivs = derivs.musse, phy = phy, ml.default = "subplex", 
			mcmc.lowerzero = TRUE, doc = NULL, reference = c("FitzJohn (submitted)"))
	} # END make.info.musse <- function(k, phy) 


# diversitree:::make.info.classe
make.info.classe <- function(k, phy) 
	{
	list(name = "classe", name.pretty = "ClaSSE", np = as.integer((k + 
		3) * k * k/2 + k), argnames = default.argnames.classe(k), 
		ny = as.integer(2 * k), k = as.integer(k), idx.e = as.integer(1:k), 
		idx.d = as.integer((k + 1):(2 * k)), derivs = derivs.classe, 
		phy = phy, ml.default = "subplex", mcmc.lowerzero = TRUE, 
		doc = NULL, reference = c("Goldberg (submitted)"))
	} # END make.info.classe <- function(k, phy) 



# diversitree:::default.argnames.musse
default.argnames.musse <- function(k) 
	{
	fmt <- sprintf("%%0%dd", ceiling(log10(k + 0.5)))
	str <- sprintf(fmt, 1:k)
	c(sprintf("lambda%s", str), sprintf("mu%s", str), sprintf("q%s%s", 
			rep(str, each = k - 1), unlist(lapply(1:k, function(i) str[-i]))))
	} # END default.argnames.musse <- function(k) 



# diversitree:::default.argnames.classe
default.argnames.classe <- function(k) 
	{
	fmt <- sprintf("%%0%dd", ceiling(log10(k + 0.5)))
	sstr <- sprintf(fmt, 1:k)
	lambda.names <- sprintf("lambda%s%s%s", rep(sstr, each = k * 
			(k + 1)/2), rep(rep(sstr, times = seq(k, 1, -1)), k), 
			unlist(lapply(1:k, function(i) sstr[i:k])))
	mu.names <- sprintf("mu%s", sstr)
	q.names <- sprintf("q%s%s", rep(sstr, each = k - 1), unlist(lapply(1:k, 
			function(i) sstr[-i])))
	c(lambda.names, mu.names, q.names)
	} # END default.argnames.classe <- function(k) 


# diversitree:::derivs.musse
derivs.musse <- function(t, y, pars) 
	{
	k <- length(y)/2L
	i1 <- seq_len(k)
	i2 <- i1 + k
	i3 <- (2 * k + 1L):(k * (k + 2L))
	E <- y[i1]
	D <- y[i2]
	lambda <- pars[i1]
	mu <- pars[i2]
	Q <- matrix(pars[i3], k, k)
	c(mu - (lambda + mu) * E + lambda * E * E + Q %*% E, -(lambda + 
			mu) * D + 2 * lambda * E * D + Q %*% D)
	} # END derivs.musse <- function(t, y, pars) 



# diversitree:::derivs.classe
derivs.classe <- function(t, y, pars) 
	{
	stop("Not yet possible")
	}


# diversitree:::check.sampling.f
check.sampling.f <- function (sampling.f, n) 
	{
	if (is.null(sampling.f)) 
			sampling.f <- rep(1, n)
	else sampling.f <- check.par.length(sampling.f, n)
	if (max(sampling.f) > 1 || min(sampling.f) <= 0) 
			stop("sampling.f must be on range (0,1]")
	sampling.f
	} # END check.sampling.f <- function (sampling.f, n) 
	

# diversitree:::check.par.length
check.par.length <- function(x, length) 
	{
	if (length(x) == 1) 
			rep(x, length)
	else if (length(x) == length) 
			x
	else stop(sprintf("'%s' of incorrect length", deparse(substitute(x))))
	} # END check.par.length <- function(x, length) 

# diversitree:::initial.tip.xxsse
initial.tip.xxsse <- function(cache, base.zero = FALSE) 
	{
	k <- cache$info$k
	f <- cache$sampling.f
	y <- matrix(rep(c(1 - f, rep(0, k)), k + 1), k + 1, 2 * k, 
			TRUE)
	y[k + 1, (k + 1):(2 * k)] <- diag(y[1:k, (k + 1):(2 * k)]) <- f
	y <- diversitree:::matrix.to.list(y)
	y.i <- cache$states
	if (base.zero) 
			y.i <- y.i + 1L
	y.i[is.na(y.i)] <- k + 1
	if (!is.null(multistate <- attr(cache$states, "multistate")))
		{
		y.multi <- unique(multistate$states)
		y.i.multi <- match(multistate$states, y.multi)
		y <- c(y, lapply(y.multi, function(x) c(1 - f, x)))
		y.i[multistate$i] <- y.i.multi + k + 1
		}
	dt.tips.grouped(y, y.i, cache)
	} # END initial.tip.xxsse <- function(cache, base.zero = FALSE) 


# diversitree:::matrix.to.list
# This one has to refer to diversitree:::matrix.to.list
#matrix.to.list <- function(m) 
#	{
#	storage.mode(m) <- "double"
#	.Call(r_matrix_to_list, m)
#	}




# diversitree:::dt.tips.grouped
dt.tips.grouped <- function(y, y.i, cache) 
	{
	tips <- cache$tips
	t <- cache$len[tips]
	if (!is.list(y)) 
		stop("'y' must be a list of initial conditions")
	if (max(y.i) > length(y) || min(y.i) < 1) 
		stop("'y.i' must be integers on 1..", length(y))
	if (length(y.i) != length(tips)) 
		stop("y must be same length as tips")
	if (length(y.i) != length(t)) 
		stop("y must be the same length as t")
	if (any(is.na(y.i)))
		{
		k <- cache$info$k
		if (!is.null(k) && !is.na(k) && length(y) == k + 1) 
				y.i[is.na(y.i)] <- k + 1
		else stop("Unhandled NA values in state vector")
		}
	if (max(abs(cache$depth[tips])) > .Machine$double.eps^0.5) 
		stop("This currently only works for ultrametric trees")
	types <- sort(unique(y.i))
	res <- vector("list", length(types))
	for (i in seq_along(types))
		{
		type <- types[i]
		j <- which(y.i == type)
		ord <- order(t[j])
		res[[i]] <- list(y = y[[type]], y.i=i, target=tips[j][ord], t=t[j][ord], type="GROUPED")
		}
	res
	} # END dt.tips.grouped <- function(y, y.i, cache) 





make.classe <- function(tree, states, k, sampling.f=NULL, strict=TRUE, control=list())
	{
	## Note that this uses MuSSE's cache...
	cache <- diversitree:::make.cache.classe(tree, states, k, sampling.f, strict)
	initial.conditions <- diversitree:::make.initial.conditions.classe(k)
	all.branches <- diversitree:::make.all.branches.dtlik(cache, control, initial.conditions)
	rootfunc <- diversitree:::rootfunc.classe
	f.pars <- diversitree:::make.pars.classe(k)

	ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL, intermediates=FALSE)
		{
		pars2 <- f.pars(pars)
		ans <- all.branches(pars2, intermediates)
		ans$branchLnL = sum(ans$lq)
		## TODO: This is different to other functions, as the
		## stationary.freq function assumes the non-expanded case.
		## However, it would be straightforward to modify stationary.freq
		## classe to use the expanded case.  At worst, it could just strip
		## off the extra parameters, but I think that it builds these
		## anyway.
		##
		## This will be an issue for creating a split or time version.
		rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
		}
	
	# 
	
	class(ll) <- c("classe", "dtlik", "function")
	ll
	}




# diversitree:::rootfunc.classe
rootfunc.classe <- function(res, pars, condition.surv, root, root.p, intermediates) 
	{
	vals <- res$vals
	lq <- res$lq
	k <- length(vals)/2
	i <- seq_len(k)
	d.root <- vals[-i]
	root.equi <- function(pars)
		{
		stationary.freq.classe(pars, k)
		}
	root.p <- root.p.calc(d.root, pars, root, root.p, root.equi)
	if (condition.surv)
		{
			nsum <- k * (k + 1)/2
			lambda <- colSums(matrix(pars[1:(nsum * k)], nrow = nsum))
			e.root <- vals[i]
			d.root <- d.root/sum(root.p * lambda * (1 - e.root)^2)
		}
	if (root == ROOT.ALL) 
			{
			loglik <- log(d.root) + sum(lq)
			} else 
			{
			loglik <- log(sum(root.p * d.root)) + sum(lq)
			}
	if (intermediates)
		{
			res$root.p <- root.p
			attr(loglik, "intermediates") <- res
			attr(loglik, "vals") <- vals
		}
	loglik
	}






# diversitree:::stationary.freq.classe
stationary.freq.classe <- function(pars, k) 
	{
    if (k == 2)
    	{
        g <- (sum(pars[1:3]) - pars[7]) - (sum(pars[4:6]) - pars[8])
        eps <- sum(pars[1:8]) * 1e-14
        ss1 <- pars[9] + 2 * pars[3] + pars[2]
        ss2 <- pars[10] + 2 * pars[4] + pars[5]
        if (abs(g) < eps)
        	{
            if (ss1 + ss2 == 0) 
            	{
                eqfreq <- 0.5
	            } else {
	            eqfreq <- ss2/(ss1 + ss2)
	            }
            eqfreq <- c(eqfreq, 1 - eqfreq)
	        } else {
            roots <- quadratic.roots(g, ss2 + ss1 - g, -ss2)
            eqfreq <- roots[roots >= 0 & roots <= 1]
            if (length(eqfreq) > 1) 
            	{
                eqfreq <- NA
                } else {
                eqfreq <- c(eqfreq, 1 - eqfreq)
                }
            } # END if (abs(g) < eps)
    	} else {
        eqfreq <- stationary.freq.classe.ev(pars, k)
    	} # END if (k == 2)
    eqfreq
	}


# diversitree:::stationary.freq.classe.ev
stationary.freq.classe.ev <- function(pars, k) 
	{
	A <- projection.matrix.classe(pars, k)
	evA <- eigen(A)
 		# https://www-sciencedirect-com.ezproxy.auckland.ac.nz/topics/mathematics/dominant-eigenvalue
		# Predicting Population Growth: Modeling with Projection Matrices
		# Janet Steven, James Kirkwood, in Mathematical Concepts and Methods in Modern Biology, 2013"
		# 7.8.4 Finding the Stable Distribution
		#
		# Suppose that A is a projection matrix that meets the assumptions of the Perron-Frobenius 
		# theorem and that va is any vector.
		# 
		# ...so the equilibrium state is the normalized eigenvector for the dominant eigenvalue.

	# OLD: i <- which(evA$values == max(evA$values))
	# NEW: 
	i <- which(abs(evA$values) == max(abs(evA$values)))
	evA$vectors[, i]/sum(evA$vectors[, i])
	
	# Translation: the stationary frequences result
	# from the instantaneous rate matrix, A
	# taking its eigenvalues and eigenvectors
	# Finding the largest eigenvector, here i=3
	# Taking the 3rd column of eigenvectors, and normalizing
	# to get the base frequencies
	}


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
	
	# For 4 states, nsum = 10 descendant combinations possible
	# sum(pars.lam[seq((i-1) * nsum + 1, i * nsum)]) means:
	# sum lambdas 1-10, 11-20, etc.
	diag(A) <- -colSums(A) + unlist(lapply(kseq, function(i) sum(pars.lam[seq((i-1) * nsum + 1, i * nsum)]) - pars.mu[i]))
	A
	}









#all.branches <- diversitree:::make.all.branches.dtlik(cache, control,
#                                          initial.conditions)

# diversitree:::make.all.branches.dtlik
make.all.branches.dtlik <- function(cache, control, initial.conditions) 
	{
    branches <- make.branches.dtlik(cache$info, control)
    # Return this function
    function(pars, intermediates, preset = NULL) all.branches.matrix(pars, 
        cache, initial.conditions, branches, preset)
	}


# diversitree:::make.branches.classe
make.branches.classe <- function(cache, control) 
	{
	make.branches.dtlik(cache$info, control)
	}
	

# diversitree:::make.branches.dtlik

# Make the function to calculate the likelihood at the 
# time-steps (dt) along a branch
#info = cache$info
control = list()
make.branches.dtlik <- function(info, control) 
	{
    info <- diversitree:::check.info.ode(info, control)
    comp.idx <- info$idx.d
    eps <- control$eps
    ode <- diversitree:::make.ode(info, control)
    branches <- function(y, len, pars, t0, idx) 
    	{
    	ode(vars=y, times=c(t0, t0 + len), pars=pars2)
    	}
    make.branches.comp(branches, comp.idx, eps)
	}


#diversitree:::make.ode
make.ode <- function(info, control) 
	{
    control <- diversitree:::check.control.ode(control)
    info <- diversitree:::check.info.ode(info, control)
    backend <- control$backend
    if (backend == "gslode") 
        ode <- diversitree:::make.ode.gslode(info, control)
    else if (backend == "deSolve") 
        ode <- diversitree:::make.ode.deSolve(info, control)
    else stop("Invalid backend", backend)
    ode
	}



#diversitree:::make.ode.gslode
make.ode.gslode <- function(info, control) 
	{
    n.var <- info$ny
    n.par <- info$np
    rtol <- atol <- control$tol
    stepper <- control$gsl.stepper
    if (length(rtol) != 1) 
        stop("rtol must (currently) be scalar")
    if (length(atol) != 1) 
        stop("atol must (currently) be scalar")
    time.varying <- isTRUE(info$time.varying)
    if (time.varying) 
        tm <- info$tm
    if (control$compiled)
    	{
        model <- info$name.ode
        dll <- info$dll
        derivs <- sprintf("derivs_%s_gslode", model)
        derivs <- getNativeSymbolInfo(derivs, PACKAGE = dll)$address
        if (time.varying)
        	{ 
            ode2 <- new(diversitree:::GslOdeTime, derivs, n.var, tm)
            } else {
            ode2 <- new(diversitree:::GslOdeCompiled, derivs, n.var)
            }
    	} else {
        derivs <- info$derivs
        ode2b <- new(diversitree:::GslOdeR, derivs, environment(derivs), n.var)
	    }
    ode$set_control(list(atol = atol, rtol = rtol, algorithm = stepper, 
        hini = 1e-04))
    do.set.tm <- time.varying && !control$compiled
    
    # Return this function
    function(vars, times, pars)
    	{
        if (length(pars) != n.par) 
            stop("Incorrect parameter length")
        if (length(vars) != n.var) 
            stop("Incorrect variable length")
        if (length(times) <= 1) 
            stop("Need >= 2 times")
        if (do.set.tm) 
            tm$set(pars)
        ode$run(times, vars, pars)
    	}
	}

# Via:
# ode2 <- new(diversitree:::GslOdeCompiled, derivs, n.var)
# ode2
# C++ object <0x7ffe8a62ff90> of class 'GslOdeCompiled' <0x7ffe9c2183c0>
# ode2$run
# Class method definition for method run()
'
ode2_run <- function(...) 
	{
    " Rcpp::Matrix<14, Rcpp::PreserveStorage> run(std::__1::vector<double, std::__1::allocator<double> >, std::__1::vector<double, std::__1::allocator<double> >, SEXP)  \n   "
    .External(list(name = "CppMethod__invoke_notvoid", address = <pointer: 0x7ffe9c40c4a0>, 
        dll = list(name = "Rcpp", path = "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/Rcpp/libs/Rcpp.so", 
            dynamicLookup = TRUE, handle = <pointer: 0x7ffe9ff190b0>, 
            info = <pointer: 0x10e980240>), numParameters = -1L), 
        <pointer: 0x7ffe9c2183c0>, <pointer: 0x7ffe9c20dd50>, 
        .pointer, ...)
	}
'

# ode2 <- new(diversitree:::GslOdeCompiled, derivs, n.var)
#ode2$run(times, vars, pars)
# > ode2$run(times, vars, pars)
#           [,1]      [,2]      [,3]
# [1,] 0.0000000 0.0000000 0.0000000
# [2,] 0.0000000 0.0000000 0.0000000
# [3,] 0.0000000 0.0000000 0.0000000
# [4,] 0.0000000 0.0000000 0.0000000
# [5,] 0.0000000 0.0000000 0.0000000
# [6,] 0.7927699 0.7927699 0.6284842
# [7,] 0.0000000 0.0000000 0.0000000
# [8,] 0.0000000 0.0000000 0.0000000


# diversitree:::GslOdeTime
# C++ class 'GslOdeTime' <0x7ffe9c2184c0>
# Constructors:
#     GslOdeTime(SEXP, int, TimeMachine)
# 
# Fields: 
#     int size [readonly]
# 
# Methods: 
#      std::__1::vector<double, std::__1::allocator<double> > derivs(double, std::__1::vector<double, std::__1::allocator<double> >, SEXP)  
#            
#      Rcpp::Matrix<14, Rcpp::PreserveStorage> run(std::__1::vector<double, std::__1::allocator<double> >, std::__1::vector<double, std::__1::allocator<double> >, SEXP)  
#            
#      void set_control(Rcpp::List)  



# diversitree:::make.ode.deSolve
#diversitree:::make.ode.deSolve(info, control)
make.ode.deSolve <- function(vars, times, pars) 
function(info, control) 
	{
    if (!is.function(info$derivs)) 
        stop("info$derivs must be a function")
    derivs <- derivs.for.deSolve(info$derivs)
    rtol <- atol <- control$tol
    if (isTRUE(info$time.varying))
    	{
        tm <- info$tm
        function(vars, times, pars)
        	{
            tm$set(pars)
            lsoda.trim(vars, times, derivs, pars, rtol = rtol, 
                atol = atol)
        	}
    	} else {
        function(vars, times, pars) {
            lsoda.trim(vars, times, derivs, pars, rtol = rtol, 
                atol = atol)
        	}
    	}
	}


diversitree:::derivs.for.deSolve
derivs.for.deSolve <- function(f)
	{
	function(...) list(f(...))
	}


# diversitree:::lsoda.trim
lsoda.trim <- function(...) 
	{
    ret <- t(lsoda(...)[-1, -1, drop = FALSE])
    dimnames(ret) <- NULL
    ret
	}





# diversitree:::make.branches.comp
make.branches.comp <- function(branches, comp.idx, eps = 0) 
	{
    if (length(comp.idx) > 0)
    	{
    	# Function you are returning
        function(y, len, pars, t0, idx)
        	{
            ret <- branches(y, len, pars, t0, idx)
            q <- colSums(ret[comp.idx, , drop = FALSE])
            if (all(q >= eps))
            	{
                i <- q > 0
                ret[comp.idx, i] <- ret[comp.idx, i]/rep(q[i], 
                  each = length(comp.idx))
                lq <- q
                lq[i] <- log(q[i])
                list(lq, ret)
            	} else {
                ti <- len[length(len)]/2
                len1 <- c(len[len <= ti], ti)
                len2 <- len[len > ti] - ti
                n1 <- length(len1)
                # Recall = recursive call of the parent function
                ret1 <- Recall(y, len1, pars, t0)
                ret2 <- Recall(ret1[[2]][, n1], len2, pars, t0 + 
                  ti)
                ret2[[1]] <- ret2[[1]] + ret1[[1]][n1]
                list(c(ret1[[1]][-n1], ret2[[1]]), cbind(ret1[[2]][, 
                  -n1], ret2[[2]]))
            	} # END if (all(q >= eps))
            } # END function(y, len, pars, t0, idx)
        } else {
    	# Function you are returning
        function(y, len, pars, t0, idx)
        	{
        	list( rep.int(0, length(len)), branches(y, len, pars, t0, idx) )
        	} # END function you are returning
        } # END if (length(comp.idx) > 0) 
	}





# f.pars <- make.pars.classe(k)

# diversitree:::make.pars.classe
# k = # of states

# Make a function that returns the 
# lambda and mu parameters, and the Q matrix
# All parsing based on k (number of states)
make.pars.classe <- function(k)
	{
	# Number of parameters (all)
    np0 <- as.integer((k + 3) * k * k/2)
    
    # Add another 4 parameters?
    np <- np0 + k
    
    # Qmat: numstates by numstates
    qmat <- matrix(0, k, k)
    
    # row, column indexes of non-diagonal elements
    idx.qmat <- cbind(rep(1:k, each = k - 1), unlist(lapply(1:k, 
        function(i) (1:k)[-i])))
    idx.qmat
    
    # Number of lambda & mu (lm) parameters
    x <- k * k * (k + 1)/2 + k
    x
    idx.lm <- seq_len(x)
    
    # Extract the q parameters
    idx.q <- seq(x + 1, np0)
    
    # Return this function
    function(pars)
    	{
        check.pars.classe(pars, k)
        qmat[idx.qmat] <- pars[idx.q]
        diag(qmat) <- -rowSums(qmat)
        
        # Return the lambda and mu parameters, and the Q matrix
        c(pars[idx.lm], qmat)
    	}
	}



# diversitree:::check.pars.classe

check.pars.classe <- function(pars, k) 
	{
	check.pars.nonnegative(pars, (k + 3) * k * k/2)
	}


# diversitree:::check.pars.nonnegative
check.pars.nonnegative <- function(pars, npar) 
	{
    if (length(pars) != npar)
    	{
        stop(sprintf("Incorrect parameter length: expected %d, got %d", npar, length(pars)))
        }
    if (any(!is.finite(pars)) || any(pars < 0)) 
    	{
        stop("Parameters must be non-negative and finite")
        }
    pars
	} # check.pars.nonnegative


