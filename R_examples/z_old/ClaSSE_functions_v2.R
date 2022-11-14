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



# Made with make.bisse
# bisse_2areas = make.bisse(tree=tr, states=states, sampling.f=sampling.f, strict=FALSE)
# printed to screen with dput(bisse_2areas)


bisse_2areas_default <- function(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL, intermediates=FALSE) 
	{
	check.pars.bisse(pars)
	preset <- branches.unresolved.bisse(pars, unresolved)
	ans <- all.branches(pars, intermediates, preset)
	rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
	}




make_bd <- function(tree, sampling.f=NULL, unresolved=NULL, times=NULL, control = list()) 
	{
	ex='
	sampling.f=NULL
	unresolved=NULL
	times=NULL
	control = list()
	'

	control <- check.control.bd(control, times)
	cache <- make.cache.bd(tree, sampling.f, unresolved, times, control)
	const <- cache$const
	if (control$method == "nee")
		{
		all.branches <- make.all.branches.bd.nee(cache, control)
		rootfunc <- rootfunc.bd.nee
		} else {
		all.branches <- make.all.branches.dtlik(cache, control, initial.conditions.bd.ode)
		rootfunc <- rootfunc.bd.ode
		}
	ll <- function(pars, condition.surv = TRUE, intermediates = FALSE)
		{
		check.pars.nonnegative(pars, 2)
		ans <- all.branches(pars, intermediates)
		rootfunc(ans, pars, condition.surv, intermediates, const)
		}
	class(ll) <- c("bd", "dtlik", "function")
	return(ll)
	}


lik_bd_default <- function(pars, condition.surv=TRUE, intermediates=FALSE) 
	{
	check.pars.nonnegative(pars, 2)
	ans <- all.branches(pars, intermediates)
	rootfunc(ans, pars, condition.surv, intermediates, const)
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
	cache <- make.cache.bisse(tree, states, unresolved, sampling.f, nt.extra, strict)
	unresolved <- cache$unresolved
	all.branches <- make.all.branches.dtlik(cache, control, initial.conditions.bisse)
	rootfunc <- rootfunc.musse
	ll <- function(pars, condition.surv = TRUE, root = ROOT.OBS, root.p = NULL, intermediates = FALSE)
		{
		check.pars.bisse(pars)
		preset <- branches.unresolved.bisse(pars, unresolved)
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

# Birth-death root function
#diversitree:::rootfunc.bd.ode
rootfunc.bd.ode <- function(res, pars, condition.surv, intermediates, const) 
	{
	vals <- res$vals
	lq <- res$lq
	d.root <- vals[2]
	if (condition.surv) {
			e.root <- vals[[1]]
			lambda <- pars[[1]]
			d.root <- d.root/(lambda * (1 - e.root)^2)
	}
	loglik <- log(d.root) + sum(lq) + const
	names(loglik) <- NULL
	if (intermediates) {
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
	i <- which(evA$values == max(evA$values))
	evA$vectors[, i]/sum(evA$vectors[, i])
	}




# diversitree:::projection.matrix.classe
projection.matrix.classe <- function(pars, k) 
	{
	A <- matrix(0, nrow = k, ncol = k)
	nsum <- k * (k + 1)/2
	kseq <- seq_len(k)
	pars.lam <- pars[seq(1, nsum * k)]
	pars.mu <- pars[seq(nsum * k + 1, (nsum + 1) * k)]
	pars.q <- pars[seq((nsum + 1) * k + 1, length(pars))]
	idx.lam <- cbind(rep(kseq, each = nsum), rep(rep(kseq, times = seq(k, 1, -1)), k), unlist(lapply(kseq, function(i) i:k)))
	idx.q <- cbind(unlist(lapply(kseq, function(i) (kseq)[-i])), rep(kseq, each = k - 1))
	for (n in seq_len(nsum * k))
		{
		r <- idx.lam[n, ]
		A[r[2], r[1]] <- A[r[2], r[1]] + pars.lam[n]
		A[r[3], r[1]] <- A[r[3], r[1]] + pars.lam[n]
		}
	A[idx.q] <- A[idx.q] + pars.q
	diag(A) <- 0
	diag(A) <- -colSums(A) + unlist(lapply(kseq, function(i) sum(pars.lam[seq((i-1) * nsum + 1, i * nsum)]) - pars.mu[i]))
	A
	}









all.branches <- make.all.branches.dtlik(cache, control,
                                          initial.conditions)

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
info = cache$info
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


# ode2 <- new(diversitree:::GslOdeCompiled, derivs, n.var)
ode2$run(times, vars, pars)
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
diversitree:::make.ode.deSolve(info, control)
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


