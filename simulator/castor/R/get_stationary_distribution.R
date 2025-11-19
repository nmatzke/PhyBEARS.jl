# calculate the stationary distribution of a continuous-time Markov model
# This is the probability vector p that satisfies p * Q = 0.
get_stationary_distribution = function(Q){
    temp_objective_function = function(x){
        p = c(x,1-sum(x)); # add one more element to x, making sure x is a probability vector
        return(sum((p %*% Q)^2));
    }
    Nstates = nrow(Q);
    if(Nstates>2){ 
        fit = stats::optim(rep(1.0/Nstates, times=Nstates-1), temp_objective_function, control=list(factr=1e4), lower=rep(0,times=Nstates-1), upper=rep(1,times=Nstates-1), method="L-BFGS-B")
        p 	= c(fit$par,1-sum(fit$par))	# add one more element to optimum, to turn it into a probability vector of size Nstates
    }else{
        fit = stats::optimize(temp_objective_function, interval=c(0,1))
        p 	= c(fit$minimum,1-sum(fit$minimum))	# add one more element to optimum, to turn it into a probability vector of size Nstates
    }
	return(p);
}