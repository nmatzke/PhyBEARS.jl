# given a piecewise polynomial (natural splines) function f(x), defined as a time series on some x-grid, determine its piecewise polynomial coefficients
# this function is most efficient when the requested target x-values are monotonically increasing or decreasing
spline_coefficients = function(	Xgrid,				# numeric vector of size NG, listing x-values in ascending order
								Ygrid,				# numeric vector of size NG, listing y-values along Xgrid
								splines_degree){	# integer, either 0,1,2 or 3, specifying the splines degree assumed for Y. E.g. 0 means piecewise constant, 1 means piecewise linear etc
	# basic error checking
	if(!(splines_degree %in% c(0,1,2,3))) stop("splines_degree must be 0, 1, 2 or 3");
	if(any(diff(Xgrid)<=0)) stop("Xgrid must be in ascending order")

	Ycoeff = get_spline_CPP(Xgrid=Xgrid, Ygrid=Ygrid, splines_degree=splines_degree)
	Ycoeff = matrix(Ycoeff, ncol=splines_degree+1, byrow=TRUE)
	
	return(Ycoeff)
}