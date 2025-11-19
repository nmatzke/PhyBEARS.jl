# given a piecewise polynomial (natural splines) function f(x), defined as a time series on some x-grid, calculate its values on an arbitrary set of x-values
# this function is most efficient when the requested target x-values are monotonically increasing or decreasing
evaluate_spline = function(	Xgrid,					# numeric vector of size NG, listing x-values in ascending order
							Ygrid,					# numeric vector of size NG, listing y-values along Xgrid
							splines_degree,			# integer, either 0,1,2 or 3, specifying the splines degree assumed for Y. E.g. 0 means piecewise constant, 1 means piecewise linear etc
							Xtarget,				# numeric vector of size N, specifying the target x values on which to evaluate the function Y. The function is most efficient if Xtarget are in ascending or descending order.
							extrapolate="const",	# character, specifying how to extrapolate Y beyond Xgrid if needed. Available options are "const" or "splines" (i.e. use the polynomial coefficients from the nearest grid point)
							derivative=0){			# integer, one of 0,1,2, specifying which derivative to return.
	# basic error checking
	if(!(splines_degree %in% c(0,1,2,3))) stop("splines_degree must be 0, 1, 2 or 3");
	if(!(derivative %in% c(0,1,2))) stop("derivative must be 0, 1 or 2");
	if(Xgrid[1]>tail(Xgrid,1)) stop("Xgrid must be in ascending order")

	Ytarget = evaluate_spline_CPP(Xgrid=Xgrid, Ygrid=Ygrid, splines_degree=splines_degree, Xtarget=Xtarget, extrapolate=extrapolate, derivative=derivative);
	return(Ytarget)
}