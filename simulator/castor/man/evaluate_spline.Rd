\name{evaluate_spline}
\alias{evaluate_spline}
\title{
Evaluate a scalar spline at arbitrary locations.
}
\description{
Given a natural spline function \eqn{Y:\R\to\R}, defined as a series of Y values on a discrete X grid, evaluate its values (or derivative) at arbitrary X points. Supported splines degrees are 0 (Y is piecewise constant), 1 (piecewise linear), 2 (piecewise quadratic) and 3 (piecewise cubic).
}
\usage{
evaluate_spline(Xgrid, 
                Ygrid,
                splines_degree,
                Xtarget,
                extrapolate = "const",
                derivative  = 0)
}
\arguments{
\item{Xgrid}{
Numeric vector, listing x-values in ascending order.
}
\item{Ygrid}{
Numeric vector of the same length as \code{Xgrid}, listing the values of Y on \code{Xgrid}.
}
\item{splines_degree}{
Integer, either 0, 1, 2 or 3, specifying the polynomial degree of the spline curve Y between grid points. For example, 0 means Y is piecewise constant, 1 means Y is piecewise linear and so on.
}
\item{Xtarget}{
Numeric vector, listing arbitrary X values on which to evaluate Y.
}
\item{extrapolate}{
Character, specifying how to extrapolate Y beyond \code{Xgrid} if needed. Available options are "const" (i.e. use the value of Y on the nearest \code{Xgrid} point) or "splines" (i.e. use the polynomial coefficients from the nearest grid point).
}
\item{derivative}{
Integer, specifying which derivative to return. To return the spline's value, set \code{derivative=0}. Currently only the options 0,1,2 are supported.
}
}


\details{
Spline functions are returned by some of castor's fitting routines, so \code{evaluate_spline} is meant to aid with the evaluation and plotting of such functions. A spline function of degree \eqn{D\geq1} has continuous derivatives up to degree \eqn{D-1}. The function \code{evaluate_spline} is much more efficient if \code{Xtarget} is monotonically increasing or decreasing.

This function is used to evaluate the spline's values at arbitrary points. To obtain the spline's polynomial coefficients, use \code{\link{spline_coefficients}}.
}


\value{
A numeric vector of the same length as \code{Xtarget}, listing the values (or derivatives, if \code{derivative>0}) of Y on \code{Xtarget}.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{spline_coefficients}}
}

\examples{
# specify Y on a coarse X grid
Xgrid = seq(from=0,to=10,length.out=10)
Ygrid = sin(Xgrid)

# define a fine grid of target X values
Xtarget = seq(from=0,to=10,length.out=1000)

# evaluate Y on Xtarget, either as piecewise linear or piecewise cubic function
Ytarget_lin = evaluate_spline(Xgrid,Ygrid,splines_degree=1,Xtarget=Xtarget)
Ytarget_cub = evaluate_spline(Xgrid,Ygrid,splines_degree=3,Xtarget=Xtarget)

# plot both the piecewise linear and piecewise cubic curves
plot(x=Xtarget, y=Ytarget_cub, type='l', col='red', xlab='X', ylab='Y')
lines(x=Xtarget, y=Ytarget_lin, type='l', col='blue', xlab='X', ylab='Y')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{random}
