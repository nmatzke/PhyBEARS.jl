\name{spline_coefficients}
\alias{spline_coefficients}
\title{
Get the polynomial coefficients of a spline.
}
\description{
Given a natural spline function \eqn{Y:\R\to\R}, defined as a series of Y values on a discrete X grid, obtain its corresponding piecewise polynomial coefficients. Supported splines degrees are 0 (Y is piecewise constant), 1 (piecewise linear), 2 (piecewise quadratic) and 3 (piecewise cubic).
}
\usage{
spline_coefficients(Xgrid, 
                    Ygrid,
                    splines_degree)
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
}


\details{
Spline functions are returned by some of castor's fitting routines, so \code{spline_coefficients} is meant to aid with the further analysis of such functions. A spline function of degree \eqn{D\geq1} has continuous derivatives up to degree \eqn{D-1}.
}


\value{
A numeric matrix of size NR x NC, where NR (number of rows) is equal to the length of \code{Xgrid} and NC (number of columns) is equal to \code{splines_degree+1}. The r-th row lists the polynomial coefficients (order 0, 1 etc) of the spline within the interval [Xgrid[r],Xgrid[r+1]]. For exampe, for a spline of order 2, the value at X=0.5*(Xgrid[1]+Xgrid[2]) will be equal to C[1,1]+C[1,2]*X+C[1,3]*X*X, where C is the matrix of coefficients.
}

\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{evaluate_spline}}
}


\examples{
# The following code defines a quadratic spline on 20 grid points
# The curve's polynomial coefficients are then determined
# and used to evaluate the spline on a fine grid for plotting.
Xgrid = seq(0,10,length.out=20)
Ygrid = sin(Xgrid)
splines_degree = 2

Ycoeff = castor::spline_coefficients(Xgrid, Ygrid, splines_degree)

plot(Xgrid, Ygrid, type='p')

for(g in seq_len(length(Xgrid)-1)){
  Xtarget = seq(Xgrid[g], Xgrid[g+1], length.out=100)
  Ytarget = rep(Ycoeff[g,1], length(Xtarget))
  for(p in seq_len(splines_degree)){
    Ytarget = Ytarget + (Xtarget^p) * Ycoeff[g,p+1];
  }
  lines(Xtarget, Ytarget, type='l', col='red')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{random}
