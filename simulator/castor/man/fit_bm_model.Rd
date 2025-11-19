\name{fit_bm_model}
\alias{fit_bm_model}
\title{
Fit a Brownian Motion model for multivariate trait evolution.
}
\description{
Given a rooted phylogenetic tree and states of one or more continuous (numeric) traits on the tree's tips, fit a multivariate Brownian motion model of correlated co-evolution of these traits. This estimates a single diffusivity matrix, which describes the variance-covariance structure of each trait's random walk. The model assumes a fixed diffusivity matrix on the entire tree. Optionally, multiple trees can be used as input, under the assumption that the trait evolved on each tree according to the same BM model.
}
\usage{
fit_bm_model(   trees, 
                tip_states, 
                isotropic   = FALSE,
                Nbootstraps = 0,
                check_input = TRUE)
}
\arguments{
\item{trees}{
Either a single rooted tree or a list of rooted trees, of class "phylo". The root of each tree is assumed to be the unique node with no incoming edge. Edge lengths are assumed to represent time intervals or a similarly interpretable phylogenetic distance.
}
\item{tip_states}{
Numeric state of each trait at each tip in each tree. If \code{trees} was a single tree, then \code{tip_states} must either be a numeric vector of size Ntips or a 2D numeric matrix of size Ntips x Ntraits, listing the trait states for each tip in the tree.
If \code{trees} is a list of Ntrees trees, then \code{tip_states} must be a list of length Ntrees, each element of which lists the trait states for the corresponding tree (as a vector or 2D matrix, similarly to the single-tree case).
}
\item{isotropic}{
Logical, specifying whether diffusion should be assumed to be isotropic (i.e., independent of the direction). Hence, if \code{isotropic=TRUE}, then the diffusivity matrix is forced to be diagonal, with all entries being equal. If \code{isotropic=FALSE}, an arbitrary diffusivity matrix is fitted (i.e., the diffusivity matrix is only constrained to be symmetric and non-negative definite).
}
\item{Nbootstraps}{
Integer, specifying the number of parametric bootstraps to perform for calculating the confidence intervals. If <=0, no bootstrapping is performed.
}
\item{check_input}{
Logical, specifying whether to perform some basic checks on the validity of the input data. If you are certain that your input data are valid, you can set this to \code{FALSE} to reduce computation.
}
}


\details{
The BM model is defined by the stochastic differential equation
\deqn{
dX = \sigma \cdot dW
}
where \eqn{W} is a multidimensional Wiener process with Ndegrees independent components and \eqn{\sigma} is a matrix of size Ntraits x Ndegrees. Alternatively, the same model can be defined as a Fokker-Planck equation for the probability density \eqn{\rho}:
\deqn{
\frac{\partial \rho}{\partial t} = \sum_{i,j}D_{ij}\frac{\partial^2\rho}{\partial x_i\partial x_j}.
}
The matrix \eqn{D} is referred to as the diffusivity matrix (or diffusion tensor), and \eqn{2D=\sigma\cdot\sigma^T}. Note that \eqn{\sigma} can be obtained from \eqn{D} by means of a Cholesky decomposition.

The function uses phylogenetic independent contrasts (Felsenstein, 1985) to retrieve independent increments of the multivariate random walk. The diffusivity matrix \eqn{D} is then fitted using maximum-likelihood on the intrinsic geometry of positive-definite matrices. If multiple trees are provided as input, then independent contrasts are extracted from all trees and combined into a single set of independent contrasts (i.e., as if they had been extracted from a single tree).

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. The tree may include multifurcations (i.e. nodes with more than 2 children) as well as monofurcations (i.e. nodes with only one child). Note that multifurcations are internally expanded to bifurcations, prior to model fitting.
}


\value{
A list with the following elements:
\item{success}{
Logical, indicating whether the fitting was successful. If \code{FALSE}, then an additional return variable, \code{error}, will contain a description of the error; in that case all other return variables may be undefined.
}
\item{diffusivity}{
Either a single non-negative number (if \code{tip_states} was a vector) or a 2D quadratic non-negative-definite matrix (if \code{tip_states} was a 2D matrix). The fitted diffusivity matrix of the multivariate Brownian motion model.
}
\item{loglikelihood}{
The log-likelihood of the fitted model, given the provided tip states data.
}
\item{Ncontrasts}{
Integer, number of independent contrasts used to estimate the diffusivity. This corresponds to the number of independent data points used.
}
\item{standard_errors}{
Either a single numeric or a 2D numeric matrix of size Ntraits x Ntraits, listing the estimated standard errors of the estimated diffusivity, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50lower}{
Either a single numeric or a 2D numeric matrix of size Ntraits x Ntraits, listing the lower bounds of the 50\% confidence interval for the estimated diffusivity (25-75\% percentile), based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI50upper}{
Either a single numeric or a 2D numeric matrix of size Ntraits x Ntraits, listing the upper bound of the 50\% confidence interval for the estimated diffusivity, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95lower}{
Either a single numeric or a 2D numeric matrix of size Ntraits x Ntraits, listing the lower bound of the 95\% confidence interval for the estimated diffusivity (2.5-97.5\% percentile), based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{CI95upper}{
Either a single numeric or a 2D numeric matrix of size Ntraits x Ntraits, listing the upper bound of the 95\% confidence interval for the estimated diffusivity, based on parametric bootstrapping. Only returned if \code{Nbootstraps>0}.
}
\item{consistency}{
Numeric between 0 and 1, estimated consistency of the data with the fitted model. If \eqn{L} denotes the loglikelihood of new data generated by the fitted model (under the same model) and \eqn{M} denotes the expectation of \eqn{L}, then \code{consistency} is the probability that \eqn{|L-M|} will be greater or equal to \eqn{|X-M|}, where \eqn{X} is the loglikelihood of the original data under the fitted model.
Only returned if \code{Nbootstraps>0}. A low consistency (e.g., <0.05) indicates that the fitted model is a poor description of the data. See Lindholm et al. (2019) for background.
}
}

\author{Stilianos Louca}

\references{
J. Felsenstein (1985). Phylogenies and the Comparative Method. The American Naturalist. 125:1-15.

A. Lindholm, D. Zachariah, P. Stoica, T. B. Schoen (2019). Data consistency approach to model validation. IEEE Access. 7:59788-59796.
}

\seealso{
\code{\link{simulate_bm_model}},
\code{\link{get_independent_contrasts}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1), 10000)$tree


# Example 1: Scalar case
# - - - - - - - - - - - - - -
# simulate scalar continuous trait on the tree
D = 1
tip_states = simulate_bm_model(tree, diffusivity=D)$tip_states

# estimate original diffusivity from the generated data
fit = fit_bm_model(tree, tip_states)
cat(sprintf("True D=\%g, fitted D=\%g\n",D,fit$diffusivity))


# Example 2: Multivariate case
# - - - - - - - - - - - - - - -
# simulate vector-valued continuous trait on the tree
D = get_random_diffusivity_matrix(Ntraits=5)
tip_states = simulate_bm_model(tree, diffusivity=D)$tip_states

# estimate original diffusivity matrix from the generated data
fit = fit_bm_model(tree, tip_states)

# compare true and fitted diffusivity matrices
cat("True D:\n"); print(D)
cat("Fitted D:\n"); print(fit$diffusivity)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BM model}
