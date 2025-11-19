\name{simulate_sbm}
\alias{simulate_sbm}
\title{
Simulate Spherical Brownian Motion on a tree.
}
\description{
Given a rooted phylogenetic tree and a Spherical Brownian Motion (SBM) model for the evolution of the geographical location of a lineage on a sphere, simulate random outcomes of the model on all nodes and/or tips of the tree. The function traverses nodes from root to tips and randomly assigns a geographical location to each node or tip based on its parent's previously assigned location and the specified model parameters. The generated states have joint distributions consistent with the SBM model (Perrin 1928; Brillinger 2012). This function generalizes the simple SBM model to support time-dependent diffusivities.
}
\usage{
simulate_sbm(tree, 
             radius, 
             diffusivity,
             time_grid      = NULL,
             splines_degree = 1,
             root_latitude  = NULL, 
             root_longitude = NULL)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge. Edge lengths are assumed to represent time intervals or a similarly interpretable phylogenetic distance.
}
\item{radius}{
Strictly positive numeric, specifying the radius of the sphere. For Earth, the mean radius is 6371 km.
}
\item{diffusivity}{
Either a single numeric, or a numeric vector of length equal to that of \code{time_grid}. Diffusivity ("\eqn{D}") of the SBM model (in units distance^2/time). If \code{time_grid} is \code{NULL}, then \code{diffusivity} should be a single number specifying the time-independent diffusivity. Otherwise \code{diffusivity} specifies the diffusivity at each time point listed in \code{time_grid}.

Under a planar approximation the squared geographical distance of a node from the root will have expectation \eqn{4LD}, where \eqn{L} is the node's phylogenetic distance from the root. Note that distance is measured in the same units as the \code{radius} (e.g., km if the radius is given in km), and time is measured in the same units as the tree's edge lengths (e.g., Myr if edge lengths are given in Myr).
}
\item{time_grid}{
Numeric vector of the same length as \code{diffusivity} and listing times since the root in ascending order, or \code{NULL}. This can be used to specify a time-variable diffusivity (see details below). If \code{NULL}, the diffusivity is assumed to be constant over time and equal to \code{diffusivity} (which should be a single numeric). Time is measured in the same units as edge lengths, with root having time 0.
}
\item{splines_degree}{
Integer, either 0,1,2 or 3, specifying the polynomial degree of the provided \code{diffusivity} between grid points in \code{time_grid}. For example, if \code{splines_degree==1}, then the provided \code{diffusivity} is interpreted as a piecewise-linear curve; if \code{splines_degree==2} it is interpreted as a quadratic spline; if \code{splines_degree==3} it is interpreted as a cubic spline. The \code{splines_degree} influences the analytical properties of the curve, e.g. \code{splines_degree==1} guarantees a continuous curve, \code{splines_degree==2} guarantees a continuous curve and continuous derivative, and so on.
}
\item{root_latitude}{
The latitude of the tree's root, in decimal degrees, between -90 and 90. If NULL, the root latitude is chosen randomly according to the stationary probability distribution of the SBM.
}
\item{root_longitude}{
The longitude of the tree's root, in decimal degrees, between -180 and 180. If NULL, the root longitude is chosen randomly according to the stationary probability distribution of the SBM.
}
}


\details{
For short expected transition distances this function uses the approximation formula by Ghosh et al. (2012). For longer expected transition distances the function uses a truncated approximation of the series representation of SBM transition densities (Perrin 1928).

The pair \code{time_grid} and \code{diffusivity} can be used to define a time-dependent diffusivity, with time counted from the root to the tips (i.e. root has time 0) in the same units as edge lengths. For example, to define a diffusivity that varies linearly with time, you only need to specify the diffusivity at two time points (one at 0, and one at the time of the youngest tip), i.e. \code{time_grid} and \code{diffusivity} would each have length 2. Note that \code{time_grid} should cover the full time range of the tree; otherwise, \code{diffusivity} will be extrapolated as a constant when needed.

If \code{tree$edge.length} is missing, each edge in the tree is assumed to have length 1. The tree may include multifurcations as well as monofurcations.
}


\value{
A list with the following elements:
\item{success}{
Logical, specifying whether the simulation was successful. If \code{FALSE}, then an additional return variable \code{error} will contain a brief description of the error that occurred, and all other return variables may be undefined.
}
\item{tip_latitudes}{
Numeric vector of length Ntips, listing simulated decimal latitudes for each tip in the tree.
}
\item{tip_longitudes}{
Numeric vector of length Ntips, listing simulated decimal longitudes for each tip in the tree.
}
\item{node_latitudes}{
Numeric vector of length Nnodes, listing simulated decimal latitudes for each internal node in the tree.
}
\item{node_longitudes}{
Numeric vector of length Nnodes, listing simulated decimal longitudes for each internal node in the tree.
}
}

\author{Stilianos Louca}

\references{
F. Perrin (1928). Etude mathematique du mouvement Brownien de rotation. 45:1-51.

D. R. Brillinger (2012). A particle migrating randomly on a sphere. in Selected Works of David Brillinger. Springer.

A. Ghosh, J. Samuel, S. Sinha (2012). A Gaussian for diffusion on the sphere. Europhysics Letters. 98:30003.

S. Louca (2021). Phylogeographic estimation and simulation of global diffusive dispersal. Systematic Biology. 70:340-359.
}

\seealso{
\code{\link{simulate_ou_model}}, 
\code{\link{simulate_rou_model}},
\code{\link{simulate_bm_model}},
\code{\link{fit_sbm_const}}
}

\examples{
\dontrun{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=100)$tree

# simulate SBM on the tree
simulation = simulate_sbm(tree, radius=6371, diffusivity=1e4,
                          root_latitude=0, root_longitude=0)

# plot latitudes and longitudes of the tips
plot(simulation$tip_latitudes,simulation$tip_longitudes)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{SBM model}
\keyword{random}
\keyword{simulation}
