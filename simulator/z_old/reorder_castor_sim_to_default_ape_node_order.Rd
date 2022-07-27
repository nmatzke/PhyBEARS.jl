\name{reorder_castor_sim_to_default_ape_node_order}
\alias{reorder_castor_sim_to_default_ape_node_order}
\title{
Reorder a castor simulation to default APE order
}
\description{
Reorder a castor simulation from e.g. simulate_tdsse2() to the default APE node order
}
\usage{
reorder_castor_sim_to_default_ape_node_order(simulation)
}
\arguments{
\item{simulation}{
A simulation object from e.g. simulate_tdsse2().
}
}


\details{
Procedure: 

1. prints simulation$tree to a newick file, then reads it back in

2. Reorders the simulation$tip_states and simulation$node_states to match 

3. Adds simulation$states, to put all the states in one place
}


\value{
A numeric vector of length Nnodes, listing the RED of each node in the tree. The REDs of tips are not included, since these would always be equal to 1.
}

\author{Stilianos Louca}

\references{
D. H. Parks, M. Chuvochina et al. (2018). A proposal for a standardized bacterial taxonomy based on genome phylogeny. bioRxiv 256800. DOI:10.1101/256800
}

%\seealso{
%}

\examples{
# generate a random tree
params = list(birth_rate_intercept=1, death_rate_intercept=0.8)
tree = generate_random_tree(params, max_time=100, coalescent=FALSE)$tree

# calculate and print REDs
REDs = reorder_castor_sim_to_default_ape_node_order(tree)
print(REDs)
}
\keyword{relative evolutionary divergence}
