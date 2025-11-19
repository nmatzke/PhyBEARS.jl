\name{count_transitions_between_clades}
\alias{count_transitions_between_clades}
\title{
Count the number of state transitions between tips or nodes.
}
\description{
Given a rooted phylogenetic tree, one or more pairs of tips and/or nodes, and the state of some discrete trait at each tip and node, calculate the number of state transitions along the shortest path between each pair of tips/nodes.
}
\usage{
count_transitions_between_clades(tree, A, B, states, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{A}{
An integer vector or character vector of size Npairs, specifying the first of the two members of each pair of tips/nodes. If an integer vector, it must list indices of tips (from 1 to Ntips) and/or nodes (from Ntips+1 to Ntips+Nnodes). If a character vector, it must list tip and/or node names.
}
\item{B}{
An integer vector or character vector of size Npairs, specifying the second of the two members of each pair of tips/nodes. If an integer vector, it must list indices of tips (from 1 to Ntips) and/or nodes (from Ntips+1 to Ntips+Nnodes). If a character vector, it must list tip and/or node names.
}
\item{states}{
Integer vector of length Ntips+Nnodes, listing the discrete state of each tip and node in the tree. The order of entries must match the order of tips and nodes in the tree; this requirement is only verified if \code{states} has names and \code{check_input==TRUE}.
}
\item{check_input}{
Logical, whether to perform basic validations of the input data. If you know for certain that your input is valid, you can set this to \code{FALSE} to reduce computation time.
}
}
\details{
The discrete state must be represented by integers (both negatives and positives are allowed); characters and other data types are not allowed. If tip/node states are originally encoded as characters rather than integers, you can use \code{map_to_state_space} to convert these to integers (for example ``\code{male}'' & ``\code{female}'' may be represented as 1 & 2). Also note that a state must be provided for each tip and ancestral node, not just for the tips. If you only know the states of tips, you can use an ancestral state reconstruction tool to estimate ancestral states first.

The tree may include multi-furcations as well as mono-furcations (i.e. nodes with only one child). If \code{A} and/or \code{B} is a character vector, then \code{tree$tip.label} must exist. If node names are included in \code{A} and/or \code{B}, then \code{tree$node.label} must also exist.
}
\value{
An integer vector of size Npairs, with the i-th element being the number of state transitions between tips/nodes \code{A[i]} and \code{B[i]} (along their shortest connecting path).
}


\author{Stilianos Louca}

%\references{
%}

%\seealso{
%}

\examples{
# generate a random tree
Ntips = 100
tree = generate_random_tree(list(birth_rate_intercept=1),Ntips)$tree

# pick 3 random pairs of tips or nodes
Npairs = 3
A = sample.int(n=(Ntips+tree$Nnode), size=Npairs, replace=FALSE)
B = sample.int(n=(Ntips+tree$Nnode), size=Npairs, replace=FALSE)

# assign a random state to each tip & node in the tree
# consider a binary trait
states = sample.int(n=2, size=Ntips+tree$Nnode, replace=TRUE)

# calculate number of transitions for each tip pair
Ntransitions = count_transitions_between_clades(tree, A, B, states=states)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Mk model}
