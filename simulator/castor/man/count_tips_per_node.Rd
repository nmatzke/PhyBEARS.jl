\name{count_tips_per_node}
\alias{count_tips_per_node}
\title{
Count descending tips.
}
\description{
Given a rooted phylogenetic tree, count the number of tips descending (directy or indirectly) from each node.
}
\usage{
count_tips_per_node(tree)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
}
\details{
The asymptotic time complexity of this function is O(Nedges), where Nedges is the number of edges.
}
\value{
An integer vector of size Nnodes, with the i-th entry being the number of tips descending (directly or indirectly) from the i-th node.
}


\author{Stilianos Louca}


\seealso{
\code{\link{get_subtree_at_node}}
}

\examples{
# generate a tree using a simple speciation model 
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=1000)$tree

# count number of tips descending from each node
tips_per_node = count_tips_per_node(tree);

# plot histogram of tips-per-node
barplot(table(tips_per_node[tips_per_node<10]), xlab="# tips", ylab="# nodes")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{MRCA}
