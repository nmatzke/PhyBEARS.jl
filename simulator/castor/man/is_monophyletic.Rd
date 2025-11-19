\name{is_monophyletic}
\alias{is_monophyletic}
\title{
Determine if a set of tips is monophyletic.
}
\description{
Given a rooted phylogenetic tree and a set of focal tips, this function determines whether the tips form a monophyletic group.
}
\usage{
is_monophyletic(tree, focal_tips, check_input=TRUE)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo". The root is assumed to be the unique node with no incoming edge.
}
\item{focal_tips}{
Either an integer vector or a character vector, listing the tips to be checked for monophyly. If an integer vector, it should list tip indices (i.e. from 1 to Ntips). If a character vector, it should list tip names; in that case \code{tree$tip.label} must exist.
}
\item{check_input}{
Logical, whether to perform basic validations of the input data. If you know for certain that your input is valid, you can set this to \code{FALSE} to reduce computation time.
}
}
\details{
This function first finds the most recent common ancestor (MRCA) of the focal tips, and then checks if all tips descending from that MRCA fall within the focal tip set.
}
\value{
A logical, indicating whether the focal tips form a monophyletic set.
}


\author{Stilianos Louca}

%\references{
%}

\seealso{
\code{\link{get_mrca_of_set}}
}

\examples{
# generate random tree
Ntips = 100
tree  = generate_random_tree(list(birth_rate_intercept=1),max_tips=Ntips)$tree

# pick a random subset of focal tips
focal_tips = which(sample.int(2,size=Ntips,replace=TRUE)==1)

# check if focal tips form a monophyletic group
is_monophyletic(tree, focal_tips)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MRCA}
