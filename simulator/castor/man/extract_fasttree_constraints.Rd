\name{extract_fasttree_constraints}
\alias{extract_fasttree_constraints}
\title{
Extract tree constraints in FastTree alignment format.
}
\description{
Given a rooted phylogenetic tree, extract binary constraints in FastTree alignment format. Every internal bifurcating node with more than 2 descending tips will constitute an separate constraint.
}
\usage{
extract_fasttree_constraints(tree)
}
\arguments{
\item{tree}{
A rooted tree of class "phylo".
}
}
\details{
This function can be used to define constraints based on a backbone subtree, to be used to generate a larger tree using FastTree (as of v2.1.11). Only bifurcating nodes with at least 3 descending tips are used as constraints. 

The constraints are returned as a 2D matrix; the actual fasta file with the constraint alignments can be written easily from this matrix. For more details on FastTree constraints see the original FastTree documentation.
}
\value{
A list with the following elements:
\item{Nconstraints}{
Integer, specifying the number of constraints extracted.
}
\item{constraints}{
2D character matrix of size Ntips x Nconstraints, with values '0', '1' or '-', specifying which side ("left" or "right") of a constraint (node) each tip is found on.
}
\item{constraint2node}{
Integer vector of size Nconstraints, with values in 1,..,Nnodes, specifying the original node index used to define each constraint.
}
}


\author{Stilianos Louca}

%\seealso{
%}

\examples{
# generate a simple rooted tree, with tip names tip.1, tip.2, ...
Ntips = 10
tree  = generate_random_tree(list(birth_rate_intercept=1),
                             max_tips=Ntips,
                             tip_basename="tip.")$tree

# extract constraints
constraints = castor::extract_fasttree_constraints(tree)$constraints

# print constraints to screen in fasta format
cat(paste(sapply(1:Ntips, 
    FUN=function(tip) sprintf(">\%s\n\%s\n",tree$tip.label[tip],
    paste(as.character(constraints[tip,]),collapse=""))),collapse=""))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{MRCA}
