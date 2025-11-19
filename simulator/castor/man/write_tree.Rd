\name{write_tree}
\alias{write_tree}
\title{
Write a tree in Newick (parenthetic) format.
}
\description{
Write a phylogenetic tree to a file or a string, in Newick (parenthetic) format. If the tree is unrooted, it is first rooted internally at the first node.
}
\usage{
write_tree (tree, 
            file                 = "",
            append               = FALSE,
            digits               = 10,
            quoting              = 0,
            include_edge_labels  = FALSE,
            include_edge_numbers = FALSE)
}

\arguments{
\item{tree}{
A tree of class "phylo".
}
\item{file}{
An optional path to a file, to which the tree should be written. The file may be overwritten without warning. If left empty (default), then a string is returned representing the tree.
}
\item{append}{
Logical, specifying whether the tree should be appended at the end of the file, rather than replacing the entire file (if it exists).
}
\item{digits}{
Integer, number of significant digits for writing edge lengths.
}
\item{quoting}{
Integer, specifying whether and how to quote tip/node/edge names, as follows: 0:no quoting at all, 1:always use single quotes, 2:always use double quotes, -1:only quote when needed and prefer single quotes if possible, -2:only quote when needed and prefer double quotes if possible.
}
\item{include_edge_labels}{
Logical, specifying whether to include edge labels (if available) in the output tree, inside square brackets. Note that this is an extension (Matsen et al. 2012) to the standard Newick format, as, and edge labels in square brackets may not be supported by all Newick readers.
}
\item{include_edge_numbers}{
Logical, specifying whether to include edge numbers (if available) in the output tree, inside curly braces. Note that this is an extension (Matsen et al. 2012) to the standard Newick format, and edge numbers in curly braces may not be supported by all Newick readers.
}
}


\details{
If your tip and/or node and/or edge labels contain special characters (round brackets, commas, colons or quotes) then you should set \code{quoting} to non-zero, as appropriate. 

If the \code{tree} contains edge labels (as a character vector named \code{edge.label}) and \code{include_edge_labels==TRUE}, then edge labels are written in square brackets (Matsen et al. 2012). If \code{tree} contains edge numbers (as an integer vector named \code{edge.number}) and \code{include_edge_numbers==TRUE}, then edge numbers are written in curly braces (Matsen et al. 2012).

This function is comparable to (but typically much faster than) the \code{ape} function \code{write.tree}.
}

\value{
If \code{file==""}, then a string is returned containing the Newick representation of the tree. Otherwise, the tree is directly written to the file and no value is returned.
}


\author{Stilianos Louca}

\references{
Frederick A. Matsen et al. (2012). A format for  phylogenetic placements. PLOS One. 7:e31009
}

\seealso{
\code{\link{read_tree}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=100)$tree

# obtain a string representation of the tree in Newick format
Newick_string = write_tree(tree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }
