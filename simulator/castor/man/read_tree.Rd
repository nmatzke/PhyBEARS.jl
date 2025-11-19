\name{read_tree}
\alias{read_tree}
\title{
Load a tree from a string or file in Newick (parenthetic) format.
}
\description{
Load a phylogenetic tree from a file or a string, in Newick (parenthetic) format. Any valid Newick format is acceptable. Extended variants including edge labels and edge numbers are also supported.
}
\usage{
read_tree(  string  = "", 
            file    = "", 
            edge_order              = "cladewise",
            include_edge_lengths    = TRUE, 
            look_for_edge_labels    = FALSE, 
            look_for_edge_numbers   = FALSE, 
            include_node_labels     = TRUE, 
            underscores_as_blanks   = FALSE, 
            check_label_uniqueness  = FALSE,
            interpret_quotes        = FALSE,
            trim_white              = TRUE)
}

\arguments{
\item{string}{
A character containing a single tree in Newick format. Can be used alternatively to \code{file}.
}
\item{file}{
Character, a path to an input text file containing a single tree in Newick format. Can be used alternatively to \code{string}.
}
\item{edge_order}{
Character, one of ``cladewise'' or ``pruningwise'', specifying the order in which edges should be listed in the returned tree. This does not influence the topology of the tree or the tip/node labeling, it only affects the way edges are numbered internally.
}
\item{include_edge_lengths}{
Logical, specifying whether edge lengths (if available) should be included in the returned tree.
}
\item{look_for_edge_labels}{
Logical, specifying whether edge labels may be present in the input tree. If edge labels are found, they are included in the returned tree as a character vector \code{edge.label}. Edge labels are sought inside square brackets, which are not part of the standard Newick format but used by some tree creation software (Matsen 2012). If \code{look_for_edge_labels==FALSE}, square brackets are read verbatim just like any other character.
}
\item{look_for_edge_numbers}{
Logical, specifying whether edge numbers (non-negative integers) may be present in the input tree. If edge numbers are found, they are included in the returned tree as an integer vector \code{edge.number}. Edge numbers are sought inside curly braces, which are not part of the standard Newick format but used by some tree creation software (Matsen 2012). If \code{look_for_edge_numbers==FALSE}, curly braces are read verbatim just like any other character.
}
\item{include_node_labels}{
Logical, specifying whether node labels (if available) should be included in the returned tree.
}
\item{underscores_as_blanks}{
Logical, specifying whether underscores ("_") in tip and node labels should be replaced by spaces (" "). This is common behavior in other tree parsers. In any case, tip, node and edge labels (if available) are also allowed to contain explicit whitespace (except for newline characters).
}
\item{check_label_uniqueness}{
Logical, specifying whether to check if all tip labels are unique.
}
\item{interpret_quotes}{
Logical, specifying whether to interpret quotes as delimiters of tip/node/edge labels. If \code{FALSE}, then quotes are read verbatim just like any other character.
}
\item{trim_white}{
Logical, specifying whether to trim flanking whitespace from tip, node and edge labels.
}
}


\details{
This function is comparable to (but typically much faster than) the \code{ape} function \code{read.tree}. The function supports trees with monofurcations and multifurcations, trees with or without tip/node labels, and trees with or without edge lengths. The time complexity is linear in the number of edges in the tree.

Either \code{file} or \code{string} must be specified, but not both. The tree may be arbitrarily split across multiple lines, but no other non-whitespace text is permitted in \code{string} or in the input file. Flanking whitespace (space, tab, newlines) is ignored.
}

\value{
A single rooted phylogenetic tree in ``phylo'' format.
}


\author{Stilianos Louca}

\references{
Frederick A. Matsen et al. (2012). A format for  phylogenetic placements. PLOS One. 7:e31009
}

\seealso{
\code{\link{write_tree}}
}

\examples{
# generate a random tree
tree = generate_random_tree(list(birth_rate_intercept=1),max_tips=100)$tree

# obtain a string representation of the tree in Newick format
Newick_string = write_tree(tree)

# re-parse tree from string
parsed_tree = read_tree(Newick_string)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Pruning }
