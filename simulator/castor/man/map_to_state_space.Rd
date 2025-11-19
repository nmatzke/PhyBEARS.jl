\name{map_to_state_space}
\alias{map_to_state_space}
\title{
Map states of a discrete trait to integers.
}
\description{
Given a list of states (e.g., for each tip in a tree), map the unique states to integers 1,..,Nstates, where Nstates is the number of possible states. This function can be used to translate states that are originally represented by characters or factors, into integer states as required by ancestral state reconstruction and hidden state prediction functions in this package.
}
\usage{
map_to_state_space(raw_states, fill_gaps=FALSE, 
                   sort_order="natural", include_state_values=FALSE)
}
\arguments{
\item{raw_states}{
A vector of values (states), each of which can be converted to a different character. This list can include the same value multiple times, for example if values represent the trait's states for tips in a tree.
}
\item{fill_gaps}{
Logical. If \code{TRUE}, then states are converted to integers using \code{as.integer(as.character())}, and then all missing intermediate integer values are included as additional possible states. For example, if \code{raw_states} contained the values 2,4,6, then 3 and 5 are assumed to also be possible states.
}
\item{sort_order}{
Character, specifying the order in which raw_states should be mapped to ascending integers. Either "natural" or "alphabetical". If "natural", numerical parts of characters are sorted numerically, e.g. as in "3"<"a2"<"a12"<"b1".
}
\item{include_state_values}{
Logical, specifying whether to also return a numerical version of the unique states. For example, the states "3","a2","4.5" will be mapped to the numeric values 3, NA, 4.5.
}
}

\details{
Several ancestral state reconstruction and hidden state prediction algorithms in the \code{castor} package (e.g., \code{asr_max_parsimony}) require that the focal trait's states are represented by integer indices within 1,..,Nstates. These indices are then associated, afor example, with column and row indices in the transition cost matrix (in the case of maximum parsimony reconstruction) or with column indices in the returned matrix containing marginal ancestral state probabilities (e.g., in \code{\link{asr_mk_model}}). The function \code{map_to_state_space} can be used to conveniently convert a set of discrete states into integers, for use with the aforementioned algorithms.
}

\value{
A list with the following elements:
\item{Nstates}{
Integer. Number of possible states for the trait, based on the unique values encountered in \code{raw_states} (after conversion to characters). This may be larger than the number of unique values in \code{raw_states}, if \code{fill_gaps} was set to \code{TRUE}.
}
\item{state_names}{
Character vector of size Nstates, storing the original name (character version) of each state. For example, if \code{raw_states} was \code{c("b1","3","a12","a2","b1","a2")} and \code{sort_order=="natural"}, then \code{Nstates} will be 4 and \code{state_names} will be \code{c("3","a2","a12","b1")}.
}
\item{state_values}{
Optional, only included if \code{include_state_values==TRUE}. A numeric vector of size \code{Nstates}, providing the numerical value for each unique state.
}
\item{mapped_states}{
Integer vector of size equal to \code{length(raw_states)}, listing the integer representation of each value in \code{raw_states}.
}
\item{name2index}{
An integer vector of size Nstates, with \code{names(name2index)} set to \code{state_names}. This vector can be used to map any new list of states (in character format) to their integer representation. In particular, \code{name2index[as.character(raw_states)]} is equal to \code{mapped_states}.
}
}

\author{Stilianos Louca}

%\references{
%}

%\seealso{
%\code{\link{ps_ouss}}
%}

\examples{
# generate a sequence of random states
unique_states = c("b","c","a")
raw_states = unique_states[sample.int(3,size=10,replace=TRUE)]

# map to integer state space
mapping = map_to_state_space(raw_states)

cat(sprintf("Checking that original unique states is the same as the one inferred:\n"))
print(unique_states)
print(mapping$state_names)

cat(sprintf("Checking reversibility of mapping:\n"))
print(raw_states)
print(mapping$state_names[mapping$mapped_states])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{algebra}
\keyword{ancestral state reconstruction}
