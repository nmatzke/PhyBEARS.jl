\name{generate_gene_tree_msc}
\alias{generate_gene_tree_msc}
\title{
Generate a gene tree based on the multi-species coalescent model.
}
\description{
Generate a random gene tree within a given species timetree, based on the multi-species coalescent (MSC) model. In this implementation of the MSC, every branch of the species tree has a specific effective population size (Ne) and a specific generation time (T), and gene alleles coalesce backward in time according to the Wright-Fisher model. This model does not account for gene duplication/loss, nor for hybridization or horizontal gene transfer. It is only meant to model "incomplete lineage sorting", otherwise known as "deep coalescence", which is one of the many mechanisms that can cause discordance between gene trees and species trees.
}
\usage{
generate_gene_tree_msc( species_tree,
                        allele_counts              = 1, 
                        population_sizes           = 1,
                        generation_times           = 1,
                        mutation_rates             = 1,
                        gene_edge_unit             = "time",
                        Nsites                     = 1,
                        bottleneck_at_speciation   = FALSE,
                        force_coalescence_at_root  = FALSE,
                        ploidy                     = 1,
                        gene_tip_labels            = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{species_tree}{
Rooted timetree of class "phylo". The tree can include multifurcations and monofurcations. The tree need not necessarily be ultrametric, i.e. it may include extinct species. Edge lengths are assumed to be in time units.
}
\item{allele_counts}{
Integer vector, listing the number of alleles sampled per species. Either \code{NULL} (1 allele per species), or a single integer (same number of alleles per species), or a vector of length Ntips listing the numbers of alleles sampled per species. In the latter case, the total number of tips in the returned gene tree will be equal to the sum of entries in \code{allele_counts}. Some entries in \code{allele_counts} may be zero (no alleles sampled from those species).
}
\item{population_sizes}{
Integer vector, listing the effective population size on the edge leading into each tip/node in the species tree. Either \code{NULL} (all population sizes are 1), or a single integer (same population sizes for all edges), or a vector of length Ntips+Nnodes, listing population sizes for each clade's incoming edge (including the root). The population size for the root's incoming edge corresponds to the population size at the tree's stem (only relevant if \code{force_coalescence_at_root=FALSE}).
}
\item{generation_times}{
Numeric vector, listing the generation time along the edge leading into each clade. Either \code{NULL} (all generation times are 1), or a single integer (same generation time for all edges) or a vector of length Ntips+Nnodes, listing generation times for each clade's incoming edge (including the root). The generation time for the root's incoming edge corresponds to the generation time at the tree's stem (only relevant if \code{force_coalescence_at_root=FALSE}).
}
\item{mutation_rates}{
Numeric vector, listing the mutation rate (per site and per generation) along the edge leading into each clade. Either \code{NULL} (all mutation rates are 1), or a single integer (same mutation rate for all edges) or a vector of length Ntips+Nnodes, listing mutation rates for each clade's incoming edge (including the root). The mutation rate for the root's incoming edge corresponds to the mutation rate at the tree's stem (only relevant if \code{force_coalescence_at_root=FALSE}). The value of \code{mutation_rates} is only relevant if \code{gene_edge_unit} is "mutations_expected" or "mutations_random". Mutation rates represent probabilities, and so they must be between 0 and 1.
}
\item{gene_edge_unit}{
Character, either "time", "generations", "mutations_expected" (expected mean number of mutations per site), or "mutations_random" (randomly generated mean number of mutations per site), specifying how edge lengths in the gene tree should be measured. By default, gene-tree edges are measured in time, as is the case for the input species tree.
}
\item{Nsites}{
Integer, specifying the number of sites (nucleotides) in the gene. Only relevant when generating edge lengths in terms of random mutation counts, i.e. if \code{gene_edge_unit=="mutations_random"}.
}
\item{bottleneck_at_speciation}{
Logical. If \code{TRUE}, then all but one children at each node are assumed to have emerged from a single mutant individual, and thus all gene lineages within these bottlenecked species lineages must coalesce at a younger or equal age as the speciation event. Only the first child at each node is excluded from this assumption, corresponding to the "resident population" during the speciation event. This option deviates from the classical MSC model, and is experimental.
}
\item{force_coalescence_at_root}{
Logical. If \code{TRUE}, all remaining orphan gene lineages that haven't coalesced before reaching the species-tree's root, will be combined at the root (via multiple adjacent bifurcations). If \code{FALSE}, coalescence events may extend beyond the species-tree's root into the stem lineage, as long as it takes until all gene lineages have coalesced.
}
\item{ploidy}{
Integer, specifying the assumed genetic ploidy, i.e. number of gene copies per individual. Typically 1 for haploids, or 2 for diploids.
}
\item{gene_tip_labels}{
Character vector specifying tip labels for the gene tree (i.e., for each of the sampled alleles) in the order of the corresponding species tips. Can also be \code{NULL}, in which case gene tips will be set to <species_tip_label>.<allele index>.
}
}


\details{
This function assumes that Kingman's coalescent assumption is met, i.e. that the effective population size is much larger than the number of allele lineages coalescing within any given branch.

The function assumes that the species tree is a time tree, i.e. with edge lengths given in actual time units. To simulate gene trees in coalescence time units, choose \code{population_sizes} and \code{generation_times} accordingly (this only makes sense if the product of \code{population_sizes} \eqn{\times} \code{generation_times} is the same everywhere). If \code{species_tree} is ultrametric and \code{gene_edge_unit=="time"}, then the gene tree will be ultrametric as well.

If \code{gene_edge_unit} is "mutations_random", then the number of generations elapsed along each time segment is translated into a randomly distributed number of accumulated mutations, according to a binomial distribution where the probability of success is equal to the mutation rate and the number of trials is equal to the number of generations multiplied by \code{Nsites}; this number of mutations is averaged across all sites, i.e. the edge lengths in the returned gene tree always refer to the mean number of mutations per site. In cases where the mutation rate varies across the species tree and a single gene edge spans multiple species edges, the gene edge length will be a sum of multiple binomially distributed mutation counts (again, divided by the number of sites), corresponding to the times spent in each species edge.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the gene tree was successfully generated. If \code{FALSE}, the only other value returned is \code{error}.
}
\item{tree}{
The generated gene tree, of class "phylo". This tree will be rooted and bifurcating. It is only guaranteed to be ultrametric if \code{species_tree} was ultrametric.
}
\item{gene_tip2species_tip}{
Integer vector of length NGtips (where NGtips is the number of tips in the gene tree), mapping gene-tree tips to species-tree tips.
}
\item{gene_node2species_edge}{
Integer vector of length NGnodes (where NGnodes is the number of internal nodes in the gene tree), mapping gene-tree nodes (=coalescence events) to the species-tree edges where the coalescences took place.
}
\item{gene_clade_times}{
Numeric vector of size NGtips+NGnodes, listing the time (total temporal distance from species root) of each tip and node in the gene tree. The units will be the same as the time units assumed for the species tree. Note that this may include negative values, if some gene lineages coalesce at a greater age than the root.
}
\item{error}{
Character, containing an explanation of the error that occurred. Only included if \code{success==FALSE}.
}
}

\author{Stilianos Louca}

\references{
J. H. Degnan, N. A. Rosenberg (2009). Gene tree discordance, phylogenetic inference and the multispecies coalescent. Trends in Ecology & Evolution. 24:332-340.

B. Rannala, Z. Yang (2003). Bayes estimation of species divergence times and ancestral population sizes using DNA sequences from multiple loci. Genetics. 164:1645-1656.
}

\seealso{
\code{\link{generate_random_tree}},
\code{\link{generate_gene_tree_msc_hgt_dl}}
}

\examples{
# Simulate a simple species tree
parameters   = list(birth_rate_factor=1)
Nspecies     = 10
species_tree = generate_random_tree(parameters,max_tips=Nspecies)$tree

# Simulate a haploid gene tree within the species tree
# Assume the same population size and generation time everywhere
# Assume the number of alleles samples per species is poisson-distributed
results = generate_gene_tree_msc(species_tree, 
                                 allele_counts      = rpois(Nspecies,3),
                                 population_sizes   = 1000,
                                 generation_times   = 1,
                                 ploidy             = 1);
if(!results$success){
    # simulation failed
    cat(sprintf("  ERROR: \%s\n",results$error))
}else{
    # simulation succeeded
    gene_tree = results$tree
    cat(sprintf("  Gene tree has \%d tips\n",length(gene_tree$tip.label)))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Only 1 keyword per line
\keyword{random}
\keyword{tree model}
\keyword{simulation}
\keyword{gene tree}
