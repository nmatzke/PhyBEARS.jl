\name{generate_gene_tree_msc_hgt_dl}
\alias{generate_gene_tree_msc_hgt_dl}
\title{
Generate gene trees based on the multi-species coalescent, horizontal gene transfers and duplications/losses.
}
\description{
Generate a random gene tree within a given species timetree, based on an extension of the multi-species coalescent (MSC) model that includes horizontal gene transfers (HGT, incorporation of non-homologous genes as new loci), gene duplication and gene loss. The simulation consists of two phases. In the first phase a random "locus tree" is generated in forward time, according to random HGT, duplication and loss events. In the 2nd phase, alleles picked randomly from each locus are coalesced in backward time according to the multispecies coalescent, an extension of the Wright-Fisher model to multiple species. This function does not account for hybridization.
}
\usage{
generate_gene_tree_msc_hgt_dl( species_tree,
                        allele_counts               = 1, 
                        population_sizes            = 1,
                        generation_times            = 1,
                        mutation_rates              = 1,
                        HGT_rates                   = 0,
                        duplication_rates           = 0,
                        loss_rates                  = 0,
                        gene_edge_unit              = "time",
                        Nsites                      = 1,
                        bottleneck_at_speciation    = FALSE,
                        force_coalescence_at_root   = FALSE,
                        ploidy                      = 1,
                        HGT_source_by_locus         = FALSE,
                        HGT_only_to_empty_clades    = FALSE,
                        no_loss_before_time         = 0,
                        max_runtime                 = NULL,
                        include_event_times         = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{species_tree}{
Rooted timetree of class "phylo". The tree can include multifurcations and monofurcations. The tree need not necessarily be ultrametric, i.e. it may include extinct species. Edge lengths are assumed to be in time units.
}
\item{allele_counts}{
Integer vector, listing the number of alleles sampled per species and per locus. This can be interpreted as the number if individual organisms surveyed from each species, assuming that all loci are included once from each individual. The number of tips in the generated gene tree will be equal to the sum of allele counts across all species.
\code{allele_counts} can either be \code{NULL} (1 allele per species), or a single integer (same number of alleles per species), or a vector of length Ntips listing the numbers of alleles sampled per species. In the latter case, the total number of tips in the returned gene tree will be equal to the sum of entries in \code{allele_counts}. Some entries in \code{allele_counts} may be zero (no alleles sampled from those species).
}
\item{population_sizes}{
Integer vector, listing the effective population size on the edge leading into each tip/node in the species tree. Either \code{NULL} (all population sizes are 1), or a single integer (same population sizes for all edges), or a vector of length Ntips+Nnodes, listing population sizes for each clade's incoming edge (including the root). The population size for the root's incoming edge corresponds to the population size at the tree's stem (only relevant if \code{force_coalescence_at_root=FALSE}).
}
\item{generation_times}{
Numeric vector, listing the generation time along the edge leading into each clade. Either \code{NULL} (all generation times are 1), or a single integer (same generation time for all edges) or a vector of length Ntips+Nnodes, listing generation times for each clade's incoming edge (including the root). The generation time for the root's incoming edge corresponds to the generation time at the tree's stem (only relevant if \code{force_coalescence_at_root=FALSE}).
}
\item{mutation_rates}{
Numeric vector, listing the probability of mutation per site and per generation along the edge leading into each clade. Either \code{NULL} (all mutation rates are 1), or a single integer (same mutation rate for all edges) or a vector of length Ntips+Nnodes, listing mutation rates for each clade's incoming edge (including the root). The mutation rate for the root's incoming edge corresponds to the mutation rate at the tree's stem (only relevant if \code{force_coalescence_at_root=FALSE}). The value of \code{mutation_rates} is only relevant if \code{gene_edge_unit} is "mutations_expected" or "mutations_random". Mutation rates represent probabilities, and so they must be between 0 and 1.
}
\item{HGT_rates}{
Numeric vector, listing horizontal gene transfer rates per lineage per time, along the edge leading into each clade. Either \code{NULL} (all HGT rates are 0) or a single integer (same HGT rate for all edges) or a vector of length Ntips+Nnodes, listing HGT rates for each clade's incoming edge (including the root).
}
\item{duplication_rates}{
Numeric vector, listing gene duplication rates per locus per lineage per time, along the edge leading into each clade. Either \code{NULL} (all duplication rates are 0) or a single integer (same duplication rate for all edges) or a vector of length Ntips+Nnodes listing duplication rates for each clade's incoming edge (including the root).
}
\item{loss_rates}{
Numeric vector, listing gene loss rates per locus per lineage per time, along the edge leading into each clade. Either \code{NULL} (all loss rates are 0) or a single integer (same loss rate for all edges) or a vector of length Ntips+Nnodes listing loss rates for each clade's incoming edge (including the root).
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
\item{HGT_source_by_locus}{
Logical. If \code{TRUE}, then at any HGT event, every extant locus is chosen as source locus with the same probability (hence the probability of a lineage to be a source is proportional to the number of current loci in it). If \code{FALSE}, source lineages are chosen with the same probability (regardless of the number of current loci in them) and the source locus within the source lineage is chosen randomly.
}
\item{HGT_only_to_empty_clades}{
Logical, specifying whether HGT transfers only occur into clades with no current loci.
}
\item{no_loss_before_time}{
Numeric, optional time since the root during which no gene losses shall occur (even if \code{loss_rate>0}). This option can be used to reduce the probability of an early extinction of the entire gene tree, by giving the gene tree some "startup time" to spread into various species lineages. If zero, gene losses are possible right from the start of the simulation.
}
\item{max_runtime}{
Numeric, optional maximum computation time (in seconds) to allow for the simulation. Use this to avoid occasional explosions of runtimes, for example due to very large generated trees. Aborted simulations will return with the flag \code{success=FALSE} (i.e., no tree is returned at all).
}
\item{include_event_times}{
Logical, specifying whether the times of HGT, duplication and loss events should be returned as well. If these are not needed, then set \code{include_event_times=FALSE} for efficiency.
}
}


\details{
This function assumes that the species tree is a time tree, i.e. with edge lengths given in actual time units. If \code{species_tree} is ultrametric and \code{gene_edge_unit=="time"}, then the gene tree (but not necessarily the locus tree) will be ultrametric as well. The root of the locus and gene tree coincides with the root of the species tree.

The meaning of \code{gene_edge_unit} is the same as for the function \code{\link{generate_gene_tree_msc}}.
}


\value{
A named list with the following elements:
\item{success}{
Logical, indicating whether the gene tree was successfully generated. If \code{FALSE}, the only other value returned is \code{error}.
}
\item{locus_tree}{
The generated locus timetree, of class "phylo". The locus tree describes the genealogy of loci due to HGT, duplication and loss events. Each tip and node of the locus tree is embedded within a specific species edge. For example, tips of the locus tree either coincide with tips of the species tree (if the locus persisted until the species went extinct or until the present) or they correspond to gene loss events. In the absence of any HGT, duplication and loss events, the locus tree will resemble the species tree.
}
\item{locus_type}{
Character vector of length NLtips + NLnodes (where NLtips and NLnodes are the number of tips and nodes in the locus tree, respectively), specifying the type/origin of each tip and node in the locus tree. For nodes, type 'h' corresponds to an HGT event, type 'd' to a duplication event, and type 's' to a speciation event. For tips, type 'l' represents a loss event, and type 't' a terminal locus (i.e., coinciding with a species tip). For example, if the input species tree was an ultrametric tree representing only extant species, then the locus tree tips of type 't' are the loci that could potentially be sampled from those extant species.
}
\item{locus2clade}{
Integer vector of length NLtips + NLnodes, with values in NStips+NSnodes, specifying for every locus tip or node the correspondng "host" species tip or node.
}
\item{HGT_times}{
Numeric vector, listing HGT event times (counted since the root) in ascending order. Only included if \code{include_event_times==TRUE}.
}
\item{HGT_source_clades}{
Integer vector of the same length as \code{HGT_times} and with values in 1,..,Ntips+Nnodes, listing the "source" species tip/node of each HGT event (in order of occurrence). The source tip/node is the tip/node from whose incoming edge a locus originated at the time of the transfer. Only included if \code{include_event_times==TRUE}.
}
\item{HGT_target_clades}{
Integer vector of the same length as \code{HGT_times} and with values in 1,..,Ntips+Nnodes, listing the "target" species tip/node of each HGT event (in order of occurrence). The target (aka. recipient) tip/node is the tip/node within whose incoming edge a locus was created by the transfer. Only included if \code{include_event_times==TRUE}.
}
\item{duplication_times}{
Numeric vector, listing gene duplication event times (counted since the root) in ascending order. Only included if \code{include_event_times==TRUE}.
}
\item{duplication_clades}{
Integer vector of the same length as \code{duplication_times} and with values in 1,..,Ntips+Nnodes, listing the species tip/node in whose incoming edge each duplication event occurred (in order of occurrence). Only included if \code{include_event_times==TRUE}.
}
\item{loss_times}{
Numeric vector, listing gene loss event times (counted since the root) in ascending order. Only included if \code{include_event_times==TRUE}.
}
\item{loss_clades}{
Integer vector of the same length as \code{loss_times} and with values in 1,..,Ntips+Nnodes, listing the species tip/node in whose incoming edge each loss event occurred (in order of occurrence). Only included if \code{include_event_times==TRUE}.
}
\item{gene_tree}{
The generated gene tree, of type "phylo".
}
\item{gene_tip2species_tip}{
Integer vector of length NGtips (where NGtips is the number of tips in the gene tree) with values in 1,..,Ntips+Nnodes, mapping gene-tree tips to species-tree tips.
}
\item{gene_tip2locus_tip}{
Integer vector of length NGtips with values in 1,..,NLtips, mapping gene-tree tips to locus-tree tips.
}
\item{gene_node2locus_edge}{
Integer vector of length NGnodes with values in 1,..,NLedges, mapping gene-tree nodes to locus-tree edges.
}
\item{gene_clade_times}{
Numeric vector of size NGtips+NGnodes, listing the time (temporal distance from species root) of each tip and node in the gene tree. The units will be the same as the time units of the species tree. Note that this may include negative values, if some gene lineages coalesce at a greater age than the root.
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
\code{\link{generate_gene_tree_msc}}
}

\examples{
# Simulate a simple species tree
parameters   = list(birth_rate_factor=1)
Nspecies     = 10
species_tree = generate_random_tree(parameters,max_tips=Nspecies)$tree

# Simulate a haploid gene tree within the species tree, including HGTs and gene loss
# Assume the same population size and generation time everywhere
# Assume the number of alleles samples per species is poisson-distributed
results = generate_gene_tree_msc_hgt_dl(species_tree, 
                                        allele_counts      = rpois(Nspecies,3),
                                        population_sizes   = 1000,
                                        generation_times   = 1,
                                        ploidy             = 1,
                                        HGT_rates          = 0.1,
                                        loss_rates         = 0.05);
if(!results$success){
    # simulation failed
    cat(sprintf("  ERROR: \%s\n",results$error))
}else{
    # simulation succeeded
    gene_tree = results$gene_tree
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
