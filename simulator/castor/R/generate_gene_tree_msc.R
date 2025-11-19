# generate a random gene tree, based on the multi-species coalescent (MSC) model
# in the MSC, every branch of the species tree has an effective population size (Ne) and generation time (T),
#   and gene alleles coalesce backward in time according to the Wright-Fisher model 
#   and depending on the branch's population size and generation time.
# This model does not account for gene duplication/loss, nor for hybridization or horizontal gene transfer
# It is only meant to model "incomplete lineage sorting", otherwise known as "deep coalescence",
#   which is one of the many mechanisms that can cause discordance between gene trees and species trees.
# Main input:
#   Rooted time-calibrated species tree (not necessarily ultrametric, i.e. extinct lineages are allowed)
#   Effective population sizes and generation times for each edge in the tree
#   Number of alleles sampled at each tip
generate_gene_tree_msc = function(	species_tree,							# rooted timetree, with NStips tips and NSedges edges
									allele_counts				= 1,		# (integer vector) number of alleles sampled per species. Either NULL (1 allele per species) or a single integer (same number of alleles per species) or a vector of length NStips listing the numbers of alleles sampled per species. The total number of tips in the returned gene tree will be NGtips=sum(allele_counts).
									population_sizes 			= 1,		# (numeric vector) effective population size along the edge leading into each clade. Either NULL (all population sizes are 1) or a single integer (same population sizes for all edges) or a vector of length NSclades listing population sizes for each clade's incoming edge (including the root).
									generation_times			= 1,		# (numeric vector) generation time along the edge leading into each clade. Either NULL (all generation times are 1) or a single integer (same generation time for all edges) or a vector of length NSclades listing generation times for each clade's incoming edge (including the root).
									mutation_rates				= 1,		# (numeric vector) mutation probabilities per generation per site, along the edge leading into each clade. Either NULL (all mutation rates are 1) or a single integer (same mutation rate for all edges) or a vector of length NSclades listing mutation rates for each clade's incoming edge (including the root). Mutation rates should be within 0 and 1.
									gene_edge_unit				= "time",	# (character) one of "time", "generations" or "mutations_expected" or "mutations_random", specifying how edge lengths in the gene-tree are to be measured.
									Nsites						= 1,		# (integer) number of sites (nucleotides) in the gene. Only relevant when generating edge lengths in terms of random mutation counts, i.e. if gene_edge_unit=="mutations_random".
									bottleneck_at_speciation	= FALSE,	# (boolean) if TRUE, all but one children at each node is assumed to have emerged from a single mutant individual, and thus all gene lineages within these bottlenecked species lineages must coalesce at a younger age than the speciation event. Only the first child at each node is excluded from this assumption, corresponding to the "resident population".
									force_coalescence_at_root	= FALSE,	# (boolean) if TRUE, all remaining orphan gene lineages that haven't coalesced before reaching the root, will be combined at the root (via multiple adjacent bifurcations). If false, coalescence events may extend beyond the root into the stem lineage, as long as it takes until all gene lineages have coalesced.
									ploidy						= 1,		# (integer) genetic ploidy, i.e. number of gene copies per individual. Typically 1 (haploids) or 2 (diploids)
									gene_tip_labels				= NULL){	# (character vector) Tip labels for the gene tree. Either NULL (gene tips will be <species_tip_label>.<allele index>) or a character vector of length NGtips (with genes listed in the order of the corresponding species tips)
	NStips   = length(species_tree$tip.label)
	NSnodes  = species_tree$Nnode
	NSclades = NStips + NSnodes
	
	# basic error checking
	if(NStips<=1) return(list(success=FALSE, error="Species tree is too small"));
	if(ploidy<=0) return(list(success=FALSE, error="Ploidy must be a strictly positive integer"));
	if(!is.null(population_sizes)){
		if((length(population_sizes)!=1) && (length(population_sizes)!=NSclades)){
			return(list(success=FALSE, error=sprintf("Wrong number of population sizes specified; expecting either 1 or %d, but got %d sizes",NSclades,length(population_sizes))));
		}
	}
	if(!is.null(generation_times)){
		if((length(generation_times)!=1) && (length(generation_times)!=NSclades)){
			return(list(success=FALSE, error=sprintf("Wrong number of generation times specified; expecting either 1 or %d, but got %d sizes",NSclades,length(generation_times))));
		}
	}
	if(!is.null(mutation_rates)){
		if((length(mutation_rates)!=1) && (length(mutation_rates)!=NSclades)){
			return(list(success=FALSE, error=sprintf("Wrong number of mutation rates specified; expecting either 1 or %d, but got %d sizes",NSclades,length(mutation_rates))));
		}
	}
	if(!is.null(allele_counts)){
		if(length(allele_counts)==1){
			allele_counts = rep(allele_counts[[1]], NStips);
		}else if(length(allele_counts)!=NStips){
			return(list(success=FALSE, error=sprintf("Wrong number of allele counts specified; expecting either 1 or %d, but got %d counts",NStips,length(allele_counts))));
		}
	}else{
		allele_counts = rep(1, NStips);
	}
	
	# simulate gene tree
	results = generate_gene_tree_in_species_tree_MSC_CPP(	NStips						= NStips,
															NSnodes						= NSnodes,
															NSedges						= nrow(species_tree$edge),
															tree_edge					= as.vector(t(species_tree$edge))-1,	# flatten in row-major format and make indices 0-based
															edge_length					= (if(is.null(species_tree$edge.length)) numeric() else species_tree$edge.length),
															population_sizes			= (if(is.null(population_sizes)) numeric() else population_sizes),
															generation_times			= (if(is.null(generation_times)) numeric() else generation_times),
															mutation_rates				= (if(is.null(mutation_rates)) numeric() else mutation_rates),
															allele_counts				= (if(is.null(allele_counts)) numeric() else allele_counts),
															gene_edge_unit				= gene_edge_unit,
															Nsites						= Nsites,
															bottleneck_at_speciation 	= bottleneck_at_speciation,
															force_coalescence_at_root	= force_coalescence_at_root,
															ploidy						= ploidy);
	if(!results$success) return(list(success=FALSE, error=results$error)); # something went wrong
	NGtips	= results$NGtips
	NGnodes = results$NGnodes
	NGedges = results$NGedges
	if((!is.null(gene_tip_labels)) && (NGtips!=length(gene_tip_labels))){
		return(list(success=FALSE, error=sprintf("Wrong number of gene tip labels specified; expecting either %d, but got %d labels",NGtips,length(gene_tip_labels))));
	}else{
		gene_tip_labels = unlist(sapply(1:NStips, FUN=function(tip) (if(allele_counts[tip]==0) NULL else (if(allele_counts[tip]==1) species_tree$tip.label[tip] else sprintf("%s.%d",species_tree$tip.label[tip],c(1:allele_counts[tip]))))))
	}
	gene_tree = list(	Nnode 		= NGnodes,
						tip.label 	= gene_tip_labels,
						node.label 	= NULL,
						edge 		= matrix(results$gene_tree_edge,ncol=2,byrow=TRUE) + 1,
						edge.length = results$gene_edge_length,
						root 		= results$gene_root+1)
	class(gene_tree) = "phylo";
	attr(gene_tree,"order") = NULL
		
	return(list(success					= TRUE,
				tree					= gene_tree,
				gene_tip2species_tip	= results$gene_tip2species_tip+1,
				gene_node2species_edge	= results$gene_node2species_edge+1,
				gene_clade_times		= results$gene_clade2time));	
}