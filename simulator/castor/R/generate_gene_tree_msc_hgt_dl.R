# generate a random gene tree, based on an extension of the multi-species coalescent (MSC) model to account for horizontal gene transfer and gene duplicatio/loss
# Main input:
#   Rooted time-calibrated species tree (not necessarily ultrametric, i.e. extinct lineages are allowed)
#   Effective population sizes and generation times for each clade (incoming edge) in the species tree. Note that only the product (generation_time*effective_population_size = expected coalescence time for 2 habloid lineages) matters
#   HGT, D/L rates for each clade (incoming edge) in the species tree
#   Number of alleles (individuals) sampled at each tip
generate_gene_tree_msc_hgt_dl = function(	species_tree,							# rooted timetree, with NStips tips and NSedges edges. May (and typically should) include extinct lineages.
											allele_counts				= 1,		# (integer vector) number of alleles sampled per locus and per species (equivalently, the number of individual organisms sampled per species). Either NULL (1 allele per locus per species) or a single integer (same number of alleles per locus per species) or a vector of length NStips listing the numbers of alleles sampled per locus per species. To only not include alleles from a subset of species (e.g. extinct species) set their allele counts to zero.
											population_sizes 			= 1,		# (numeric vector) effective population size along the edge leading into each clade. Either NULL (all population sizes are 1) or a single integer (same population sizes for all edges) or a vector of length NSclades listing population sizes for each clade's incoming edge (including the root).
											generation_times			= 1,		# (numeric vector) generation time along the edge leading into each clade. Either NULL (all generation times are 1) or a single integer (same generation time for all edges) or a vector of length NSclades listing generation times for each clade's incoming edge (including the root).
											mutation_rates				= 1,		# (numeric vector) mutation rates per generation per site, along the edge leading into each clade. Either NULL (all mutation rates are 1) or a single integer (same mutation rate for all edges) or a vector of length NSclades listing mutation rates for each clade's incoming edge (including the root).
											HGT_rates					= 0,		# (numeric vector) horizontal gene transfer rates per lineage per time, along the edge leading into each clade. Either NULL (all HGT rates are 0) or a single integer (same HGT rate for all edges) or a vector of length NSclades listing HGT rates for each clade's incoming edge (including the root).
											duplication_rates			= 0,		# (numeric vector) gene duplication rates per locus per lineage per time, along the edge leading into each clade. Either NULL (all duplication rates are 0) or a single integer (same duplication rate for all edges) or a vector of length NSclades listing duplication rates for each clade's incoming edge (including the root).
											loss_rates					= 0,		# (numeric vector) gene loss rates per locus per lineage per time, along the edge leading into each clade. Either NULL (all loss rates are 0) or a single integer (same loss rate for all edges) or a vector of length NSclades listing loss rates for each clade's incoming edge (including the root).
											gene_edge_unit				= "time",	# (character) one of "time", "generations" or "mutations_expected" or "mutations_random", specifying how edge lengths in the gene-tree are to be measured.
											Nsites						= 1,		# (integer) number of sites (nucleotides) in the gene. Only relevant when generating edge lengths in terms of random mutation counts, i.e. if gene_edge_unit=="mutations_random".
											bottleneck_at_speciation	= FALSE,	# (boolean) if TRUE, all but one children at each node is assumed to have emerged from a single mutant individual, and thus all gene lineages within these bottlenecked species lineages must coalesce at a younger age than the speciation event. Only the first child at each node is excluded from this assumption, corresponding to the "resident population".
											force_coalescence_at_root	= FALSE,	# (boolean) if TRUE, all remaining orphan gene lineages that haven't coalesced before reaching the root, will be combined at the root (via multiple adjacent bifurcations). If false, coalescence events may extend beyond the root into the stem lineage, as long as it takes until all gene lineages have coalesced.
											ploidy						= 1,		# (integer) genetic ploidy, i.e. number of gene copies per individual. Typically 1 (haploids) or 2 (diploids)
											HGT_source_by_locus			= FALSE,	# (boolean) if TRUE, at any HGT event, every extant locus is chosen as source locus with the same probability (hence the probability of a lineage to be a source is proportional to the number of current loci in it). If FALSE, source lineages are chosen with the same probability (regardless of their number of loci) and the source locus within the source lineage is chosen randomly.
											HGT_only_to_empty_clades	= FALSE,	# (boolean) if TRUE, HGT transfers are only done to clades with no current loci
											no_loss_before_time			= 0,		# (numeric) optional time since the root during which no gene losses shall occur (even if loss_rate>0). This is to reduce the probability of an early extinction of the entire locus tree, i.e. give the locus some "startup time" to spread into various species lineages. The default value should be 0.
											max_runtime					= NULL,		# maximum time (in seconds) to allow for the computation; if the computation roughly exceeds this threshold, it is aborted. Use this as protection against badly parameterized models. If NULL or <=0, this option is ignored.
											include_event_times			= TRUE){	# (boolean) if TRUE, then HGT & DL event times and associated clades will be returned
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
	if(!is.null(HGT_rates)){
		if((length(HGT_rates)!=1) && (length(HGT_rates)!=NSclades)){
			return(list(success=FALSE, error=sprintf("Wrong number of HGT rates specified; expecting either 1 or %d, but got %d sizes",NSclades,length(HGT_rates))));
		}
	}
	if(!is.null(duplication_rates)){
		if((length(duplication_rates)!=1) && (length(duplication_rates)!=NSclades)){
			return(list(success=FALSE, error=sprintf("Wrong number of gene duplication rates specified; expecting either 1 or %d, but got %d sizes",NSclades,length(duplication_rates))));
		}
	}
	if(!is.null(loss_rates)){
		if((length(loss_rates)!=1) && (length(loss_rates)!=NSclades)){
			return(list(success=FALSE, error=sprintf("Wrong number of gene loss rates specified; expecting either 1 or %d, but got %d sizes",NSclades,length(loss_rates))));
		}
	}
	if(!is.null(allele_counts)){
		if(length(allele_counts)==1){
			allele_counts = rep(allele_counts[[1]], NStips);
		}else if(length(allele_counts)!=NStips){
			return(list(success=FALSE, error=sprintf("Wrong number of allele counts specified; expecting either 1 or %d, but got %d sizes",NStips,length(allele_counts))));
		}
	}else{
		allele_counts = rep(1, NStips);
	}
	if(is.null(max_runtime)) max_runtime = 0
	
	# simulate gene tree
	results = generate_gene_tree_in_species_tree_MSC_HGT_DL_CPP(NStips						= NStips,
																NSnodes						= NSnodes,
																NSedges						= nrow(species_tree$edge),
																tree_edge					= as.vector(t(species_tree$edge))-1,	# flatten in row-major format and make indices 0-based
																edge_length					= (if(is.null(species_tree$edge.length)) numeric() else species_tree$edge.length),
																population_sizes			= (if(is.null(population_sizes)) numeric() else population_sizes),
																generation_times			= (if(is.null(generation_times)) numeric() else generation_times),
																mutation_rates				= (if(is.null(mutation_rates)) numeric() else mutation_rates),
																HGT_rates					= (if(is.null(HGT_rates)) numeric() else HGT_rates),
																duplication_rates			= (if(is.null(duplication_rates)) numeric() else duplication_rates),
																loss_rates					= (if(is.null(loss_rates)) numeric() else loss_rates),
																allele_counts				= (if(is.null(allele_counts)) numeric() else allele_counts),
																gene_edge_unit				= gene_edge_unit,
																Nsites						= Nsites,
																bottleneck_at_speciation 	= bottleneck_at_speciation,
																force_coalescence_at_root	= force_coalescence_at_root,
																ploidy						= ploidy,
																HGT_source_by_locus			= HGT_source_by_locus,
																HGT_only_to_empty_clades	= HGT_only_to_empty_clades,
																no_loss_before_time			= no_loss_before_time,
																runtime_out_seconds			= max_runtime,
																include_event_times			= include_event_times);
	if(!results$success) return(list(success=FALSE, error=results$error)); # something went wrong
	
	# extract locus tree
	NLtips	= results$locus_tree$NLtips
	NLnodes = results$locus_tree$NLnodes
	NLedges = results$locus_tree$NLedges
	locus2clade = results$locus_tree$locus2clade+1;
	locus_tree = list(	Nnode 		= NLnodes,
						tip.label 	= unlist(sapply(1:NLtips, FUN=function(locus) sprintf("%s.%d",(if(locus2clade[locus]<=NStips) species_tree$tip.label[locus2clade[locus]] else "locus"),locus))),
						node.label 	= NULL,
						edge 		= matrix(results$locus_tree$locus_tree_edge,ncol=2,byrow=TRUE) + 1,
						edge.length = results$locus_tree$locus_edge_length,
						root 		= results$locus_tree$locus_root+1)
	class(locus_tree) = "phylo";
	attr(locus_tree,"order") = NULL

	# extract gene (coalescent) tree, simulated according to the MSC along the locus tree
	NGtips	= results$gene_tree$NGtips
	NGnodes = results$gene_tree$NGnodes
	NGedges = results$gene_tree$NGedges
	gene_tip2species_tip = results$gene_tree$gene_tip2species_tip+1;
	gene_tree = list(	Nnode 		= NGnodes,
						tip.label 	= unlist(sapply(1:NGtips, FUN=function(gtip) sprintf("%s.%d",species_tree$tip.label[gene_tip2species_tip[gtip]],gtip))),
						node.label 	= NULL,
						edge 		= matrix(results$gene_tree$gene_tree_edge,ncol=2,byrow=TRUE) + 1,
						edge.length = results$gene_tree$gene_edge_length,
						root 		= results$gene_tree$gene_root+1)
	class(gene_tree) = "phylo";
	attr(gene_tree,"order") = NULL
		
	return(list(success					= TRUE,
				locus_tree				= locus_tree,
				locus_type				= strsplit(intToUtf8(results$locus_tree$locus_type),"")[[1]],
				locus2clade				= locus2clade,
				HGT_times				= (if(include_event_times) results$locus_tree$HGT_times else NULL),
				HGT_source_clades		= (if(include_event_times) results$locus_tree$HGT_source_clades+1 else NULL),
				HGT_target_clades		= (if(include_event_times) results$locus_tree$HGT_target_clades+1 else NULL),
				duplication_times		= (if(include_event_times) results$locus_tree$duplication_times else NULL),
				duplication_clades		= (if(include_event_times) results$locus_tree$duplication_clades+1 else NULL),
				loss_times				= (if(include_event_times) results$locus_tree$loss_times else NULL),
				loss_clades				= (if(include_event_times) results$locus_tree$loss_clades+1 else NULL),
				gene_tree				= gene_tree,
				gene_tip2species_tip	= gene_tip2species_tip,
				gene_tip2locus_tip		= results$gene_tree$gene_tip2locus_tip+1,
				gene_node2locus_edge	= results$gene_tree$gene_node2locus_edge+1,
				gene_clade_times		= results$geneTree$gene_clade2time));
}