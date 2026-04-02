

library(ape)
library(BioGeoBEARS)

wd = "~/GitHub/PhyBEARS.jl/sims/mvd_ant/"
setwd(wd)

trfn = "output_mcctree_rev_Aug28.newick"
tr = read.tree(trfn)

tipnames = tr$tip.label
length(tipnames)

geogfn = "geog_5areas.data"
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
dim(tipranges@df)

tipnames_140 = row.names(tipranges@df)
tipnames_140

TF = tipnames_140 %in% tipnames
sum(TF)

tr2 = keep.tip(phy=tr, tip=tipnames_140[TF])
length(tr2$tip.label)

plot(tr2)


TF = tipnames_140 %in% tr2$tip.label
sum(TF)
tipranges_131 = tipranges
tipranges_131@df = tipranges@df[TF,]
dim(tipranges_131@df)

write.tree(tr2, file="tree_131sp.newick")
write.tree(tr2, file="tree.newick")

save_tipranges_to_LagrangePHYLIP(tipranges_131, lgdata_fn="geog_131sp.data")
save_tipranges_to_LagrangePHYLIP(tipranges_131, lgdata_fn="geog.data")

