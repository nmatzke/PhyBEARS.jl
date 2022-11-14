library(ape)
library(BioGeoBEARS)

trfn = "/GitHub/PhyBEARS.jl/notes/amph_shl_dates_frogs.tre"
geogfn = "/GitHub/PhyBEARS.jl/notes/frog_ecoregions_reduced.data"
tr = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
tips_in_geog = row.names(tipranges@df)
tips_in_geog

tips_in_tr = tr$tip.label

tips_to_dropTF = (tips_in_tr %in% tips_in_geog) == FALSE
sum(tips_to_dropTF)

tr2 = ladderize(tr) 
tr2$node.label = NULL
write.tree(tr2, file="/GitHub/PhyBEARS.jl/notes/amph_shl_dates_frogs2.tre")
tr2$node.label
