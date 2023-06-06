

library(BioGeoBEARS)
library(ape)

wd = "/GitHub/PhyBEARS.jl/simulator/examples/ex1/"
setwd(wd)


# Load simulated tree
trfn="living_tree_noNodeLabels.newick"
tr=read.tree(trfn)

# Load simulated states
statesfn="simstates_living.txt"
states=read.table(statesfn)

# Labels and Colors (4 states)
tmplabels = c("_", "A", "B", "AB")
tmpcols= c("black","blue","green","white")

# Internal node numbers
intnums = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode)

# Tip node numbers
tipnums = 1:length(tr$tip.label)


# Make the plot
plot(tr, label.offset=1)

# Color the background (bg) of the boxes by state
nodelabels(text=tmplabels[c(states[intnums,])], node=intnums, bg=tmpcols[c(states[intnums,])])

# Add tip.labels
tiplabels(text=tmplabels[c(states[tipnums,])], tip=tipnums, bg=tmpcols[c(states[tipnums,])])