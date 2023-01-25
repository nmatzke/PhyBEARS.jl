

#######################################################
# SETUP INSTRUCTIONS FOR BIOGEOBEARS
#######################################################

#######################################################
# INSTALL R
#######################################################
# Install R from:
# https://cran.r-project.org/

# Unless you really love purely-command-line-R, install either:
# RStudio for Windows / Mac / Linux: https://posit.co/downloads/ -- RStudio Desktop, Free version
#          ...OR...
# R.app for Mac (if you like), from: https://cran.r-project.org/ --> Download R for macOS

# Alternatively, you can run R online with a free RStudio Cloud account:
# https://posit.cloud/
#######################################################


#######################################################
# INSTALLING BIOGEOBEARS
#######################################################
# Run the install.packages commands ONCE

# Please run them one at a time.

# Installing devtools
install.packages("devtools", type="binary", repos="https://cloud.r-project.org")

# 
# IF YOU GET A MESSAGE LIKE THIS, TYPE "n" FOR NO:
#  There are binary versions available but the source versions are later:
#           binary source needs_compilation
# jsonlite   1.8.3  1.8.4              TRUE
# htmltools  0.5.3  0.5.4              TRUE
# 
# Do you want to install from sources the packages which need compilation? (Yes/no/cancel)
# 

install.packages("ape", type="binary", repos="https://cloud.r-project.org")
install.packages("Rcpp", type="binary", repos="https://cloud.r-project.org")
install.packages("ape", type="binary", repos="https://cloud.r-project.org")
install.packages("FD", type="binary", repos="https://cloud.r-project.org")
install.packages("snow", type="binary", repos="https://cloud.r-project.org")
install.packages("phytools", type="binary", repos="https://cloud.r-project.org")
install.packages("phangorn", type="binary", repos="https://cloud.r-project.org")
install.packages("phylobase", type="binary", repos="https://cloud.r-project.org")
install.packages("optimx", type="binary", repos="https://cloud.r-project.org")
install.packages("GenSA", type="binary", repos="https://cloud.r-project.org")

# R packages by Nick Matzke -- dependencies of BioGeoBEARS
install.packages("rexpokit", type="binary", repos="https://cloud.r-project.org")
install.packages("cladoRcpp", type="binary", repos="https://cloud.r-project.org")


# Install BioGeoBEARS from GitHub
# (BioGeoBEARS is pure R, so installation is easy *if* the above 
#  packages have been installed)
library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS", INSTALL_opts="--byte-compile", upgrade="never")
#######################################################


# Check that your BioGeoBEARS installation loads
library(ape)
library(phytools)
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
                  # seems to sometimes fail on simple problems (2-3 parameters)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions 
library(parallel)	# prefer library(parallel), try either)


#######################################################
# 2018-10-10 update: I have been putting the 
# updates on CRAN/GitHub
# You should use:
# rexpokit version 0.26.6 from CRAN
# cladoRcpp version 0.15 from CRAN
# BioGeoBEARS version 1.1 from GitHub, install with:
# library(devtools)
# devtools::install_github(repo="nmatzke/BioGeoBEARS")
#######################################################
library(rexpokit)
library(cladoRcpp)

# IGNORE WARNINGS - these are just about some example files & documentation
library(BioGeoBEARS)
# WHEN YOU SEE: "DONE (BioGeoBEARS)", you are done installing


# Try some basic commands: rexpokit
Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504,
0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)

# Make a series of t values
tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)

# Exponentiate each with EXPOKIT's dgpadm (good for small dense matrices)
for (t in tvals)
	{
	Pmat = expokit_dgpadm_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
	cat("\n\nTime=", t, "\n", sep="")
	print(Pmat)
	}


# Try some basic commands: cladoRcpp
numstates_from_numareas(numareas=4, maxareas=4, include_null_range=TRUE)
numstates_from_numareas(numareas=4, maxareas=3, include_null_range=TRUE)
numstates_from_numareas(numareas=4, maxareas=2, include_null_range=TRUE)
numstates_from_numareas(numareas=4, maxareas=1, include_null_range=TRUE)

numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=3, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=2, include_null_range=FALSE)
numstates_from_numareas(numareas=4, maxareas=1, include_null_range=FALSE)




# Try some basic commands: Tree reading with APE
newick_string = "((((human:2.5,Lucy:0.5):3.5,chimp:6.0):1.0,gorilla:7.0):5.0,orang:12.0);"
tr = read.tree(file="", text=newick_string)
tr
plot(tr)
title("Example phylogeny: great apes")
axisPhylo() # plots timescale
mtext(text="Mega-annum (Ma)", side=1, line=2)


# Try some basic commands: BioGeoBEARS

# The function "prt": prints the tree to a table
# This can be very handy for understanding R node numbers,
# the APE phylo object tree structure, etc.
trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE, fossils_older_than=0.001)
trtable

# Compare to:
names(tr)
tr$Nnode
tr$tip.label
tr$edge.length
tr$edge


# If you like, you can continue on to the basic tutorial script:
# http://phylo.wikidot.com/biogeobears#script





