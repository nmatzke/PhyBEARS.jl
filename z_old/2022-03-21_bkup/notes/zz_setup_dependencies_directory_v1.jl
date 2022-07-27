

#######################################################
# Set up a dependencies directory
# (this is a Julia environment, where you hold your
#  dependencies, i.e. the packages you need)
#######################################################

# 
# https://medium.com/@Jernfrost/my-new-workflow-with-julia-1-0-99711103d97c

# Terminal
cd /GitHub
mkdir /GitHub/BioGeoDependencies


# Julia
cd("/GitHub/BioGeoDependencies")
Pkg.activate(".")
# Activating environment at `/GitHub/BioGeoDependencies/Project.toml`


# Terminal
cat Project.toml
cat Manifest.toml