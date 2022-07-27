# PhyBEARS.jl
BioGeoBEARS but faster and new stuff with Julia

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaLang.github.io/PhyBEARS.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaLang.github.io/PhyBEARS.jl/dev)

Linux and macOS: [![Build Status](https://travis-ci.org/JuliaLang/PhyBEARS.jl.svg?branch=master)](https://travis-ci.org/JuliaLang/PhyBEARS.jl)

Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/JuliaLang/PhyBEARS.jl?branch=master&svg=true)](https://ci.appveyor.com/project/tkelman/example-jl/branch/master)

[![Coverage Status](https://coveralls.io/repos/JuliaLang/PhyBEARS.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaLang/PhyBEARS.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaLang/PhyBEARS.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaLang/PhyBEARS.jl?branch=master)



# Local add instructions:
Pkg.add(PackageSpec(path="/GitHub/PhyBEARS.jl"))
using PhyBEARS

# Get the list of installed packages:
x = Pkg.installed()
list_of_installed_packages = collect(x);
println.(list_of_installed_packages)	# Messy

# Or:
x["PhyBEARS"]
# v"0.0.1"

