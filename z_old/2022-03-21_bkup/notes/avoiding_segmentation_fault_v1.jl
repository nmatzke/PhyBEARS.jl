#######################################################
# Julia setup from scratch
# (e.g. after a new Julia version) 
# Nick Matzke 2020-08-30
#######################################################


#######################################################
# SHORT VERSION:
#######################################################

# 1. Delete entire ~/.julia directory and old julia versions

# 2. Install latest julia

# 3. Add the packages in 
# /GitHub/BioGeoJulia.jl/notes/Fezzik_precompile_packages_v1.jl

# 4. Close, reopen Julia, follow instructions to do (long) Fezzik 
#    precompile, then close/open as instructed.


# 5. In your startup.jl file, do this to load BioGeoJulia:

print("Loading BioGeoJulia...\n")

# DON'T DO THIS:
#Pkg.rm("BioGeoJulia")
#Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
#Pkg.add(PackageSpec(url="https://github.com/nmatzke/BioGeoJulia.jl"))

# DON'T ADD BioGeoJulia, INSTEAD DO:
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

using BioGeoJulia

using BioGeoJulia.TrUtils
using BioGeoJulia.StateSpace
using BioGeoJulia.TreePass
using BioGeoJulia.SSEs

# 6. Run the example BioGeoJulia functions here:
#   include("/GitHub/BioGeoJulia.jl/notes/startup_basic_tests_v1.jl")
#
# ...rather than in startup (speeds it up).

# 7. If startup.jl ever gives problems, start with:

# julia --startup=no



#######################################################
# Other tricks below
# (none worked, but may be handy for other issues)
#######################################################




# Setting up Fezzik:
/GitHub/BioGeoJulia.jl/notes/Fezzik_precompile_packages_v1.jl

# (Set up Fezzik, close Julia, reopen, install packages)


Things to try to remove this error:

Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))

signal (11): Segmentation fault: 11
in expression starting at /Users/nickm/.julia/config/startup.jl:131
unknown function (ip: 0x0)
Allocations: 55354997 (Pool: 55342213; Big: 12784); GC: 32
Segmentation fault: 11


Pkg.add(PackageSpec(url="https://github.com/nmatzke/BioGeoJulia.jl"))
    Cloning git-repo `https://github.com/nmatzke/BioGeoJulia.jl`

signal (11): Segmentation fault: 11
in expression starting at REPL[7]:1
unknown function (ip: 0x0)
Allocations: 62471383 (Pool: 62455325; Big: 16058); GC: 40



Pkg.add(PackageSpec(url="https://github.com/JuliaLang/Example.jl"))
    Cloning git-repo `https://github.com/JuliaLang/Example.jl`

signal (11): Segmentation fault: 11
in expression starting at REPL[2]:1
unknown function (ip: 0x0)
Allocations: 2571649 (Pool: 2570883; Big: 766); GC: 3


Pkg.add(url="https://github.com/JuliaLang/Example.jl")

signal (11): Segmentation fault: 11
in expression starting at REPL[4]:1
unknown function (ip: 0x0)
Allocations: 2571653 (Pool: 2570884; Big: 769); GC: 3
Segmentation fault: 11





#######################################################
# Things to try deleting
#######################################################

Search ~/.julia/ on BioGeoJulia and delete old version

cd ~/.julia/
open .


# Registry
~/.julia/registries/General

# Compiled versions:




#######################################################
# Find the uuid (universal ID for the package)(
#######################################################
using Pkg
Pkg.PackageSpec("BioGeoJulia").uuid

using DataFrames
Pkg.PackageSpec("DataFrames").uuid
Pkg.PackageSpec("DataFrames")

Removing packages from Julia:

# 1. Remove from the current Julia version
#    Manifest.toml	Project.toml
#    ...with in-Julia commands:
using Pkg
]remove BioGeoJulia
# (backspace to exit the "]" Pkg commands)

# 2. But, you may have other manifests hiding:
Pkg.gc()

more ~/.julia/environments/v1.4/Manifest.toml

# Remove any old version:
cd ~/.julia/environments/v1.4/
rm *.toml
cd ..
rm -r v1.4


# Re-do the General registry:
registry remove General
   Removing registry `General` from ~/.julia/registries/General


Pkg.status()


registry remove General
registry add General

JULIA_PKG_SERVER=pkg.julialang.org julia

using Pkg
Pkg.rm("BioGeoJulia.jl")
Pkg.gc()

using LibGit2

LibGit2.version()



# Check if any packages can be installed:
using Pkg
Pkg.add(PackageSpec(url="https://github.com/JuliaLang/Example.jl"))



# If you:
cd /Users/nickm/.julia/registries
rm -r General

# Then...
Pkg.add(PackageSpec(url="https://github.com/JuliaLang/Example.jl"))

# You get:
ERROR: GitError(Code:ENOTFOUND, Class:Repository, could not find repository from '/Users/nickm/.julia/clones/15037770950778237535')


# If you do:
JULIA_PKG_SERVER=pkg.julialang.org julia --startup=no
]
registry remove General

# You get:
# registry `General` not found.


# But this still works:
registry remove General
registry add General

######################################################################## 100.0%
      Added registry `General` to `~/.julia/registries/General`

(@v1.5) pkg> 


