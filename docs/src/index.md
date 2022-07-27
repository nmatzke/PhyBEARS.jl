# Example

Example Julia package repo.

```@autodocs
Modules = [Example, BioGeoJulia]
```


``` NOTE: If a module is missing, this causes e.g.:

Pkg.add(PackageSpec(path="/GitHub/BioGeoJulia.jl"))
  Updating git-repo `/GitHub/BioGeoJulia.jl`
ERROR: expected the file `src/BioGeoJulia.jl` to exist for package BioGeoJulia at /var/folders/_l/ph0n0x6d1l7_8ywbfqrmpdtm0000gp/T/jl_PgtjWR
Stacktrace:
 [1] pkgerror(::String) at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/Types.jl:113
 [2] read_package(::String) at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/Types.jl:480
 [3] parse_package!(::Pkg.Types.Context, ::Pkg.Types.PackageSpec, ::String) at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/Types.jl:752
 [4] resolve_repo_add!(::Pkg.Types.Context, ::Pkg.Types.PackageSpec) at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/Types.jl:711
 [5] handle_repo_add!(::Pkg.Types.Context, ::Pkg.Types.PackageSpec) at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/Types.jl:724
 [6] handle_repos_add!(::Pkg.Types.Context, ::Array{Pkg.Types.PackageSpec,1}) at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/Types.jl:743
 [7] #add#25(::Bool, ::Pkg.BinaryPlatforms.MacOS, ::Base.Iterators.Pairs{Union{},Union{},Tuple{},NamedTuple{(),Tuple{}}}, ::typeof(Pkg.API.add), ::Pkg.Types.Context, ::Array{Pkg.Types.PackageSpec,1}) at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/API.jl:90
 [8] add(::Pkg.Types.Context, ::Array{Pkg.Types.PackageSpec,1}) at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/API.jl:72
 [9] #add#20 at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/API.jl:69 [inlined]
 [10] add(::Pkg.Types.PackageSpec) at /Users/sabae/buildbot/worker/package_macos64/build/usr/share/julia/stdlib/v1.3/Pkg/src/API.jl:66
 [11] top-level scope at REPL[16]:1

```
