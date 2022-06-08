using Pkg
Pkg.activate("documentation")

# Remove old `docs/`
rm("docs"; recursive = true, force = true)

# Build new docs
include("documentation/make.jl")

# Copy new docs to replace old docs.
cp("documentation/build/", "docs/"; force = true)
