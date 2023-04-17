#=
This adds and precompiles the Finch and Plots packages
=#

using Pkg

Pkg.add("https://github.com/paralab/Finch")
Pkg.add("Plots")

using Finch
using Plots
