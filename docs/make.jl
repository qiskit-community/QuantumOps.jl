using Documenter
using QuantumOps

makedocs(
    sitename = "QuantumOps",
    format = Documenter.HTML(),
    modules = [QuantumOps]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
