using Documenter, Nclusion
using Pkg
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()
makedocs(sitename="Nclusion.jl",
         modules=[Nclusion],
         pages=[
            "Home" => "index.md",
            "API Reference" => [
                "Math Functions" => "api/math_utils.md",
                "Processing Functions" => "api/processing.md",
                "CAVI Functions" => "api/cavi.md",
                "Custom Types" => "api/custom_types.md",
                "ELBO Functions" => "api/elbo_calculations.md",
                "Expectations Functions" => "api/expectations.md",
                "Initialization Functions" => "api/initialization.md",
                "Variational Update Functions" => "api/variational_updates.md",
            ]
         ],
         checkdocs = :none,

)