# push!(LOAD_PATH,"../src/")
using Nclusion
using Documenter
makedocs(
         sitename = "Nclusion.jl",
         modules  = [Nclusion],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/microsoft/Nclusion.jl",
)
