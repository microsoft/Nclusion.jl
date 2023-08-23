push!(LOAD_PATH,"../src/")
using nclusion
using Documenter
makedocs(
         sitename = "nclusion",
         modules  = [nclusion],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/microsoft/nclusion",
)