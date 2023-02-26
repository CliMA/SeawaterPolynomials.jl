using
  Documenter,
  SeawaterPolynomials

format = Documenter.HTML(
    collapselevel = 2,
       prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://clima.github.io/SeawaterPolynomials/stable/",
       mathengine = MathJax3()
)

makedocs(
    modules = [SeawaterPolynomials],
    doctest = false,
      clean = true,
   checkdocs = :all,
     format = format,
    authors = "Climate Modeling Alliance and contributors",
   sitename = "SeawaterPolynomials.jl",

      pages = Any[
              "Home" => "index.md",
              "API" => "API.md",
                 ]
)

deploydocs(
          repo = "github.com/CliMA/SeawaterPolynomials.jl.git",
      versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
     forcepush = true,
  push_preview = true,
     devbranch = "main"
)
