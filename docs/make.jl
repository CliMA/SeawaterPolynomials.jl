using
  Documenter,
  SeawaterPolynomials

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://clima.github.io/SeawaterPolynomials/dev/"
)

makedocs(
    modules = [SeawaterPolynomials],
    doctest = false,
      clean = true,
   checkdocs = :all,
     format = format,
    authors = "Gregory L. Wagner and Ali Ramadhan",
   sitename = "SeawaterPolynomials.jl",
  
      pages = Any[
              "Home" => "index.md",
              "DocStrings" => Any[
                  "man/types.md",
                  "man/functions.md"]
                 ]
)

deploydocs(
  repo = "github.com/CliMA/SeawaterPolynomials.jl.git",
)
