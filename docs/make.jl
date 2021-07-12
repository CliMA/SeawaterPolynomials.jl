using
  Documenter,
  SeawaterPolynomials

format = Documenter.HTML(
    collapselevel = 1,
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
    authors = "Gregory L. Wagner and Ali Ramadhan",
   sitename = "SeawaterPolynomials.jl",
  
      pages = Any[
              "Home" => "index.md",
                 ]
)

deploydocs(
        repo = "github.com/CliMA/SeawaterPolynomials.jl.git",
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
push_preview = true
)
