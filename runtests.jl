if !haskey(Pkg.installed(), "SeaMass")
  Pkg.clone("https://github.com/biospi/SeaMass.jl")
else
  Pkg.checkout("SeaMass")
end

Pkg.test("SeaMass")
