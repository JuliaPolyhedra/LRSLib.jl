using GLPK
const polyhedra_test = joinpath(dirname(dirname(pathof(Polyhedra))), "test")

include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "polyhedra.jl"))
lpsolver = tuple()
@testset "Polyhedra tests" begin
    polyhedratest(LRSLib.Library(with_optimizer(GLPK.Optimizer)),
                  ["empty", "cubedecompose", "largedecompose"])
end
