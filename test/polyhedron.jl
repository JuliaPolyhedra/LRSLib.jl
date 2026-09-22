using Test
using Polyhedra
using LRSLib
using GLPK
const polyhedra_test = joinpath(dirname(dirname(pathof(Polyhedra))), "test")

include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "polyhedra.jl"))
lpsolver = tuple()
@testset "Polyhedra tests" begin
    polyhedratest(LRSLib.Library(GLPK.Optimizer),
                  ["empty", "cubedecompose", "largedecompose", "recipe"])
end

@testset "Representation recovery" begin
    @testset "H-representation from a fresh matrix" begin
        A = Rational{BigInt}[1 -1 0; 0 1 0; 0 0 1]
        p = polyhedron(LiftedHRepresentation(A, BitSet([3])), LRSLib.Library())
        matrix = LRSLib.getinem(p, :Fresh)
        p.ine = nothing
        @test !hrepiscomputed(p)
        @test LRSLib.checkfreshness(matrix, :Fresh)

        rep = hrep(p)
        @test rep.A == A[[3, 1, 2], :]
        @test rep.linset == BitSet([1])
        @test hrepiscomputed(p)
        @test hrep(p) === rep
        @test !vrepiscomputed(p)
    end

    @testset "V-representation from inequalities" begin
        p = polyhedron(hrep(Rational{BigInt}[1 0; -1 0; 0 1; 0 -1],
                            Rational{BigInt}[1, 0, 1, 0]), LRSLib.Library())
        @test !vrepiscomputed(p)
        @test p.extm === nothing
        @test !p.vlinearitydetected

        rep = vrep(p)
        @test Set(points(rep)) == Set([[0, 0], [0, 1], [1, 0], [1, 1]])
        @test isempty(collect(rays(rep)))
        @test isempty(collect(lines(rep)))
        @test p.extm === nothing
        @test p.vlinearitydetected
        @test vrepiscomputed(p)
        @test vrep(p) === rep
    end
end
