using LRSLib: initmatrix, HMatrix, extractrow, getd, getfirstbasis


@testset "Tests for matrix.jl" begin

    @testset "Test for issue #23" begin
        M = Rational{BigInt}[
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            0 -1 0 -4 1
            0 0 -2 -3 1
            -1 1 1 1 0
        ]
        m, d = size(M, 1), size(M, 2) - 1
        linearity = BitSet(m)
        P, Q = initmatrix(M, linearity, true)
        hmatrix = HMatrix(d, P, Q)

        row = extractrow(hmatrix, 1)
        @test length(row[1]) == getd(hmatrix)

        getfirstbasis(hmatrix)
        row = extractrow(hmatrix, 1)
        @test length(row[1]) == getd(hmatrix)
    end

end
