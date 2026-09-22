using LRSLib: nashsolve, readgame


@testset "Tests for nash.jl" begin

    @testset "Tests for nashsolve" begin
        A = [3 3; 2 5; 0 6]
        B = [3 2; 2 6; 3 1]
        NEs = [
            ([4//5, 1//5, 0//1], [2//3, 1//3]),
            ([0//1, 1//3, 2//3], [1//3, 2//3]),
            ([1//1, 0//1, 0//1], [1//1, 0//1])
        ]

        NEs_computed = @inferred(nashsolve(A, B))
        @test sort(NEs_computed) == sort(NEs)

        A_neg = A .- 10
        B_neg = B .- 10
        NEs_computed = @inferred(nashsolve(A_neg, B_neg))
        @test sort(NEs_computed) == sort(NEs)

        B_rat = B .// 7
        NEs_computed = @inferred(nashsolve(A, B_rat))
        @test sort(NEs_computed) == sort(NEs)

        @testset "$T" for T in [BigInt, Int32, Int8, UInt8, UInt64,
                                Rational{BigInt}, Rational{Int32},
                                Rational{UInt8}]
            NEs_computed = @inferred(nashsolve(T.(A), T.(B)))
            @test sort(NEs_computed) == sort(NEs)
        end
    end

    @testset "Tests for nashsolve with degenerate game" begin
        A = [3 3; 2 5; 0 6]
        B = [3 3; 2 6; 3 1]
        NEs = [
            ([1//1, 0//1, 0//1], [1//1, 0//1]),
            ([1//1, 0//1, 0//1], [2//3, 1//3]),
            ([0//1, 1//3, 2//3], [1//3, 2//3])
        ]

        NEs_computed = @inferred(nashsolve(A, B))
        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "Tests for nashsolve with standard format" begin
        filename = "game"
        NEs = [
            ([2//3, 1//3, 0//1], [1//3, 2//3]),
            ([0//1, 1//3, 2//3], [2//3, 1//3]),
            ([0//1, 0//1, 1//1], [1//1, 0//1])
        ]

        NEs_computed = @inferred(nashsolve(filename))
        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "Tests for nashsolve with legacy format" begin
        filename1 = "game1"
        filename2 = "game2"
        NEs = [
            ([2//3, 1//3, 0//1], [1//3, 2//3]),
            ([0//1, 1//3, 2//3], [2//3, 1//3]),
            ([0//1, 0//1, 1//1], [1//1, 0//1])
        ]

        NEs_computed = @inferred(nashsolve(filename1, filename2))
        @test sort(NEs_computed) == sort(NEs)
    end

    @testset "Tests for readgame" begin
        filename = "game"
        A = Rational{Int}[0 6
                          2 5
                          3 3]
        B = Rational{Int}[1 0
                          0 2
                          4 3]
        @inferred(readgame(filename)) == (A, B)
    end

end
