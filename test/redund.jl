@testset "Tests for redundancy removal" begin

    @testset "Test for issue #52" begin
        A = [0 1 0; -1 0 0; 0 -1 0; 0 0 -1; 1 0 0]
        b = [2, -1, -1, -1, 2]
        linset = BitSet(5)
        hr = hrep(A, b, linset)

        V = [2 1 1;
             2 2 1]
        R = [0 0 1]
        exp = vrep(V, R)

        @testset "removevredundancy!" begin
            p = polyhedron(hr, LRSLib.Library())
            removevredundancy!(p)
            @test nrays(p) == nrays(exp)
        end

        @testset "detectvlinearity! and then removevredundancy!" begin
            p = polyhedron(hr, LRSLib.Library())
            detectvlinearity!(p)
            removevredundancy!(p)
            @test nrays(p) == nrays(exp)
        end
    end

end
