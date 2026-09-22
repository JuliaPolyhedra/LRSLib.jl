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

    # The V-representation computed by LRS has the rays interleaved with the
    # vertices and contains duplicated rays
    @testset "Unbounded polyhedron with duplicated rays" begin
        A = [-1 0 1 0; -1 -1 0 1; 1 0 -1 1; 1 1 -1 1; -1 1 1 -1]
        b = [2, 0, 1, 2, 1]
        hr = hrep(A, b)
        exp_points = [[1, 3//2, 3, 5//2], [2, 1, 4, 3]]
        exp_rays = [[0, -1, 0, -1], [1, 0, 1, 0], [-1, 0, -2, -1], [1, -1, 1, 0]]

        @testset "removevredundancy!" begin
            p = polyhedron(hr, LRSLib.Library())
            removevredundancy!(p)
            @test npoints(p) == 2
            @test nrays(p) == 4
            @test Set(collect(points(p))) == Set(exp_points)
            @test Set(coord.(rays(p))) == Set(exp_rays)
        end

        @testset "isredundant" begin
            p = polyhedron(hr, LRSLib.Library())
            vrep(p)
            @test npoints(p) == 2
            @test nrays(p) == 6
            for idx in eachindex(points(p))
                @test !Polyhedra.isredundant(p, idx)
            end
            for r in exp_rays
                idxs = [idx for idx in eachindex(rays(p)) if coord(get(p, idx)) == r]
                # At most one copy of a duplicated ray is redundant
                @test count(idx -> Polyhedra.isredundant(p, idx), idxs) == length(idxs) - 1
            end
        end
    end

    @testset "Rows of the LRS matrices are in the order of the lifted representations" begin
        T = Rational{BigInt}
        rows(m, N, offset) = [LRSLib.extractrow(unsafe_load(m.P), unsafe_load(m.Q), N, i, offset) for i in 1:length(m)]
        # ray, point, ray, point, line
        R = T[0 0 1; 1 0 0; 0 1 1; 1 1 0; 0 1 0]
        p = polyhedron(LiftedVRepresentation{T}(R, BitSet([5])), LRSLib.Library())
        extm = LRSLib.getextm(p)
        @test rows(extm, 2, 1) == [R[i, :] for i in 1:5]
        @test LRSLib.linset(extm) == BitSet([5])
        # ray, point, ray, point, ray: only the ray (1, 1) and the point (1, 1) are redundant
        R = T[0 1 0; 1 0 0; 0 1 1; 1 1 1; 0 0 1]
        p = polyhedron(LiftedVRepresentation{T}(R), LRSLib.Library())
        @test [idx.value for idx in eachindex(rays(p)) if Polyhedra.isredundant(p, idx)] == [3]
        @test [idx.value for idx in eachindex(points(p)) if Polyhedra.isredundant(p, idx)] == [4]
        # hyperplane in the second row
        A = T[1 1 0; 0 0 1; 1 0 1]
        p = polyhedron(LiftedHRepresentation{T}(A, BitSet([2])), LRSLib.Library())
        inem = LRSLib.getinem(p)
        @test rows(inem, 2, 0) == [A[i, :] for i in 1:3]
        @test LRSLib.linset(inem) == BitSet([2])
    end

    @testset "Test for removehredundancy!" begin
        A = [0 0; -1 0; 0 -1]
        b = [1, 0, 0]
        hr = hrep(A, b)
        p = polyhedron(hr, LRSLib.Library())
        @test !([-1, 1] in p)
        removehredundancy!(p)
        @test !([-1, 1] in p)
    end

    @testset "removehredundancy! with hyperplane" begin
        A = [1 0 0; 1 1 0; 0 1 0; -1 0 0; 0 -1 0; 0 0 1]
        b = [1, 5, 1, 0, 0, 0]
        p = polyhedron(hrep(A, b, BitSet(6)), LRSLib.Library())
        removehredundancy!(p)
        @test nhyperplanes(p) == 1
        @test nhalfspaces(p) == 4
        @test !([2, 2, 0] in p)
        @test !([1//2, 1//2, 7] in p)
        @test [1//2, 1//2, 0] in p
    end

end
