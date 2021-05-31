#=
Julia wrapper/translation of lrs_solve_nash

=#
const RatOrInt = Union{Rational,Int}

"""
    nashsolve(A::AbstractMatrix, B::AbstractMatrix)

Compute all Nash equilibria (NE) for a two person noncooperative game with two
interleaved reverse search vertex enumeration steps. The inputs for the problem
are two m by n matrices `A`, `B` of integers or rationals. The first player is
the row player, the second is the column player. If row `i` and column `j` are
played, player 1 receives `A[i, j]` and player 2 receives `B[i, j]`.

# Examples

```julia-repl
julia> using LRSLib: nashsolve

julia> A = [0 6; 2 5; 3 3];

julia> B = [1 0; 0 2; 4 3];

julia> nashsolve(A, B)
3-element Array{Tuple{Array{Rational{BigInt},1},Array{Rational{BigInt},1}},1}:
 ([2//3, 1//3, 0//1], [1//3, 2//3])
 ([0//1, 1//3, 2//3], [2//3, 1//3])
 ([0//1, 0//1, 1//1], [1//1, 0//1])
```

The returned array contains tuples of NE. For example, in the first tuple,
player 1 uses row probabilities 2/3, 2/3, and 0, while player 2 uses column
probabilities 1/3 and 2/3.
"""
function nashsolve(A::AbstractMatrix{<:RatOrInt}, B::AbstractMatrix{<:RatOrInt})
    size(A) == size(B) || throw(ArgumentError(
        "A and B must have same size (got $(size(A)) and $(size(B)))"
    ))
    hr1 = buildrep(1, B')
    hr2 = buildrep(2, A)
    return solve_nash(hr1, hr2)
end

# Julia version of lrs_solve_nash in lrsnashlib.c
function solve_nash(hr1::HMatrix, hr2::HMatrix)
    NEs = NTuple{2,Vector{Rational{BigInt}}}[]

    # Step 1
    ccall((:resetNashSolver, liblrsnash), Cvoid, ())

    unsafe_field_store!(hr1.Q, :nash, Clrs_true)
    unsafe_field_store!(hr2.Q, :nash, Clrs_true)
    A2orig = unsafe_load(hr2.P).A

    linindex = zeros(Clong, unsafe_load(hr2.P).m+unsafe_load(hr2.P).d+2)

    # Step 2
    getfirstbasis(hr1)

    if unsafe_load(hr1.Q).dualdeg == Clrs_true
        @warn("Dual degenerate, ouput may be incomplete")
    end
    if unsafe_load(hr1.Q).unbounded == Clrs_true
        @warn("Unbounded starting dictionary for p1, output may be incomplete")
    end

    startcol = 0
    if (unsafe_load(hr1.Q).homogeneous == Clrs_true) &&
       (unsafe_load(hr1.Q).hull == Clrs_true)
        startcol += 1
    end

    col = startcol
    for i in startcol:(unsafe_load(hr1.Q).nredundcol-1)
        col += 1
    end

    # Step 3
    while true
        prune = @lrs_ccall(
            checkbound, Clong, (Ptr{Clrs_dic}, Ptr{Clrs_dat}), hr1.P, hr1.Q
        )
        output1 = getsolution(hr1, col)
        if prune == Clrs_false && output1 != nothing
            outputs2 = nash2_main(hr1, hr2, linindex)
            if !isempty(outputs2)
                for output2 in outputs2
                    NE = (output1[2:end-1], output2[2:end-1])
                    push!(NEs, NE)
                end
            end
        end
        getnextbasis(hr1) || break
    end

    unsafe_field_store!(hr2.Q, :Qhead, hr2.P)
    unsafe_field_store!(hr2.P, :A, A2orig)

    return NEs
end

# Julia version of nash2_main in lrsnashlib.c
function nash2_main(hr1::HMatrix, hr2::HMatrix, linindex::Vector{Clong})
    P1, Q1 = hr1.P, hr1.Q
    P2orig, Q2 = hr2.P, hr2.Q
    outputs = Vector{Rational{BigInt}}[]

    # Step 1
    P2 = @lrs_ccall(getdic, Ptr{Clrs_dic}, (Ptr{Clrs_dat},), Q2)
    P2ptr = Ref{Ptr{Clrs_dic}}(P2)
    @lrs_ccall2(
        copy_dict, Cvoid, (Ptr{Clrs_dat}, Ptr{Clrs_dic}, Ptr{Clrs_dic}),
        Q2, P2ptr[], P2orig
    )

    linearity = unsafe_load(Q2).linearity
    nlinearity = Clong(0)
    A1 = unsafe_load(P1).A
    for i in (unsafe_load(Q1).lastdv+1):(unsafe_load(P1).m)
        Row_i = unsafe_load(unsafe_load(P1).Row, i+1)
        if !iszero(extractbigintat(unsafe_load(A1, Row_i+1), 1))
            j = unsafe_load(
                unsafe_load(Q1).inequality,
                unsafe_load(unsafe_load(P1).B, i+1)-unsafe_load(Q1).lastdv+1
            )
            if ((unsafe_load(Q1).nlinearity == 0) ||
                (j < unsafe_load(unsafe_load(Q1).linearity, 1)))
                unsafe_store!(linearity, j, nlinearity+1)
                nlinearity += 1
            end
        end
    end

    if unsafe_load(Q1).nlinearity > 0
        unsafe_store!(
            linearity, unsafe_load(unsafe_load(Q1).linearity, 1), nlinearity+1
        )
        nlinearity += 1
    end

    for i in 1:nlinearity-1
        @lrs_ccall2(reorder, Cvoid, (Ptr{Clong}, Clong), linearity, nlinearity)
    end

    unsafe_field_store!(Q2, :nlinearity, nlinearity)
    unsafe_field_store!(Q2, :polytope, Clrs_false)

    # Step 2
    Lin2 = Ref{Clrs_mp_matrix}(C_NULL)
    ccall(
        (:lrs_getfirstbasis2, liblrsnash), Clong,
        (Ptr{Ptr{Clrs_dic}}, Ptr{Clrs_dat}, Ptr{Clrs_dic},
         Ptr{Clrs_mp_matrix}, Clong, Ptr{Clong}),
        P2ptr, Q2, P2orig, Lin2, Clrs_true, linindex
    ) == Clrs_true || @goto sayonara

    if unsafe_load(Q2).dualdeg == Clrs_true
        @warn("Dual degenerate, ouput may be incomplete")
    end

    if unsafe_load(Q2).unbounded == Clrs_true
        @warn("Unbounded starting dictionary for p2, output may be incomplete")
    end

    # startcol not in use
    # startcol = 0
    # if (unsafe_load(Q2).homogeneous == Clrs_true) &&
    #    (unsafe_load(Q2).hull == Clrs_true)
    #     startcol += 1
    # end

    # Step 3
    while true
        prune = @lrs_ccall(
            checkbound, Clong, (Ptr{Clrs_dic}, Ptr{Clrs_dat}), P2ptr[], Q2
        )
        col = 0
        output2 = getsolution(P2ptr[], Q2, col)
        if prune == Clrs_false && output2 != nothing
            if !iszero(view(output2, 2:length(output2)))
                push!(outputs, output2)
            end
        end
        @lrs_ccall(
            getnextbasis, Clong, (Ptr{Ptr{Clrs_dic}}, Ptr{Clrs_dat}, Clong),
            P2ptr, Q2, prune
        ) == Clrs_true || break
    end

    @label sayonara
    @lrs_ccall(free_dic, Cvoid, (Ptr{Clrs_dic}, Ptr{Clrs_dat}), P2ptr[], Q2)
    return outputs
end

# for player_idx==i, opponent_payoff_matrix is that of player 3-i
function buildrep(player_idx::Integer,
                  opponent_payoff_matrix::AbstractMatrix{<:RatOrInt})
    sz = size(opponent_payoff_matrix)
    n = sz[2] + 2
    m = sum(sz) + 1
    M = zeros(Rational{BigInt}, m, n)

    if player_idx == 1
        nonnegativity_subarray = view(M, 1:sz[2], 2:sz[2]+1)
        constraint_subarray = view(M, sz[2]+1:sum(sz), 2:n)
    else
        constraint_subarray = view(M, 1:sz[1], 2:n)
        nonnegativity_subarray = view(M, sz[1]+1:sum(sz), 2:sz[2]+1)
    end

    # FillNonnegativityRows
    nonnegativity_subarray[diagind(nonnegativity_subarray)] .= 1

    # FillConstraintRows
    constraint_subarray[:, 1:end-1] = -opponent_payoff_matrix
    constraint_subarray[:, end] .= 1

    # FillLinearityRow
    M[end, 1] = 1
    M[end, 2:end-1] .= -1
    linearity = BitSet(m)

    P, Q = initmatrix(M, linearity, true)
    return HMatrix(sz[2]+1, P, Q)
end


# Reading games

# Function to convert String to Rational
function tl_readrat!(v::AbstractArray{<:Rational,0}, s::AbstractString)
    div = split(s, '/')
    num, den = 0, 1
    if length(div) == 1
        num = parse(Int, div[1])
    elseif length(div) == 2
        num, den = parse.(Int, div)
    else
        return false
    end
    v[1] = num//den
    return true
end

# readGame in lrsnash.c
function readgame(filename::AbstractString)
    entries = split(read(filename, String))
    num_entries = length(entries)
    read_error_msg = "Premature end of input file '$(filename)'."
    num_entries < 2 && error(read_error_msg)

    nr::Int, nc::Int = parse.(Int, entries[1:2])
    num_entries < 2 + (nr*nc)*2 && error(read_error_msg)
    payoff_matrices::NTuple{2,Vector{Rational{BigInt}}} =
        ntuple(i -> Vector{Rational{BigInt}}(undef, nr*nc), Val(2))

    i = 3
    for p in 1:2
        for j in 1:nr*nc
            tl_readrat!(view(payoff_matrices[p], j), entries[i]) ||
                @warn(String("String '$(entries[i])' is not a rational number",
                             "in file $(filename)."))
            i += 1
        end
    end

    num_entries > 2 + (nr*nc)*2 && @warn("Excess data in file $(filename).")
    return ntuple(i -> transpose(reshape(payoff_matrices[i], nc, nr)), Val(2))
end

# For standard format
nashsolve(filename::AbstractString) = nashsolve(readgame(filename)...)

# For legacy format
function nashsolve(filename1::AbstractString, filename2::AbstractString)
    hrs = HMatrix.([filename1, filename2])
    Q2 = unsafe_load(hrs[2].Q)

    # From lrs_solve_nash_legacy
    if Q2.nlinearity > 0
        ccall((:free, liblrs), Cvoid, (Ptr{Cvoid},), Q2.linearity)
    end
    ptr = @lrs_ccall2(
        xcalloc, Ptr{Cvoid}, (Clong, Clong, Clong, Cstring),
        Q2.m+2, sizeof(Clong), @__LINE__, @__FILE__
    )
    unsafe_field_store!(hrs[2].Q, :linearity, ptr)

    return solve_nash(hrs...)
end
