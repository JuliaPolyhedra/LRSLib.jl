export LRSMatrix, LRSInequalityMatrix, LRSGeneratorMatrix, setdebug
import Base.eltype, Base.copy

function lrs_alloc_dat()
    @lrs_ccall alloc_dat Ptr{Clrs_dat} (Ptr{Cchar},) C_NULL
end

function lrs_alloc_dic(Q::Ptr{Clrs_dat})
    @lrs_ccall alloc_dic Ptr{Clrs_dic} (Ptr{Clrs_dat},) Q
end

function initmatrix(filename::AbstractString)
    Q = lrs_alloc_dat()
    # The first element does not matter
    argv = ["", filename]
    ok = Clrs_true == @lrs_ccall read_dat Clong (Ptr{Clrs_dat}, Cint, Ptr{Ptr{Cchar}}) Q length(argv) argv
    if !ok
        error("Invalid file $filename")
    end
    P = lrs_alloc_dic(Q)
    ok = Clrs_true == @lrs_ccall read_dic Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}) P Q
    if !ok
        error("Invalid file $filename")
    end
    (P,Q)
end

# FIXME still needed ?
function initmatrix(M::Matrix{Rational{BigInt}}, linset, Hrep::Bool)
    m = Clong(size(M, 1))
    n = Clong(size(M, 2))
    Q = lrs_alloc_dat()
    @lrs_ccall init_dat Void (Ptr{Clrs_dat}, Clong, Clong, Clong) Q m n Clong(Hrep ? 0 : 1)
    P = lrs_alloc_dic(Q)
    #Q->getvolume= TRUE; # compute the volume # TODO cheap do it
    for i in 1:m
        #   num = map(x -> GMPInteger(x.num.alloc, x.num.size, x.num.d), M[i,:])
        #   den = map(x -> GMPInteger(x.den.alloc, x.den.size, x.den.d), M[i,:])
        ineq = !(i in linset)
        setrow(P, Q, i, M[i,:], !(i in linset))
        #   @lrs_ccall set_row_mp Void (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong, Clrs_mp_vector, Clrs_mp_vector, Clong) P Q i num den Clong(ineq)
    end
    # This is the objective. If I have no objective LRS might fail
    if Hrep
        setrow(P, Q, 0, ones(Rational{BigInt}, n), true)
    end
    (P, Q)
end


function fillmatrix(inequality::Bool, P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}, itr1, offset=0)
    for (i, item) in enumerate(itr1)
        a = vec(coord(lift(item)))
        setrow(P, Q, offset+i, inequality ? -a : a, !islin(item))
    end
end

function initmatrix(inequality, itr1, itr2=nothing)
    n = fulldim(itr1)+1
    m = length(itr1)
    if !(itr2 === nothing)
        m += length(itr2)
    end
    Q = lrs_alloc_dat()
    @lrs_ccall init_dat Void (Ptr{Clrs_dat}, Clong, Clong, Clong) Q m n Clong(inequality ? 0 : 1)
    P = lrs_alloc_dic(Q)
    #Q->getvolume= TRUE; # compute the volume # TODO cheap do it
    fillmatrix(inequality, P, Q, itr1)
    if !(itr2 === nothing)
        fillmatrix(inequality, P, Q, itr2, length(itr1))
    end
    # This is the objective. If I have no objective LRS might fail
    if inequality
        setrow(P, Q, 0, ones(Rational{BigInt}, n), true)
    end
    (P, Q)
end


# Representation

type LRSLinearitySpace{N}
    Lin::Clrs_mp_matrix
    nlin::Int
    n::Int
    hull::Bool
    homogeneous::Bool

    function LRSLinearitySpace{N}(Lin::Clrs_mp_matrix, nlin, n, hull, homogeneous) where {N}
        m = new{N}(Lin, nlin, n, hull, homogeneous)
        finalizer(m, myfree)
        m
    end
end

type LRSInequalityMatrix{N} <: HRepresentation{N, Rational{BigInt}}
    P::Ptr{Clrs_dic}
    Q::Ptr{Clrs_dat}
    status::Symbol
    lin::Nullable{LRSLinearitySpace{N}}
    function LRSInequalityMatrix{N}(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}) where {N}
        m = new{N}(P, Q, :AtNoBasis, nothing)
        finalizer(m, myfree)
        m
    end
end
changefulldim{N}(::Type{LRSInequalityMatrix{N}}, NewN) = LRSInequalityMatrix{NewN}
decomposedfast(ine::LRSInequalityMatrix) = false
eltype{N}(::Type{LRSInequalityMatrix{N}}) = Rational{BigInt}
eltype(::LRSInequalityMatrix) = Rational{BigInt}

type LRSGeneratorMatrix{N} <: VRepresentation{N, Rational{BigInt}}
    P::Ptr{Clrs_dic}
    Q::Ptr{Clrs_dat}
    status::Symbol
    lin::Nullable{LRSLinearitySpace{N}}
    function LRSGeneratorMatrix{N}(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}) where {N}
        m = new{N}(P, Q, :AtNoBasis, nothing)
        finalizer(m, myfree)
        m
    end
end
changefulldim{N}(::Type{LRSGeneratorMatrix{N}}, NewN) = LRSGeneratorMatrix{NewN}
decomposedfast(ine::LRSGeneratorMatrix) = false
eltype{N}(::Type{LRSGeneratorMatrix{N}}) = Rational{BigInt}
eltype(::LRSGeneratorMatrix) = Rational{BigInt}

const LRSMatrix{N} = Union{LRSInequalityMatrix{N}, LRSGeneratorMatrix{N}}

function linset(matrix::LRSMatrix)
    extractinputlinset(unsafe_load(matrix.Q))
end
function Base.length(matrix::LRSMatrix)
    unsafe_load(matrix.P).m
end

LRSMatrix(hrep::HRepresentation) = LRSInequalityMatrix(hrep)
LRSMatrix(vrep::VRepresentation) = LRSGeneratorMatrix(vrep)

function checkfreshness(m::LRSMatrix, fresh::Symbol)
    fresh == :AnyFreshNess ||
    (fresh == :Fresh && m.status in [:AtNoBasis, :AtFirstBasis, :Empty]) ||
    (fresh == :AlmostFresh && m.status in [:AtNoBasis, :AtFirstBasis, :Empty, :RedundancyChecked])
end

function myfree(l::LRSLinearitySpace)
    if l.nlin > 0
        @lrs_ccall clear_mp_matrix Void (Clrs_mp_matrix, Clong, Clong) l.Lin l.nlin l.n
    end
end

function myfree(m::LRSMatrix)
    @lrs_ccall free_dic_and_dat Void (Ptr{Clrs_dic}, Ptr{Clrs_dat}) m.P m.Q
end

# H-representation

#LRSInequalityMatrix{N,T}(rep::Rep{N,T}) = LRSInequalityMatrix{N,polytypefor(T), mytypefor(T)}(rep) # TODO finish this line

function LRSInequalityMatrix(filename::AbstractString)
    P, Q = initmatrix(filename)
    LRSInequalityMatrix{unsafe_load(P).d}(P, Q)
end
LRSInequalityMatrix{N}(rep::HRep{N}) = LRSInequalityMatrix{N}(rep)

copy{N}(ine::LRSInequalityMatrix{N}) = LRSInequalityMatrix{N}(hreps(ine))

function (::Type{LRSInequalityMatrix{N}}){N}(it::HRepIterator{N, Rational{BigInt}})
    P, Q = initmatrix(true, it)
    LRSInequalityMatrix{N}(P, Q)
end
function (::Type{LRSInequalityMatrix{N}}){N}(; eqs=nothing, ineqs=nothing)
    P, Q = initmatrix(true, eqs, ineqs)
    LRSInequalityMatrix{N}(P, Q)
end

nhreps(matrix::LRSInequalityMatrix) = length(matrix)
neqs(matrix::LRSInequalityMatrix) = unsafe_load(matrix.Q).nlinearity
nineqs(matrix::LRSInequalityMatrix) = length(matrix) - neqs(matrix)

starthrep(ine::LRSInequalityMatrix) = 1
donehrep(ine::LRSInequalityMatrix, state) = state > length(ine)
nexthrep(ine::LRSInequalityMatrix, state) = (extractrow(ine, state), state+1)

# V-representation

function LRSGeneratorMatrix(filename::AbstractString)
    d, P, Q = initmatrix(filename)
    LRSGeneratorMatrix{unsafe_load(P).d-1}(P, Q)
end
LRSGeneratorMatrix{N}(rep::VRep{N}) = LRSGeneratorMatrix{N}(rep)

copy{N}(ext::LRSGeneratorMatrix{N}) = LRSGeneratorMatrix{N}(vreps(ext))

function (::Type{LRSGeneratorMatrix{N}}){N}(it::VRepIterator{N, Rational{BigInt}})
    P, Q = initmatrix(false, it)
    LRSGeneratorMatrix{N}(P, Q)
end
function (::Type{LRSGeneratorMatrix{N}}){N}(; rays=nothing, points=nothing)
    P, Q = initmatrix(false, rays, points)
    LRSGeneratorMatrix{N}(P, Q)
end

nvreps(ext::LRSGeneratorMatrix) = length(ext)

startvrep(ext::LRSGeneratorMatrix) = 1
donevrep(ext::LRSGeneratorMatrix, state) = state > length(ext)
nextvrep(ext::LRSGeneratorMatrix, state) = (extractrow(ext, state), state+1)

#I should also remove linearity (should I remove one if hull && homogeneous ?)
#getd{N}(m::LRSInequalityMatrix{N}) = N
#getd{N}(m::LRSGeneratorMatrix{N}) = N+1
#Let's do it the easy way
getd{N}(m::LRSInequalityMatrix{N}) = unsafe_load(m.P).d
getd{N}(m::LRSGeneratorMatrix{N}) = unsafe_load(m.P).d

function setrow(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}, i::Int, row::Vector{Rational{BigInt}}, ineq::Bool)
    num = map(x -> GMPInteger(x.num.alloc, x.num.size, x.num.d), row)
    den = map(x -> GMPInteger(x.den.alloc, x.den.size, x.den.d), row)
    @lrs_ccall set_row_mp Void (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong, Clrs_mp_vector, Clrs_mp_vector, Clong) P Q i num den Clong(ineq)
end

function setdebug(m::LRSMatrix, debug::Bool)
    @lrs_ccall setdebug Void (Ptr{Clrs_dat}, Clong) m.Q (debug ? Clrs_true : Clrs_false)
end

function extractrow(P::Clrs_dic, Q::Clrs_dat, N, i, offset)
    #d = Q.n-offset-1 # FIXME when it is modified...
    a = Vector{Rational{BigInt}}(N+1)
    gcd = extractbigintat(Q.Gcd, 1+i) # first row is the objective
    lcm = extractbigintat(Q.Lcm, 1+i)
    row = unsafe_load(P.A, 1+i)
    extractthisrow(i::Int) = (extractbigintat(row, offset+i) * gcd) // lcm
    for j in 1:N+1
        a[j] = extractthisrow(j)
    end
    a
end

function warn_fresh(m::LRSMatrix)
    if !checkfreshness(m, :Fresh)
        warn("Extracting the rows of an LRS matrix after it has been used for representation conversion does not give the correct elements of the polyhedron it represents")
    end
end

function extractrow{N}(matrix::LRSInequalityMatrix{N}, i::Int)
    P = unsafe_load(matrix.P)
    Q = unsafe_load(matrix.Q)
    b = extractrow(P, Q, N, i, 0)
    β = b[1]
    a = -b[2:end]
    if isininputlinset(Q, i)
        HyperPlane(a, β)
    else
        HalfSpace(a, β)
    end
end

function extractrow{N}(matrix::LRSGeneratorMatrix{N}, i::Int)
    P = unsafe_load(matrix.P)
    Q = unsafe_load(matrix.Q)
    #d = Q.n-offset-1 # FIXME when it is modified...
    b = extractrow(P, Q, N, i, 1)
    ispoint = b[1]
    islin = isininputlinset(Q, i)
    @assert ispoint == zero(Rational{BigInt}) || ispoint == one(Rational{BigInt})
    a = b[2:end]
    if ispoint == zero(Rational{BigInt})
        if islin
            Line(a)
        else
            Ray(a)
        end
    else
        if islin
            SymPoint(a)
        else
            a
        end
    end
end

function isininputlinset(Q::Clrs_dat, j)
    for i in 1:Q.nlinearity
        if j == unsafe_load(Q.linearity, i)
            return true
        end
    end
    false
end

function extractinputlinset(Q::Clrs_dat)
    linset = IntSet([])
    for i in 1:Q.nlinearity
        push!(linset, unsafe_load(Q.linearity, i))
    end
    linset
end

function extractoutputlinset(Q::Clrs_dat)
    k = (Q.hull == Clrs_true && Q.homogeneous == Clrs_true) ? 1 : 0
    nredundcol = Q.nredundcol
    IntSet(1:(nredundcol-k))
end

# FIXME The only think that is done is that the linearities given that were redundant have been
# removed from the linearity set so that they can be marked as redundant inequalities.
# New linearities are detected but getinputlinsubset does not give them.
# I should check in redundcols
function getinputlinsubset(m::LRSMatrix)
    if m.status == :AtNoBasis
        getfirstbasis(m)
    end
    linset(m)
end
function getoutputlinset(m::LRSMatrix)
    if m.status == :AtNoBasis
        getfirstbasis(m)
    end
    extractoutputlinset(unsafe_load(m.Q))
end


function convertoutput(x::Clrs_mp_vector, n, hull)
    first = extractbigintat(x, 1)
    rest = Vector{BigInt}(n-1)
    for i = 2:n
        rest[i-1] = extractbigintat(x, i)
    end
    if hull || first == 0
        Rational{BigInt}[first; rest]
    else
        [one(Rational{BigInt}); rest // first]
    end
end

function getmat{N}(lin::LRSLinearitySpace{N})
    startcol = lin.hull && lin.homogeneous ? 2 : 1 # col zero not treated as redundant
    A = Matrix{BigInt}(lin.nlin-startcol+1, lin.n)
    for col in startcol:lin.nlin # print linearity space */
        A[col-startcol+1,:] = convertoutput(unsafe_load(lin.Lin, col), lin.n, lin.hull)
    end
    A
end

function getfirstbasis{N}(m::LRSMatrix{N})
    Lin = Ref{Clrs_mp_matrix}(C_NULL)
    Pptr = Ref{Ptr{Clrs_dic}}(m.P)
    # The "Clrs_true" at the last argument since that it should not be verbose
    found = Clrs_true == (@lrs_ccall getfirstbasis Clong (Ptr{Ptr{Clrs_dic}}, Ptr{Clrs_dat}, Ptr{Clrs_mp_matrix}, Clong) Pptr m.Q Lin Clrs_true)
    m.P = Pptr[]
    if !found
        # Note that I can have a basis found while the polyhedron is empty
        m.status = :Empty
        # FIXME in that case does redundancy checking with getindex still works ?
    else
        m.status = :AtFirstBasis
        # FIXME does this linearity also works if the first basis is not found ?
        #       I could say that there are linearities which are x_1 = 0 and x_1 = 1
        Q = unsafe_load(m.Q)
        if Q.nredundcol > 0
            # There may have been column redundancy
            # If so the linearity space is obtained and redundant
            # columns are removed. User can access linearity space
            # from lin dimensions nredundcol x d+1

            m.lin = LRSLinearitySpace{N}(Lin[], Q.nredundcol, Q.n, Q.hull == Clrs_true, Q.homogeneous == Clrs_true)
        end
    end
end

function getnextbasis(m::LRSMatrix)
    Pptr = Ref{Ptr{Clrs_dic}}(m.P)
    x = (@lrs_ccall getnextbasis Clong (Ptr{Ptr{Clrs_dic}}, Ptr{Clrs_dat}, Clong) Pptr m.Q Clrs_false)
    found = Clrs_true == x
    m.P = Pptr[]
    m.status = :AtSomeBasis
    found
end

function getsolution(m::LRSMatrix, col::Int)
    Q = unsafe_load(m.Q)
    output = @lrs_ccall alloc_mp_vector Clrs_mp_vector (Clong,) Q.n
    found = Clrs_true == (@lrs_ccall getsolution Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clrs_mp_vector, Clong) m.P m.Q output col)
    if found
        out = convertoutput(output, Q.n, Q.hull == Clrs_true)
    else
        out = nothing
    end
    @lrs_ccall clear_mp_vector Void (Clrs_mp_vector, Clong) output Q.n
    out
end

function checkindex(m::LRSMatrix, index::Int)
    if m.status == :AtNoBasis
        getfirstbasis(m)
    end
    # FIXME if it is at some basis or last basis, does this still works ?
    ret = @lrs_ccall2 checkindex Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong) m.P m.Q index
    m.status = :RedundancyChecked
    [:nonredundant, :redundant, :linearity][ret+1]
end
