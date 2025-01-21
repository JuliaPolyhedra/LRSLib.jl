function lrs_alloc_dat()
    @lrs_ccall alloc_dat Ptr{Clrs_dat} (Ptr{Cchar},) NONAME
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

function _init_dat(Q::Ptr{Clrs_dat}, m, n, hull)
    unsafe_field_store!(Q, :m, m)
    unsafe_field_store!(Q, :n, n)
    if hull
        unsafe_field_store!(Q, :hull, Clrs_true)
        unsafe_field_store!(Q, :polytope, Clrs_true)
    else
        unsafe_field_store!(Q, :hull, Clrs_false)
    end
end

# FIXME still needed ?
function initmatrix(M::Matrix{Rational{BigInt}}, linset, Hrep::Bool)
    m = Clong(size(M, 1))
    n = Clong(size(M, 2))
    Q = lrs_alloc_dat()
    _init_dat(Q, m, n, !Hrep)
    P = lrs_alloc_dic(Q)
    #Q->getvolume= TRUE; # compute the volume # TODO cheap do it
    for i in 1:m
        #   num = map(x -> GMPInteger(x.num.alloc, x.num.size, x.num.d), M[i,:])
        #   den = map(x -> GMPInteger(x.den.alloc, x.den.size, x.den.d), M[i,:])
        ineq = !(i in linset)
        setrow(P, Q, i, M[i,:], !(i in linset))
        #   @lrs_ccall set_row_mp Nothing (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong, Clrs_mp_vector, Clrs_mp_vector, Clong) P Q i num den Clong(ineq)
    end
    # This is the objective. If I have no objective LRS might fail
    if Hrep
        setrow(P, Q, 0, ones(Rational{BigInt}, n), true)
    end
    (P, Q)
end

function unsafe_load_inequality(m::Matrix, idx)
    Q = unsafe_load(m.Q)
    return unsafe_load(Q.inequality, idx - Q.lastdv + 1)
end

function fillmatrix(inequality::Bool, P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}, itr::Polyhedra.ElemIt, offset::Int)
    for (i, item) in enumerate(itr)
        a = convert(Vector{Rational{BigInt}}, vec(coord(lift(item))))
        setrow(P, Q, offset+i, inequality ? -a : a, !islin(item))
    end
end

# If ElemIt contains AbstractVector{T}, we cannot infer FullDim so we add it as argument
function initmatrix(d::Polyhedra.FullDim, inequality::Bool, itr::Polyhedra.ElemIt...)
    n = fulldim(d)+1
    cs = cumsum(collect(length.(itr))) # cumsum not defined for tuple :(
    m = cs[end]
    offset = [0; cs[1:end-1]]
    Q = lrs_alloc_dat()
    _init_dat(Q, m, n, !inequality)
    P = lrs_alloc_dic(Q)
    #Q->getvolume= TRUE; # compute the volume # TODO cheap do it
    fillmatrix.(inequality, P, Q, itr, offset)
    # This is the objective. If I have no objective LRS might fail
    if inequality
        setrow(P, Q, 0, ones(Rational{BigInt}, n), true)
    end
    (P, Q)
end


# Representation

mutable struct LRSLinearitySpace
    Lin::Clrs_mp_matrix
    nlin::Int
    n::Int
    hull::Bool
    homogeneous::Bool

    function LRSLinearitySpace(Lin::Clrs_mp_matrix, nlin, n, hull, homogeneous)
        m = new(Lin, nlin, n, hull, homogeneous)
        finalizer(myfree, m)
        m
    end
end

mutable struct HMatrix <: Polyhedra.MixedHRep{Rational{BigInt}}
    N::Int
    P::Ptr{Clrs_dic}
    Q::Ptr{Clrs_dat}
    status::Symbol
    lin::Union{Nothing, LRSLinearitySpace}
    function HMatrix(N::Integer, P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat})
        @assert N == getd(P)
        m = new(N, P, Q, :AtNoBasis, nothing)
        finalizer(myfree, m)
        m
    end
end
Polyhedra.similar_type(::Type{<:HMatrix}, ::Polyhedra.FullDim, ::Type{Rational{BigInt}}) = HMatrix

mutable struct VMatrix <: Polyhedra.MixedVRep{Rational{BigInt}}
    N::Int
    P::Ptr{Clrs_dic}
    Q::Ptr{Clrs_dat}
    status::Symbol
    lin::Union{Nothing, LRSLinearitySpace}
    cone::Bool # If true, LRS will not return any point so we need to add the origin
    function VMatrix(N::Integer, P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat})
        @assert N + 1 == getd(P)
        m = _length(P)
        cone = !iszero(m) # If there is no ray and no point, it is empty so we should not add the origin
        for i in 1:m
            if isrowpoint(P, Q, i)
                cone = false
                break
            end
        end
        m = new(N, P, Q, :AtNoBasis, nothing, cone)
        finalizer(myfree, m)
        m
    end
end
Polyhedra.similar_type(::Type{<:VMatrix}, ::Polyhedra.FullDim, ::Type{Rational{BigInt}}) = VMatrix

const RepMatrix = Union{HMatrix, VMatrix}
Polyhedra.FullDim(m::RepMatrix) = m.N
Polyhedra.hvectortype(::Union{HMatrix, Type{<:HMatrix}}) = Vector{Rational{BigInt}}
Polyhedra.vvectortype(::Union{VMatrix, Type{<:VMatrix}}) = Vector{Rational{BigInt}}
Polyhedra.coefficient_type(::Union{RepMatrix, Type{<:RepMatrix}}) = Rational{BigInt}

function linset(matrix::RepMatrix)
    extractinputlinset(unsafe_load(matrix.Q))
end
_length(P::Ptr{Clrs_dic}) = unsafe_load(P).m
Base.length(matrix::HMatrix) = _length(matrix.P)
Base.length(matrix::VMatrix) = _length(matrix.P) + matrix.cone

RepMatrix(hrep::HRepresentation) = convert(HMatrix, hrep)
RepMatrix(vrep::VRepresentation) = convert(VMatrix, vrep)

function checkfreshness(m::RepMatrix, fresh::Symbol)
    fresh == :AnyFreshNess ||
    (fresh == :Fresh && m.status in [:AtNoBasis, :AtFirstBasis, :Empty]) ||
    (fresh == :AlmostFresh && m.status in [:AtNoBasis, :AtFirstBasis, :Empty, :RedundancyChecked])
end

function myfree(l::LRSLinearitySpace)
    if l.nlin > 0
        @lrs_ccall clear_mp_matrix Nothing (Clrs_mp_matrix, Clong, Clong) l.Lin l.nlin l.n
    end
end

function myfree(m::RepMatrix)
    @lrs_ccall free_all_memory Nothing (Ptr{Clrs_dic}, Ptr{Clrs_dat}) m.P m.Q
end

_islin(rep::RepMatrix, idx::Polyhedra.Index) = isininputlinset(unsafe_load(rep.Q), idx.value)
Polyhedra.islin(hrep::HMatrix, idx::Polyhedra.HIndex{Rational{BigInt}}) = _islin(hrep, idx)
Polyhedra.islin(vrep::VMatrix, idx::Polyhedra.VIndex{Rational{BigInt}}) = !(vrep.cone && idx.value == nvreps(vrep)) && _islin(vrep, idx)

Polyhedra.done(idxs::Polyhedra.Indices{Rational{BigInt}, ElemT, <:RepMatrix}, idx::Polyhedra.Index{Rational{BigInt}, ElemT}) where {ElemT} = idx.value > length(idxs.rep)
Base.get(rep::RepMatrix, idx::Polyhedra.Index{Rational{BigInt}}) = Polyhedra.valuetype(idx)(extractrow(rep, idx.value)...)

# H-representation

#HMatrix{T}(rep::Rep{T}) = HMatrix{polytypefor(T), mytypefor(T)}(rep) # TODO finish this line

function HMatrix(filename::AbstractString)
    P, Q = initmatrix(filename)
    HMatrix(unsafe_load(P).d, P, Q)
end

Base.copy(ine::HMatrix) = HMatrix(Polyhedra.hreps(ine)...)

function HMatrix(d::Polyhedra.FullDim, hits::Polyhedra.HIt...)
    P, Q = initmatrix(d, true, hits...)
    HMatrix(fulldim(d), P, Q)
end

nhreps(matrix::HMatrix) = length(matrix)
neqs(matrix::HMatrix) = unsafe_load(matrix.Q).nlinearity
Base.length(idxs::Polyhedra.Indices{Rational{BigInt}, <:HyperPlane{Rational{BigInt}}, <:HMatrix}) = neqs(idxs.rep)
Base.length(idxs::Polyhedra.Indices{Rational{BigInt}, <:HalfSpace{Rational{BigInt}}, <:HMatrix}) = nhreps(idxs.rep) - neqs(idxs.rep)

function Base.isvalid(hrep::HMatrix, idx::Polyhedra.HIndex{Rational{BigInt}})
    0 < idx.value <= nhreps(hrep) && Polyhedra.islin(hrep, idx) == islin(idx)
end

# V-representation

function VMatrix(filename::AbstractString)
    d, P, Q = initmatrix(filename)
    VMatrix(unsafe_load(P).d-1, P, Q)
end

Base.copy(ext::VMatrix) = VMatrix(Polyhedra.vreps(ext)...)

function VMatrix(d::Polyhedra.FullDim, vits::Polyhedra.VIt...)
    P, Q = initmatrix(d, false, vits...)
    VMatrix(fulldim(d), P, Q)
end

nvreps(ext::VMatrix) = length(ext)
function Base.length(idxs::Polyhedra.PointIndices{Rational{BigInt}, <:VMatrix})
    if idxs.rep.cone
        1
    else
        Polyhedra.mixedlength(idxs)
    end
end

function Base.isvalid(vrep::VMatrix, idx::Polyhedra.VIndex{Rational{BigInt}})
    isp = isrowpoint(vrep, idx.value)
    isl = Polyhedra.islin(vrep, idx)
    0 < idx.value <= length(vrep) && isl == islin(idx) && isp == ispoint(idx)
end

getd(P::Ptr{Clrs_dic}) = unsafe_load(P).d

# `P.d` may be modified, e.g. by `getfirstbasis` so it may be different to `m.N`
getd(m::HMatrix) = getd(m.P)
# For the same reason, it might be different to `m.N + 1`
getd(m::VMatrix) = getd(m.P)

# FIXME the difference between getd and `m.N (+ 1)` is for removing linearity.
#       Wouldn't it work if I do `-1` when hull && homogeneous ?

function setrow(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}, i::Int, row::Vector{Rational{BigInt}}, ineq::Bool)
    num = map(x -> GMPInteger(x.num.alloc, x.num.size, x.num.d), row)
    den = map(x -> GMPInteger(x.den.alloc, x.den.size, x.den.d), row)
    @lrs_ccall set_row_mp Nothing (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong, Clrs_mp_vector, Clrs_mp_vector, Clong) P Q i num den Clong(ineq)
end
setrow(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}, i::Int, row::AbstractVector{Rational{BigInt}}, ineq::Bool) = setrow(P, Q, i, collect(row), ineq) # e.g. for sparse a

function setdebug(m::RepMatrix, debug::Bool)
    unsafe_field_store!(m.Q, :debug, (debug ? Clrs_true : Clrs_false))
end

function setverbose(m::RepMatrix, verbose::Bool)
    unsafe_field_store!(m.Q, :verbose, (verbose ? Clrs_true : Clrs_false))
end

function printA(m::RepMatrix)
    @lrs_ccall2 printA Cvoid (Ptr{Clrs_dic}, Ptr{Clrs_dat}) m.P m.Q
end

function isrowpoint(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}, i)
    offset = unsafe_load(Q).homogeneous == Clrs_true ? 0 : 1
    row = unsafe_load(unsafe_load(P).A, 1+i)
    !iszero(extractbigintat(row, offset+1))
end
function isrowpoint(matrix::VMatrix, i::Int)
    (matrix.cone && i == nvreps(matrix)) || isrowpoint(matrix.P, matrix.Q, i)
end

function extractrow(P::Clrs_dic, Q::Clrs_dat, N, i, offset)
    #d = Q.n-offset-1 # FIXME when it is modified...
    a = Vector{Rational{BigInt}}(undef, N+1)
    gcd = extractbigintat(Q.Gcd, 1+i) # first row is the objective
    lcm = extractbigintat(Q.Lcm, 1+i)
    row = unsafe_load(P.A, 1+i)
    extractthisrow(i::Int) = (extractbigintat(row, offset+i) * gcd) // lcm
    for j in 1:N+1
        a[j] = extractthisrow(j)
    end
    a
end

function warn_fresh(m::RepMatrix)
    if !checkfreshness(m, :Fresh)
        warn("Extracting the rows of an LRS matrix after it has been used for representation conversion does not give the correct elements of the polyhedron it represents")
    end
end

function extractrow(matrix::HMatrix, i::Int)
    P = unsafe_load(matrix.P)
    Q = unsafe_load(matrix.Q)
    b = extractrow(P, Q, getd(matrix), i, 0)
    β = b[1]
    a = -b[2:end]
    a, β
end

function extractrow(matrix::VMatrix, i::Int)
    if matrix.cone && i == nvreps(matrix)
        a = Polyhedra.origin(Polyhedra.arraytype(matrix), FullDim())
    else
        P = unsafe_load(matrix.P)
        Q = unsafe_load(matrix.Q)
        #d = Q.n-offset-1 # FIXME when it is modified...
        @assert Q.hull == Clrs_true
        offset = Q.homogeneous == Clrs_true ? 0 : 1
        b = extractrow(P, Q, getd(matrix)-1, i, offset)
        a = b[2:end]
    end
    (a,) # Needs to be a tuple, see Base.get(::RepMatrix, ...)
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
    linset = BitSet([])
    for i in 1:Q.nlinearity
        push!(linset, unsafe_load(Q.linearity, i))
    end
    linset
end

function extractoutputlinset(Q::Clrs_dat)
    k = (Q.hull == Clrs_true && Q.homogeneous == Clrs_true) ? 1 : 0
    nredundcol = Q.nredundcol
    BitSet(1:(nredundcol-k))
end

# FIXME The only think that is done is that the linearities given that were redundant have been
# removed from the linearity set so that they can be marked as redundant inequalities.
# New linearities are detected but getinputlinsubset does not give them.
# I should check in redundcols
function getinputlinsubset(m::RepMatrix)
    if m.status == :AtNoBasis
        getfirstbasis(m)
    end
    linset(m)
end
function getoutputlinset(m::RepMatrix)
    if m.status == :AtNoBasis
        getfirstbasis(m)
    end
    extractoutputlinset(unsafe_load(m.Q))
end


function convertoutput(x::Clrs_mp_vector, n, hull)
    first = extractbigintat(x, 1)
    rest = Vector{BigInt}(undef, n-1)
    for i = 2:n
        rest[i-1] = extractbigintat(x, i)
    end
    if hull || first == 0
        Rational{BigInt}[first; rest]
    else
        [one(Rational{BigInt}); rest // first]
    end
end

function getmat(lin::LRSLinearitySpace)
    startcol = lin.hull && lin.homogeneous ? 2 : 1 # col zero not treated as redundant
    A = Matrix{BigInt}(undef, lin.nlin-startcol+1, lin.n)
    for col in startcol:lin.nlin # print linearity space */
        A[col-startcol+1,:] = convertoutput(unsafe_load(lin.Lin, col), lin.n, lin.hull)
    end
    A
end

function getfirstbasis(m::RepMatrix)
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

            m.lin = LRSLinearitySpace(Lin[], Q.nredundcol, Q.n, Q.hull == Clrs_true, Q.homogeneous == Clrs_true)
        end
    end
end

function getnextbasis(m::RepMatrix)
    Pptr = Ref{Ptr{Clrs_dic}}(m.P)
    x = (@lrs_ccall getnextbasis Clong (Ptr{Ptr{Clrs_dic}}, Ptr{Clrs_dat}, Clong) Pptr m.Q Clrs_false)
    found = Clrs_true == x
    m.P = Pptr[]
    m.status = :AtSomeBasis
    found
end

function getsolution(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat}, col::Int)
    n = unsafe_load(Q).n
    output = @lrs_ccall alloc_mp_vector Clrs_mp_vector (Clong,) n
    found = Clrs_true == (@lrs_ccall getsolution Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clrs_mp_vector, Clong) P Q output col)
    if found
        out = convertoutput(output, n, unsafe_load(Q).hull == Clrs_true)
    else
        out = nothing
    end
    @lrs_ccall clear_mp_vector Nothing (Clrs_mp_vector, Clong) output n
    out
end

getsolution(m::RepMatrix, col::Int) = getsolution(m.P, m.Q, col)

function checkindex(m::RepMatrix, index::Int)
    if m.status == :AtNoBasis
        getfirstbasis(m)
    end
    # FIXME if it is at some basis or last basis, does this still works ?
    ret = @lrs_ccall2 checkindex Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clong) m.P m.Q index
    m.status = :RedundancyChecked
    [:nonredundant, :redundant, :linearity][ret+1]
end
