# It is immutable so that it is stored by value in the structures
# Clrs_dic and Clrs_dat and not by reference
struct GMPInteger
    alloc::Cint
    size::Cint
    data::Ptr{UInt32}
end

const Clrs_true = Clong(1)
const Clrs_false = Clong(0)

const Clrs_mp = GMPInteger
const Clrs_mp_vector = Ptr{Clrs_mp}
const Clrs_mp_matrix = Ptr{Ptr{Clrs_mp}}
primitive type Clrs_fname 32768 end  # 4096 * 8

function extractbigintat(array::Clrs_mp_vector, i::Int)
    tmp = BigInt()
    ccall((:__gmpz_set, :libgmp), Nothing, (Ptr{BigInt}, Ptr{BigInt}), pointer_from_objref(tmp), array + (i-1) * sizeof(GMPInteger))
    tmp
end

function extractbigint(a::Clrs_mp)
    tmp = BigInt()
    ccall((:__gmpz_set, :libgmp), Cvoid, (Ptr{BigInt}, Ref{Clrs_mp}), pointer_from_objref(tmp), a)
    tmp
end

mutable struct Clrs_dic  # dynamic dictionary data
    A::Clrs_mp_matrix
    m::Clong           # A has m+1 rows, row 0 is cost row
    m_A::Clong         # =m or m-d if nonnegative flag set
    d::Clong           # A has d+1 columns, col 0 is b-vector
    d_orig::Clong      # value of d as A was allocated  (E.G.)
    lexflag::Clong     # true if lexmin basis for this vertex
    depth::Clong       # depth of basis/vertex in reverse search tree
    i::Clong           # last pivot row pivot indices
    j::Clong           # last pivot column pivot indices
    det::Clrs_mp       # current determinant of basis
    objnum::Clrs_mp    # objective numerator value
    objden::Clrs_mp    # objective denominator value
    B::Ptr{Clong}      # basis location indices
    Row::Ptr{Clong}    # row location indices
    C::Ptr{Clong}      # cobasis location indices
    Col::Ptr{Clong}    # column location indices
    prev::Ptr{Clrs_dic}
    next::Ptr{Clrs_dic}
end

mutable struct Clrs_dat      # global problem data
    redineq::Clong      # holds indices of rows for redundancy check
    Ain::Clrs_mp_matrix # used only in redund to hold input matrix
    Gcd::Clrs_mp_vector     # Gcd of each row of numerators
    Lcm::Clrs_mp_vector     # Lcm for each row of input denominators
    output::Clrs_mp_vector     # One line of output dimensioned to n

    sumdet::Clrs_mp    # sum of determinants
    Nvolume::Clrs_mp    # volume numerator
    Dvolume::Clrs_mp    # volume denominator
    boundn::Clrs_mp    # objective bound numerator
    boundd::Clrs_mp    # objective bound denominator
    unbounded::Clong    # lp unbounded
    fname::Clrs_fname    # input file name from line 1 of input

    # initially holds order used to find starting
    # basis, default: m,m-1,...,2,1
    facet::Ptr{Clong}    # cobasic indices for restart in needed
    redundcol::Ptr{Clong}    # holds columns which are redundant
    inequality::Ptr{Clong}    # indices of inequalities corr. to cobasic ind
    linearity::Ptr{Clong}    # holds cobasic indices of input linearities
    remain::Ptr{Clong}
    startcob::Ptr{Clong}
    minratio::Ptr{Clong}    # used for lexicographic ratio test
    temparray::Ptr{Clong}    # for sorting indices, dimensioned to d
    isave::Ptr{Clong}
    jsave::Ptr{Clong}  # arrays for estimator, malloc'ed at start
    inputd::Clong    # input dimension: n-1 for H-rep, n for V-rep

    m::Clong          # number of rows in input file
    n::Clong      # number of columns in input file
    lastdv::Clong    # index of last dec. variable after preproc
    # given by inputd-nredundcol
    count0::Clong
    count1::Clong
    count2::Clong
    count3::Clong
    count4::Clong
    count5::Clong
    count6::Clong
    count7::Clong
    count8::Clong
    count9::Clong    # count[0]=rays [1]=verts. [2]=base [3]=pivots
    # count[4]=integer vertices

    startcount0::Clong
    startcount1::Clong
    startcount2::Clong
    startcount3::Clong
    startcount4::Clong

    deepest::Clong    # max depth ever reached in search
    nredundcol::Clong    # number of redundant columns
    nlinearity::Clong    # number of input linearities
    totalnodes::Clong    # count total number of tree nodes evaluated
    runs::Clong      # probes for estimate function
    seed::Clong      # seed for random number generator
    cest0::Clong
    cest1::Clong
    cest2::Clong
    cest3::Clong
    cest4::Clong
    cest5::Clong
    cest6::Clong
    cest7::Clong
    cest8::Clong
    cest9::Clong    # ests: 0=rays,1=vert,2=bases,3=vol,4=int vert
    nextineq::Clong
    #*** flags  **********
    allbases::Clong    # TRUE if all bases should be printed
    bound::Clong                 # TRUE if upper/lower bound on objective given
    countonly::Clong             # TRUE if only count totals should be output
    debug::Clong
    dualdeg::Clong    # TRUE if start dictionary is dual degenerate
    etrace::Clong    # turn off debug at basis # strace
    extract::Clong
    frequency::Clong    # frequency to print cobasis indices
    geometric::Clong    # TRUE if incident vertex prints after each ray
    getvolume::Clong    # do volume calculation
    givenstart::Clong    # TRUE if a starting cobasis is given
    giveoutput::Clong
    verifyredund::Clong
    noredundcheck::Clong
    homogeneous::Clong    # TRUE if all entries in column one are zero
    hull::Clong      # do convex hull computation if TRUE
    incidence::Clong             # print all tight inequalities (vertices/rays)
    lponly::Clong    # true if only lp solution wanted
    maxdepth::Clonglong    # max depth to search to in treee
    maximize::Clong    # flag for LP maximization
    maxoutput::Clong       # if positive, maximum number of output lines
    maxcobases::Clong       # if positive, after maxcobasis unexplored subtrees reported
    messages::Clong
    minimize::Clong    # flag for LP minimization
    mindepth::Clonglong    # do not backtrack above mindepth
    mplrs::Clong
    nash::Clong                  # TRUE for computing nash equilibria
    nonnegative::Clong    # TRUE if last d constraints are nonnegativity
    polytope::Clong    # TRUE for facet computation of a polytope
    printcobasis::Clong    # TRUE if all cobasis should be printed
    printslack::Clong    # TRUE if indices of slack inequal. printed
    truncate::Clong              # TRUE: truncate tree when moving from opt vert
    verbose::Clong               # FALSE for minimalist output
    restart::Clong    # TRUE if restarting from some cobasis
    strace::Clong    # turn on  debug at basis # strace
    voronoi::Clong    # compute voronoi vertices by transformation
    subtreesize::Clonglong   # in estimate mode, iterates if cob_est >= subtreesize
    triangulation::Clong     # TRUE: the cobases printed triangulate the polytope
    newstart::Clong          # TRUE: lrs is restarted with new arithmetic

    # Variables for saving/restoring cobasis,  db

    id::Clong      # numbered sequentially
    name::Cstring      # passed by user

    saved_count0::Clong
    saved_count1::Clong
    saved_count2::Clong
    saved_count3::Clong
    saved_count4::Clong  # How often to print out current cobasis
    saved_C::Ptr{Clong}
    saved_det::Clrs_mp
    saved_sumdet::Clrs_mp
    saved_depth::Clong
    saved_d::Clong

    saved_flag::Clong    # There is something in the saved cobasis

    # Variables for cacheing dictionaries, db
    Qhead::Ptr{Clrs_dic}
    Qtail::Ptr{Clrs_dic}
end


# Utility function

function unsafe_field_store!(p::Ptr{T}, name::Symbol, x, err::Bool=true) where T
    i = Base.fieldindex(T, name, err)
    offset = fieldoffset(T, i)
    unsafe_store!(convert(Ptr{UInt}, p + offset), x)
    return p
end

