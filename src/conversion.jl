export enumtomat, generatorproducer

function generatorproducer(m::RepMatrix)
    Channel(; ctype = Vector{Rational{BigInt}}) do c
        # code from here is borrowed from lrs_main
        if !(m.status in [:AtNoBasis, :AtFirstBasis, :Empty])
            error("I am not at first basis")
        end

        # Pivot to a starting dictionary
        if m.status == :AtNoBasis
            getfirstbasis(m)
        end
        if m.lin !== nothing # FIXME should I do that if m.status is :Empty ?
            L = getmat(m.lin)
            for i in 1:size(L, 1)
                put!(c, L[i, :])
            end
        end

        if m.status != :Empty
            # We initiate reverse search from this dictionary
            # getting new dictionaries until the search is complete
            # User can access each output line from output which is
            # vertex/ray/facet from the lrs_mp_vector output

            while true
                for col in 0:getd(m)
                    output = getsolution(m, col)
                    if output !== nothing
                        put!(c, output)
                    end
                end
                if !getnextbasis(m)
                    break
                end
            end
        end
    end
end

function enumtomat(m::RepMatrix)
    rows = Vector{Vector{Rational{BigInt}}}(undef, 0)
    sizehint!(rows, fulldim(m)) # heuristic size hint
    for row in generatorproducer(m)
        push!(rows, row)
    end
    M = Matrix{Rational{BigInt}}(undef, length(rows), fulldim(m)+1)
    for i in eachindex(rows)
        M[i,:] .= rows[i]
    end
    return M
end

function Base.convert(::Type{LiftedHRepresentation{Rational{BigInt}}}, m::VMatrix)
    linset = getoutputlinset(m)
    A = enumtomat(m)
    LiftedHRepresentation{Rational{BigInt}}(A, linset)
end
HRepresentation(m::VMatrix) = Base.convert(LiftedHRepresentation{Rational{BigInt}}, m)
LiftedHRepresentation(m::VMatrix) = Base.convert(LiftedHRepresentation{Rational{BigInt}}, m)

function Base.convert(::Type{LiftedVRepresentation{Rational{BigInt}}}, m::HMatrix)
    linset = getoutputlinset(m)
    R = enumtomat(m)
    LiftedVRepresentation{Rational{BigInt}}(R, linset)
end
VRepresentation(m::HMatrix) = Base.convert(LiftedVRepresentation{Rational{BigInt}}, m)
LiftedVRepresentation(m::HMatrix) = Base.convert(LiftedVRepresentation{Rational{BigInt}}, m)

function Base.convert(::Type{HMatrix}, m::VMatrix)
    linset = getoutputlinset(m)
    M = enumtomat(m)
    (P, Q) = initmatrix(M, linset, true)
    HMatrix(fulldim(m), P, Q)
end
HMatrix(m::VMatrix) = HMatrix(m)
function Base.convert(::Type{VMatrix}, m::HMatrix)
    linset = getoutputlinset(m)
    M = enumtomat(m)
    (P, Q) = initmatrix(M, linset, false)
    VMatrix(fulldim(m), P, Q)
end
VMatrix(m::HMatrix) = VMatrix(m)
