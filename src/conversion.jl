export enumtomat, generatorproducer
function generatorproducer(m::LRSMatrix)
    Channel() do c
        # code from here is borrowed from lrs_main
        if !(m.status in [:AtNoBasis, :AtFirstBasis, :Empty])
            error("I am not at first basis")
        end

        # Pivot to a starting dictionary
        if m.status == :AtNoBasis
            getfirstbasis(m)
        end
        if !isnull(m.lin) # FIXME should I do that if m.status is :Empty ?
            L = getmat(get(m.lin))
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
function enumtomat{N}(m::LRSMatrix{N})
    M = Matrix{Rational{BigInt}}(0, N+1)
    for output in generatorproducer(m)
        M = [M; output']
    end
    M
end

function Base.convert{N}(::Type{LiftedHRepresentation{N, Rational{BigInt}}}, m::LRSGeneratorMatrix{N})
    linset = getoutputlinset(m)
    A = enumtomat(m)
    LiftedHRepresentation{N, Rational{BigInt}}(A, linset)
end
HRepresentation{N}(m::LRSGeneratorMatrix{N}) = Base.convert(LiftedHRepresentation{N, Rational{BigInt}}, m)
LiftedHRepresentation{N}(m::LRSGeneratorMatrix{N}) = Base.convert(LiftedHRepresentation{N, Rational{BigInt}}, m)

function Base.convert{N}(::Type{LiftedVRepresentation{N, Rational{BigInt}}}, m::LRSInequalityMatrix{N})
    linset = getoutputlinset(m)
    R = enumtomat(m)
    LiftedVRepresentation{N, Rational{BigInt}}(R, linset)
end
VRepresentation{N}(m::LRSInequalityMatrix{N}) = Base.convert(LiftedVRepresentation{N, Rational{BigInt}}, m)
LiftedVRepresentation{N}(m::LRSInequalityMatrix{N}) = Base.convert(LiftedVRepresentation{N, Rational{BigInt}}, m)

function Base.convert{N}(::Type{LRSInequalityMatrix{N}}, m::LRSGeneratorMatrix{N})
    linset = getoutputlinset(m)
    M = enumtomat(m)
    (P, Q) = initmatrix(M, linset, true)
    LRSInequalityMatrix{N}(P, Q)
end
LRSInequalityMatrix{N}(m::LRSGeneratorMatrix{N}) = LRSInequalityMatrix{N}(m)
function Base.convert{N}(::Type{LRSGeneratorMatrix{N}}, m::LRSInequalityMatrix{N})
    linset = getoutputlinset(m)
    M = enumtomat(m)
    (P, Q) = initmatrix(M, linset, false)
    LRSGeneratorMatrix{N}(P, Q)
end
LRSGeneratorMatrix{N}(m::LRSInequalityMatrix{N}) = LRSGeneratorMatrix{N}(m)
