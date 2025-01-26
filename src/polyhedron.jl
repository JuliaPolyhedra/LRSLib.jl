struct Library <: Polyhedra.Library
    solver
    function Library(solver=nothing)
        new(solver)
    end
end
Polyhedra.similar_library(::Library, ::Polyhedra.FullDim, ::Type{T}) where T<:Union{Integer,Rational} = Library()
Polyhedra.similar_library(::Library, d::Polyhedra.FullDim, ::Type{T}) where T = Polyhedra.default_library(d, T)

mutable struct Polyhedron <: Polyhedra.Polyhedron{Rational{BigInt}}
    ine::Union{Nothing, LiftedHRepresentation{Rational{BigInt}, Matrix{Rational{BigInt}}}}
    inem::Union{Nothing, HMatrix}
    ext::Union{Nothing, LiftedVRepresentation{Rational{BigInt}, Matrix{Rational{BigInt}}}}
    extm::Union{Nothing, VMatrix}
    hlinearitydetected::Bool
    vlinearitydetected::Bool
    noredundantinequality::Bool
    noredundantgenerator::Bool
    solver

    function Polyhedron(ine::HRepresentation{Rational{BigInt}}, ext::VRepresentation{Rational{BigInt}}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool, solver)
        new(ine, nothing, ext, nothing, hld, vld, nri, nrg, solver)
    end
    function Polyhedron(ine::HRepresentation{Rational{BigInt}}, ::Nothing, hld::Bool, vld::Bool, nri::Bool, nrg::Bool, solver)
        new(ine, nothing, nothing, nothing, hld, vld, nri, nrg, solver)
    end
    function Polyhedron(::Nothing, ext::VRepresentation{Rational{BigInt}}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool, solver)
        new(nothing, nothing, ext, nothing, hld, vld, nri, nrg, solver)
    end
    function Polyhedron(ine::HRepresentation{Rational{BigInt}}, solver)
        new(ine, nothing, nothing, nothing, false, false, false, false, solver)
    end
    function Polyhedron(ext::VRepresentation{Rational{BigInt}}, solver)
        new(nothing, nothing, ext, nothing, false, false, false, false, solver)
    end
end
Polyhedron(h::HRepresentation, solver) = Polyhedron(HRepresentation{Rational{BigInt}}(h), solver)
Polyhedron(v::VRepresentation, solver) = Polyhedron(VRepresentation{Rational{BigInt}}(v), solver)

Polyhedra.FullDim(p::Polyhedron) = Polyhedra.FullDim_rep(p.ine, p.inem, p.ext, p.extm)
Polyhedra.library(::Polyhedron) = Library()
Polyhedra.default_solver(p::Polyhedron; T=nothing) = p.solver
Polyhedra.supportssolver(::Type{<:Polyhedron}) = true

Polyhedra.hvectortype(::Union{Polyhedron, Type{Polyhedron}}) = Polyhedra.hvectortype(LiftedHRepresentation{Rational{BigInt}, Matrix{Rational{BigInt}}})
Polyhedra.vvectortype(::Union{Polyhedron, Type{Polyhedron}}) = Polyhedra.vvectortype(LiftedVRepresentation{Rational{BigInt}, Matrix{Rational{BigInt}}})
Polyhedra.similar_type(::Type{<:Polyhedron}, ::Polyhedra.FullDim, ::Type{Rational{BigInt}}) = Polyhedron
Polyhedra.similar_type(::Type{<:Polyhedron}, d::Polyhedra.FullDim, ::Type{T}) where T = Polyhedra.default_type(d, T)

# Helpers
function getine(p::Polyhedron)
    if p.ine === nothing
        if p.inem !== nothing && checkfreshness(p.inem, :Fresh)
            p.ine = p.inem
        else
            p.ine = LiftedHRepresentation(getextm(p, :Fresh))
            p.inem = nothing
            p.hlinearitydetected = true
        end
    end
    return p.ine
end
function getinem(p::Polyhedron, fresh::Symbol=:AnyFreshness)
    if p.inem === nothing || !checkfreshness(p.inem, fresh)
        p.inem = RepMatrix(getine(p))
    end
    return p.inem
end
function getext(p::Polyhedron)
    if p.ext === nothing
        if p.extm !== nothing && checkfreshness(p.extm, :Fresh)
            p.ext = p.extm
        else
            p.ext = LiftedVRepresentation(getinem(p, :Fresh))
            p.extm = nothing
            p.vlinearitydetected = true
        end
    end
    return p.ext
end
function getextm(p::Polyhedron, fresh::Symbol=:AnyFreshness)
    if p.extm === nothing || !checkfreshness(p.extm, fresh)
        p.extm = RepMatrix(getext(p))
    end
    return p.extm
end

function clearfield!(p::Polyhedron)
    p.ine = nothing
    p.inem = nothing
    p.ext = nothing
    p.extm = nothing
    hlinearitydetected = false
    vlinearitydetected = false
    noredundantinequality = false
    noredundantgenerator = false
end

function Polyhedra.sethrep!(p::Polyhedron, h::HRepresentation{Rational{BigInt}})
    p.ine = h
end
function Polyhedra.setvrep!(p::Polyhedron, v::VRepresentation{Rational{BigInt}})
    p.ext = v
end
function Polyhedra.resethrep!(p::Polyhedron, h::HRepresentation{Rational{BigInt}})
    clearfield!(p)
    p.ine = h
end
function Polyhedra.resetvrep!(p::Polyhedron, v::VRepresentation{Rational{BigInt}})
    clearfield!(p)
    p.ext = v
end


# Implementation of Polyhedron's mandatory interface
Polyhedra.polyhedron(rep::Representation, lib::Library) = Polyhedron(rep, lib.solver)

function Polyhedron(d::Polyhedra.FullDim, hits::Polyhedra.HIt...; solver=nothing)
    Polyhedron(HMatrix(d, hits...), solver)
end
function Polyhedron(d::Polyhedra.FullDim, vits::Polyhedra.VIt...; solver=nothing)
    Polyhedron(VMatrix(d, vits...), solver)
end

function Base.copy(p::Polyhedron)
    ine = nothing
    if p.ine !== nothing
        ine = copy(p.ine)
    end
    ext = nothing
    if p.ext !== nothing
        ext = copy(p.ext)
    end
    Polyhedron(ine, ext, p.hlinearitydetected, p.vlinearitydetected, p.noredundantinequality, p.noredundantgenerator, p.solver)
end
function Polyhedra.hrepiscomputed(p::Polyhedron)
    p.ine !== nothing
end
function Polyhedra.hrep(p::Polyhedron)
    getine(p)
end
function Polyhedra.vrepiscomputed(p::Polyhedron)
    p.ext !== nothing
end
function Polyhedra.vrep(p::Polyhedron)
    getext(p)
end
#eliminate(p::Polyhedron, delset::BitSet)                     = error("not implemented")
function Polyhedra.detecthlinearity!(p::Polyhedron)
    if !p.hlinearitydetected
        getext(p)
        p.inem = nothing
        p.ine = nothing
        getine(p)
        # getine sets hlinearity as detected
    end
end
function Polyhedra.detectvlinearity!(p::Polyhedron)
    if !p.vlinearitydetected
        getine(p)
        p.extm = nothing
        p.ext = nothing
        getext(p)
        # getext sets vlinearity as detected
    end
end
function Polyhedra.removehredundancy!(p::Polyhedron)
    #if !p.noredundantinequality
    ine = getine(p)
    inem = getinem(p, :AlmostFresh) # FIXME does it need to be fresh ?
    linset = getinputlinsubset(inem)
    redset = redund(inem)
    nonred = setdiff(BitSet(1:size(ine.A, 1)), redset)
    nonred = collect(setdiff(nonred, linset))
    lin = collect(linset)
    ine.A = [ine.A[lin,:]; ine.A[nonred,:]]
    ine.linset = BitSet(1:length(linset))
    p.noredundantinequality = true
    #end
end
function Polyhedra.removevredundancy!(p::Polyhedron)
    if !p.noredundantgenerator
        detectvlinearity!(p)
        ext = getext(p)
        extm = getextm(p, :AlmostFresh) # FIXME does it need to be fresh ?
        redset = BitSet(redund(extm) .+ 1)
        nonred = setdiff(BitSet(1:size(ext.R, 1)), redset)
        nonred = collect(setdiff(nonred, ext.linset))
        lin = collect(ext.linset)
        ext.R = [ext.R[lin,:]; ext.R[nonred,:]]
        ext.linset = BitSet(1:length(ext.linset))
        p.noredundantgenerator = true
    end
end
#function getredundantinequalities(p::Polyhedron)
#  redund(getinem(p, :AlmostFresh))
#end
_getrepfor(p::Polyhedron, ::Polyhedra.HIndex, status::Symbol) = getinem(p, status)
_getrepfor(p::Polyhedron, ::Polyhedra.VIndex, status::Symbol) = getextm(p, status)
function Polyhedra.isredundant(p::Polyhedron, idx::Polyhedra.Index; strongly=false, cert=false, solver=nothing)
    @assert !strongly && !cert
    Polyhedra.isredundant(_getrepfor(p, idx, :AlmostFresh), idx) # FIXME does it need to be fresh ?
end
# Optional interface
function Polyhedra.loadpolyhedron!(p::Polyhedron, filename::AbstractString, ::Type{Val{:ext}})
    clearfield!(p)
    p.extm = VMatrix(string(filename, ".ext"))
end
