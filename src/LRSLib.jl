module LRSLib

using BinDeps
using Polyhedra
using LinearAlgebra
using Markdown

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
    else
    error("LRSLib not properly installed. Please run Pkg.build(\"LRSLib\")")
end

macro lrs_ccall(f, args...)
    quote
        ret = ccall(($"lrs_$(f)_gmp", liblrs), $(map(esc,args)...))
        ret
    end
end
macro lrs_ccall2(f, args...)
    quote
        ret = ccall(($"$(f)_gmp", liblrs), $(map(esc,args)...))
        ret
    end
end


include("lrstypes.jl")

function __init__()
    # lrslib segfault if the name is `C_NULL`.
    if Clrs_false == (@lrs_ccall init Clong (Ptr{Cchar},) "LRSLib Julia wrapper")
        error("Initialization of LRS failed")
    end
end

include("matrix.jl")
include("conversion.jl")
include("redund.jl")
include("polyhedron.jl")
include("lp.jl")
include("nash.jl")

end # module
