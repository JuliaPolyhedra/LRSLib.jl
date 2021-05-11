module LRSLib

using Polyhedra
using LinearAlgebra
using lrslib_jll
using Markdown

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
    if Clrs_false == (@lrs_ccall init Clong (Ptr{Cchar},) C_NULL)
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
