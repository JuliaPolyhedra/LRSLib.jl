#=
Interface to lrs_solve_lp

=#


mutable struct LRSLinprogSolution
    status::Symbol  # :Optimal, :Infeasible, :Unbounded, :UserLimit, or :Error
    objval::Union{Nothing,Rational{BigInt}}
    sol::Vector{Rational{BigInt}}
    attrs::Dict
end


function setobj(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat},
                c::Vector{<:Union{Rational,Integer}}, maximize::Bool)
    field = maximize ? :maximize : :minimize
    row = Vector{Rational{BigInt}}(undef, unsafe_load(P).d+1)
    row[1] = zero(eltype(row))
    row[2:end] = maximize ? c : -c
    unsafe_field_store!(Q, field, Clrs_true)
    setrow(P, Q, 0, row, true)
end

setobj(m::HMatrix, c::Vector{<:Real}, maximize::Bool) =
    setobj(m.P, m.Q, c, maximize)

@doc """
    setobj(m::HMatrix, c::Vector{<:Union{Rational,Integer}}, maximize::Bool)

Set the objective function vector `c` for `m`. Set `maximize` to true (false,
resp.) if the objective function is to be maximized (minimized, resp.).
""" setobj


# ccall to `lrs_solve_lp`
function lrs_solve_lp(P::Ptr{Clrs_dic}, Q::Ptr{Clrs_dat})
    found = (
        Clrs_true ==
        @lrs_ccall solve_lp Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}) P Q
    )
    found
end

objnum(m::HMatrix) = extractbigint(unsafe_load(m.P).objnum)
objden(m::HMatrix) = extractbigint(unsafe_load(m.P).objden)

@doc doc"""
    lpsolve(m::HMatrix)

Solve the linear program given by `m`. Return an instance of the type:

    mutable struct LRSLinprogSolution
        status::Symbol
        objval::Rational{BigInt}
        sol::Vector{Rational{BigInt}}
        attrs::Dict
    end

# Examples

Suppose we want to solve the problem,
``\max x_1`` subject to ``4 x_1 + 2 x_2 \leq 3, x_1 \geq 0, x_2 \geq 0``:

```julia-repl
julia> using LRSLib: RepMatrix, setobj, lpsolve

julia> using Polyhedra

julia> c = [1, 0];

julia> A = [4 2; -1 0; 0 -1];

julia> b = [3, 0, 0];

julia> maximize = true;

julia> hr = hrep(A, b);

julia> m = RepMatrix(hr);

julia> setobj(m, c, maximize);  # Set the objective

julia> sol = lpsolve(m);  # Solve the problem
```

Status:

```julia-repl
julia> sol.status
:Optimal
```

Optimal value:

```julia-repl
julia> sol.objval
3//4
```

Optimal solution:

```julia-repl
julia> sol.sol
2-element Array{Rational{BigInt},1}:
 3//4
 0//1
```
"""
function lpsolve(m::HMatrix)
    found = lrs_solve_lp(m.P, m.Q)
    unbounded = Bool(unsafe_load(m.Q).unbounded)
    attr = Dict()

    if unbounded || !found
        status = unbounded ? :Unbounded : :Infeasible
        return LRSLinprogSolution(status, nothing, Rational{BigInt}[], attr)
    end

    # found
    status = :Optimal
    objval = objnum(m) // objden(m)
    sol = Vector{Rational{BigInt}}(undef, unsafe_load(m.P).d_orig)
    out = getsolution(m, 0)
    copyto!(sol, out[2:end])
    return LRSLinprogSolution(status, objval, sol, attr)
end


# lrs_lpoutput from lrslib.c
# For debugging
function print_lpoutput(m::HMatrix, print_solution::Bool=false)
    P, Q = unsafe_load(m.P), unsafe_load(m.Q)
    output = @lrs_ccall alloc_mp_vector Clrs_mp_vector (Clong,) Q.n

    for col in 0:(Q.n)-1
        found = (
            Clrs_true ==
            @lrs_ccall getsolution Clong (Ptr{Clrs_dic}, Ptr{Clrs_dat}, Clrs_mp_vector, Clong) m.P m.Q output col
        )
        if print_solution
            if found
                println("col=$(col): found")
                @lrs_ccall printoutput Cvoid (Ptr{Clrs_dat}, Clrs_mp_vector) m.Q output
                println()
            else
                println("col=$(col): not found")
            end
        end
    end

    num = extractbigint(P.objnum)
    den = extractbigint(P.objden)
    println("*Objective function has value $(num)/$(den)")

    print("*Primal:")
    for i in 1:Q.n-1
        print(" x_$i=")
        num = extractbigint(unsafe_load(output, i+1))
        den = extractbigint(unsafe_load(output, 1))
        print(num, "/", den)
    end
    println()

    @lrs_ccall clear_mp_vector Nothing (Clrs_mp_vector, Clong) output Q.n

    if Q.nlinearity > 0
        println("*Linearities in input file - partial dual solution only")
    end

    print("*Dual:")
    for i in 0:P.d-1
        C_i = unsafe_load(P.C, i+1)
        Col_i = unsafe_load(P.Col, i+1)
        idx = _unsafe_load_inequality(m, C_i)
        print(" y_$(idx)=")
        temp1 = extractbigint(unsafe_load(Q.Lcm, Col_i+1)) * (-1)
        temp1 *= extractbigint(unsafe_load(unsafe_load(P.A, 1), Col_i+1))
        temp2 = extractbigint(unsafe_load(Q.Gcd, Col_i+1))
        temp2 *= extractbigint(P.det)
        print(temp1, "/", temp2)
    end
    println()
end
