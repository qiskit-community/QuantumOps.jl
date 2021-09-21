## Might be more elegant and concise to have Fermi terms and Pauli terms be
## OpTerm{FermiOp} and OpTerm{Pauli}, instead of FermiTerm, etc.
## We implement this idea here.

"""
    struct OpTerm{W<:AbstractOp, T<:AbstractVector{W}, CoeffT} <: AbstractTerm{W, CoeffT}

The operator term
"""
struct OpTerm{W<:AbstractOp, T<:AbstractVector{W}, CoeffT} <: AbstractTerm{W, CoeffT}
    ops::T
    coeff::CoeffT
end

"""
    op_string(t::OpTerm)

Return the string of operators in `t`, omitting the coefficient.

# Examples
```jldocstring
julia> op_string(PauliTerm("XY", 2.0))
2-element Vector{Pauli}:
 Pauli: X
 Pauli: Y
```
"""
op_string(t::OpTerm) = t.ops

"""
    term_type(::Type{T}) where T <: AbstractOp

Return the concrete type `OpTerm{T}` associated with `T`.
This is used in generic code to construct an `OpTerm` when `T` is known.
"""
term_type(::Type{T}) where {T <: AbstractOp} = OpTerm{T}

"""
    strip_typeof(::OpTerm{W, T, CoeffT}) where {W, T, CoeffT}

Return the type of the argument with concrete `T` and `CoeffT` replaced
by variables. `strip_typeof` is useful when a constructor is available only
for the simpler type returned by `strip_typeof`.
# Examples
```jldoctest
julia> typeof(PauliTerm("XX"))
PauliTerm{Vector{Pauli}, Complex{Int64}} (alias for OpTerm{Pauli, Array{Pauli, 1}, Complex{Int64}})

julia> strip_typeof(PauliTerm("XX"))
PauliTerm (alias for OpTerm{Pauli, T} where T<:AbstractArray{Pauli, 1})

julia> typeof(OpTerm{PauliI}("XX"))
OpTerm{PauliI, Vector{PauliI}, Complex{Int64}}

julia> QuantumOps.strip_typeof(OpTerm{PauliI}("XX"))
OpTerm{PauliI, T} where T<:AbstractVector{PauliI}
```

The following attempt to construct an operator will throw an exception because the requested
method does not exist.
```julia
julia> typeof(OpTerm{PauliI}("XX"))("YY")
ERROR: MethodError: no method matching OpTerm{PauliI, Vector{PauliI}, Complex{Int64}}(::String)
```

But, this works as intended.
```jldoctest
julia> QuantumOps.strip_typeof(OpTerm{PauliI}("XX"))("YY")
2-factor OpTerm{PauliI, Vector{PauliI}, Complex{Int64}}:
YY * (1 + 0im)
```
"""
strip_typeof(::OpTerm{W, T, CoeffT}) where {W, T, CoeffT} = OpTerm{W}

"""
    DenseOpTerm

An alias for `OpTerm` with dense storage.
"""
const DenseOpTerm = OpTerm{W, T} where {W<:AbstractOp, T<:Union{Vector{W}, StaticArrays.SVector}}

## TODO: Commented out and tests pass. Is this necessary?
# OpTerm(term::AbstractVector, coeff=_DEFAULT_COEFF) = OpTerm(term, coeff)

"""
    OpTerm{T}(term::V, coeff=_DEFAULT_COEFF) where {T, V<:AbstractVector{T}}

Construct an `OpTerm{T}` from a vector of `T`.

# Examples
```jldoctest
julia> OpTerm{Pauli}([X, Y, Z])
3-factor PauliTerm{Vector{Pauli}, Complex{Int64}}:
XYZ * (1 + 0im)
```
"""
OpTerm{T}(term::V, coeff=_DEFAULT_COEFF) where {T, V<:AbstractVector{T}} = OpTerm(term, coeff)

OpTerm(term::V, coeff=_DEFAULT_COEFF) where {T, V<:AbstractVector{T}} = OpTerm(term, coeff)

"""
    OpTerm{T}(s::AbstractString, coeff=_DEFAULT_COEFF)

Construct `OpTerm{T}` from a string representation.

# Examples
```jldoctest
julia> OpTerm{Pauli}("XYZI")
4-factor PauliTerm{Vector{Pauli}, Complex{Int64}}:
XYZI * (1 + 0im)

julia> OpTerm{FermiOp}("++--")
4-factor FermiTerm{Vector{FermiOp}, Complex{Int64}}:
++-- * (1 + 0im)
```

Recall that `FermiTerm` is an alias for `OpTerm{FermiOp}`.

```jldoctest
julia> FermiTerm("++--", 2.0)
4-factor FermiTerm{Vector{FermiOp}, Complex{Int64}}:
++-- * 2.0
```
"""
OpTerm{T}(s::AbstractString, coeff=_DEFAULT_COEFF) where T = OpTerm(Vector{T}(s), coeff)

## TODO: document and give examples for all signatures below
"""
    OpTerm{OpT}(inds::AbstractVector{<:Integer}, coeff=_DEFAULT_COEFF) where OpT <: AbstractOp
    OpTerm{OpT}(ops::AbstractOp...; coeff=_DEFAULT_COEFF)
    OpTerm(ops::AbstractOp...; coeff=_DEFAULT_COEFF)
    OpTerm{OpT}(s::Symbol, coeff=_DEFAULT_COEFF) where OpT

Construct an `OpTerm{OpT}` from an array of integer indices (codes) `inds` for `OpT`,
from arguments `ops`, etc.

# Examples

```jldoctest
julia> PauliTerm([0,1,2,3])
4-factor PauliTerm{Vector{Pauli}, Complex{Int64}}:
IXYZ * (1 + 0im)

julia> FermiTerm([0,1,2,3,4,5])
6-factor FermiTerm{Vector{FermiOp}, Complex{Int64}}:
INE+-0 * (1 + 0im)

julia> OpTerm{PauliI}([3, 2, 1, 0])
4-factor OpTerm{PauliI, Vector{PauliI}, Complex{Int64}}:
ZYXI * (1 + 0im)

julia> PauliTerm(0:3)
4-factor PauliTerm{Vector{Pauli}, Complex{Int64}}:
IXYZ * (1 + 0im)
```
"""
function OpTerm{OpT}(inds::AbstractVector{<:Integer}, coeff=_DEFAULT_COEFF) where OpT <: AbstractOp
    return OpTerm(OpT.(inds), coeff)
end

## Documented elsewhere (above, when this comment was written)
function OpTerm{OpT}(ops::AbstractOp...; coeff=_DEFAULT_COEFF) where OpT
    if isempty(ops)
        return OpTerm(OpT[], coeff)
    end
    return OpTerm([ops...], coeff)
end

## Documented elsewhere (above, when this comment was written)
function OpTerm(ops::AbstractOp...; coeff=_DEFAULT_COEFF)
    return OpTerm([ops...], coeff)
end

function OpTerm{OpT}(s::Symbol, coeff=_DEFAULT_COEFF) where OpT
    return OpTerm{OpT}(Vector{OpT}(String(s)), coeff)
end

function _op_term_macro_helper(::Type{T}, str) where T
    v = split(str, "*")
    if length(v) == 1
        return :($T($str))
    end
    if length(v) == 2
        num = parse(ComplexF64, v[1])
        return :($T($(strip(v[2])), $num))
    end
    throw(ErrorException("bad"))
end

"""
    rand_op_term(::Type{OpT}, n::Integer; coeff=_DEFAULT_COEFF) where {OpT <: AbstractOp}

Return a random `OpTerm` of `n` tensor factors.
"""
rand_op_term(::Type{OpT}, n::Integer; coeff=_DEFAULT_COEFF) where {OpT <: AbstractOp} =
    term_type(OpT)(rand(OpT, n), coeff)

####
#### Algebra
####

function multiply_keeping_phase!(target::AbstractArray{<:AbstractOp},
                                s1::AbstractArray{<:AbstractOp}, s2::AbstractArray{<:AbstractOp})
    length(s1) != length(s2) && throw(DimensionMismatch())
    _phase_data = AbstractOps.phase_data(eltype(s1))
    @inbounds for i in eachindex(s1)
        AbstractOps.accumulate_phase!(_phase_data, s1[i], s2[i])
        _prod = s1[i] * s2[i]
        if iszero(_prod)  # If any factor is zero, the entire product is zero
            fill!(target, Zero)
            return (target, 1)
        end
        target[i] = _prod
    end
    cum_phase = AbstractOps.compute_phase(_phase_data, s1, s2)
    return(target, cum_phase)
end

function mul!(target::AbstractArray{<:AbstractOp}, ps1::DenseOpTerm, ps2::DenseOpTerm)
    s_new, phase = multiply_keeping_phase!(target, op_string(ps1), op_string(ps2))
    if iszero(s_new)
        return strip_typeof(ps1)(s_new, zero(ps1.coeff))
    end
    return strip_typeof(ps1)(s_new, ps1.coeff * ps2.coeff * phase)
end

function Base.:*(ps1::DenseOpTerm, ps2::DenseOpTerm)
    return mul!(similar(op_string(ps1)), ps1, ps2)
end

function Base.zero(term::OpTerm{T}) where T
    new_string = similar(op_string(term))
    fill!(new_string, zero(T))
    return strip_typeof(term)(new_string, zero(term.coeff))
end
