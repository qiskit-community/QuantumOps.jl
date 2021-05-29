## Might be more elegant and concise to have Fermi terms and Pauli terms be
## OpTerm{FermiOp} and OpTerm{Pauli}, instead of FermiTerm, etc.
## We implement this idea here.

struct OpTerm{W<:AbstractOp, T<:AbstractVector{W}, CoeffT} <: AbstractTerm{W, CoeffT}
    ops::T
    coeff::CoeffT
end

op_string(t::OpTerm) = t.ops
term_type(::Type{T}) where T<:AbstractOp = OpTerm{T}
strip_typeof(::OpTerm{W, T, CoeffT}) where {W, T, CoeffT} = OpTerm{W}

OpTerm(term::AbstractVector, coeff=_DEFAULT_COEFF) = OpTerm(term, coeff)
OpTerm{T}(term::V, coeff=_DEFAULT_COEFF) where {T, V<:AbstractVector{T}} = OpTerm(term, coeff)
OpTerm{T}(s::AbstractString, coeff=_DEFAULT_COEFF) where T = OpTerm(Vector{T}(s), coeff)

function OpTerm{OpT}(inds::AbstractVector{<:Integer}, coeff=_DEFAULT_COEFF) where OpT <: AbstractOp
    return OpTerm(OpT.(inds), coeff)
end

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FIXME: empty ops causes overflow here
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function OpTerm{OpT}(ops::AbstractOp...; coeff=_DEFAULT_COEFF) where OpT
    return OpTerm([ops...], coeff)
end

function OpTerm(ops::AbstractOp...; coeff=_DEFAULT_COEFF)
    return OpTerm([ops...], coeff)
end

function OpTerm{OpT}(s::Symbol, coeff=_DEFAULT_COEFF) where OpT
    return OpTerm{OpT}(Vector{OpT}(String(s)), coeff)
end

"""
    rand_op_term(::Type{OpT}, n::Integer; coeff=_DEFAULT_COEFF) where {OpT <: AbstractOp}

Return a random `OpTerm` of `n` tensor factors.
"""
rand_op_term(::Type{OpT}, n::Integer; coeff=_DEFAULT_COEFF) where {OpT <: AbstractOp} =
    term_type(OpT)(rand(OpT, n), coeff)

####
#### OpSum
####

struct OpSum{OpT<:AbstractOp, StringT, CoeffT} <: AbstractSum{OpT, StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function OpSum(strings, coeffs; already_sorted=false)
        _abstract_sum_inner_constructor_helper!(strings, coeffs; already_sorted=already_sorted)
        return new{eltype(eltype(strings)), typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

strip_typeof(::OpSum{W, T, CoeffT}) where {W, T, CoeffT} = OpSum{W}

#term_type(::Type{T}) where {V, T <: OpSum{V}} = OpTerm{V}
term_type(::Type{<:OpSum{T}}) where T = OpTerm{T}

sum_type(::Type{<:OpTerm{T}}) where T = OpSum{T}

#OpSum{T}(strings, coeffs; already_sorted=false) where T = OpSum(strings, coeffs; already_sorted=already_sorted)
#OpSum{T}(args...) where T = OpSum(args...)
OpSum{T}(args...; kwargs...) where T = OpSum(args...; kwargs...)

function OpSum{T}(strings::AbstractVector{<:AbstractString}, coeffs) where T
    return OpSum(Vector{T}.(strings), coeffs)
end

# function OpSum(v::AbstractVector{<:OpTerm}; already_sorted=false)
#     strings = [x.ops for x in v]
#     coeffs = [x.coeff for x in v]
#     return OpSum(strings, coeffs; already_sorted=already_sorted)
# end

function OpSum(v::AbstractVector{<:OpTerm}; already_sorted=false)
    strings = [op_string(x) for x in v]
    coeffs = [x.coeff for x in v]
    return OpSum(strings, coeffs; already_sorted=already_sorted)
end

## Enable this again. why broken
#term_type(::Type{T}) where {T <: AbstractOp} = OpTerm{T}

"""
    OpSum(v::AbstractMatrix{<:AbstractOp}, coeffs=fill(_DEFAULT_COEFF, size(v, 1)))

Construct a sum from a matrix of single-particle operators. If `size(v) == (m, n)`, then
the the sum has `m` terms with `n` operators in each string.
"""
function OpSum(v::AbstractMatrix{<:AbstractOp}, coeffs=fill(_DEFAULT_COEFF, size(v, 1)))
    strings = @inbounds [v[i,:] for i in 1:size(v,1)]
    return OpSum(strings, coeffs)
end

OpSum{OpT}() where OpT <: AbstractOp = OpSum(Vector{OpT}[], Complex{Float64}[])

"""
    rand_op_sum(::Type{OpT}, n_factors::Integer, n_terms::Integer; coeff_func=nothing) where {OpT <: AbstractOp}

Return an `OpSum{OpT}` of `n_terms` terms of `n_factors` factors each.

If `coeff_func` is `nothing`, then the coefficients are all equal to one. Otherwise `coeff_func` must be a function that
takes one argument `n_terms`, and returns `n_terms` coefficients.

# Examples
```julia
julia> rand_op_sum(Pauli, 4, 3)
(1 + 0im) * IIYX
(1 + 0im) * IIZY
(1 + 0im) * ZXIX

julia> rand_op_sum(FermiOp, 3, 4; coeff_func=randn)
N-I * -1.125005096910286
EI- * 0.8400663971914751
EE- * -0.11901005957119838
-E- * 0.4507511196805674
```
"""
function rand_op_sum(::Type{OpT}, n_factors::Integer, n_terms::Integer; coeff_func=nothing) where {OpT <: AbstractOp}
    paulis = [rand(OpT, n_factors) for i in 1:n_terms]
    if coeff_func == nothing
        coeffs = fill(_DEFAULT_COEFF, n_terms)
    else
        coeffs = coeff_func(n_terms)
    end
    return OpSum{OpT}(paulis, coeffs)
end

#####
##### Conversion
#####

"""
    Matrix(ps::PauliSumA)

Convert `ps` to a dense `Matrix`.

# Examples

We convert a matrix to a `PauliSumA` and then back to a matrix.
```jldoctest
julia> m = [0.1 0.2; 0.3 0.4];

julia> PauliSumA(m)
(0.25 + 0.0im) * I
(0.25 + 0.0im) * X
(0.0 - 0.04999999999999999im) * Y
(-0.15000000000000002 + 0.0im) * Z

julia> Matrix(PauliSumA(m)) ≈ m
true
```
"""
Base.Matrix(ps::OpSum) = Matrix(SparseArrays.sparse(ps))

## ThreadsX helps enormously for large sums. 22x faster for 4^8 terms
"""
    sparse(ps::PauliSum)

Convert `ps` to a sparse matrix.
"""
SparseArrays.sparse(ps::OpSum) = ThreadsX.sum(SparseArrays.sparse(ps[i]) for i in eachindex(ps))

####
#### Algebra / mathematical operations
####

# We use lmul! because that's how LinearAlgebra offers "scaling" of a Matrix (or rmul!)
"""
    lmul!(psum::PauliSum, n)

Left multiplies the coefficient of `psum` by `n` in place.
"""
function LinearAlgebra.lmul!(opsum::OpSum, n)
    @. opsum.coeffs = n * opsum.coeffs
    return opsum
end

function Base.zero(term::OpTerm{T}) where T
    new_string = similar(op_string(term))
    fill!(new_string, zero(T))
    return strip_typeof(term)(new_string, zero(term.coeff))
end
