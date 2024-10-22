####
#### OpSum
####

struct OpSum{OpT<:AbstractOp, StringT, CoeffT} <: AbstractSum{OpT, StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function OpSum(strings, coeffs; already_sorted=false)
        if length(strings) != length(coeffs)
            throw(DimensionMismatch("bad dims"))
        end
        if ! isempty(strings)
            if ! already_sorted  # Slightly dangerous to only do length check if not sorted
                n = length(first(strings))
                if ! all(x -> length(x) == n, strings)
                    throw(DimensionMismatch("Operator strings are of differing lengths."))
                end
                sort_and_sum_duplicates!(strings, coeffs)
            end
        end
        return new{eltype(eltype(strings)), typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

strip_typeof(::OpSum{W, T, CoeffT}) where {W, T, CoeffT} = OpSum{W}

"""
    term_type(::Type{<:OpSum{T}})

Return the the concrete `OpTerm` type associated with `OpSum{T}`.
This is used in generic code to construct terms given the type of a sum.
"""
term_type(::Type{<:OpSum{T}}) where T = OpTerm{T}

"""
    sum_type(::Type{<:OpTerm{T}})

Return the the concrete `OpSum` type associated with `OpTerm{T}`.
This is used in generic code to construct sums given the type of a term.
"""
sum_type(::Type{<:OpTerm{T}}) where T = OpSum{T}

#OpSum{T}(strings, coeffs; already_sorted=false) where T = OpSum(strings, coeffs; already_sorted=already_sorted)
#OpSum{T}(args...) where T = OpSum(args...)
OpSum{T}(args...; kwargs...) where T = OpSum(args...; kwargs...)

function OpSum{T}(strings::AbstractVector{<:AbstractString}, coeffs) where T
    return OpSum(Vector{T}.(strings), coeffs)
end

function OpSum{T}(strings::AbstractVector{<:AbstractString}) where T
    coeffs = fill(_DEFAULT_COEFF, length(strings))
    return OpSum{T}(strings, coeffs)
end

function OpSum(v::AbstractVector{<:OpTerm}; already_sorted=false)
    strings = [op_string(x) for x in v]
    coeffs = [x.coeff for x in v]
    return OpSum(strings, coeffs; already_sorted=already_sorted)
end

"""
    empty(s::OpSum)

Create an empty `OpSum`, that is with no terms, similar to `s`.
"""
Base.empty(s::OpSum) = strip_typeof(s)(empty(s.strings), empty(s.coeffs))

Base.isempty(s::OpSum) = isempty(s.strings) && isempty(s.coeffs)  # Do we need both. Should we check ?

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
    term_string = [rand(OpT, n_factors) for i in 1:n_terms]
    if coeff_func === nothing
        coeffs = fill(_DEFAULT_COEFF, n_terms)
    else
        coeffs = coeff_func(n_terms)
    end
    return OpSum{OpT}(term_string, coeffs)
end

#####
##### Conversion
#####

"""
    Matrix(ps::OpSum)

Convert `ps` to a dense `Matrix`.

# Examples

We convert a matrix to a `OpSum` and then back to a matrix.
```jldoctest
julia> m = [0.1 0.2; 0.3 0.4];

julia> PauliSum(m)
(0.25 + 0.0im) * I
(0.25 + 0.0im) * X
(0.0 - 0.04999999999999999im) * Y
(-0.15000000000000002 + 0.0im) * Z

julia> Matrix(PauliSum(m)) ≈ m
true
```
"""
Base.Matrix(ps::OpSum) = Matrix(SparseArrays.sparse(ps))

## ThreadsX helps enormously for large sums. 22x faster for 4^8 terms
## But, ThreadsX is five times slower for 1-q sum
"""
    sparse(ps::OpSum)

Convert `ps` to a sparse matrix.
"""
function SparseArrays.sparse(ps::OpSum)
    if length(ps) >= 64
        return ThreadsX.sum(SparseArrays.sparse(ps[i]) for i in eachindex(ps))
    else
        return sum(SparseArrays.sparse(ps[i]) for i in eachindex(ps))
    end
end

####
#### Compare / predicates
####

## TODO: improve this. Check for coefficients, etc.
function Base.iszero(opsum::OpSum)
    return iszero(length(opsum)) || all(iszero, opsum.coeffs)
end

Base.zero(opsum::OpSum) = empty(opsum)

function ZChop.zchop(opsum::OpSum, zeps::Real=ZChop.ZEPS)
    newopsum = empty(opsum)
    newstrings = newopsum.strings
    newcoeffs = newopsum.coeffs
    for i in eachindex(opsum.coeffs)
        c = opsum.coeffs[i]
        if ! iszero(ZChop.zchop(c, zeps))
            push!(newstrings, opsum.strings[i])
            push!(newcoeffs, c)
        end
    end
    return newopsum
end


function ZChop.zchop!(opsum::OpSum, zeps::Real=ZChop.ZEPS)
    inds = Int[]
    for i in eachindex(opsum.coeffs)
        if iszero(ZChop.zchop(opsum.coeffs[i], zeps))
            push!(inds, i)
        end
    end
    deleteat!(opsum.coeffs, inds)
    deleteat!(opsum.strings, inds)
    return opsum
end


####
#### Algebra / mathematical operations
####

# We use lmul! because that's how LinearAlgebra offers "scaling" of a Matrix (or rmul!)
"""
    lmul!(psum::OpSum, n)

Left multiplies the coefficient of `psum` by `n` in place.
"""
function LinearAlgebra.lmul!(opsum::OpSum, n)
    coeffs = opsum.coeffs
    # FastBroadcast.@.. opsum.coeffs = n * opsum.coeffs
    # @. opsum.coeffs = n * opsum.coeffs
    @inbounds for i in eachindex(coeffs)
        coeffs[i] = coeffs[i] * n
    end
    return opsum
end
