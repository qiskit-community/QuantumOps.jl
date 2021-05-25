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

#abstract type AbstractSumA{OpT, StringT, CoeffT} end

struct OpSum{OpT<:AbstractOp, StringT, CoeffT} <: AbstractSum{OpT, StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function OpSum(strings, coeffs; already_sorted=false)
        _abstract_sum_inner_constructor_helper!(strings, coeffs; already_sorted=already_sorted)
        return new{eltype(eltype(strings)), typeof(strings), typeof(coeffs)}(strings, coeffs)
    end

    # function OpSum{OpT}(strings, coeffs; already_sorted=false) where {OpT}
    #     _abstract_sum_inner_constructor_helper!(strings, coeffs; already_sorted=already_sorted)
    #     return new{eltype(eltype(strings)), typeof(strings), typeof(coeffs)}(strings, coeffs)
    # end
end

OpSum{T}(strings, coeffs; already_sorted=false) where T = OpSum(strings, coeffs; already_sorted=already_sorted)

strip_typeof(::OpSum{W, T, CoeffT}) where {W, T, CoeffT} = OpSum{W}
term_type(::Type{T}) where {V, T <: OpSum{V}} = OpTerm{V}

OpSum{T}(args...) where T = OpSum(args...)

function OpSum{T}(strings::AbstractVector{<:AbstractString}, coeffs) where T
    return OpSum(Vector{T}.(strings), coeffs)
end

function OpSum(v::AbstractVector{<:OpTerm}; already_sorted=false)
    strings = [x.ops for x in v]
    coeffs = [x.coeff for x in v]
    return OpSum(strings, coeffs; already_sorted=already_sorted)
end

sum_type(::Type{<:OpTerm{T}}) where T = OpSum{T}
term_type(::Type{<:OpSum{T}}) where T = OpTerm{T}

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

function OpSum(v::AbstractVector{<:OpTerm}; already_sorted=false)
    strings = [op_string(x) for x in v]
    coeffs = [x.coeff for x in v]
    return OpSum(strings, coeffs; already_sorted=already_sorted)
end
