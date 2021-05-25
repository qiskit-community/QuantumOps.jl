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
end

strip_typeof(::OpSum{W, T, CoeffT}) where {W, T, CoeffT} = OpSum{W}
term_type(::Type{T}) where {V, T <: OpSum{V}} = OpTerm{V}

function OpSum{T}(strings::AbstractVector{<:AbstractString}, coeffs) where T
    return OpSum(Vector{T}.(strings), coeffs)
end

OpSum{T}(strings, coeffs) where T = OpSum(strings, coeffs)

function OpSum(v::AbstractVector{<:OpTerm}; already_sorted=false)
    strings = [x.ops for x in v]
    coeffs = [x.coeff for x in v]
    return OpSum(strings, coeffs; already_sorted=already_sorted)
end
