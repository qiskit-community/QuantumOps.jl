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

OpTerm{T}(term::V, coeff=_DEFAULT_COEFF) where {T, V<:AbstractVector{T}} = OpTerm(term, coeff)

OpTerm{T}(s::AbstractString, coeff=_DEFAULT_COEFF) where T = OpTerm(Vector{T}(s), coeff)

####
#### OpTerm{AbstractFermiOp}
####

function OpTerm{T}(inds::NTuple{N, Int}, coeff, n_modes::Integer) where {T<:AbstractFermiOp, N}
    (ops, phase) = index_to_ops_phase(inds)
    (factors, new_coeff) = _fermi_term(ops, phase, coeff, n_modes)
    return OpTerm{T}(factors, new_coeff)
end

####
#### AbstractSumA
####

abstract type AbstractSumA{OpT, StringT, CoeffT} end

struct OpSum{OpT<:AbstractOp, StringT, CoeffT} <: AbstractSumA{OpT, StringT, CoeffT}
#struct OpSum{OpT, StringT, CoeffT} <: AbstractSumA{OpT, StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function OpSum(strings, coeffs; already_sorted=false)
        _abstract_sum_inner_constructor_helper!(strings, coeffs; already_sorted=already_sorted)
        return new{eltype(eltype(strings)), typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

strip_typeof(::OpSum{W, T, CoeffT}) where {W, T, CoeffT} = OpSum{W}

function OpSum{T}(strings::AbstractVector{<:AbstractString}, coeffs) where T
    return OpSum(Vector{T}.(strings), coeffs)
end
