####
#### FermiOp
####

module FermiOps

import Random
import LinearAlgebra

export FermiOp, FermiDefault
export id_op, number_op, empty_op, raise_op, lower_op, zero_op,
    count_raise_lower

import .._AbstractOp, ..op_symbols, ..AbstractOp, ..AbstractFermiOp

struct FermiOp <: AbstractFermiOp
    ind::UInt8
end

const FermiDefault = FermiOp

op_index(op::FermiOp) = op.ind

Base.one(op::FermiOp) = id_op
Base.one(::Type{FermiOp}) = id_op

Base.zero(op::FermiOp) = zero_op
Base.zero(::Type{FermiOp}) = zero_op

Base.isless(op1::FermiOp, op2::FermiOp) = isless(op_index(op1), op_index(op2))

const id_op = FermiOp(0)
const number_op = FermiOp(1)
const empty_op = FermiOp(2)
const raise_op = FermiOp(3)
const lower_op = FermiOp(4)
const zero_op = FermiOp(5)
const Z_op = FermiOp(6)
const no_op = FermiOp(7)

const _fermi_chars = ('I', 'N', 'E', '+', '-', '0', 'Z', 'X')
const _fermi_string_chars = string.(_fermi_chars)
const _fermi_symbol_chars = Symbol.(_fermi_chars)

op_symbols(::Type{<:AbstractFermiOp}, ::Type{<:AbstractChar}) = _fermi_chars
op_symbols(::Type{<:AbstractFermiOp}, ::Type{<:AbstractString}) = _fermi_string_chars
op_symbols(::Type{<:AbstractFermiOp}, ::Type{Symbol}) = _fermi_symbol_chars

FermiOp(s::Union{AbstractChar, AbstractString, Symbol}) = _AbstractOp(FermiOp, s)

function Base.show(io::IO, op::FermiOp)
    print(io, _fermi_chars[op_index(op) + 1])
end

"""
    ferm_op_mult

Multiplication for single-fermion operators. We include `Z=N-E`.
Multiplication does not keep track of phase incurred when `Z` is
an operand
"""
const ferm_op_mult =
    #   I       N          E        +         -         0       Z
    ((id_op, number_op, empty_op, raise_op, lower_op, zero_op, Z_op),      # I
     (number_op, number_op, zero_op, raise_op, zero_op, zero_op, Z_op),    # N
     (empty_op, zero_op, empty_op, zero_op, lower_op, zero_op, empty_op),  # E
     (raise_op, zero_op, raise_op, zero_op, number_op, zero_op, raise_op), # +
     (lower_op, lower_op, zero_op, empty_op, zero_op, zero_op, lower_op),  # -
     (zero_op, zero_op, zero_op, zero_op, zero_op, zero_op, zero_op),      # 0
     (Z_op, number_op, empty_op, raise_op, lower_op, zero_op, id_op))      # Z

# const ferm_op_mult =
#     #   I       N          E        +         -         0
#     ((id_op, number_op, empty_op, raise_op, lower_op, zero_op),   # I
#      (number_op, number_op, zero_op, raise_op, zero_op, zero_op), # N
#      (empty_op, zero_op, empty_op, zero_op, lower_op, zero_op),   # E
#      (raise_op, zero_op, raise_op, zero_op, number_op, zero_op),  # +
#      (lower_op, lower_op, zero_op, empty_op, zero_op, zero_op),   # -
#      (zero_op, zero_op, zero_op, zero_op, zero_op, zero_op))      # 0

Base.:*(op1::FermiOp, op2::FermiOp) = ferm_op_mult[op_index(op1)+1][op_index(op2)+1]

function Base.inv(op::FermiOp)
    if op === id_op
        return id_op
    else
        throw(DomainError(op, "Operator has no inverse"))
    end
end

function Base.:^(op::FermiOp, n::Integer)
    n < 0 && throw(DomainError(n))
    n == 0 && return id_op
    if op === id_op || op === number_op || op === empty_op
        return op
    end
    return zero_op
end

function Base.adjoint(op::FermiOp)
    if op === id_op || op === empty_op || op === number_op || op === zero_op
        return op
    end
    if op === raise_op
        return lower_op
    end
    if op === lower_op
        return raise_op
    end
    throw(DomainError(op))
end

function LinearAlgebra.ishermitian(op::FermiOp)
    if op === id_op || op === empty_op || op === number_op || op === zero_op
        return true
    else
        return false
    end
end

function count_raise_lower(op::FermiOp)
    if op === raise_op || op === lower_op
        return 1
    elseif op === number_op || op === empty_op
        return 2
    else
        return 0
    end
end

const _id_mat = [1. 0.; 0. 1.]
const _zero_mat = [0. 0.; 0. 0.]
const _number_mat = [0. 0.; 0. 1.]
const _empty_mat = [1. 0.; 0. 0.]
const _raise_mat = [0. 1.; 0. 0.]
const _lower_mat = [0. 0.; 1. 0.]

const _fermi_mats = (_id_mat, _number_mat, _empty_mat, _raise_mat, _lower_mat, _zero_mat)

function Base.Matrix(op::FermiOp)
    return _fermi_mats[op_index(op) + 1]
end

end # module FermiOps
