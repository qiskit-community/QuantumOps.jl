<####
#### FermiOp
####

module FermiOps

import Random

export FermiOp, FermiDefault
export id_op, number_op, empty_op, raise_op, lower_op, zero_op

import .._AbstractOp, ..op_symbols, ..AbstractOp, ..AbstractFermiOp

struct FermiOp <: AbstractFermiOp
    ind::UInt8
end

const FermiDefault = FermiOp

op_index(op::FermiOp) = op.ind

Base.one(op::FermiOp) = id_op
Base.one(::Type{FermiOp}) = id_op

Base.isless(op1::FermiOp, op2::FermiOp) = isless(op_index(op1), op_index(op2))

const id_op = FermiOp(0)
const number_op = FermiOp(1)
const empty_op = FermiOp(2)
const raise_op = FermiOp(3)
const lower_op = FermiOp(4)
const zero_op = FermiOp(5)
const no_op = FermiOp(6)

const _fermi_chars = ('I', 'N', 'E', '+', '-', '0', 'X')
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

Multiplication for single-fermion operators.
"""
const ferm_op_mult =
    ((id_op, number_op, empty_op, raise_op, lower_op, zero_op),
     (number_op, number_op, zero_op, raise_op, zero_op, zero_op),
     (empty_op, zero_op, empty_op, zero_op, lower_op, zero_op),
     (raise_op, zero_op, raise_op, zero_op, number_op, zero_op),
     (lower_op, lower_op, zero_op, empty_op, zero_op, zero_op),
     (zero_op, zero_op, zero_op, zero_op, zero_op, zero_op))

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

end # module FermiOps
