<####
#### FermiOp
####

module FermiOps

export FermiOp, AbstractFermiOp

import .._AbstractOp, ..op_symbols, ..AbstractOp

abstract type AbstractFermiOp <: AbstractOp end

struct FermiOp <: AbstractFermiOp
    ind::UInt8
end

op_index(op::FermiOp) = op.ind

Base.isless(op1::FermiOp, op2::FermiOp) = isless(op_index(op1), op_index(op2))

const id_op = FermiOp(0)
const number_op = FermiOp(1)
const empty_op = FermiOp(2)
const raise_op = FermiOp(3)
const lower_op = FermiOp(4)
const no_op = FermiOp(5)

const _fermi_chars = ('I', 'N', 'E', '+', '-', 'X')
const _fermi_string_chars = string.(_fermi_chars)
const _fermi_symbol_chars = Symbol.(_fermi_chars)

op_symbols(::Type{<:AbstractFermiOp}, ::Type{<:AbstractChar}) = _fermi_chars
op_symbols(::Type{<:AbstractFermiOp}, ::Type{<:AbstractString}) = _fermi_string_chars
op_symbols(::Type{<:AbstractFermiOp}, ::Type{Symbol}) = _fermi_symbol_chars

FermiOp(s::Union{AbstractChar, AbstractString, Symbol}) = _AbstractOp(FermiOp, s)

function Base.show(io::IO, op::FermiOp)
    print(io, _fermi_chars[op_index(op) + 1])
end

function Base.show(io::IO, ps::AbstractArray{<:AbstractFermiOp})
    for p in ps
        show(io, p)
    end
end

end # module FermiOps
