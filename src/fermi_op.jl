####
#### FermiOp
####

module FermiOps

import Random
import IsApprox

export FermiOp, FermiDefault
export I, NumberOp, Empty, Raise, Lower, Zero, count_raise_lower, is_raise_lower

import ..AbstractOps: _AbstractOp, op_symbols, AbstractOp, _show_op_plain, op_index
import ..AbstractFermiOp

struct FermiOp <: AbstractFermiOp
    ind::Int  # Int is often faster than UInt8
#    ind::UInt8
end

const FermiDefault = FermiOp

op_index(op::FermiOp) = op.ind

Base.one(op::FermiOp) = I
Base.one(::Type{FermiOp}) = I

Base.zero(op::FermiOp) = Zero
Base.zero(::Type{FermiOp}) = Zero

"""
    isless(op1::FermiOp, op2::FermiOp)

Canonical (lexical) order of `FermiOp` is I < N < (I - N) < Raise < Lower < Zero < (Z = N - E)
"""
Base.isless(op1::FermiOp, op2::FermiOp) = isless(op_index(op1), op_index(op2))

const _fermi_chars =
     ('I',   'N',    'E',   '+',   '-',  '0', 'Z', 'X')
#      0      1       2      3      4     5    6    7
const (I, NumberOp, Empty, Raise, Lower, Zero, Z,  NoOp) = FermiOp.(0:7)

const _fermi_string_chars = string.(_fermi_chars)
const _fermi_symbol_chars = Symbol.(_fermi_chars)

op_symbols(::Type{<:AbstractFermiOp}, ::Type{<:AbstractChar}) = _fermi_chars
op_symbols(::Type{<:AbstractFermiOp}, ::Type{<:AbstractString}) = _fermi_string_chars
op_symbols(::Type{<:AbstractFermiOp}, ::Type{Symbol}) = _fermi_symbol_chars

FermiOp(s::Union{AbstractChar, AbstractString, Symbol}) = _AbstractOp(FermiOp, s)

function _show_op_plain(io::IO, op::FermiOp)
    i = op_index(op) + 1
    if i < 1 || i > 8
        i = 8  # uninitialized data is coerced to 8 which means NoOp
    end
    print(io, _fermi_chars[i])
end

function Base.show(io::IO, op::FermiOp)
    print(io, typeof(op), ": ")
    _show_op_plain(io, op)
end

"""
    ferm_op_mult

Multiplication for single-fermion operators. We include `Z=N-E`. Multiplication does not
keep track of phase incurred when `Z` is an operand
"""
const ferm_op_mult =
    #   I       N          E        +         -         0       Z
    ((I, NumberOp, Empty, Raise, Lower, Zero, Z),      # I
     (NumberOp, NumberOp, Zero, Raise, Zero, Zero, Z),    # N
     (Empty, Zero, Empty, Zero, Lower, Zero, Empty),  # E
     (Raise, Zero, Raise, Zero, NumberOp, Zero, Raise), # +
     (Lower, Lower, Zero, Empty, Zero, Zero, Lower),  # -
     (Zero, Zero, Zero, Zero, Zero, Zero, Zero),      # 0
     (Z, NumberOp, Empty, Raise, Lower, Zero, I))      # Z

# const ferm_op_mult =
#     #   I       N          E        +         -         0
#     ((I, NumberOp, Empty, Raise, Lower, Zero),   # I
#      (NumberOp, NumberOp, Zero, Raise, Zero, Zero), # N
#      (Empty, Zero, Empty, Zero, Lower, Zero),   # E
#      (Raise, Zero, Raise, Zero, NumberOp, Zero),  # +
#      (Lower, Lower, Zero, Empty, Zero, Zero),   # -
#      (Zero, Zero, Zero, Zero, Zero, Zero))      # 0

Base.:*(op1::FermiOp, op2::FermiOp) = ferm_op_mult[op_index(op1)+1][op_index(op2)+1]

function Base.inv(op::FermiOp)
    if op === I
        return I
    else
        throw(DomainError(op, "Operator has no inverse"))
    end
end

function Base.:^(op::FermiOp, n::Integer)
    n < 0 && throw(DomainError(n))
    n == 0 && return I
    if op === I || op === NumberOp || op === Empty
        return op
    end
    return Zero
end

function Base.adjoint(op::FermiOp)
    if op === I || op === Empty || op === NumberOp || op === Zero
        return op
    end
    if op === Raise
        return Lower
    end
    if op === Lower
        return Raise
    end
    throw(DomainError(op))
end

"""
    is_raise_lower(x)

Return `true` if `x` is a raising or lowering operator.
"""
is_raise_lower(x) = x === FermiOps.Raise || x === FermiOps.Lower

IsApprox.ishermitian(op::FermiOp) = !is_raise_lower(op)

"""
    count_raise_lower(op::FermiOp)

Return the sum of the numbers of raising and lowering operators in `op`.
Number and complement of number operator each count as two.
"""
function count_raise_lower(op::FermiOp)
    if is_raise_lower(op)
        return 1
    elseif op === NumberOp || op === Empty
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

# const _fermi_mats = (_id_mat, _number_mat, _empty_mat, _raise_mat, _lower_mat, _zero_mat)

function Base.Matrix(op::FermiOp)
    _fermi_mats = (_id_mat, _number_mat, _empty_mat, _raise_mat, _lower_mat, _zero_mat)
    _fermi_mats[op_index(op) + 1]
end

end # module FermiOps
