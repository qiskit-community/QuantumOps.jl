module Paulis

import ..AbstractOps
import ..AbstractOps: op_index, unsafe_op, symplectic, _AbstractOp
using ..AbstractPaulis: AbstractPauli

export Pauli, unsafe_op
import IsApprox

"""
    struct Pauli <: AbstractPauli

An implementation of `AbstractPauli`.
"""
struct Pauli <: AbstractPauli{Complex{Int}}
    hi::Bool
    lo::Bool
end

####
#### Constructors
####

const I = Pauli(0, 0)
const X = Pauli(0, 1)
const Y = Pauli(1, 0)
const Z = Pauli(1, 1)

# Indexing into static Tuple seems faster
"""
    Pauli(ind::Union{Integer, Symbol, AbstractString, AbstractChar})

Return a `Pauli` indexed by `[0, 3]` or a representation of `[I, X, Y, Z]`.
"""
Pauli(ind::Integer)::Pauli = (I, X, Y, Z)[ind + 1]

## The unsafe_op, unsafe_pauli thing does not seem to be more efficient.

@inline AbstractOps.unsafe_op(::Type{Pauli}, ind::Integer) = unsafe_pauli(ind)

@inline function unsafe_pauli(ind::Integer)::Pauli
    return  @inbounds (I, X, Y, Z)[ind + 1]
end

Pauli(s::Union{Symbol, AbstractString, AbstractChar}) = _AbstractOp(Pauli, s)

Base.copy(p::Pauli) = p

op_index(p::Pauli) = 2 * p.hi + p.lo

function symplectic(p::Pauli)
    if p.hi == 1
        if p.lo == 1
            return (false, true)
        else
            return (true, true)
        end
    else
        if p.lo == 1
            return (true, false)
        else
            return (false, false)
        end
    end
end

####
#### Compare / predicates
####

Base.:(==)(p1::Pauli, p2::Pauli) = p1.hi == p2.hi && p1.lo == p2.lo

## This is maybe not really necessary, but this is stricter than the fallback method.
## Also, this would have prevented a perf regression bug.
Base.isone(p::Pauli) = p === one(Pauli)

####
#### Algebra / mathematical operations
####

"""
    *(p1::Pauli, p2::Pauli)

Multiplication that returns a `Pauli`, but ignores the phase.
"""
@inline Base.:*(p1::Pauli, p2::Pauli) = Pauli(p1.hi ⊻ p2.hi, p1.lo ⊻ p2.lo)

end # module Paulis
