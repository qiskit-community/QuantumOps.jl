module Paulis

import ..AbstractOps: op_index
using ..AbstractPaulis: AbstractPauli, _AbstractPauli

export Pauli

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

# TODO: This is ok, we don't currently have a convenient way to use it.
# unsafe_pauli(ind::Integer)::Pauli @inbounds (I, X, Y, Z)[ind + 1]

Pauli(s::Union{Symbol, AbstractString, AbstractChar}) = _AbstractPauli(Pauli, s)

Base.copy(p::Pauli) = p

op_index(p::Pauli) = 2 * p.hi + p.lo

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
