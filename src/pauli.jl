module Paulis

import ..AbstractPauli
import .._AbstractPauli
import ..op_index

export Pauli

# TODO: Make this a `module` to create a namespace

"""
    struct Pauli <: AbstractPauli

This is the only implementation of `AbstractPauli`
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

const PAULIS = (I, X, Y, Z)

# Indexing into static Tuple seems faster
"""
    Pauli(ind::Union{Integer, Symbol, AbstractString, AbstractChar})

Return a `Pauli` indexed by `[0, 3]` or a representation of `[I, X, Y, Z]`.
"""
function Pauli(ind::Integer)::Pauli
    return PAULIS[ind + 1]
end

Pauli(s::Union{Symbol, AbstractString, AbstractChar}) = _AbstractPauli(Pauli, s)

Base.copy(p::Pauli) = p

#op_index(p::Pauli) = (p.hi << 1) + p.lo
op_index(p::Pauli) = 2 * p.hi + p.lo

####
#### Compare / predicates
####

Base.:(==)(p1::Pauli, p2::Pauli) = p1.hi == p2.hi  && p1.lo == p2.lo

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
