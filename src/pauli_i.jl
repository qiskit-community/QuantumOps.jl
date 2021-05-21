"""
    PaulisI

This module contains `PauliI <: AbstractPauli` and supporting code.
"""
module PaulisI

import ..AbstractPauli
import .._AbstractPauli
import ..op_index

export PauliI

# TODO: Make this a `module` to create a namespace

"""
    struct PauliI <: AbstractPauli

Another implementation of `AbstractPauli`.

This implementation encodes the Pauli operators via the index `0-3`
in a single byte.
"""
struct PauliI <: AbstractPauli{Complex{Int}}
    ind::Int8
end

####
#### Constructors
####

const I = PauliI(0)
const X = PauliI(1)
const Y = PauliI(2)
const Z = PauliI(3)

const PAULIS = (I, X, Y, Z)

# Indexing into static Tuple seems faster
"""
    Pauli(ind::Union{Integer, Symbol, AbstractString, AbstractChar})

Return a `Pauli` indexed by `[0, 3]` or a representation of `[I, X, Y, Z]`.
"""
PauliI(s::Union{Symbol, AbstractString, AbstractChar}) = _AbstractPauli(PauliI, s)

Base.copy(p::PauliI) = p

op_index(p::PauliI) = p.ind

####
#### Compare / predicates
####

Base.:(==)(p1::PauliI, p2::PauliI) = p1.ind == p2.ind

## This is maybe not really necessary, but this is stricter than the fallback method.
## Also, this would have prevented a perf regression bug.
Base.isone(p::PauliI) = p === one(PauliI)

####
#### Algebra / mathematical operations
####

## We try two versions of the multiplication table.
## Vector and matrix versions have similar performance.

const pauli_i_mult =
    ((0, 1, 2, 3),
     (1, 0, 3, 2),
     (2, 3, 0, 1),
     (3, 2, 1, 0))

const pauli_i_mult_1d =
    (0, 1, 2, 3,
     1, 0, 3, 2,
     2, 3, 0, 1,
     3, 2, 1, 0)


using ..Paulis

"""
    *(p1::PauliI, p2::PauliO)

Multiplication that returns a `PauliI`, but ignores the phase.
"""
Base.:*(p1::PauliI, p2::PauliI) = PauliI(@inbounds pauli_i_mult_1d[p1.ind * 4 + p2.ind + 1])
#Base.:*(p1::PauliI, p2::PauliI) = PauliI(pauli_i_mult[p1.ind+1][p2.ind+1])
# This one allocates for some reason
#Base.:*(p1::PauliI, p2::PauliI) = PauliI(op_index(Pauli(p1.ind) * Pauli(p2.ind)))

end # module PaulisI
