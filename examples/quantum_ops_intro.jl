# Activate the package from its development location
# import Pkg
# Pkg.activate("/home/lapeyre/.julia/gjl-quantum/QuantumOps")

# First, import some identifiers for types
import QuantumOps ## Import the identifier `QuantumOps` for using fully qualified identifiers
using QuantumOps: AbstractOp, AbstractFermiOp, FermiOp, AbstractPauli, Pauli, PauliI

# # QuantumOps
# ## Overview
# `QuantumOps` is built mainly around three levels of operator types
# * `AbstractOp` representing single-particle Fermionic or Pauli operators
# * `OpTerm` a parametric type containing a string of operators as an `AbstractVector{<:AbstractOp}`
#   and a coefficient. `OpTerm` represents a multi-particle term. The string type may be a dense
#  `Vector` or a `SparseVector` (with identities omitted). The coefficient may be numeric or symbolic.
# * `OpSum` a parametric type representing a sum of `OpTerm`s.
#
# Some features are
# * `QuantumOps` tries to make the structures and methods efficient, but there is little done to optimize
#    them. Still, `QuantumOps` is often one to three orders of magnitude faster
#    than Qiskit.
# * `QuantumOps` works together `ElectronicStructure.jl` to represent electonic Hamiltonians.
# *  The Jordan-Wigner transform is implemented as an example.
#
# ## Simple operators -- `AbstractOp`
# Types for Fermionic and Pauli operators on a single mode or qubit have `AbstractOp` as an ancestor.
# We have
FermiOp <: AbstractFermiOp <: AbstractOp
# and
Pauli <: AbstractPauli <: AbstractOp
# There is an alternative encoding for Pauli operators
PauliI <: AbstractPauli <: AbstractOp
# `PauliI` may be more efficient in some circumstances.
# None (almost) of `QuantumOps` is written explicitly against a subtype of `AbstractPauli`,
# so `Pauli` and `PauliI` may be used interchangeably. Eventually probably only one of the two will be retained.
#
# All subtypes of `AbstractPauliOp` may instantiated by an integer index like this
(Pauli(0), Pauli(1), Pauli(2), Pauli(3))

# Recall that this expression may also be written using broadcasting like this
Pauli.((0, 1, 2, 3))

# The alternate encoding works the same way
PauliI.((0, 1, 2, 3))

# Seven (or maybe six) Fermi operators are supported
FermiOp.((0:6...,))

# The symbols mean the following


# (I = identity, N = number, E = complement of number (empty), + = raising, - = lowering, 0 = zero, Z = N - E)

# For convenience, you can use variables that name the operators. For example
using QuantumOps.Paulis: I, X, Y, Z
(I, X, Y, Z)

# The same variables are defined for `PauliI`
QuantumOps.PaulisI.X

# The details of `Pauli`, `PauliI`, and `FermiOp` are hidden. For example, even `OpTerm` and supporting code
# knows nothing about them. But here is a look indside
dump(QuantumOps.Paulis.X)
#-
dump(QuantumOps.PaulisI.X)
#-
dump(QuantumOps.FermiOps.NumberOp)

# ## Multi-qubit/mode operators -- `OpTerm`
# Operators on multiple qubits or degrees of freedom, together with a coefficient, are represented by `OpTerm`.
using QuantumOps: OpTerm, PauliTerm, FermiTerm

# `PauliTerm` and `FermiTerm` are aliases
OpTerm{Pauli}
#-
OpTerm{FermiOp}
#-
PauliTerm === OpTerm{Pauli}
# `Pauli` is the default encoding, as there is no alias for `OpTerm{PauliI}`
OpTerm{PauliI}

# One way to instantiate an `OpTerm` is from a string, like this
OpTerm{Pauli}("IXIYIZ", 1.0 + 1.0im)
# Or, using `PauliI`
OpTerm{PauliI}("IXIYIZ", 1.0 + 1.0im)
# Or, for Fermionic operators
OpTerm{FermiOp}("++I--NE", 1.0)
# Alternatively, we can use the aliases for convenience
PauliTerm("IXIYIZ", 1.0 + 1.0im)
FermiTerm("++I--NE")
# To be clear, note that
FermiTerm("+-INE")
# means ``a_0^\dagger a_1 a^\dagger_3 a_3 a_4 a^\dagger_4``.
#
#
# `OpTerm` is a parametric type with three parameters. The aliases `PauliTerm` and `FermiTerm` each "eat"
# the first parameter.
# The first parameter is an operator type `<:AbstractOp`.
# The second parameter specifies the storage for the string of operators.
# By default, as above, it is a dense `Vector`. For the concrete operator types that we
# have implemented: `Pauli`, `PauliI`, and `FermiOp`, it will be an effcient, packed array of bitstype objects.
# The third type parameter is the type of the coefficient, which can be anything, but should support multiplication
# and addition. This may be, for example, `ComplexF64`, or a symbolic type.

# ## Sparse operators
# Terms with sparse storage of operator strings are supported by using a sparse array type with `OpTerm`.
# We use a generalized (at no runtime cost) version of the standard Julia `SparseVector` that allows one to
# specify that the neutral element of `AbstractOp` in the sparse vector should be the identity rather than zero.
# For example
using QuantumOps: sparse_op, dense_op
sparse_op(OpTerm{Pauli}("IIIIIXIYIZ", 1.0 + 1.0im))
# You can convert back like this
dense_op(sparse_op(OpTerm{Pauli}("IIIIIXIYIZ", 1.0 + 1.0im)))

## Sums of multi-qubit/mode operators -- `OpSum`

# ## Arithemetic on operators

# Multiplication is defined between simple operators, but types are not promoted to
# types capable of representing phase, so the phase is not tracked.
# For example
X * Y

import QuantumOps.FermiOps as FOps
FOps.NumberOp * FOps.NumberOp
#-
FOps.Raise * FOps.Lower

# Multiplying terms does correctly preserve the phase.
t1 = PauliTerm("XIYIZ")
#-
t2 = PauliTerm("YXIZZ")
#-
t1 * t2
# Sparse terms also support `*`
sparse_op(t1) * sparse_op(t2)
# Multiplication between sparse and dense terms is currently not supported.

# ## `OpTerm` features

# You can create a generator of the ``n``-qubit Pauli basis operators like this
collect(QuantumOps.pauli_basis(2))


nothing;
