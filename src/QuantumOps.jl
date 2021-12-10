"""
    QuantumOps

`QuantumOps` provides types and functions that represent quantum operators.
Currently fermionic and Pauli operators are supported. They are organized
in three levels: single particle/mode operators; multi-particle terms with a coefficient;
sums of multi-particle terms.
"""
module QuantumOps

import Requires
import IsApprox
import IsApprox: isunitary
import Random
import LinearAlgebra, Random
using StaticArrays
using DocStringExtensions
import ThreadsX
#import FLoops
import ILog2
import ZChop
import FastBroadcast

import SparseArraysN
import SparseArraysN: neutral, isneutral, neutrals
export neutral, isneutral
export commutes, anticommutes

#import SparseArrays
const SparseArrays = SparseArraysN

export FermiOp, FermiTerm, FermiSum, AbstractPauli, Pauli, PauliI, PauliTerm, PauliSum,
    OpTerm, OpSum # , FermiTermA  #, PauliTermA  , PauliSumA
export unsafe_op
#export APauliTerm
export op_string, op_strings
export count_bodies, jordan_wigner, jordan_wigner_fermi
export Z4Group0, Z4Group, AbstractZ4Group

#export SparseVec
export sparse_op, dense_op

export isunitary
export pauli_basis, mul!, rand_op_term, rand_op_sum
export @pauli_str, @fermi_str
export op_index, phase, weight, pauli_vector
export add!, lmul!, numeric_function
export group_paulis, property_graph
export z4group0
export z4group
export kron_alt

export ⊗

const ⊗ = kron

"""
    _DEFAULT_COEFF = Complex(1, 0)

The default coefficient for `PauliTerm`.
"""
const _DEFAULT_COEFF = Complex(1, 0)

# This is sort of a hack because in Symbolics `==` does not always
# return Bool.
@inline function isbool_and_equal(x, y)
    test = x == y
    return isa(test, Bool) && test
end

function __init__()
    # Add convenience functions if PyCall is loaded.
    # sympy can be used in any case, but this makes it more convenient
    Requires.@require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
        include("pycall.jl")
    end
    Requires.@require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" begin
        include("symbolics.jl")
    end
end

# include("sparse_vec.jl")
# using .SparseVecs

include("util.jl")
include("z4group.jl")

using .Z4Groups

include("z4group0.jl")

using .Z4Group0s

include("abstract_op.jl")

using .AbstractOps

include("abstract_pauli.jl")

include("pauli.jl")
using .Paulis

include("pauli_i.jl")
using .PaulisI

const PauliDefault = Pauli

"""
    PauliDefault

An alias for the default implementation of `AbstractPauli`. Functions that that take
an optional argument that is a subtype of `AbstractPauli` default to `PauliDefault`
if this argument is omitted. There are currently two implementations of `AbstractPauli`:
`Pauli` and `PauliI`.
"""
PauliDefault

include("abstract_term.jl")
include("abstract_sum.jl")
include("op_term.jl")
include("op_sum.jl")
include("abstract_fermi_op.jl")
include("fermi_op.jl")

using .FermiOps

include("op_term_fermi.jl")
include("op_term_pauli.jl")

include("sparse.jl")

include("jordan_wigner.jl")

using .JordanWigner

include("from_interaction_op.jl")

# We are abandoning this
# include("abstract_stabilizers.jl")
# include("stabilizers.jl")

end # module
