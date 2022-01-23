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
using IsApprox: isunitary, ishermitian, commutes,
    isposdef, isdiag, issymmetric, isunitary
using IsApprox: AbstractApprox, Equal, Approx
import Random
import LinearAlgebra, Random
using StaticArrays
import ThreadsX
import ILog2
import ZChop

import SparseArraysN
import SparseArraysN: neutral, isneutral, neutrals
export neutral, isneutral
export commutes
# TODO: find if and how to export following
# export anticommutes

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
export op_index, weight, pauli_vector
export add!, numeric_function
# TODO: export following ?
# export lmul!
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
using .AbstractOps: AbstractOps # trying to satisfy JET

include("abstract_pauli.jl")
using .AbstractPaulis

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

# The following does *not* make QuantumOps load faster.
# In fact in makes the following twice as slow on the first call:
# PauliSum(Matrix(rand_op_sum(Pauli, 3, 4)))
# let
#     while true
#         group_paulis(rand_op_sum(Pauli, 10, 20))
#         rand_op_sum(FermiOp, 10, 20)
#         rand_op_term(Pauli, 10)
#         rand_op_term(FermiOp, 10)

#         h2_hamiltonian = FermiSum([FermiTerm("IIII", 0.7137539936876182),FermiTerm("IIIN", -0.47594871522096355),
#                                    FermiTerm("IINI", -0.47594871522096355),FermiTerm("IINN", 0.6973937674230275),
#                                    FermiTerm("INII", -1.2524635735648981),FermiTerm("ININ", 0.482179288212072),
#                                    FermiTerm("INNI", 0.663468096423568),FermiTerm("NIII", -1.2524635735648981),
#                                    FermiTerm("NIIN", 0.663468096423568),FermiTerm("NINI", 0.482179288212072),
#                                    FermiTerm("NNII", 0.6744887663568375),FermiTerm("++--", -0.181288808211496),
#                                    FermiTerm("+--+", 0.181288808211496),FermiTerm("-++-", 0.181288808211496),
#                                    FermiTerm("--++", -0.181288808211496)])
#         jordan_wigner(h2_hamiltonian)
#         PauliSum(Matrix(rand_op_sum(Pauli, 3, 4)))
#         println("Done chores")
#         break
#     end
# end

# We are abandoning this
# include("abstract_stabilizers.jl")
# include("stabilizers.jl")

end # module
