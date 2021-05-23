module QuantumOps

import Requires
import IsApprox
import IsApprox: isunitary
import Random
import LinearAlgebra, Random
using StaticArrays
import ThreadsX
#import FLoops
import ILog2
import SparseArrays

export FermiOp, FermiTerm, FermiSum, AbstractPauli, Pauli, PauliI, PauliTerm, PauliSum
export count_bodies
export Z4Group0, Z4Group, AbstractZ4Group

export isunitary
export pauli_basis, mul!, rand_pauli_term, rand_pauli_sum, rand_fermi_term
export op_index, phase, weight, pauli_vector
export add!, lmul!, numeric_function
export z4group0
export z4group

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

include("util.jl")
include("z4group.jl")

using .Z4Groups

include("z4group0.jl")

using .Z4Group0s

include("abstract_op.jl")
include("abstract_pauli.jl")

include("pauli.jl")
using .Paulis

include("pauli_i.jl")
using .PaulisI

include("abstract_term.jl")
include("pauli_term.jl")
include("abstract_sum.jl")
include("pauli_sum.jl")

include("default_abstract_pauli.jl")

include("abstract_fermi_op.jl")
include("fermi_op.jl")

using .FermiOps

include("fermi_term.jl")
include("fermi_sum.jl")

end # module
