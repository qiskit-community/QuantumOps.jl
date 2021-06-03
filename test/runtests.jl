using QuantumOps
using Test
import SparseArraysN
const SparseArrays = SparseArraysN
using QuantumOps: Paulis, FermiOps

include("./test_sparse.jl")
include("./test_jordan_wigner.jl")
# include("./test_ind_op.jl")
include("./test_op_term_sum.jl")
include("./test_pauli.jl")
include("./test_fermi.jl")


# For testing. Sometimes twice slower than constructing all at once
# Testing that this gives the same result as standard construction
# is a good idea
# TODO: move this to test suite
# function make_pauli_sum(strings)
#     s = PauliSum([first(strings)])
#     for i in 2:length(strings)
#         add!(s, strings[i])
#     end
#     return s
# end
