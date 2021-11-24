using QuantumOps
using Test
using StaticArrays
import SparseArraysN
import LinearAlgebra
const SparseArrays = SparseArraysN
using QuantumOps: Paulis, FermiOps

include("./test_pauli.jl")
include("./test_jordan_wigner.jl")
include("./test_op_term_sum.jl")
include("./test_fermi.jl")
include("./test_sparse.jl")
include("./test_z4group.jl")

# DISABLE this because it confuses PackageCompiler.
# @testset "quantum_ops_intro" begin
#     include("../examples/quantum_ops_intro.jl")
# end

# @testset "jordan_wigner_example" begin
#     include("../examples/jordan_wigner_example.jl")
# end
