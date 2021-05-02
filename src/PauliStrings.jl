module PauliStrings

export ⊗
const ⊗ = kron

include("abstract_pauli.jl")
include("pauli_string.jl")
include("pauli.jl")

end # module
