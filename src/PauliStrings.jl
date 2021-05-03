module PauliStrings

using Requires: @require

export ⊗
const ⊗ = kron

function __init__()
    # Add a convenience function if PyCall is loaded before PauliStrings.
    # sympy can be used in any case, but this makes it more convenient
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
        Base.:*(z::PyCall.PyObject, ps::PauliTerm) = PauliTerm(ps.s, ps.coeff * z)
        Base.:*(ps::PauliTerm, z::PyCall.PyObject) = *(z, ps)
    end
end

include("abstract_pauli.jl")
include("pauli_term.jl")
include("pauli_sum.jl")
include("pauli.jl")
include("default_types.jl")

end # module
