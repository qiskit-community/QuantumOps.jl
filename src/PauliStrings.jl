module PauliStrings

using Requires: @require

export ⊗
const ⊗ = kron

function __init__()
    # Add convenience functions if PyCall is loaded.
    # sympy can be used in any case, but this makes it more convenient
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
        include("pycall.jl")
    end
    @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" begin
        include("symbolics.jl")
    end
end

include("abstract_pauli.jl")
include("pauli_term.jl")
include("pauli_sum.jl")
include("pauli.jl")
include("default_types.jl")
include("iphase.jl")

end # module
