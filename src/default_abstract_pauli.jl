####
#### Convenience methods making one subtype of `AbstractPauli` the default.
####

## Choose one of the subtypes of `AbstractPauli`; `Pauli`, or `PauliI`.
## Reducing a large array with `*` is 5x faster with Pauli than PauliI

const PauliDefault = Pauli

####
#### PauliTerm
####

PauliTerm(s::AbstractString, coeff=_DEFAULT_COEFF) = PauliTerm(PauliDefault, s, coeff)
PauliTerm(s::Symbol, coeff=_DEFAULT_COEFF) = PauliTerm(PauliDefault, s, coeff)

"""
    PauliTerm(inds::AbstractVector{<:Integer}, coeff=_DEFAULT_COEFF)

Initialize a `PauliTerm` from a vector of integers in `[0, 3]`.
"""
PauliTerm(inds::AbstractVector{<:Integer}, coeff=_DEFAULT_COEFF) = PauliTerm(PauliDefault.(inds), coeff)

PauliTerm(index::Integer, n_paulis::Integer, coeff=_DEFAULT_COEFF) = PauliTerm(PauliDefault, index, n_paulis, coeff)

Base.rand(::Type{PauliTerm}, n::Integer; coeff=_DEFAULT_COEFF) = rand(PauliTerm{PauliDefault}, n; coeff=coeff)

## This works. But, see notes in pauli_term.jl for why it is not satisfactory
# function PauliTermSampler(npaulis::Integer, coeff_func=() -> _DEFAULT_COEFF)
#     return PauliTermSampler{PauliDefault,typeof(coeff_func)}(npaulis, coeff_func)
# end

####
#### PauliSum
####

PauliSum(strings::AbstractVector{<:AbstractString}, coeffs) = PauliSum(Vector{PauliDefault}.(strings), coeffs)

"""
    PauliSum()

Return an empty `PauliSum` with Paulis of type `PauliDefault`
and `Complex{Float64}` coefficients.
"""
PauliSum() = PauliSum(PauliDefault)

PauliSum(m::AbstractMatrix{<:Number}) = PauliSum(PauliDefault, m)

function pauli_vector(pauli_index::Integer, n_qubits::Integer, indices=Vector{Int}(undef, n_qubits))
    return pauli_vector(PauliDefault, pauli_index, n_qubits, indices)
end

"""
    pauli_basis(n_qubits; coeff=_DEFAULT_COEFF)

Return an iterator over all `PauliTerm{PauliDefault}`s of `n_qubits`.
"""
pauli_basis(n_qubits::Integer; coeff=_DEFAULT_COEFF) = pauli_basis(PauliDefault, n_qubits; coeff=coeff)
