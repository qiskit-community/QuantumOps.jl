####
#### Convenience methods making one subtype of `AbstractPauli` the default.
####

## Choose one of the subtypes of `AbstractPauli`; `Pauli`, or `PauliI`.
## Reducing a large array with `*` is 5x faster with Pauli than PauliI

const PauliDefault = Pauli

"""
    PauliDefault

An alias for the default implementation of `AbstractPauli`. Functions that that take
an optional argument that is a subtype of `AbstractPauli` default to `PauliDefault`
if this argument is omitted. There are currently two implementations of `AbstractPauli`:
`Pauli` and `PauliI`.
"""
PauliDefault

####
#### PauliTermA
####

PauliTermA(s::AbstractString, coeff=_DEFAULT_COEFF) = PauliTermA(PauliDefault, s, coeff)
PauliTermA(s::Symbol, coeff=_DEFAULT_COEFF) = PauliTermA(PauliDefault, s, coeff)

"""
    PauliTermA(inds::AbstractVector{<:Integer}, coeff=_DEFAULT_COEFF)

Construct a `PauliTermA` from a vector of integers in `[0, 3]`.
"""
PauliTermA(inds::AbstractVector{<:Integer}, coeff=_DEFAULT_COEFF) = PauliTermA(PauliDefault, inds, coeff)

PauliTermA(index::Integer, n_paulis::Integer, coeff=_DEFAULT_COEFF) = PauliTermA(PauliDefault, index, n_paulis, coeff)

Base.rand(::Type{PauliTermA}, n::Integer; coeff=_DEFAULT_COEFF) = rand(PauliTermA{PauliDefault}, n; coeff=coeff)

## This works. But, see notes in pauli_term.jl for why it is not satisfactory
# function PauliTermASampler(npaulis::Integer, coeff_func=() -> _DEFAULT_COEFF)
#     return PauliTermASampler{PauliDefault,typeof(coeff_func)}(npaulis, coeff_func)
# end

####
#### PauliSumA
####

PauliSumA(strings::AbstractVector{<:AbstractString}, coeffs) = PauliSumA(Vector{PauliDefault}.(strings), coeffs)

"""
    PauliSumA()

Return an empty `PauliSumA` with Paulis of type `PauliDefault`
and `Complex{Float64}` coefficients.
"""
PauliSumA() = PauliSumA(PauliDefault)

PauliSumA(m::AbstractMatrix{<:Number}) = PauliSumA(PauliDefault, m)

function pauli_vector(pauli_index::Integer, n_qubits::Integer, indices=Vector{Int}(undef, n_qubits))
    return pauli_vector(PauliDefault, pauli_index, n_qubits, indices)
end

"""
    pauli_basis(n_qubits; coeff=_DEFAULT_COEFF)

Return an iterator over all `PauliTermA{PauliDefault}`s of `n_qubits`.
"""
pauli_basis(n_qubits::Integer; coeff=_DEFAULT_COEFF) = pauli_basis(PauliDefault, n_qubits; coeff=coeff)
