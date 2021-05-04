####
#### Convenience methods making `Pauli` the default `AbstractPauli
####

PauliTerm(s::AbstractString, coeff=_default_coeff) = PauliTerm(Pauli, s, coeff)
PauliTerm(s::Symbol, coeff=_default_coeff) = PauliTerm(Pauli, s, coeff)

"""
    PauliTerm(inds::AbstractVector{<:Integer}, coeff=_default_coeff)

Initialize a `PauliTerm` from a vector of integers in `[0, 3]`.
"""
PauliTerm(inds::AbstractVector{<:Integer}, coeff=_default_coeff) = PauliTerm(Pauli.(inds), coeff)

Base.rand(::Type{PauliTerm}, n::Integer) = rand(PauliTerm{Pauli}, n)

"""
    PauliSum()

Return an empty `PauliSum` with Paulis of type `Pauli`
and `Complex{Float64}` coefficients.
"""
PauliSum() = PauliSum(Vector{Vector{Pauli}}[], Complex{Float64}[])
