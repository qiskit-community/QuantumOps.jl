####
#### Convenience methods making `Pauli` the default `AbstractPauli
####

PauliTerm(s::AbstractString, coeff=Complex(1, 0)) = PauliTerm(Pauli, s, coeff)
Base.rand(::Type{PauliTerm}, n::Integer) = rand(PauliTerm{Pauli}, n)

"""
    PauliSum()

Return an empty `PauliSum` with Paulis of type `Pauli`
and `Complex{Float64}` coefficients.
"""
PauliSum() = PauliSum(Vector{Vector{Pauli}}[], Complex{Float64}[])
