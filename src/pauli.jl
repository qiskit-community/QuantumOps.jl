export Pauli

####
#### Pauli
####

struct Pauli <: AbstractPauli
    hi::Bool
    lo::Bool
end

pauli_index(p::Pauli) = 2 * p.hi + p.lo

# For convenience. Must be explicitly imported
const I = Pauli(0, 0)
const X = Pauli(0, 1)
const Y = Pauli(1, 0)
const Z = Pauli(1, 1)

function Pauli(ind::Integer)
    if ind == 0
        return Pauli(0, 0)
    elseif ind == 1
        return Pauli(0, 1)
    elseif ind == 2
        return Pauli(1, 0)
    elseif ind == 3
        return Pauli(1, 1)
    else
        throw(ArgumentError("Invalid Pauli index $ind"))
    end
end

Pauli(s::Union{Symbol, AbstractString, AbstractChar}) = _pauli(Pauli, s)

"""
    p1::Pauli * p2::Pauli

Multiply single-qubit operators ignoring phase.
"""
Base.:*(p1::Pauli, p2::Pauli) = Pauli(p1.hi ⊻ p2.hi, p1.lo ⊻ p2.lo)

function phase(p1::Pauli, p2::Pauli)
    if isone(p1) || isone(p2)
        has_imag_unit = false
        has_sign_flip = false
    else
        d = pauli_index(p1) - pauli_index(p2)
        if d == 0
            has_imag_unit = false
            has_sign_flip = false
        elseif d == 1 || d == -2
            has_imag_unit = true
            has_sign_flip = true
        else
            has_imag_unit = true
            has_sign_flip = false
        end
    end
    return (has_sign_flip=has_sign_flip, has_imag_unit=has_imag_unit)
end

####
#### Convenience methods making `Pauli` the default `AbstractPauli
####

PauliString(s::AbstractString, coeff=Complex(1, 0)) = PauliString(Pauli, s, coeff)
Base.rand(::Type{PauliString}, n::Integer) = rand(PauliString{Pauli}, n)
