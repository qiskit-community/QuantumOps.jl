export Pauli

# TODO: Make this a `module` to create a namespace

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

Pauli(s::Union{Symbol, AbstractString, AbstractChar}) = _AbstractPauli(Pauli, s)

"""
    *(p1::Pauli, p2::Pauli) 

Multiplication that returns a `Pauli`, but ignores the phase.
"""
Base.:*(p1::Pauli, p2::Pauli) = Pauli(p1.hi ⊻ p2.hi, p1.lo ⊻ p2.lo)
