export Pauli, PauliString, pauli_index, phase


####
#### AbstractPauli
####

abstract type AbstractPauli end

function _pauli(::Type{PauliT}, s::Symbol) where PauliT
    if s == :I
        return PauliT(0)
    elseif s == :X
        return PauliT(1)
    elseif s == :Y
        return PauliT(2)
    elseif s == :Z
        return PauliT(3)
    else
        throw(ArgumentError("Invalid Pauli symbol $s"))
    end
end

function _pauli(::Type{PauliT}, s::AbstractString) where PauliT
    if s == "I"
        return PauliT(0)
    elseif s == "X"
        return PauliT(1)
    elseif s == "Y"
        return PauliT(2)
    elseif s == "Z"
        return PauliT(3)
    else
        throw(ArgumentError("Invalid Pauli symbol \"$s\""))
    end
end

function _pauli(::Type{PauliT}, s::AbstractChar) where PauliT
    if s == 'I'
        return PauliT(0)
    elseif s == 'X'
        return PauliT(1)
    elseif s == 'Y'
        return PauliT(2)
    elseif s == 'Z'
        return PauliT(3)
    else
        throw(ArgumentError("Invalid Pauli Char '$s'"))
    end
end

Base.one(::Type{PauliT}) where PauliT <: AbstractPauli = PauliT(0)
Base.one(::PauliT) where PauliT <: AbstractPauli = PauliT(0)

Vector{T}(ps::AbstractString) where {T <: AbstractPauli} = [T(s) for s in ps]

const _pauli_chars = ('I', 'X', 'Y', 'Z')

Base.show(io::IO, p::AbstractPauli) = print(io, _pauli_chars[pauli_index(p) + 1])

# May want to define `display` as well
function Base.show(io::IO, ps::AbstractArray{<:AbstractPauli})
    for p in ps
        show(io, p)
    end
end

Base.:^(p::AbstractPauli, n::Integer) = iseven(n) ? typeof(p)(:I) : p


####
#### PauliString
####

struct PauliString{W<:AbstractPauli, T<:AbstractVector{W}, V}
    s::T
    coeff::V
end

PauliString(s) = PauliString(s, Complex(1, 0))
function PauliString(::Type{T}, s::AbstractString, coeff=Complex(1, 0)) where T <: AbstractPauli
    return PauliString(Vector{T}(s), coeff)
end


# Forward some functions to the array containing the Pauli string
for func in (:length, :size, :eachindex, :insert!, :push!, :popat!, :splice!, :eltype, :getindex, :setindex!)
    @eval begin
        function Base.$func(ps::PauliString, args...)
            return $func(ps.s, args...)
        end
    end
end

Base.:(==)(ps1::PauliString, ps2::PauliString) = ps1.coeff == ps2.coeff && ps1.s == ps2.s

# TODO:
# isless

function Base.show(io::IO, ps::PauliString)
    print(io, ps.coeff, " * ")
    print(io, ps.s)
end

function Base.:*(ps1::PauliString, ps2::PauliString)
    s_new, phase = ps1.s * ps2.s
    return PauliString(s_new, ps1.coeff * ps2.coeff * phase)
end

Base.:*(z::Number, ps::PauliString) = PauliString(ps.s, ps.coeff * z)
Base.:*(ps::PauliString, z::Number) = z * ps

Base.kron(ps1::PauliString, ps2::PauliString) = PauliString(vcat(ps1.s, ps2.s), ps1.coeff * ps2.coeff)

Base.one(ps::PauliString{T}) where T = PauliString(fill(one(T), length(ps)))

function Base.rand(::Type{<:PauliString{T}}, n::Integer) where {T <: AbstractPauli}
    return PauliString([T(i) for i in rand(0:3, n)])
end

####
#### Pauli
####

struct Pauli <: AbstractPauli
    hi::Bool
    lo::Bool
end

"""
    pauli_index(p::Pauli)::Int

Return the Pauli index in `[0,3]` of `p`.
"""
pauli_index(p::Pauli) = 2 * p.hi + p.lo

const _I = Pauli(0, 0)
const _X = Pauli(0, 1)
const _Y = Pauli(1, 0)
const _Z = Pauli(1, 1)

function Pauli(ind::Integer)
    if ind == 0
        return _I
    elseif ind == 1
        return _X
    elseif ind == 2
        return _Y
    elseif ind == 3
        return _Z
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
    if p1 == one(p1) || p2 == one(p1)
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

"""
    mul(p1::Pauli, p2::Pauli)

Multiply Paulis returning a `Tuple` of the bare product and the phase.
See `PauliString.phase`.
"""
mul(p1::Pauli, p2::Pauli) = (p1 * p2, phase(p1, p2))

function Base.:*(s1::AbstractArray{<:Pauli}, s2::AbstractArray{<:Pauli})
    length(s1) != length(s2) && throw(DimensionMismatch())
    s_out = similar(s1)
    cum_sign_flips = 0
    cum_imag_units = 0
    @inbounds for i in eachindex(s1)
        (product, phase) = mul(s1[i], s2[i])
        cum_sign_flips += phase[:has_sign_flip]
        cum_imag_units += phase[:has_imag_unit]
        s_out[i] = product
    end
    return(s_out, (-1) ^ cum_sign_flips * im ^ cum_imag_units)
end

####
#### Convenience methods making `Pauli` the default `AbstractPauli
####

PauliString(s::AbstractString, coeff=Complex(1, 0)) = PauliString(Pauli, s, coeff)
Base.rand(::Type{PauliString}, n::Integer) = rand(PauliString{Pauli}, n)
