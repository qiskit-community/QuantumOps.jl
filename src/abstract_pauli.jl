####
#### AbstractPauli
####

export AbstractPauli, pauli_index, phase

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
Base.one(::PauliT) where PauliT <: AbstractPauli = one(PauliT)

Base.isless(p1::AbstractPauli, p2::AbstractPauli) = pauli_index(p1) < pauli_index(p2)

"""
    Vector{<:AbstractPauli}(ps::AbstractString)

Construct a `Vector{<:AbstractPauli}` by parsing `ps` containing
characters I, X, Y, Z.
"""
Vector{T}(ps::AbstractString) where {T <: AbstractPauli} = [T(s) for s in ps]

const _pauli_chars = ('I', 'X', 'Y', 'Z')

Base.show(io::IO, p::AbstractPauli) = print(io, _pauli_chars[pauli_index(p) + 1])

# May want to define `display` as well
function Base.show(io::IO, ps::AbstractArray{<:AbstractPauli})
    for p in ps
        show(io, p)
    end
end

Base.:^(p::AbstractPauli, n::Integer) = iseven(n) ? one(p) : p

"""
    mul(p1::AbstractPauli, p2::AbstractPauli)

Multiply Paulis returning a `Tuple` of the bare product and the phase.
See `PauliString.phase`.
"""
mul(p1::AbstractPauli, p2::AbstractPauli) = (p1 * p2, phase(p1, p2))

function Base.:*(s1::AbstractArray{<:AbstractPauli}, s2::AbstractArray{<:AbstractPauli})
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

#### Subtypes of AbstractPauli must implement these

"""
    pauli_index(p::AbstractPauli)::Int

Return the Pauli index in `[0,3]` of `p`.
"""
function pauli_index end

"""
    AbstractPauli(ind::Integer)::AbstractPauli

This function must be called as a concrete subtype of `AbstractPauli`.

Create an instance of (a subtype of) `AbstractPauli` given
the Pauli index in `[0, 3]`.

    AbstractPauli(s::Union{Symbol, AbstractString, AbstractChar})::AbstractPauli

Create an instance of (a subtype of) `AbstractPauli` from the `Symbol`s `:I, :X, :Y, :Z`,
or the `AbstractString`s `"I", "X", "Y", "Z"` or the `AbstractChar`s `'I', 'X', 'Y', 'Z'`.
In general creating `AbstractPauli`s with `String`s is slower than with the other types.
"""
function AbstractPauli end

# The second method above can be implemented like this.
# Pauli(s::Union{Symbol, AbstractString, AbstractChar}) = _pauli(Pauli, s)
