####
#### AbstractPauli
####

export AbstractPauli, pauli_index, phase
import LinearAlgebra, Random

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
LinearAlgebra.ishermitian(::PauliT) where PauliT <: AbstractPauli = true
LinearAlgebra.isposdef(::PauliT) where PauliT <: AbstractPauli = p === PauliT(0)
LinearAlgebra.isdiag(p::PauliT) where PauliT <: AbstractPauli = p === PauliT(0)
LinearAlgebra.issymmetric(p::PauliT) where PauliT <: AbstractPauli = p != PauliT(2)

# TODO: This function is not defined in Julia Base and stdlibs
# We need to decide which package should own it.
isunitary(::PauliT) where PauliT <: AbstractPauli = true

"""
    isless(p1::AbstractPauli, p2::AbstractPauli)

Canonical (lexical) order of `AbstractPauli` is I < X < Y < Z.
"""
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

Multiply Paulis returning a  `NamedTuple` with members
`:bare_pauli`, `:has_sign_flip`, `:has_imag_unit`.
See `PauliString.phase`.
"""
function mul(p1::AbstractPauli, p2::AbstractPauli)
    bare_pauli = p1 * p2
    _phase = phase(p1, p2)
    return((bare_pauli=bare_pauli, has_sign_flip=_phase[:has_sign_flip], has_imag_unit=_phase[:has_imag_unit]))
end

function multiply_keeping_phase(s1::AbstractArray{<:AbstractPauli}, s2::AbstractArray{<:AbstractPauli})
    length(s1) != length(s2) && throw(DimensionMismatch())
    s_out = similar(s1)
    cum_sign_flips = 0
    cum_imag_units = 0
    @inbounds for i in eachindex(s1)
        product = mul(s1[i], s2[i])
        cum_sign_flips += product[:has_sign_flip]
        cum_imag_units += product[:has_imag_unit]
        s_out[i] = product[:bare_pauli]
    end
    return(s_out, (-1) ^ cum_sign_flips * im ^ cum_imag_units)
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{PauliT}) where {PauliT <: AbstractPauli}
    return PauliT(rand(rng, 0:3))
end

####
#### Subtypes of AbstractPauli must implement these
####

# The function `*` must return only the bare Pauli. That is the phase +-im must
# be discarded

"""
    p1::AbstractPauli * p2::AbstractPauli

Multiply single-qubit operators ignoring phase.
"""
function Base.:*(p1::AbstractPauli, p2::AbstractPauli)
    throw(MethodError(*, (p1, p2)))
end

"""
    phase(p1::AbstractPauli, p2::AbstractPauli)

Return a `NamedTuple` of two `Bool`s representing the phase
of the product of `p1` and `p2`. The elements are `:has_sign_flip`
and `:has_imag_unit`. The phase can be reconstructed as
`im^has_imag_unit * (-1)^has_sign_flip`.
"""
function phase end

"""
    pauli_index(p::AbstractPauli)::Int

Return the Pauli index in `[0,3]` of `p`.
"""
function pauli_index end

"""
    abstract type AbstractPauli

Represents single-qubit Pauli matrices and their algebra.

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
