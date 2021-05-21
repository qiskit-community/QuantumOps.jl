#abstract type AbstractPauli{T} <: AbstractMatrix{T} end
abstract type AbstractPauli{T} end

####
#### Constructors
####

## These are not really constructors, since they don't exist for abstract types.
## But, they are the closest we can get, and are called by the concrete subtypes.

## The following methods are called by concrete subtypes like this:
## Pauli(s::Union{Symbol, AbstractString, AbstractChar}) = _AbstractPauli(Pauli, s)

function _AbstractPauli(::Type{PauliT}, s::Symbol) where PauliT
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

function _AbstractPauli(::Type{PauliT}, s::AbstractString) where PauliT
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

function _AbstractPauli(::Type{PauliT}, s::AbstractChar) where PauliT
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

"""
    Vector{<:AbstractPauli}(ps::AbstractString)

Construct a `Vector{<:AbstractPauli}` by parsing `ps` containing
characters I, X, Y, Z.
"""
Vector{T}(ps::AbstractString) where {T <: AbstractPauli} = [T(s) for s in ps]

"""
    pauli_vector(pauli_index, n_qubits, indices=Vector{Int}(undef, n_qubits))

Return a `Vector` of Pauli's representing a multi-qubit Pauli indexed by the base-4
representation of `pauli_index`. A temporary array `indices` may be passed to avoid
allocation.
"""
function pauli_vector(::Type{PauliT}, pauli_index::Integer, n_qubits::Integer,
                     indices=Vector{Int}(undef, n_qubits)) where PauliT
    digits!(indices, pauli_index; base=4)
    reverse!(indices)  # agree with sort order
    return PauliT.(indices)
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{PauliT}) where {PauliT <: AbstractPauli}
    return PauliT(rand(rng, 0:3))
end

####
#### IO
####

function Base.show(io::IO, p::AbstractPauli)
    _pauli_chars = ('I', 'X', 'Y', 'Z')
    print(io, _pauli_chars[pauli_index(p) + 1])
end

## Use this if `AbstractPauli <: AbstractMatrix`
## We need both definitions.
Base.display(io::IO, p::AbstractPauli) = show(io, p)
Base.display(p::AbstractPauli) = show(stdout, p)

## We may want to define `display` as well.
function Base.show(io::IO, ps::AbstractArray{<:AbstractPauli})
    for p in ps
        show(io, p)
    end
end

####
#### Conversion
####

## Use static matrices when constructing large matrices in loops.
## But, I don't seem to get any performance gain.
const _Imat = float.([1 0; 0 1])
const _Xmat = float.([0 1; 1 0])
const _Ymat = float.([0 -im; im 0])
const _Zmat = float.([1 0; 0 -1])

const _Isparse = SparseArrays.sparse(_Imat)
const _Xsparse = SparseArrays.sparse(_Xmat)
const _Ysparse = SparseArrays.sparse(_Ymat)
const _Zsparse = SparseArrays.sparse(_Zmat)

const _SImat = @SMatrix [1.0 0; 0 1]
const _SXmat = @SMatrix [0.0 1; 1 0]
const _SYmat = @SMatrix [0.0 -im; im 0]
const _SZmat = @SMatrix [1.0 0; 0 -1]

function Base.Matrix(p::AbstractPauli)
    matrices = (_Imat, _Xmat, _Ymat, _Zmat)
    return matrices[pauli_index(p) + 1]
end

function SparseArrays.sparse(p::AbstractPauli)
    matrices = (_Isparse, _Xsparse, _Ysparse, _Zsparse)
    return matrices[pauli_index(p) + 1]
end

####
#### Container interface
####

Base.size(::AbstractPauli) = (2, 2)

function smatrix(p::AbstractPauli)
    static = (_SImat, _SXmat, _SYmat, _SZmat)
    return static[pauli_index(p) + 1]
end

Base.getindex(p::AbstractPauli, i, j) = smatrix(p)[i, j]

####
#### Compare / predicates
####

Base.one(::Type{PauliT}) where PauliT <: AbstractPauli = PauliT(:I)

Base.one(::PauliT) where PauliT <: AbstractPauli = one(PauliT)

Base.:(==)(p1::AbstractPauli, p2::AbstractPauli) = pauli_index(p1) == pauli_index(p2)

"""
    isless(p1::AbstractPauli, p2::AbstractPauli)

Canonical (lexical) order of `AbstractPauli` is I < X < Y < Z.
"""
Base.isless(p1::AbstractPauli, p2::AbstractPauli) = pauli_index(p1) < pauli_index(p2)

LinearAlgebra.ishermitian(::AbstractPauli) = true

LinearAlgebra.isposdef(::PauliT) where PauliT <: AbstractPauli = p === PauliT(:I)

LinearAlgebra.isdiag(p::PauliT) where PauliT <: AbstractPauli = (p === PauliT(:I) || p === PauliT(:Z))

LinearAlgebra.issymmetric(p::PauliT) where PauliT <: AbstractPauli = p != PauliT(:Y)

IsApprox.isunitary(::AbstractPauli) = true

####
#### Algebra / mathematical operations
####

Base.:*(p::AbstractPauli, m::AbstractMatrix) = smatrix(p) * m

Base.:*(m::AbstractMatrix, p::AbstractPauli) = m * smatrix(p)

Base.:^(p::AbstractPauli, n::Integer) = iseven(n) ? one(p) : p

Base.inv(p::AbstractPauli) = p

"""
    mul(p1::AbstractPauli, p2::AbstractPauli)

Multiply Paulis returning a  `NamedTuple` with members
`:bare_pauli`, `:has_sign_flip`, `:has_imag_unit`.

See `QuantumOps.phase`.
"""
function mul(p1::AbstractPauli, p2::AbstractPauli)
    bare_pauli = p1 * p2
    _phase = phase(p1, p2)
    return((bare_pauli=bare_pauli, has_sign_flip=_phase[:has_sign_flip], has_imag_unit=_phase[:has_imag_unit]))
end

function multiply_keeping_phase(s1::AbstractArray{<:AbstractPauli}, s2::AbstractArray{<:AbstractPauli})
    return multiply_keeping_phase!(similar(s1), s1, s2)
end

function multiply_keeping_phase!(target::AbstractArray{<:AbstractPauli},
                                s1::AbstractArray{<:AbstractPauli}, s2::AbstractArray{<:AbstractPauli})
    length(s1) != length(s2) && throw(DimensionMismatch())
    cum_sign_flips = 0
    cum_imag_units = 0
    @inbounds for i in eachindex(s1)
        product = mul(s1[i], s2[i])
        cum_sign_flips += product[:has_sign_flip]
        cum_imag_units += product[:has_imag_unit]
        target[i] = product[:bare_pauli]
    end
    sign = iseven(cum_sign_flips) ? 1 : -1
    cum_imag_units = mod(cum_imag_units, 4)
    return(target, sign * im ^ cum_imag_units)
end

"""
    phase(p1::AbstractPauli, p2::AbstractPauli)

Return a `NamedTuple` of two `Bool`s representing the phase of the product of `p1` and
`p2`. The elements are `:has_sign_flip` and `:has_imag_unit`. The phase can be reconstructed
as `im^has_imag_unit * (-1)^has_sign_flip`.
"""
function phase(p1::AbstractPauli, p2::AbstractPauli)
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

### Linear algebra

LinearAlgebra.tr(::AbstractPauli) = 0

## This would otherwise be computed by LinearAlgebra calling getindex
LinearAlgebra.eigvals(::AbstractPauli) = [-1.0, 1.0]

####
#### Other
####

"""
    weight(v::AbstractArray{<:AbstractPauli})
    weight(ps::PauliTerm)

Count the number of Paulis in the string that are not the identity.
"""
weight(v::AbstractArray{<:AbstractPauli}) = count(pauli -> pauli_index(pauli) != 0, v)

####
#### Interface for subtypes
####

## Subtypes of AbstractPauli must implement these.

## The function `*` must return only the bare Pauli. That is, the phase +-im must be discarded.

"""
    p1::AbstractPauli * p2::AbstractPauli

Multiply single-qubit operators ignoring phase.
"""
function Base.:*(p1::AbstractPauli, p2::AbstractPauli)
    throw(MethodError(*, (p1, p2)))
end

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
