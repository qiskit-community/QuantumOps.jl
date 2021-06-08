module AbstractPaulis

using ..AbstractOps
import ..AbstractOps: _show_op_plain
import Random
import SparseArraysN
const SparseArrays = SparseArraysN
using StaticArrays
import LinearAlgebra
import IsApprox
import ..pow_of_minus_one

export AbstractPauli

# """
#     AbstractPauli{T}
# Represents a single-qubit Pauli operator, one
# of `I`, `X`, `Y`, `Z`.
# """
abstract type AbstractPauli{T} <: AbstractOp end

const (Iop, Xop, Yop, Zop) = (0, 1, 2, 3)

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
    pauli_vector(op_index, n_qubits, indices=Vector{Int}(undef, n_qubits))

Return a `Vector` of Pauli's representing a multi-qubit Pauli indexed by the base-4
representation of `op_index`. A temporary array `indices` may be passed to avoid
allocation.
"""
function pauli_vector(::Type{PauliT}, op_index::Integer, n_qubits::Integer,
                     indices=Vector{Int}(undef, n_qubits)) where PauliT
    digits!(indices, op_index; base=4)
    reverse!(indices)  # agree with sort order
    return PauliT.(indices)
end

AbstractOps.rand_ind_range(::Type{<:AbstractPauli}) = 0:3

# function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{PauliT}) where {PauliT <: AbstractPauli}
#     return PauliT(rand(rng, 0:3))
# end

####
#### IO
####

const _pauli_chars = ('I', 'X', 'Y', 'Z')
op_chars(::Type{AbstractPauli}) = _pauli_chars

## I think all of this was due to a bug elsewhere. Perhaps it can be reverted.
## Bizarre, some Int prints as 132, is converted to float -123.0
## and is not > 4, nor < 1
## So, we enumerate the values of valid indices.
function Base.show(io::IO, p::AbstractPauli)
    print(io, typeof(p), ": ")
    _show_op_plain(io, p)
end

function _show_op_plain(io::IO, p::AbstractPauli)
    i::Int = op_index(p) + 1
    if i != 1 && i != 2 && i != 3 && i != 4
        char = '0'
    else
        char = _pauli_chars[i]
    end
    print(io, char)
end

## Use this if `AbstractPauli <: AbstractMatrix`
## We need both definitions.
Base.display(io::IO, p::AbstractPauli) = show(io, p)
Base.display(p::AbstractPauli) = show(stdout, p)

## We may want to define `display` as well.
function Base.show(io::IO, ps::AbstractArray{<:AbstractPauli})
    for p in ps
        _show_op_plain(io, p)
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
    return matrices[op_index(p) + 1]
end

function SparseArrays.sparse(p::AbstractPauli)
    matrices = (_Isparse, _Xsparse, _Ysparse, _Zsparse)
    return matrices[op_index(p) + 1]
end

####
#### Container interface
####

Base.size(::AbstractPauli) = (2, 2)

function smatrix(p::AbstractPauli)
    static = (_SImat, _SXmat, _SYmat, _SZmat)
    return static[op_index(p) + 1]
end

Base.getindex(p::AbstractPauli, i, j) = smatrix(p)[i, j]

####
#### Compare / predicates
####

Base.one(::Type{PauliT}) where PauliT <: AbstractPauli = PauliT(:I)

Base.one(::PauliT) where PauliT <: AbstractPauli = one(PauliT)

Base.:(==)(p1::AbstractPauli, p2::AbstractPauli) = op_index(p1) == op_index(p2)

"""
    isless(p1::AbstractPauli, p2::AbstractPauli)

Canonical (lexical) order of `AbstractPauli` is I < X < Y < Z.
"""
Base.isless(p1::AbstractPauli, p2::AbstractPauli) = op_index(p1) < op_index(p2)

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
        d = op_index(p1) - op_index(p2)
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
    return PauliPhase(has_sign_flip, has_imag_unit)
end

### Linear algebra

LinearAlgebra.tr(::AbstractPauli) = 0

## This would otherwise be computed by LinearAlgebra calling getindex
LinearAlgebra.eigvals(p::AbstractPauli) = isone(p) ? [1.0, 1.0] : [-1.0, 1.0]

Base.length(::AbstractPauli) = 1
Base.iterate(p::AbstractPauli) = (p, nothing)
Base.iterate(::AbstractPauli, ::Any) = nothing
Base.isempty(::AbstractPauli) = false

Base.iszero(::AbstractPauli) = false

abstract type AbstractPhaseData end

struct PauliPhase
    sign_flip::Bool
    imag_unit::Bool
end

mutable struct PauliPhaseCounts <: AbstractPhaseData
    n_sign_flips::Int
    n_imag_units::Int
end

PauliPhaseCounts() = PauliPhaseCounts(0, 0)

AbstractOps.phase_data(::Type{<:AbstractPauli}) = PauliPhaseCounts()

function AbstractOps.accumulate_phase!(p_counts::PauliPhaseCounts, op1::T, op2::T) where {T <: AbstractPauli}
    pauli_phase = phase(op1, op2)
    return AbstractOps.accumulate_phase!(p_counts, pauli_phase)
end

function AbstractOps.accumulate_phase!(p_counts::PauliPhaseCounts, p_phase::PauliPhase)
    p_counts.n_sign_flips += p_phase.sign_flip
    p_counts.n_imag_units += p_phase.imag_unit
    return p_counts
end

function AbstractOps.accumulate_phase(old_phase_data, op1::T, op2::T) where {T <: AbstractPauli}
    phase_info = phase(op1, op2)
    new_phase_data = (old_phase_data[1] + phase_info.sign_flip, old_phase_data[2] + phase_info.imag_unit)
    return new_phase_data
end

AbstractOps.compute_phase(phase_counts::PauliPhaseCounts, _1, _2) =
    AbstractOps.compute_phase(AbstractPauli, (phase_counts.n_sign_flips, phase_counts.n_imag_units))

function AbstractOps.compute_phase(::Type{<:AbstractPauli}, phase_data)
    (n_sign_flips, n_imag_units) = phase_data
    _sign = pow_of_minus_one(n_sign_flips)
    n = n_imag_units % 4
    if n == 0
        cim = complex(1)
    elseif n == 1
        cim = complex(im)
    elseif n == 2
        cim = complex(-1)
    else
        cim = complex(-im)
    end
    return _sign * cim
end

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

end # module AbstractPaulis
