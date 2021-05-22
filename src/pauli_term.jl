"""
    struct PauliTerm{W<:AbstractPauli, T<:AbstractVector{W}, V}

Represents a Pauli string (tensor product of Paulis) with a coefficient.
"""
struct PauliTerm{W<:AbstractPauli, T<:AbstractVector{W}, V} <: AbstractTerm{W}
    paulis::T
    coeff::V
end

op_string(t::PauliTerm) = t.paulis

term_type(::Type{<:AbstractPauli}) = PauliTerm

####
#### Constructors
####

"""
    PauliTerm(s)

Construct a `PauliTerm` with default coefficient.
"""
PauliTerm(s) = PauliTerm(s, _DEFAULT_COEFF)

"""
    PauliTerm(::Type{T}=PauliDefault, s::AbstractString, coeff=_DEFAULT_COEFF) where T <: AbstractPauli

Construct a `PauliTerm`, where `s` is of the form "XYZ", etc. If `::Type{T}` is
ommited the default implementation of `AbstractPauli` is used.

# Examples
```jldoctest
julia> PauliTerm("XX")
(1 + 0im) * XX

julia> PauliTerm(Pauli, "IXYZ")
(1 + 0im) * IXYZ

julia> PauliTerm(Pauli, "IXYZ", 2.0)
2.0 * IXYZ
```
"""
function PauliTerm(::Type{T}, s::AbstractString, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTerm(Vector{T}(s), coeff)
end

"""
    PauliTerm(::Type{T}=PauliDefault, s::Symbol, coeff=_DEFAULT_COEFF) where T <: AbstractPauli

Construct a `PauliTerm`, where `s` is of the form `:XYZ`, etc.

# Examples
```jldoctest
julia> PauliTerm(Pauli, :IZXY)
(1 + 0im) * IZXY
```
"""
function PauliTerm(::Type{T}, s::Symbol, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTerm(Vector{T}(String(s)), coeff)
end

function PauliTerm(::Type{T}, index::Integer, n_paulis::Integer, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTerm(pauli_vector(T, index, n_paulis), coeff)
end

function PauliTerm(paulis::AbstractPauli...; coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTerm([paulis...], coeff)
end

function PauliTerm(::Type{T}, inds::AbstractVector{<:Integer}, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTerm(T.(inds), coeff)
end

"""
    rand_pauli_term(::Type{PauliT}=PauliDefault, n::Integer; coeff=_DEFAULT_COEFF) where {PauliT <: AbstractPauli}

Return a random `PauliTerm` of `n` tensor factors.
"""
function rand_pauli_term(::Type{PauliT}, n::Integer; coeff=_DEFAULT_COEFF) where {PauliT <: AbstractPauli}
    return PauliTerm(rand(PauliT, n), coeff)
end

rand_pauli_term(n::Integer; coeff=_DEFAULT_COEFF) = rand_pauli_term(PauliDefault, n, coeff=coeff)

####
#### Conversion
####

function _multiply_coefficient(coeff, matrix)
    if isbool_and_equal(coeff, one(eltype(matrix)))
        return matrix
    end
    if coeff isa eltype(matrix)
        return LinearAlgebra.lmul!(coeff, matrix)
    end
    coeff1 = convert(promote_type(typeof(coeff), eltype(matrix)), coeff)
    return coeff1 .* matrix
end

Base.Matrix(pt::PauliTerm) = Matrix(Float64, pt)

# FIXME: Not type stable.
function Base.Matrix(::Type{Z4Group0}, pt::PauliTerm)
    matrix = _kron((Z4Group0.(m) for m in Matrix.(pt.paulis))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

function Base.Matrix(::Type{Float64}, pt::PauliTerm)
    matrix = _kron(Matrix.(pt.paulis)...)
    return _multiply_coefficient(pt.coeff, matrix)
end

SparseArrays.sparse(pt::PauliTerm) = SparseArrays.sparse(Float64, pt)

function SparseArrays.sparse(::Type{Float64}, pt::PauliTerm)
    matrix = _kron(SparseArrays.sparse.(pt.paulis)...)
    return _multiply_coefficient(pt.coeff, matrix)
end

# FIXME: broken, should throw inexact error earlier rather than return wrong type
function SparseArrays.sparse(::Type{Z4Group0}, pt::PauliTerm)
    matrix = _kron((Z4Group0.(m) for m in SparseArrays.sparse.(pt.paulis))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

####
#### IO
####

function Base.show(io::IO, ps::PauliTerm{T,V,Z4Group}) where {T,V}
    print(io, ps.coeff, " ")
    print(io, ps.paulis)
end

####
#### Compare / predicates
####

"""
    isunitary(pt::PauliTerm)

Return `true` if `pt` is a unitary operator.
"""
IsApprox.isunitary(pt::PauliTerm) = IsApprox.isunitary(pt.coeff)

"""
    ishermitian(pt::PauliTerm)

Return `true` if `pt` is a Hermitian operator.
"""
LinearAlgebra.ishermitian(pt::PauliTerm) = isreal(pt.coeff)

####
#### Algebra
####

function mul!(target::AbstractArray{<:AbstractPauli}, ps1::PauliTerm, ps2::PauliTerm)
    s_new, phase = multiply_keeping_phase!(target, ps1.paulis, ps2.paulis)
    return PauliTerm(s_new, ps1.coeff * ps2.coeff * phase)
end

function Base.:*(ps1::PauliTerm, ps2::PauliTerm)
    return mul!(similar(ps1.paulis), ps1, ps2)
end

Base.inv(p::PauliTerm) = PauliTerm(p.paulis, inv(p.coeff))

function Base.:^(p::PauliTerm, n::Integer)
    new_coeff = p.coeff^n
    if iseven(n)
        return PauliTerm(fill(Pauli(:I), length(p)), new_coeff)
    else
        return PauliTerm(p.paulis, new_coeff)
    end
end

function Base.conj(pt::PauliTerm)
    num_ys = count(p -> op_index(p) == 2, pt.paulis)
    fac = iseven(num_ys) ? 1 : -1
    return PauliTerm(pt.paulis, conj(pt.coeff) * fac)
end

function Base.transpose(pt::PauliTerm)
    num_ys = count(p -> op_index(p) == 2, pt.paulis)
    fac = iseven(num_ys) ? 1 : -1
    return PauliTerm(pt.paulis, pt.coeff * fac)
end

function Base.adjoint(pt::PauliTerm)
    return PauliTerm(pt.paulis, conj(pt.coeff))
end

function LinearAlgebra.eigvals(pt::PauliTerm)
    vals = Vector{promote_type(Float64, typeof(pt.coeff))}(undef, 2 * length(pt))
    pos_eigval = pt.coeff
    neg_eigval = -pos_eigval
    if real(neg_eigval) > real(pos_eigval)
        (pos_eigval, neg_eigval) = (neg_eigval, pos_eigval)
    end
    @inbounds for i in eachindex(pt)
        vals[i] = neg_eigval  # follow the usual order
        vals[i+length(pt)] = pos_eigval
    end
    return vals
end

Base.kron(ps1::PauliTerm, ps2::PauliTerm) = PauliTerm(vcat(ps1.paulis, ps2.paulis), ps1.coeff * ps2.coeff)

Base.kron(paulis::AbstractPauli...) = PauliTerm([paulis...])

# TODO: @code_warntype shows red here
function Base.kron(ps::Union{PauliTerm, AbstractPauli}...)
    if ps[1] isa AbstractPauli
        v = typeof(ps[1])[]
    else
        v = eltype(ps[1])[]
    end
    coeffs = []
    for x in ps
        if x isa PauliTerm
            push!(coeffs, x.coeff)
            append!(v, x)
        else
            push!(v, x)
        end
    end
    if isempty(coeffs)
        return PauliTerm(v)
    else
        return PauliTerm(v, reduce(*, coeffs))
    end
end

####
#### Other ?
####

"""
    pauli_basis(::Type{PauliT}, n_qubits; coeff=_DEFAULT_COEFF)

Return a `Generator` over all `PauliTerm`s of `n_qubits`.
"""
function pauli_basis(::Type{PauliT}, n_qubits; coeff=_DEFAULT_COEFF) where PauliT
    return (PauliTerm(PauliT, i, n_qubits, coeff) for i in 0:(4^n_qubits - 1))
end

weight(ps::PauliTerm) = weight(ps.paulis)
