"""
    struct PauliTermA{W<:AbstractPauli, T<:AbstractVector{W}, V}

Represents a Pauli string (tensor product of Paulis) with a coefficient.
"""
struct PauliTermA{W<:AbstractPauli, T<:AbstractVector{W}, CoeffT} <: AbstractTerm{W, CoeffT}
    paulis::T
    coeff::CoeffT
end

const PauliTerms = Union{PauliTermA, OpTerm{<:AbstractPauli}}

op_string(t::PauliTermA) = t.paulis

# Use OpTerm{Pauli} instead
# term_type(::Type{<:AbstractPauli}) = PauliTermA

strip_typeof(::PauliTermA) = PauliTermA

####
#### Constructors
####

"""
    PauliTermA(s)

Construct a `PauliTermA` with default coefficient.
"""
PauliTermA(s) = PauliTermA(s, _DEFAULT_COEFF)

"""
    PauliTermA(::Type{T}=PauliDefault, s::AbstractString, coeff=_DEFAULT_COEFF) where T <: AbstractPauli

Construct a `PauliTermA`, where `s` is of the form "XYZ", etc. If `::Type{T}` is
ommited the default implementation of `AbstractPauli` is used.

# Examples
```jldoctest
julia> PauliTermA("XX")
(1 + 0im) * XX

julia> PauliTermA(Pauli, "IXYZ")
(1 + 0im) * IXYZ

julia> PauliTermA(Pauli, "IXYZ", 2.0)
2.0 * IXYZ
```
"""
function PauliTermA(::Type{T}, s::AbstractString, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTermA(Vector{T}(s), coeff)
end

## This is defined just for construction by copying the printed type.
function PauliTermA{W, T, V}(s::AbstractString,
                            coeff=_DEFAULT_COEFF) where {W<:AbstractPauli,T<:AbstractVector{W},V}
    return PauliTermA{W,T,V}(Vector{W}(s), coeff)
end

"""
    PauliTermA(::Type{T}=PauliDefault, s::Symbol, coeff=_DEFAULT_COEFF) where T <: AbstractPauli

Construct a `PauliTermA`, where `s` is of the form `:XYZ`, etc.

# Examples
```jldoctest
julia> PauliTermA(Pauli, :IZXY)
(1 + 0im) * IZXY
```
"""
function PauliTermA(::Type{T}, s::Symbol, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTermA(Vector{T}(String(s)), coeff)
end

function PauliTermA(::Type{T}, index::Integer, n_paulis::Integer, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTermA(pauli_vector(T, index, n_paulis), coeff)
end

function PauliTermA(paulis::AbstractPauli...; coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTermA([paulis...], coeff)
end

function PauliTermA(::Type{T}, inds::AbstractVector{<:Integer}, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTermA(T.(inds), coeff)
end

# """
#     rand_pauli_term(::Type{PauliT}=PauliDefault, n::Integer; coeff=_DEFAULT_COEFF) where {PauliT <: AbstractPauli}

# Return a random `PauliTermA` of `n` tensor factors.
# """
# function rand_pauli_term(::Type{PauliT}, n::Integer; coeff=_DEFAULT_COEFF) where {PauliT <: AbstractPauli}
#     return PauliTermA(rand(PauliT, n), coeff)
# end
# rand_pauli_term(n::Integer; coeff=_DEFAULT_COEFF) = rand_pauli_term(PauliDefault, n, coeff=coeff)

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

Base.Matrix(pt::PauliTerms) = Matrix(Float64, pt)

# FIXME: Not type stable.
function Base.Matrix(::Type{Z4Group0}, pt::PauliTerms)
    matrix = _kron((Z4Group0.(m) for m in Matrix.(op_string(pt)))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

function Base.Matrix(::Type{Float64}, pt::PauliTerms)
    matrix = _kron(Matrix.(op_string(pt))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

SparseArrays.sparse(pt::PauliTerms) = SparseArrays.sparse(Float64, pt)

function SparseArrays.sparse(::Type{Float64}, pt::PauliTerms)
    matrix = _kron(SparseArrays.sparse.(op_string(pt))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

# FIXME: broken, should throw inexact error earlier rather than return wrong type
function SparseArrays.sparse(::Type{Z4Group0}, pt::PauliTerms)
    matrix = _kron((Z4Group0.(m) for m in SparseArrays.sparse.(op_string(pt)))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

####
#### Compare / predicates
####

"""
    isunitary(pt::PauliTermA)

Return `true` if `pt` is a unitary operator.
"""
IsApprox.isunitary(pt::PauliTerms) = IsApprox.isunitary(pt.coeff)

"""
    ishermitian(pt::PauliTermA)

Return `true` if `pt` is a Hermitian operator.
"""
LinearAlgebra.ishermitian(pt::PauliTerms) = isreal(pt.coeff)

####
#### Algebra
####

function mul!(target::AbstractArray{<:AbstractPauli}, ps1::PauliTerms, ps2::PauliTerms)
    s_new, phase = multiply_keeping_phase!(target, op_string(ps1), op_string(ps2))
    return strip_typeof(ps1)(s_new, ps1.coeff * ps2.coeff * phase)
end

function Base.:*(ps1::PauliTerms, ps2::PauliTerms)
    return mul!(similar(op_string(ps1)), ps1, ps2)
end

Base.inv(p::PauliTerms) = strip_typeof(p)(op_string(p), inv(p.coeff))

function Base.:^(p::PauliTerms, n::Integer)
    new_coeff = p.coeff^n
    if iseven(n)
        return strip_typeof(p)(fill(Pauli(:I), length(p)), new_coeff)
    else
        return strip_typeof(p)(op_string(p), new_coeff)
    end
end

function Base.conj(pt::PauliTerms)
    num_ys = count(p -> op_index(p) == 2, op_string(pt))
    fac = iseven(num_ys) ? 1 : -1
    return strip_typeof(pt)(op_string(pt), conj(pt.coeff) * fac)
end

function Base.transpose(pt::PauliTerms)
    num_ys = count(p -> op_index(p) == 2, op_string(pt))
    fac = iseven(num_ys) ? 1 : -1
    return strip_typeof(pt)(op_string(pt), pt.coeff * fac)
end

function Base.adjoint(pt::PauliTerms)
    return strip_typeof(pt)(op_string(pt), conj(pt.coeff))
end

function LinearAlgebra.eigvals(pt::PauliTerms)
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

Base.kron(ps1::PauliTerms, ps2::PauliTerms) = strip_typeof(ps1)(vcat(op_string(ps1), op_string(ps2)), ps1.coeff * ps2.coeff)

Base.kron(paulis::AbstractPauli...) = PauliTerm([paulis...])

# TODO: @code_warntype shows red here
function Base.kron(ps::Union{PauliTerms, AbstractPauli}...)
    if ps[1] isa AbstractPauli
        v = typeof(ps[1])[]
    else
        v = eltype(ps[1])[]
    end
    coeffs = []
    for x in ps
        if x isa PauliTerms
            push!(coeffs, x.coeff)
            append!(v, x)
        else
            push!(v, x)
        end
    end
    if isempty(coeffs)
        return PauliTermA(v)
    else
        return PauliTermA(v, reduce(*, coeffs))
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
    return (PauliTermA(PauliT, i, n_qubits, coeff) for i in 0:(4^n_qubits - 1))
end

weight(ps::PauliTerms) = weight(op_string(ps))
