# TODO: Allow Tuple as well as Vector
# struct PauliTerm{W<:AbstractPauli, T, V}
struct PauliTerm{W<:AbstractPauli, T<:AbstractVector{W}, V}
    paulis::T
    coeff::V
end

####
#### Constructors
####

"""
    PauliTerm(s)

Construct a `PauliTerm` with default coefficient.
"""
PauliTerm(s) = PauliTerm(s, _DEFAULT_COEFF)


"""
    PauliTerm(::Type{T}, s::AbstractString, coeff=_DEFAULT_COEFF)

Construct a `PauliTerm`, where `s` is of the form "XYZ", etc.
"""
function PauliTerm(::Type{T}, s::AbstractString, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTerm(Vector{T}(s), coeff)
end

function PauliTerm(::Type{T}, s::Symbol, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTerm(Vector{T}(String(s)), coeff)
end

function PauliTerm(::Type{T}, index::Integer, n_paulis::Integer, coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTerm(pauli_vector(T, index, n_paulis), coeff)
end

function PauliTerm(paulis::AbstractPauli...; coeff=_DEFAULT_COEFF) where T <: AbstractPauli
    return PauliTerm([paulis...], coeff)
end

Base.copy(pt::PauliTerm) = PauliTerm(copy(pt.paulis), copy(pt.coeff))

### This is simple, but not flexible
function Base.rand(::Type{<:PauliTerm{PauliT}}, n::Integer; coeff=_DEFAULT_COEFF) where {PauliT <: AbstractPauli}
    return PauliTerm(rand(PauliT, n), coeff)
end

###
### The following is partially broken. It works, except the a vector of random PauliTerms will
### have return type Any. This is explained the docs for samplers. But, it is a PITA in the case.
### The example in the Julia docs is `Die` which returns an `Int`, so this is easy to fix.
### We return a rather complicated parametric type. I think the doc is asking us to compute this
### type.
### I don't see any easy way around it.
###

# struct PauliTermSampler{PauliT, V}
#     npaulis::Int # number of sides
#     coeff_func::V
# end

# function Base.eltype(::Type{PauliTermSampler{PauliT,V}})
#     return
# end

# function Random.rand(rng::Random.AbstractRNG, d::Random.SamplerTrivial{PauliTermSampler{PauliT,V}}) where {PauliT, V}
#     v = rand(rng, PauliT, d[].npaulis)
#     return PauliTerm(v, d[].coeff_func())
# end

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
        matrix = _kron(SparseArrays.sparse.(Matrix.(pt.paulis))...)  # TODO: precompute sparse arrays.
    return _multiply_coefficient(pt.coeff, matrix)
end

# FIXME: broken, should throw inexact error earlier rather than return wrong type
function SparseArrays.sparse(::Type{Z4Group0}, pt::PauliTerm)
    matrix = _kron((Z4Group0.(m) for m in SparseArrays.sparse.(Matrix.(pt.paulis)))...)
    return _multiply_coefficient(pt.coeff, matrix)
end

####
#### IO
####

## type params only to get correct dispatch. There must be a better way
function Base.show(io::IO, ps::PauliTerm{T,V,CoeffT}) where {T,V,CoeffT}
    if ps.coeff isa Real  # could use CoeffT here.
        print(io, ps.coeff)
    else
        print(io, "(", ps.coeff, ")")
    end
    print(io, " * ")
    print(io, ps.paulis)
end

function Base.show(io::IO, ps::PauliTerm{T,V,Z4Group}) where {T,V}
    print(io, ps.coeff, " ")
    print(io, ps.paulis)
end

####
#### Container interface
####

# :popat!
for func in (:length, :size, :eltype, :eachindex, :axes, :splice!, :getindex,
             :setindex!, :iterate, :pop!, :popfirst!)
    @eval begin
        Base.$func(ps::PauliTerm, args...) = $func(ps.paulis, args...)
    end
end

for func in (:push!, :pushfirst!, :insert!)
    @eval begin
        Base.$func(ps::PauliTerm, args...) = ($func(ps.paulis, args...); ps)
    end
end

"""
    reverse!(pt::PauliTerm)

Reverse the order of the factors in `pt` in place.
See `reverse`.
"""
Base.reverse!(pt::PauliTerm) = (reverse!(pt.paulis); pt)

"""
    reverse(pt::PauliTerm)

Reverse the order of the factors in `pt`.
See `reverse!`.
"""
Base.reverse(pt::PauliTerm) = reverse!(copy(pt))

####
#### Compare / predicates
####

Base.one(ps::PauliTerm{W}) where {W} = PauliTerm(fill(one(W), length(ps)), one(ps.coeff))

Base.:(==)(ps1::PauliTerm, ps2::PauliTerm) = ps1.coeff == ps2.coeff && ps1.paulis == ps2.paulis

function Base.isless(ps1::PauliTerm, ps2::PauliTerm)
    return isless(ps1.paulis, ps1.paulis) || isless(ps1.paulis, ps2.paulis)
end

# Fails until new IsApprox is registered
IsApprox.isunitary(pt::PauliTerm) = IsApprox.isunitary(pt.coeff)

LinearAlgebra.ishermitian(pt::PauliTerm) = isreal(pt.coeff)

## TODO: Symptectic has a faster way, I think
## This is generic. We can use IsApprox.commutes instead.
## function commute(ps1::PauliTerm, ps2::PauliTerm)
#     p1 = ps1 * ps2
#     p2 = ps2 * ps1
#     return p1 == p2
# end

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

Base.:*(z::Number, ps::PauliTerm) = PauliTerm(ps.paulis, ps.coeff * z)
Base.:*(ps::PauliTerm, z::Number) = z * ps

Base.:*(z::Number, p::AbstractPauli) = PauliTerm([p], z)
Base.:*(p::AbstractPauli, z::Number) = z * p

Base.:/(ps::PauliTerm, z::Number) = PauliTerm(ps.paulis, ps.coeff / z)
Base.:/(p::AbstractPauli, z::Number) = PauliTerm([p], inv(z))

Base.:-(p::PauliTerm) = PauliTerm(p.paulis, -p.coeff)
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
    num_ys = count(p -> pauli_index(p) == 2, pt.paulis)
    fac = iseven(num_ys) ? 1 : -1
    return PauliTerm(pt.paulis, conj(pt.coeff) * fac)
end

function Base.transpose(pt::PauliTerm)
    num_ys = count(p -> pauli_index(p) == 2, pt.paulis)
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
