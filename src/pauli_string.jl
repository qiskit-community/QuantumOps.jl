export PauliString

####
#### PauliString
####

struct PauliString{W<:AbstractPauli, T<:AbstractVector{W}, V}
    paulis::T
    coeff::V
end

const _default_coeff = Complex(1, 0)

PauliString(s) = PauliString(s, _default_coeff)
function PauliString(::Type{T}, s::AbstractString, coeff=_default_coeff) where T <: AbstractPauli
    return PauliString(Vector{T}(s), coeff)
end

PauliString(ps::Tuple, coeff=_default_coeff) = PauliString([ps...], coeff)

# Forward some functions to the array containing the Pauli string
for func in (:length, :size, :eachindex, :insert!, :push!, :popat!, :splice!, :eltype, :getindex, :setindex!,
             :iterate)
    @eval begin
        function Base.$func(ps::PauliString, args...)
            return $func(ps.paulis, args...)
        end
    end
end

# Huh? The follow destroys proper printing of psvec
# convert psvec[i, j] to psvec[i][j]
# Base.getindex(psvec::Vector{<:PauliString}, ind1::Int, ind2::Int) = psvec[ind1][ind2]
# Base.getindex(psvec::Vector{<:PauliString}, ind1, ind2) = psvec[ind1][ind2]
# Base.setindex!(psvec::Vector{<:PauliString}, val, ind1::Int, ind2::Int) = (psvec[ind1][ind2] = val)
# Base.setindex!(psvec::Vector{<:PauliString}, val, ind1, ind2) = psvec[ind1][ind2]

Base.:(==)(ps1::PauliString, ps2::PauliString) = ps1.coeff == ps2.coeff && ps1.paulis == ps2.paulis
Base.isless(ps1::PauliString, ps2::PauliString) = isless(ps1.paulis, ps1.paulis) || isless(ps1.paulis, ps2.paulis)

function Base.show(io::IO, ps::PauliString)
    print(io, ps.coeff, " * ")
    print(io, ps.paulis)
end

function Base.:*(ps1::PauliString, ps2::PauliString)
    s_new, phase = multiply_keeping_phase(ps1.paulis, ps2.paulis)
    return PauliString(s_new, ps1.coeff * ps2.coeff * phase)
end

Base.:*(z::Number, ps::PauliString) = PauliString(ps.paulis, ps.coeff * z)
Base.:*(ps::PauliString, z::Number) = z * ps

Base.kron(ps1::PauliString, ps2::PauliString) = PauliString(vcat(ps1.paulis, ps2.paulis), ps1.coeff * ps2.coeff)

Base.one(ps::PauliString{T}) where T = PauliString(fill(one(T), length(ps)))

# TODO: Probably get rid of this
function Base.rand(::Type{<:PauliString{T}}, n::Integer) where {T <: AbstractPauli}
    return PauliString(rand(T, n))
end

function Base.kron(paulis::AbstractPauli...)
    return PauliString([paulis...])
end

function Base.kron(ps::Union{PauliString, AbstractPauli}...)
    v = eltype(ps[1])[]
    coeffs = []
    for x in ps
        if x isa PauliString
            push!(coeffs, x.coeff)
            append!(v, x)
        else
            push!(v, x)
        end
    end
    if isempty(coeffs)
        return PauliString(v)
    else
        return PauliString(v, *(coeffs...))
    end
end
