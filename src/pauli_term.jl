export PauliTerm

struct PauliTerm{W<:AbstractPauli, T<:AbstractVector{W}, V}
    paulis::T
    coeff::V
end

const _default_coeff = Complex(1, 0)

PauliTerm(s) = PauliTerm(s, _default_coeff)
function PauliTerm(::Type{T}, s::AbstractString, coeff=_default_coeff) where T <: AbstractPauli
    return PauliTerm(Vector{T}(s), coeff)
end

PauliTerm(ps::Tuple, coeff=_default_coeff) = PauliTerm([ps...], coeff)

# Forward some functions to the array containing the Pauli string
for func in (:length, :size, :eachindex, :insert!, :push!, :popat!, :splice!, :eltype, :getindex, :setindex!,
             :iterate)
    @eval begin
        function Base.$func(ps::PauliTerm, args...)
            return $func(ps.paulis, args...)
        end
    end
end

Base.:(==)(ps1::PauliTerm, ps2::PauliTerm) = ps1.coeff == ps2.coeff && ps1.paulis == ps2.paulis
Base.isless(ps1::PauliTerm, ps2::PauliTerm) = isless(ps1.paulis, ps1.paulis) || isless(ps1.paulis, ps2.paulis)

weight(ps::PauliTerm) = weight(ps.paulis)

function Base.show(io::IO, ps::PauliTerm)
    print(io, ps.coeff, " * ")
    print(io, ps.paulis)
end

function Base.:*(ps1::PauliTerm, ps2::PauliTerm)
    s_new, phase = multiply_keeping_phase(ps1.paulis, ps2.paulis)
    return PauliTerm(s_new, ps1.coeff * ps2.coeff * phase)
end

Base.:*(z::Number, ps::PauliTerm) = PauliTerm(ps.paulis, ps.coeff * z)
Base.:*(ps::PauliTerm, z::Number) = z * ps
Base.:*(z::Number, p::AbstractPauli) = PauliTerm([p], z)
Base.:*(p::AbstractPauli, z::Number) = z * p

Base.:-(p::PauliTerm) = PauliTerm(p.paulis, -p.coeff)

Base.kron(ps1::PauliTerm, ps2::PauliTerm) = PauliTerm(vcat(ps1.paulis, ps2.paulis), ps1.coeff * ps2.coeff)

Base.one(ps::PauliTerm{W}) where {W} = PauliTerm(fill(one(W), length(ps)), one(ps.coeff))

# TODO: Probably get rid of this
function Base.rand(::Type{<:PauliTerm{T}}, n::Integer) where {T <: AbstractPauli}
    return PauliTerm(rand(T, n))
end

function Base.kron(paulis::AbstractPauli...)
    return PauliTerm([paulis...])
end

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
