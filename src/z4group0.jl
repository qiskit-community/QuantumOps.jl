module Z4Group0s

export Z4Group0, z4group0
import Random

using ..Z4Groups

"""
    struct Z4Group0

Type that can represent 0, +1, -1, +im, -im
and supports sufficient algebra to compute and
represent tensor products of Pauli matrices.
In particular, addition is only supported if
at least one operand is zero. This type is closed
under multiplication and this restricted addition.
"""
struct Z4Group0
    z4::Z4Group
    zero::Bool
end

####
#### Constructors
####

"""
    z4group0(n::Integer)::Z4Group0

Perform one-based index `n` into `(1, -1, i, -i, 0)`.
"""
function z4group0(n::Integer)
    if n == 5
        return zero(Z4Group0)
    end
    return Z4Group0(z4group(n), true)
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Z4Group0})
    return z4group0(rand(rng, 1:5))
end

####
#### Conversion
####

Base.promote_type(::Type{T}, ::Type{Z4Group0}) where {T <: Real} = Complex{T}
Base.promote_type(::Type{Complex{T}}, ::Type{Z4Group0}) where {T <: Real} = Complex{T}

Z4Group0(x::ConvertableNumber) = convert(Z4Group0, x)

function Base.convert(::Type{Z4Group0}, x::ConvertableNumber)
    if iszero(x)
        return zero(Z4Group0)
    end
    if isreal(x)
        if isone(x)
            z = one(Z4Group)
        elseif isone(-x)
            z = -one(Z4Group)
        else
            error("Can't convert $x to {0, 1, -1, i, -i}.")
        end
    elseif iszero(x.re)
        if  isone(x.im)
            z = Z4Group(true, false) # can probably use constant prop here.
        elseif  isone(-x.im)
            z = Z4Group(true, true)
        else
            error("Can't convert $x to {0, 1, -1, i, -i}.")
        end
    else
        error("Can't convert $x to {0, 1, -1, i, -i}.")
    end
    return Z4Group0(z, true)
end

function Base.convert(::Type{T}, x::Z4Group0) where {T<:Number}
    if iszero(x)
        return zero(T)
    end
    return convert(T, x.z4)
end

####
#### IO
####

function Base.show(io::IO, ip::Z4Group0)
    if ! ip.zero
        print(io, "0")
        return nothing
    end
    show(io, ip.z4)
    return nothing
end

####
#### Container interface
####

Base.length(::Z4Group0) = 1
Base.iterate(n::Z4Group0) = (n, nothing)
Base.iterate(n::Z4Group0, ::Nothing) = nothing

####
#### Compare / predicates
####

Base.iszero(x::Z4Group0) = ~(x.zero)

Base.zero(::Type{Z4Group0}) = Z4Group0(one(Z4Group), false)

Base.one(::Type{Z4Group0}) = Z4Group0(one(Z4Group), true)
# Consider making Z4Group0 a Number so that the following need not be defined

Base.:(==)(x::Z4Group0, y::Number) = convert(Complex{Int}, x) == y

Base.isreal(x::Z4Group0) = x.zero && isreal(x.z4)

Base.real(x::Z4Group0) = ~x.zero ? 0 : real(x.z4)

Base.imag(x::Z4Group0) = ~x.zero ? 0 : imag(x.z4)

####
#### Algebra / mathematical operations
####

function Base.:*(xz::Z4Group0, yz::Z4Group0)
    zerox = xz.zero
    zeroy = yz.zero
    zero0 = zerox & zeroy
    if ! zero0
        return zero(Z4Group0)
    end
    x = xz.z4
    y = yz.z4
    sign = zero0 & (~(~(x.minus ⊻ y.minus) ⊻ (x.imag & y.imag)))
    imag = zero0 & (x.imag ⊻ y.imag)
    return Z4Group0(Z4Group(imag, sign), zero0)
end

Base.:*(xz::Z4Group, yz::Z4Group0) = yz * xz

function Base.:*(x::Z4Group0, y::Z4Group)
    if ! x.zero
        return x
    end
    return Z4Group0(x.z4 * y, true)
end

Base.:*(x::Z4Group0, y::Number) = *(promote(y, x)...)

Base.:*(y::Number, x::Z4Group0) = x * y

function Base.:+(x::Z4Group0, y::Z4Group0)
    if iszero(x)
        return y
    elseif iszero(y)
        return x
    else
        throw(DomainError("cant do this."))
    end
end

function Base.:-(x::Z4Group0)
    z4 = Z4Group(x.z4.imag, ~x.z4.minus)
    return Z4Group0(z4, x.zero)
end

end # module Z4Group0s
