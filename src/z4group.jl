# There is very litte written against the abstract type.
# So this is of limited use.

module Z4Groups

import Random

export Z4Group, z4group, ConvertableNumber

abstract type AbstractZ4Group end

Base.one(n::AbstractZ4Group) = one(typeof(n))

const ConvertableNumber = Union{Integer, Complex{<:Integer}, Float64, Complex{<:Float64}}

"""
    struct Z4Group

Type representing the group `{i, -1, -i, 1}`
"""
struct Z4Group <: AbstractZ4Group
    imag::Bool
    minus::Bool
end

const Z4GROUP_i = Z4Group(1, 0)
const Z4GROUP_minus_one = Z4Group(0, 1)
const Z4GROUP_minus_i = Z4Group(1, 1)
const Z4GROUP_one = Z4Group(0, 0)

const Z4GROUP_ELEMENTS = (Z4GROUP_i, Z4GROUP_minus_one, Z4GROUP_minus_i, Z4GROUP_one)

#####
##### Constructors
#####

"""
    z4group(n::Integer)::Z4Group

Perform zero-based index `n` into `(i, -1, -i, 1)`.
"""
z4group(n::Integer) = Z4GROUP_ELEMENTS[n]

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Z4Group})
    return z4group(rand(rng, 1:4))
end

#####
##### Conversion
#####

Base.promote_type(::Type{Z4Group}, ::Type{T}) where {T <: Real} = Complex{T}
Base.promote_type(::Type{T}, ::Type{Z4Group}) where {T <: Real} = Complex{T}
Base.promote_type(::Type{Complex{T}}, ::Type{Z4Group}) where {T <: Real} = Complex{T}

Z4Group(x::ConvertableNumber) = convert(Z4Group, x)
function Base.convert(::Type{Z4Group}, x::ConvertableNumber)
    if isreal(x)
        isone(x) && return one(Z4Group)
        isone(-x) && return -one(Z4Group)
    else iszero(x.re)
        isone(x.im) && return Z4Group(true, false)
        isone(-x.im) && return Z4Group(true, true)
    end
    throw(InexactError(:convert, Z4Group, x))
#    error("Can't convert $x to {1, -1, i, -i}.")
end

function Base.convert(::Type{T}, x::Z4Group) where {T<:Real}
    if ~x.imag
        if x.minus
            return -one(T)
        else
            return one(T)
        end
    else
        throw(InexactError(:convert, T, x))
    end
end

function Base.convert(::Type{Complex{T}}, x::Z4Group) where {T<:Real}
    if isreal(x)
        return Complex(convert(T, x), zero(T))
    end
    if x.minus
        return Complex(zero(T), -one(T))
    else
        return Complex(zero(T), one(T))
    end
end

####
#### IO
####

function Base.show(io::IO, ip::Z4Group)
    if ip.minus
        print(io, "-")
    else
        print(io, "+")
    end
    if ip.imag
        print(io, "i")
    else
        print(io, "1")
    end
end

####
#### Container interface
####

Base.length(::Z4Group) = 1
Base.iterate(n::Z4Group) = (n, nothing)
Base.iterate(n::Z4Group, ::Nothing) = nothing

####
#### Compare / predicates
####

Base.one(::Type{Z4Group}) = Z4Group(false, false)
Base.isreal(x::Z4Group) = ~x.imag

function Base.real(x::Z4Group)
    if x.imag
        return 0
    else
        return x.minus ? -1 : 1
    end
end

function Base.imag(x::Z4Group)
    if ~x.imag
        return 0
    else
        return x.minus ? -1 : 1
    end
end

## Convert n in (0, 1, 2, 3) to binary as a Tuple of Bools
_int_to_bin(n::Integer) = (Bool((n & 2) >> 1), Bool(n & 1))

####
#### Algebra / mathematical operations
####

Base.:*(x::ConvertableNumber, y::Z4Group) = convert(Z4Group, x) * y
Base.:*(y::Z4Group, x::ConvertableNumber) = x * y

function Base.:*(x::Z4Group, y::Z4Group)
    minus = ~(~(x.minus ⊻ y.minus) ⊻ (x.imag & y.imag))
    imag = x.imag ⊻ y.imag
    return Z4Group(imag, minus)
end

Base.:-(x::Z4Group) = Z4Group(x.imag, ~x.minus)

end # module Z4Groups

# TODO: finish this
# function Base.:^(x::Z4Group, n::Integer)
#     nr = mod(n, 5)
#     if nr == 0
#         ip = Z4Group(1)
#     elseif nr == 1
#         ip = Z4Group(im)
#     elseif nr == 2
#         ip = Z4Group(-im)
# end
