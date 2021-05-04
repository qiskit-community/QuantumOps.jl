export IGroup

"""
    struct IGroup

The group `{i, -1, -i, 1}`
"""
struct IGroup
    imag::Bool
    minus::Bool
end

Base.one(::Type{IGroup}) = IGroup(false, false)
Base.one(n::IGroup) = one(IGroup)

function Base.:*(x::IGroup, y::IGroup)
    sign = ~(~(x.minus ⊻ y.minus) ⊻ (x.imag & y.imag))
    imag = x.imag ⊻ y.imag
    return IGroup(imag, sign)
end

const ConvertableNumber = Union{Integer, Complex{<:Integer}}

Base.:*(x::ConvertableNumber, y::IGroup) = convert(IGroup, x) * y
Base.:*(y::IGroup, x::ConvertableNumber) = x * y

# TODO: finish this
# function Base.:^(x::IGroup, n::Integer)
#     nr = mod(n, 5)
#     if nr == 0
#         ip = IGroup(1)
#     elseif nr == 1
#         ip = IGroup(im)
#     elseif nr == 2
#         ip = IGroup(-im)
# end

IGroup(x::ConvertableNumber) = convert(IGroup, x)
function Base.convert(::Type{IGroup}, x::ConvertableNumber)
    if isreal(x)
        isone(x) && return one(IGroup)
        isone(-x) && return -one(IGroup)
    else iszero(x.re)
        isone(x.im) && return IGroup(true, false)
        isone(-x.im) && return IGroup(true, true)
    end
    error("Can't convert $x to {1, -1, i, -i}.")
end

function Base.:-(x::IGroup)
    return IGroup(x.imag, ~x.minus)
end

function Base.show(io::IO, ip::IGroup)
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

function Base.show(io::IO, ps::PauliTerm{T,V,IGroup}) where {T,V}
    print(io, ps.coeff, " ")
    print(io, ps.paulis)
end
