export IPhase

"""
    struct IPhase

The group `{i, -1, -i, 1}`
"""
struct IPhase
    imag::Bool
    minus::Bool
end

Base.one(::Type{IPhase}) = IPhase(false, false)
Base.one(n::IPhase) = one(IPhase)

function Base.:*(x::IPhase, y::IPhase)
    sign = ~(~(x.minus ⊻ y.minus) ⊻ (x.imag & y.imag))
    imag = x.imag ⊻ y.imag
    return IPhase(imag, sign)
end

const ConvertableNumber = Union{Integer, Complex{<:Integer}}

Base.:*(x::ConvertableNumber, y::IPhase) = convert(IPhase, x) * y
Base.:*(y::IPhase, x::ConvertableNumber) = x * y

# function Base.:^(x::IPhase, n::Integer)
#     nr = mod(n, 5)
#     if nr == 0
#         ip = IPhase(1)
#     elseif nr == 1
#         ip = IPhase(im)
#     elseif nr == 2
#         ip = IPhase(-im)
# end

IPhase(x::ConvertableNumber) = convert(IPhase, x)
function Base.convert(::Type{IPhase}, x::ConvertableNumber)
    if isreal(x)
        isone(x) && return one(IPhase)
        isone(-x) && return -one(IPhase)
    else iszero(x.re)
        isone(x.im) && return IPhase(true, false)
        isone(-x.im) && return IPhase(true, true)
    end
    error("Can't convert $x to {1, -1, i, -i}.")
end

function Base.:-(x::IPhase)
    return IPhase(x.imag, ~x.minus)
end

function Base.show(io::IO, ip::IPhase)
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
