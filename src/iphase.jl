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

IPhase(x::Union{Integer, Complex{<:Integer}}) = convert(IPhase, x)

function Base.:*(x::IPhase, y::IPhase)
    sign = ~(~(x.minus ⊻ y.minus) ⊻ (x.imag & y.imag))
    imag = x.imag ⊻ y.imag
    return IPhase(imag, sign)
end

Base.:*(x::Union{Integer, Complex{<:Integer}}, y::IPhase) = convert(IPhase, x) * y
Base.:*(y::IPhase, x::Union{Integer, Complex{<:Integer}}) = x * y

function Base.convert(::Type{IPhase}, x::Union{Integer, Complex{<:Integer}})
    if isreal(x) # || iszero(x.im)
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
