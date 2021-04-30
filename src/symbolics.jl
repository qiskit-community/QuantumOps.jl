function isapprox_zero(x::Symbolics.Num)
    return  x === zero(x)
end

# Cheat and don't convert to a Num
function Base.convert(::Type{Symbolics.Num}, z::Z4Group0)
    if isreal(z)
        return convert(Int, z)
    else
        return convert(Complex{Int}, z)
    end
end
