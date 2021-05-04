function isapprox_zero(x::Symbolics.Num)
    return  x === zero(x)
end
