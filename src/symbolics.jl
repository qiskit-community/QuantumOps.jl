# A hack because there is no method in Symbolics.jl
function isapprox_zero(x::Symbolics.Num)
    return  x === zero(x)
end
