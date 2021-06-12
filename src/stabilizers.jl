module Stabilizers

import ..AbstractStabilizers: AbstractStabilizer

export Stabilizer
export pI, mI, pX, mX, pY, mY, pZ, mZ

struct Stabilizer <: AbstractStabilizer
    x::Bool
    z::Bool
    p::Bool
end

const pI = Stabilizer(0, 0, 0)
const mI = Stabilizer(0, 0, 1)
const pX = Stabilizer(1, 0, 0)
const mX = Stabilizer(1, 0, 1)
const pZ = Stabilizer(0, 1, 0)
const mZ = Stabilizer(0, 1, 1)
const pY = Stabilizer(1, 1, 0)
const mY = Stabilizer(1, 1, 1)

function Base.show(io::IO, s::Stabilizer)
    if s.p
        print(io, "-")
    else
        print(io, "+")
    end
    if s.x && s.z
        print(io, "𝒴")
    elseif s.x
        print(io, "𝒳")
    elseif s.z
        print(io, "𝒵")
    else
        print(io, "ℐ")
    end
end

function Base.:*(o1::Stabilizer, o2::Stabilizer)
    (x1, z1) = (o1.x, o1.z)
    (x2, z2) = (o2.x, o2.z)
    _sign = (x1 & z2 & (x2 | z1)) | (~x1 & x2 & z1 & ~z2)
    phase = (o1.p + o2.p + _sign) % 2
    return Stabilizer(x1 ⊻ x2, z1 ⊻ z2, phase)
end

end # module Stabilizers
