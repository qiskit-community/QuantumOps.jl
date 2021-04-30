module Stabilizers

## I am abondoning this experiment. This is not a useful way to implement stabilizers

import Random
import ..AbstractStabilizers: AbstractStabilizer
import ..AbstractOps: op_index
import ..AbstractPaulis
import ..Paulis
import ..Paulis: Pauli

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

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Stabilizer})
    return Stabilizer(rand(rng, Bool), rand(rng, Bool), rand(rng, Bool))
end

function Base.:*(o1::Stabilizer, o2::Stabilizer)
    (x1, z1) = (o1.x, o1.z)
    (x2, z2) = (o2.x, o2.z)
    _sign = (x1 & z2 & (x2 | z1)) | (~x1 & x2 & z1 & ~z2)
    phase = (o1.p + o2.p + _sign) % 2
    return Stabilizer(x1 ⊻ x2, z1 ⊻ z2, phase)
end

## The symplectic order is X, Z, Y if we interpret bits as a binary index.
## So, we have to use another method
"""
    op_index(s::Stabilizer)

Return the index in `0--7` of `s`. The order is
`pI, pX, pY, pZ, mI, mX, mY, mZ`.
"""
function op_index(s::Stabilizer)
    if s.x
        if s.z
            res = 2
        else
            res = 1
        end
    else
        if s.z
            res = 3
        else
            res = 0
        end
    end
    return res + 4 * s.p
end

"""
    Pauli(s::Stabilizer)

Return the `Pauli` corresponding to `s`, disregarding
the sign and the factor of `i` in the `Y` stabilizer.
"""
function Paulis.Pauli(s::Stabilizer)
    if s.x
        if s.z
            return(Paulis.Y)
        else
            return(Paulis.X)
        end
    else
        if s.z
            return(Paulis.Z)
        else
            return(Paulis.I)
        end
    end
end

const _pI_mat = Matrix(Paulis.I)
const _pX_mat = Matrix(Paulis.X)
const _pZ_mat = Matrix(Paulis.Z)
const _pY_mat = [0.0 1.0; -1.0 0.0]
const _mI_mat = - _pI_mat
const _mX_mat = - _pX_mat
const _mY_mat = - _pY_mat
const _mZ_mat = - _pZ_mat

# function Base.Matrix(s::Stabilizer)
# end

end # module Stabilizers
