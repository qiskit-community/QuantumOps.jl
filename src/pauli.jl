export Pauli, PauliString, pauli_string, pauli_index, phase

struct Pauli
    x::Bool
    z::Bool
end

pauli_index(p::Pauli) = 2 * p.x + p.z

const _I = Pauli(0, 0)
const _X = Pauli(0, 1)
const _Y = Pauli(1, 0)
const _Z = Pauli(1, 1)

const _pauli_chars = ('I', 'X', 'Y', 'Z')

function Pauli(s::Symbol)
    if s == :I
        return _I
    elseif s == :X
        return _X
    elseif s == :Y
        return _Y
    elseif s == :Z
        return _Z
    else
        error("invalid Pauli symbol")
    end
end

function Pauli(s::AbstractString)
    if s == "I"
        return _I
    elseif s == "X"
        return _X
    elseif s == "Y"
        return _Y
    elseif s == "Z"
        return _Z
    else
        error("invalid Pauli symbol")
    end
end

function Pauli(s::AbstractChar)
    if s == 'I'
        return _I
    elseif s == 'X'
        return _X
    elseif s == 'Y'
        return _Y
    elseif s == 'Z'
        return _Z
    else
        error("invalid Pauli character")
    end
end

function Pauli(s::Integer)
    if s == 0
        return _I
    elseif s == 1
        return _X
    elseif s == 2
        return _Y
    elseif s == 3
        return _Z
    else
        error("invalid Pauli index")
    end
end

Base.:*(p1::Pauli, p2::Pauli) = Pauli(xor(p1.x, p2.x), xor(p1.z, p2.z))
Base.:^(p::Pauli, n::Integer) = iseven(n) ? Pauli(:I) : p

# 0 for phase 1
# 1 for phase -1, ie 1 for each sign flip
function phase(p1::Pauli, p2::Pauli)
    if p1 == _I || p2 == _I
        return (0, 0)
    end
    d = pauli_index(p1) - pauli_index(p2)
    if d == 0
        return (0, 0)
    end
    if d == 1 || d == -2
        return (1, 1)
    else
        return (0, 1)
    end
end

Base.show(io::IO, p::Pauli) = print(io, _pauli_chars[pauli_index(p) + 1])

# May want to define `display` as well
function Base.show(io::IO, ps::AbstractArray{<:Pauli})
    for p in ps
        show(io, p)
    end
end

pauli_string(ps::AbstractString) = [Pauli(s) for s in ps]

function Base.:*(s1::AbstractArray{<:Pauli}, s2::AbstractArray{<:Pauli})
    length(s1) != length(s2) && error("bad lengths")
    s_out = similar(s1)
    cum_flips = 0
    cum_i = 0
    for i in eachindex(s1)
        p1 = s1[i]
        p2 = s2[i]
        sign, inum = phase(p1, p2)
        cum_flips += sign
        cum_i += inum
        s_out[i] = p1 * p2
    end
    return(s_out, (-1)^cum_flips * im^cum_i)
end

struct PauliString{T<:AbstractVector{Pauli}, V}
    s::T
    coeff::V
end

PauliString(s) = PauliString(s, 1.0)
PauliString(s::AbstractString, coeff=1.0) = PauliString(pauli_string(s), coeff)

function Base.show(io::IO, ps::PauliString)
    print(io, ps.coeff, " * ")
    print(io, ps.s)
end

function Base.:*(ps1::PauliString, ps2::PauliString)
    s_new, phase = ps1.s * ps2.s
    return PauliString(s_new, ps1.coeff * ps2.coeff * phase)
end

Base.:*(z::Number, ps::PauliString) = PauliString(ps.s, ps.coeff * z)

Base.rand(::Type{PauliString}, n::Integer) = PauliString([Pauli(i) for i in rand(0:3, n)])
