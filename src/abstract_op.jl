module AbstractOps

import Random
export AbstractOp
export weight, op_index, op_symbols, unsafe_op

abstract type AbstractOp end

"""
    op_index(op::AbstractOp)::Int

Return the `Int` index corresponding to `op`.
"""
function op_index end

function unsafe_op(::Type{T}, ind::Integer) where T
    return T(ind)
end

"""
    op_symbols(::Type{OpT}, ::Type{SymT})

Return a list of objects of type `SymT` that represent operators of type `OpT`.
The list should be in a canonical order.
# Examples
```julia
op_symbols(::Type{<:AbstractPauli}, ::Type{Symbol}) = (:I, :X, :Y, :Z)
```
"""
function op_symbols end

"""
    _AbstractOp(::Type{T}, ind::V)

Create an operator of type `T` from value `ind`.
"""
function _AbstractOp(::Type{T}, ind::V) where {T, V}
    syms = op_symbols(T, V)
    j = 0
    @inbounds for i in eachindex(syms) # 1:length(syms)
        if ind == syms[i]
            j = i
            break
        end
    end
    if j == 0
        throw(ArgumentError("Unrecognized operator symbol"))
    end
    return T(j - 1)  # TODO: abstract this index shift
end

"""
    Vector{T}(opstring::AbstractString) where {T <: AbstractOp}

Initialize a `Vector{T}` by converting each character to type `T`.
"""
Vector{T}(opstring::AbstractString) where {T <: AbstractOp} = [T(s) for s in opstring]

## TODO: If we trust rand_ind_range, we could use unsafe_pauli (or a generic version)
"""
    docstring
"""
function rand_ind_range end

# Using unsafe_op does not seem to be faster
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{OpT}) where {OpT <: AbstractOp}
#    return unsafe_op(OpT, rand(rng, rand_ind_range(OpT)))
    return OpT(rand(rng, rand_ind_range(OpT)))
end

"""
    weight(v::AbstractArray{<:AbstractOp})
    weight(ps::OpTerm{<:AbstractOp})

Count the number of Paulis in the string that are not the identity.
"""
weight(v::AbstractArray{<:AbstractOp}) = count(!isone, v)

"""
    docstring
"""
function accumulate_phase end


"""
    docstring
"""
function accumulate_phase! end

accumulate_phase!(args...) = nothing

"""
    docstring
"""
function compute_phase end

"""
    phase_data(::Type{T}) where T

Return a phase data object for type `T`.
"""
function phase_data end

"""
    docstring
"""
struct DummyPhaseData
end

phase_data(::Type{<:AbstractOp}) = DummyPhaseData()

function _show_op_plain end

function anticommutes end

# Return x and z bits as Tuple{Bool,Bool}
function symplectic end

end # module AbstractOps
