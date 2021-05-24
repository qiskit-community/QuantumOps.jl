module AbstractOps

import Random
export AbstractOp
export weight, op_index, op_symbols

abstract type AbstractOp end

"""
    op_index(op::AbstractOp)::Int

Return the `Int` index corresponding to `op`.
"""
function op_index end

"""
    docstring
"""
function op_symbols end

function _AbstractOp(::Type{T}, ind::V) where {T, V}
    syms = op_symbols(T, V)
    j = 0
    for i in eachindex(syms) # 1:length(syms)
        if ind == syms[i]
            j = i
            break
        end
    end
    if j == 0
        throw(ArgumentError("Unrecognized operator symbol"))
    end
    return T(j - 1)  # TODO: abstract this
end

"""
    Vector{T}(opstring::AbstractString) where {T <: AbstractOp}

Initialize a `Vector{T}` by converting each character to type `T`.
"""
Vector{T}(opstring::AbstractString) where {T <: AbstractOp} = [T(s) for s in opstring]

"""
    docstring
"""
function rand_ind_range end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{OpT}) where {OpT <: AbstractOp}
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

end # module AbstractOps
