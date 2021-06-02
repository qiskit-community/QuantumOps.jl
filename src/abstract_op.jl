#import .SparseVecs

abstract type AbstractOp end

function op_index end
function op_symbols end

function _AbstractOp(::Type{T}, ind::V) where {T, V}
    syms = op_symbols(T, V)
    j = 0
    for i in 1:length(syms)
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
    weight(v::AbstractArray{<:AbstractPauli})
    weight(ps::PauliTerm)

Count the number of Paulis in the string that are not the identity.
"""
weight(v::AbstractArray{<:AbstractOp}) = count(op -> !isone(op), v)
