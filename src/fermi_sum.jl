import ElectronicStructure

struct FermiSumA{OpT, StringT, CoeffT} <: AbstractSum{OpT, StringT, CoeffT}
    strings::StringT
    coeffs::CoeffT

    function FermiSumA(strings, coeffs; already_sorted=false)
        _abstract_sum_inner_constructor_helper!(strings, coeffs; already_sorted=already_sorted)
        return new{eltype(eltype(strings)), typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

function term_type(::Type{T}) where T <: FermiSumA
    return FermiTermA
end

function sum_type(::Type{T}) where T <: FermiTermA
    return FermiSumA
end

strip_typeof(::FermiSumA) = FermiSumA

const FermiSums = Union{FermiSumA, OpSum{<:AbstractFermiOp}}

## Factor this out of here and PauliSum
## Factoring this out may be over-abstraction. Or maybe there is a good way to do it.
function FermiSumA(v::AbstractVector{<:FermiTermA}; already_sorted=false)
    strings = [x.ops for x in v]
    coeffs = [x.coeff for x in v]
    return FermiSumA(strings, coeffs; already_sorted=already_sorted)
end

#FermiSum{T,V}(v, c; already_sorted=false) where {T, V} = FermiSum(v, c; already_sorted=already_sorted)
#FermiSum{T,V}(v, c) where {T, V} = FermiSum(v, c, false)

_skip_inds(i, j, k, l) = (i == j || k == l) # skip if two raising or two lowering ops on one index
_skip_inds(i, j) = false # always one raising and one lowering, so never skip

function tensor_to_fermi_sum!(fsum, tensor)
    n_modes = first(size(tensor))
    for ind in CartesianIndices(tensor)
        if _skip_inds((ind.I)...) || iszero(tensor[ind])
            continue
        end
        add!(fsum, term_type(typeof(fsum))(ind.I, tensor[ind], n_modes))
    end
    return fsum
end

function FermiSumA(tensors::AbstractArray{<:Number}...)
    # fsum = FermiSumA{Vector{Vector{FermiOp}},
    #                 Vector{eltype(tensors[1])}}(Vector{FermiOp}[], eltype(tensors[1])[])
    fsum = FermiSumA(Vector{FermiOp}[], eltype(tensors[1])[])
    for tensor in tensors
        tensor_to_fermi_sum!(fsum, tensor)
    end
    return fsum
end

"""
    FermiSumA(iop::ElectronicStructure.InteractionOperator)

Convert `iop` to a `FermiSumA`.
"""
function FermiSumA(iop::ElectronicStructure.InteractionOperator)
    fsum = FermiSumA(iop.one_body_tensor, iop.two_body_tensor)
    add!(fsum, iop.nuclear_repulsion * one(fsum[1]))
    return fsum
end

function Base.adjoint(ft::FermiSums)
    op_strings = [adjoint.(x) for x in ft.strings]
    coeffs = adjoint.(ft.coeffs)
    return strip_typeof(ft)(op_strings, coeffs)
end
