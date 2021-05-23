import ElectronicStructure

struct FermiSum{StringT, CoeffT} <: AbstractSum
    strings::StringT
    coeffs::CoeffT

    function FermiSum(strings, coeffs, already_sorted=false)
        _abstract_sum_inner_constructor_helper!(strings, coeffs, already_sorted)
        return new{typeof(strings), typeof(coeffs)}(strings, coeffs)
    end
end

function term_type(::Type{T}) where T <: FermiSum
    return FermiTerm
end

function sum_type(::Type{T}) where T <: FermiTerm
    return FermiSum
end

## Factor this out of here and PauliSum
## Factoring this out may be over-abstraction. Or maybe there is a good way to do it.
function FermiSum(v::AbstractVector{<:FermiTerm}, already_sorted=false)
    strings = [x.ops for x in v]
    coeffs = [x.coeff for x in v]
    return FermiSum(strings, coeffs, already_sorted)
end

FermiSum{T,V}(v, c, already_sorted) where {T, V} = FermiSum(v, c, already_sorted)
FermiSum{T,V}(v, c) where {T, V} = FermiSum(v, c, false)

_skip_inds(i, j, k, l) = (i == j || k == l) # skip if two raising or two lowering ops on one index
_skip_inds(i, j) = false # always one raising and one lowering, so never skip

function tensor_to_fermi_sum!(fsum, tensor)
    n_modes = first(size(tensor))
    for ind in CartesianIndices(tensor)
        if _skip_inds((ind.I)...) || iszero(tensor[ind])
            continue
        end
        add!(fsum, FermiTerm(ind.I, tensor[ind], n_modes))
    end
    return fsum
end

function FermiSum(tensors::AbstractArray{<:Number}...)
    fsum = FermiSum{Vector{Vector{FermiOp}},
                    Vector{eltype(tensors[1])}}(Vector{FermiOp}[], eltype(tensors[1])[])
    for tensor in tensors
        tensor_to_fermi_sum!(fsum, tensor)
    end
    return fsum
end

"""
    FermiSum(iop::ElectronicStructure.InteractionOperator)

Convert `iop` to a `FermiSum`.
"""
function FermiSum(iop::ElectronicStructure.InteractionOperator)
    fsum = FermiSum(iop.one_body_tensor, iop.two_body_tensor)
    add!(fsum, iop.nuclear_repulsion * one(fsum[1]))
    return fsum
end

function Base.adjoint(ft::FermiSum)
    op_strings = [adjoint.(x) for x in ft.strings]
    coeffs = adjoint.(ft.coeffs)
    return FermiSum(op_strings, coeffs)
end
