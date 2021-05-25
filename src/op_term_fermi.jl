####
#### OpTerm{AbstractFermiOp}
####

const FermiTerm = OpTerm{FermiOp}
const FermiSum = OpSum{FermiOp}

function OpTerm{T}(inds::NTuple{N, Int}, coeff, n_modes::Integer) where {T<:AbstractFermiOp, N}
    (ops, phase) = index_to_ops_phase(inds)
    (factors, new_coeff) = _fermi_term(ops, phase, coeff, n_modes)
    return OpTerm{T}(factors, new_coeff)
end

function OpSum{OpT}(tensor::AbstractArray{<:Number}, tensors::AbstractArray{<:Number}...) where {OpT<:AbstractFermiOp}
    fsum = OpSum(Vector{OpT}[], eltype(tensor)[])
    for tensor in (tensor, tensors...)
        tensor_to_fermi_sum!(fsum, tensor)
    end
    return fsum
end

function OpSum{OpT}(iop::ElectronicStructure.InteractionOperator) where {OpT<:AbstractFermiOp}
    fsum = OpSum{OpT}(iop.one_body_tensor, iop.two_body_tensor)
    add!(fsum, iop.nuclear_repulsion * one(fsum[1]))
    return fsum
end
