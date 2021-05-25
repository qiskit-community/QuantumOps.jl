####
#### OpTerm{AbstractFermiOp}
####

function OpTerm{T}(inds::NTuple{N, Int}, coeff, n_modes::Integer) where {T<:AbstractFermiOp, N}
    (ops, phase) = index_to_ops_phase(inds)
    (factors, new_coeff) = _fermi_term(ops, phase, coeff, n_modes)
    return OpTerm{T}(factors, new_coeff)
end
