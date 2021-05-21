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

## Factor this out of here and PauliSum
## Factoring this out may be over-abstraction. Or maybe there is a good way to do it.
function FermiSum(v::AbstractVector{<:FermiTerm}, already_sorted=false)
    strings = [x.ops for x in v]
    coeffs = [x.coeff for x in v]
    return FermiSum(strings, coeffs, already_sorted)
end

# function FermiSum(two_body::AbstractArray{T,4}) where T
#     n_modes = first(size(two_body))
#     fsum = FermiSum{FermiOp, Vector{Vector{FermiOp}}, eltype(two_body)}()
#     for ind in CartesianIndices(two_body)
#         if (ind[1] == ind[2] || ind[3] == ind[4]) || two_body[ind] == 0
#             continue
#         end
#         add!(fsum, ferm_term(ind.I, two_body[ind], n_modes)
#     end
#     return terms
# end
