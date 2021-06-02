function mul(v1::SparseArraysN.SparseVector{T}, v2::SparseArraysN.SparseVector{T}) where {T <: AbstractOp}
    ops = similar(v1, 0)  # Length 0
    phase_data = (0, 0)
    i1 = 1  # index into terms in t1
    i2 = 1
    while i1 <= weight(t1) && i2 <= weight(t2)
        ind1 = v1.nzind[i1]
        ind2 = v2.nzind[i2]
        val1 = v1.nzval[i1]
        val2 = v2.nzval[i2]
        if ind1 < ind2
            push!(ops.nzind, ind1)
            push!(ops.nzval, val1)
            i1 += 1
        elseif tt1.ind > tt2.ind
            push!(ops.nzind, ind2)
            push!(ops.nzval, val2)
            i2 += 1
        else  # tt1.op and tt2.op operate on the same index (DOF)
            new_op = val1 * val2
            if iszero(new_op)  # if any factor vanishes, the term vanishes.
                return (empty!(ops), zero(t1.coeff))
            end
            i1 += 1
            i2 += 1
            if isone(new_op)  # Identity is not stored in sparse representation
                continue
            end
            phase_data = accumulate_phase(phase_data, va1, val2)
            push!(ops.nzind, ind1)
            push!(ops.nzval, new_op)
        end
    end
end

# function Base.:*(t1::OpTerm{T}, t2::OpTerm{T}) where {OpT, T <: IndOp{OpT}}
#     i1 = 1  # index into terms in t1
#     i2 = 1
#     ops = T[]
#     tcount = 0
#     tmax = weight(t1) + weight(t2) + 1
#     ## Include all factors up to the end of the shorter string
#     phase_data = (0, 0)
#     while i1 <= weight(t1) && i2 <= weight(t2)
#         tcount += 1
#         if tcount > tmax ## Remove this check after this routine is tested thoroughly
#             println("Too many terms")
#             return OpTerm(ops, t1.coeff * t2.coeff)
#         end
#         tt1 = t1[i1]
#         tt2 = t2[i2]
#         if tt1.ind < tt2.ind
#             push!(ops, tt1)
#             i1 += 1
#         elseif tt1.ind > tt2.ind
#             push!(ops, tt2)
#             i2 += 1
#         else  # tt1.op and tt2.op operate on the same index (DOF)
#             new_op = tt1.op * tt2.op
#             if iszero(new_op)  # if any factor vanishes, the term vanishes.
#                 return OpTerm(empty!(ops), zero(t1.coeff))
#             end
#             i1 += 1
#             i2 += 1
#             if isone(new_op)  # Identity is not stored in sparse representation
#                 continue
#             end
#             phase_data = accumulate_phase(phase_data, tt1.op, tt2.op)
#             push!(ops, IndOp(new_op, tt1.ind))
#         end
#     end
#     ## Include remaining factors from the longer string
#     if i1 <= weight(t1)
#         for i in i1:weight(t1)
#             push!(ops, t1[i])
#         end
#     elseif i2 <= weight(t2)
#         for i in i2:weight(t2)
#             push!(ops, t2[i])
#         end
#     end
#     return OpTerm(ops, t1.coeff * t2.coeff * compute_phase(OpT, phase_data))
# end
