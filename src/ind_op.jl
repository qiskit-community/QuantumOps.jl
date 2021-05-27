struct IndOp{OpT} <: AbstractOp
    op::OpT
    ind::Int
end

function Base.show(io::IO, op::IndOp)
    print(io, "(", op.op, op.ind, ")")
#    print(io, "(", op.op, ", ", op.ind, ")")
end

"""
    OpTerm{IndOp}(term::OpTerm{T}) where T

Convert the dense representation `term` to a sparser `OpTerm{<:IndOp}` representation.
"""
function OpTerm{IndOp}(term::OpTerm{T}) where T
    v = IndOp{T}[]
    for (i, op) in enumerate(term)
        if ! isone(op)
            push!(v, IndOp(op, i))
        end
    end
    return OpTerm(v, term.coeff)
end

weight(term::OpTerm{<:IndOp}) = length(op_string(term))
Base.length(term::OpTerm{<:IndOp}) = length(op_string(term)) == 0 ? 0 : op_string(term)[end].ind

import Base: isone, iszero
import LinearAlgebra: ishermitian
for f in (:isone, :iszero, :isunitary, :ishermitian)
    @eval $f(op::IndOp) = $f(op.op)
end
