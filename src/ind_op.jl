struct IndOp{OpT} <: AbstractOp
    op::OpT
    ind::Int
end

reg_index(op::IndOp) = op.ind

Base.copy(op::IndOp) = IndOp(copy(op.op), copy(op.ind))

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

## Don't wrap IndOp in IndOp
OpTerm{IndOp}(term::OpTerm{<:IndOp}) = term

"""
    OpSum{IndOp}(_sum::OpSum{T}) where T

Convert `_sum::OpSum` to sparser `OpSum{IndOp}`. If `_sum` is already of type `OpSum{IndOp}`,
it is returned  unchanged.
"""
function OpSum{IndOp}(_sum::OpSum{T}) where T
    isum = OpSum{IndOp{T}}()
    for term in _sum
        push!(isum, OpTerm{IndOp}(term))
    end
    return isum
end

## Don't wrap IndOp in IndOp
OpSum{IndOp}(_sum::OpSum{<:IndOp}) = _sum

####
#### Container interface
####

weight(term::OpTerm{<:IndOp}) = length(op_string(term))
Base.length(term::OpTerm{<:IndOp}) = length(op_string(term)) == 0 ? 0 : reg_index(op_string(term)[end])

####
#### Math
####

import Base: isone, iszero
import LinearAlgebra: ishermitian
for f in (:isone, :iszero, :isunitary, :ishermitian)
    @eval $f(op::IndOp) = $f(op.op)
end

import Base: one, zero, ^, adjoint

for f in (:one, :zero, :^, :adjoint, :ishermitian)
    @eval $f(op::IndOp, args...) = typeof(op)($f(op.op, args...), reg_index(op))
end

Base.isless(op1::IndOp, op2::IndOp) = Base.isless(op1.op, op2.op)

"""
    *(op1::IndOp, op2::IndOp)

Compute `op1 * op2` ignoring possible phase.

# Throws
- `ArgumentError` if `op1` and `op2` do not have the same index.
"""
function Base.:*(op1::IndOp, op2::IndOp)
    if reg_index(op1)== reg_index(op2)
        return IndOp(op1.op * op2.op, reg_index(op1))
    end
    throw(ArgumentError("Operands are not on the same index."))
end

function Base.:*(t1::OpTerm{T}, t2::OpTerm{T}) where {OpT, T <: IndOp{OpT}}
    i1 = 1  # index into terms in t1
    i2 = 1
    ops = T[]
    tcount = 0
    tmax = weight(t1) + weight(t2) + 1
    ## Include all factors up to the end of the shorter string
    phase_data = (0, 0)
    while i1 <= weight(t1) && i2 <= weight(t2)
        tcount += 1
        if tcount > tmax ## Remove this check after this routine is tested thoroughly
            println("Too many terms")
            return OpTerm(ops, t1.coeff * t2.coeff)
        end
        tt1 = t1[i1]
        tt2 = t2[i2]
        if tt1.ind < tt2.ind
            push!(ops, tt1)
            i1 += 1
        elseif tt1.ind > tt2.ind
            push!(ops, tt2)
            i2 += 1
        else  # tt1.op and tt2.op operate on the same index (DOF)
            new_op = tt1.op * tt2.op
            if iszero(new_op)  # if any factor vanishes, the term vanishes.
                return OpTerm(empty!(ops), zero(t1.coeff))
            end
            i1 += 1
            i2 += 1
            if isone(new_op)  # Identity is not stored in sparse representation
                continue
            end
            phase_data = accumulate_phase(phase_data, tt1.op, tt2.op)
            push!(ops, IndOp(new_op, tt1.ind))
        end
    end
    ## Include remaining factors from the longer string
    if i1 <= weight(t1)
        for i in i1:weight(t1)
            push!(ops, t1[i])
        end
    elseif i2 <= weight(t2)
        for i in i2:weight(t2)
            push!(ops, t2[i])
        end
    end
    return OpTerm(ops, t1.coeff * t2.coeff * compute_phase(OpT, phase_data))
end

function test_mult(n)
    t1 = rand_op_term(Pauli, n)
    t2 = rand_op_term(Pauli, n)
    ot1 = OpTerm{IndOp}(t1)
    ot2 = OpTerm{IndOp}(t2)
    return op_string(OpTerm{IndOp}(t1 * t2)) == op_string(ot1 * ot2)
#    println(op_string(OpTerm{IndOp}(t1 * t2)) == op_string(ot1 * ot2))
#    return(ot1 * ot2, t1 * t2)
end
