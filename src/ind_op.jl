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

####
#### Container interface
####

weight(term::OpTerm{<:IndOp}) = length(op_string(term))
Base.length(term::OpTerm{<:IndOp}) = length(op_string(term)) == 0 ? 0 : op_string(term)[end].ind

####
#### Math
####

import Base: isone, iszero
import LinearAlgebra: ishermitian
for f in (:isone, :iszero, :isunitary, :ishermitian)
    @eval $f(op::IndOp) = $f(op.op)
end

Base.isless(op1::IndOp, op2::IndOp) = Base.isless(op1.op, op2.op)

function Base.:*(op1::IndOp, op2::IndOp)
    if op1.ind == op2.ind
        return IndOp(op1.op * op2.op, op1.ind)
    end
    throw(ArgumentError("Operands are not on the same index."))
end

Base.:^(op::IndOp, n::Integer) = IndOp(op.op^n, op.ind)

## TODO: Use the following two functions to compute phase elsewhere as well.

function accumulate_phase(old_phase_data, op1::T, op2::T) where {T <: AbstractPauli}
    phase_info = phase(op1, op2)
    new_phase_data = (old_phase_data[1] + phase_info[1], old_phase_data[2] + phase_info[2])
    return new_phase_data
end

function compute_phase(::Type{<:AbstractPauli}, phase_data)
    (n_sign_flips, n_imag_units) = phase_data
    _sign = iseven(n_sign_flips) ? 1 : -1
    n = n_imag_units % 4
    if n == 0
        cim = complex(1)
    elseif n == 1
        cim = complex(im)
    elseif n == 2
        cim = complex(-im)
    else
        cim = complex(-1)
    end
    return _sign * cim
end

## We do not need to track phase for FermiOps (for the set of Fermi ops that we support)
accumulate_phase(old_phase_data, op1::T, op2::T) where {T <: AbstractFermiOp} = old_phase_data
compute_phase(::Type{<:AbstractFermiOp}, _) = 1

function Base.:*(t1::OpTerm{T}, t2::OpTerm{T}) where {V, T <: IndOp{V}}
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
        else
            new_op = tt1.op * tt2.op
            i1 += 1
            i2 += 1
            if isone(new_op)
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
    return OpTerm(ops, t1.coeff * t2.coeff * compute_phase(V, phase_data))
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
