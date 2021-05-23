## Jordan-Wigner helper
function fill_pauli(pad, op_ind, fill_op, end_op)
    str = Vector{Paulis.Pauli}(undef, pad)
    if pad <  op_ind
        throw(DimensionMismatch("pad is less than op_ind."))
    end
    @inbounds for i in 1:op_ind-1
        str[i] = fill_op
    end
    str[op_ind] = end_op
    @inbounds for i in op_ind+1:pad
        str[i] = Paulis.I
    end
    return str
end

## Jordan-Wigner
function jordan_wigner(op::FermiOp, op_ind::Integer, pad::Integer)
    already_sorted = true
    if op === lower_op
        strx = fill_pauli(pad, op_ind, Paulis.Z, Paulis.X)
        stry = fill_pauli(pad, op_ind, Paulis.Z, Paulis.Y)
        return PauliSum([strx, stry], [-1/2, -im/2], already_sorted)
    elseif op === raise_op
        strx = fill_pauli(pad, op_ind, Paulis.Z, Paulis.X)
        stry = fill_pauli(pad, op_ind, Paulis.Z, Paulis.Y)
        return PauliSum([strx, stry], [-1/2, im/2], already_sorted)
    elseif op === number_op
        strx = fill(Paulis.I, pad)
        stry = fill_pauli(pad, op_ind, Paulis.I, Paulis.Z)
        return PauliSum([strx, stry], complex.([1/2, -1/2]), already_sorted)
    elseif op === empty_op
        strx = fill(Paulis.I, pad)
        stry = fill_pauli(pad, op_ind, Paulis.I, Paulis.Z)
        return PauliSum([strx, stry], complex.([1/2, 1/2]), already_sorted)
    else
        raise(DomainError(op))
    end
end


function fill_fermi(pad, op_ind, fill_op, end_op)
    str = Vector{FermiOps.FermiOp}(undef, pad)
    if pad <  op_ind
        throw(DimensionMismatch("pad is less than op_ind."))
    end
    @inbounds for i in 1:op_ind-1
        str[i] = fill_op
    end
    str[op_ind] = end_op
    @inbounds for i in op_ind+1:pad
        str[i] = FermiOps.id_op
    end
    return str
end

## Jordan-Wigner
function jordan_wigner_fermi(op::FermiOp, op_ind::Integer, pad::Integer)
    already_sorted = true
    if op === lower_op
        str = fill_fermi(pad, op_ind, FermiOps.Z_op, FermiOps.lower_op)
        return FermiTerm(str, complex(1.0))
    elseif op === raise_op
        str = fill_fermi(pad, op_ind, FermiOps.Z_op, FermiOps.raise_op)
        return FermiTerm(str, complex(1.0))
    elseif op === number_op
        str = fill_fermi(pad, op_ind, FermiOps.id_op, FermiOps.number_op)
        return FermiTerm(str, complex(1.0))
    elseif op === empty_op
        str = fill_fermi(pad, op_ind, FermiOps.id_op, FermiOps.empty_op)
        return FermiTerm(str, complex(1.0))
    elseif op === id_op
        str = fill(FermiOps.id_op, pad)
        return FermiTerm(str, complex(1.0))
    else
        raise(DomainError(op))
    end
end
