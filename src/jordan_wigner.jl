module JordanWigner

# using DocStringExtensions

import LinearAlgebra
using ..AbstractOps
using ..AbstractPaulis: Iop, Xop, Yop, Zop
import ..Paulis, ..FermiTerm, ..FermiSum
using ..FermiOps
import ..PauliDefault
import ..OpSum, ..sort_and_sum_duplicates!

export jordan_wigner, jordan_wigner_fermi

## Jordan-Wigner helper
"""
    fill_pauli(pad, op_ind, fill_op::PauliT, end_op) where PauliT

Return a `Vector{PauliT}` of length `pad` with indices `1:op_ind-1` filled with `fill_op`
and index `op_ind` set to `end_op` and the remaining indices filled with `one(Pauli)`
"""
function fill_pauli(pad, op_ind, fill_op::PauliT, end_op) where PauliT
    str = Vector{PauliT}(undef, pad)
    if pad <  op_ind
        throw(DimensionMismatch("pad is less than op_ind."))
    end
    for i in 1:op_ind-1
       @inbounds str[i] = fill_op
    end
    str[op_ind] = end_op
    for i in op_ind+1:pad
        @inbounds str[i] = PauliT(Iop)
    end
    return str
end

# Future compiler improvements may make this optimization worthless.
const _c1 = [complex(-1/2), complex(-im/2)]
const _c2 = [complex(-1/2), complex(im/2)]
const _c3 = [complex(1/2), complex(-1/2)]
const _c4 = [complex(1/2), complex(1/2)]

## Jordan-Wigner
function jordan_wigner(op::FermiOp, op_ind::Integer, pad::Integer, ::Type{PauliT}=PauliDefault) where PauliT
    _fill_pauli(fill_op_ind, end_op_ind) = fill_pauli(pad, op_ind, PauliT(fill_op_ind), PauliT(end_op_ind))
    if op === Lower
        strx = _fill_pauli(Zop, Xop)
        stry = _fill_pauli(Zop, Yop)
        coeffs = _c1
    elseif op === Raise
        strx = _fill_pauli(Zop, Xop)
        stry = _fill_pauli(Zop, Yop)
        coeffs = _c2
    elseif op === NumberOp
        strx = fill(PauliT(Iop), pad)
        stry = _fill_pauli(Iop, Zop)
        coeffs = _c3
    elseif op === Empty
        strx = fill(PauliT(Iop), pad)
        stry = _fill_pauli(Iop, Zop)
        coeffs = _c4
    else
        raise(DomainError(op))
    end
    return OpSum{PauliT}([strx, stry], coeffs; already_sorted=true)
end

function jordan_wigner(term::FermiTerm, ::Type{PauliT}=PauliDefault) where PauliT
    pad = length(term)
    facs = []
    for (i, op) in enumerate(term)
        if op !== I #  op === Raise || op === Lower || op === NumberOp || op === Empty
            push!(facs, jordan_wigner(op, i, pad, PauliT))
        end
    end
    if isempty(facs)  # String is all I
        return(OpSum{PauliT}([fill(one(PauliT), length(term))], [complex(term.coeff)]))
    end
#    return term.coeff * reduce(*, facs)
    psum = length(facs) == 1 ? only(facs) : reduce(*, facs) # TODO: better performance
    LinearAlgebra.lmul!(psum, term.coeff) # inplace mult does not appear to be faster
    return psum
end

"""
    jordan_wigner(fsum::FermiSum, ::Type{PauliT}=PauliDefault) where PauliT

Using the Jordan-Wigner transform, convert the electronic Hamiltonian `fsum` to a qubit operator
of type `OpSum{PauliT}`.

    jordan_wigner(term::FermiTerm, ::Type{PauliT}=PauliDefault) where PauliT

Convert `term` to a qubit operator of type `OpSum{PauliT}`.
"""
function jordan_wigner(fsum::FermiSum, ::Type{PauliT}=PauliDefault) where PauliT
    psum = jordan_wigner(fsum[1], PauliT) # could use already sorted flag
    @inbounds for i in 2:length(fsum)
        append!(psum, jordan_wigner(fsum[i], PauliT))
    end
    return sort_and_sum_duplicates!(psum)
end

####
#### Jordan-Wigner using Fermi operators augmented by Z = N - E
####

## These are experimental, and so far do not seem very useful

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
        str[i] = FermiOps.I
    end
    return str
end

## Jordan-Wigner
function jordan_wigner_fermi(op::FermiOp, op_ind::Integer, pad::Integer)
    if op === Lower
        str = fill_fermi(pad, op_ind, FermiOps.Z, FermiOps.Lower)
        return FermiTerm(str, complex(1.0))
    elseif op === Raise
        str = fill_fermi(pad, op_ind, FermiOps.Z, FermiOps.Raise)
        return FermiTerm(str, complex(1.0))
    elseif op === NumberOp
        str = fill_fermi(pad, op_ind, FermiOps.I, FermiOps.NumberOp)
        return FermiTerm(str, complex(1.0))
    elseif op === Empty
        str = fill_fermi(pad, op_ind, FermiOps.I, FermiOps.Empty)
        return FermiTerm(str, complex(1.0))
    elseif op === I
        str = fill(FermiOps.I, pad)
        return FermiTerm(str, complex(1.0))
    else
        raise(DomainError(op))
    end
end

function jordan_wigner_fermi(term::FermiTerm)
    pad = length(term)
    facs = [] # FermiTerm{FermiOp, Vector{FermiOp}, ComplexF64}[]
    for (i, op) in enumerate(term)
        if op === Raise || op === Lower || op === NumberOp
            push!(facs, jordan_wigner_fermi(op, i, pad))
        end
    end
    if isempty(facs)
        return (1.0 + 0.0im) * term
    end
    return term.coeff * reduce(*, facs)  # TODO: performance
end

function jordan_wigner_fermi(fsum::FermiSum)
    terms = FermiTerm{FermiOp, Vector{FermiOp}, ComplexF64}[]
    for i in eachindex(fsum) # 1:length(fsum)
        push!(terms, jordan_wigner_fermi(fsum[i]))
    end
    nterms = [x for x in terms]  # nterms is never used
    return FermiSum(terms)
end

end # module JordanWigner
