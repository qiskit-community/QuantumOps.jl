####
#### OpTerm and OpSum backed by SparseVector
####

# using ..AbstractPaulis: PauliPhaseCounts

SparseArraysN.neutral(T::Type{<:AbstractOp}) = one(T)
SparseArraysN.neutral(x::AbstractOp) = one(x)
SparseArraysN.isneutral(x::AbstractOp) = isone(x)
SparseArraysN.neutrals(T::Type{<:AbstractOp}, args...) = ones(T, args...)
SparseArraysN.neutrals(T::Type{<:AbstractOp}, args::Integer...) = ones(T, args...)

const SparseOpTerm = OpTerm{T, V} where {T<:AbstractOp, V<:SparseArraysN.SparseVector{T}}

#function Base.:*(t1::OpTerm{T, V}, t2::OpTerm{T, V}) where {T<:AbstractOp, V<:SparseArraysN.SparseVector{T}}
function Base.:*(t1::SparseOpTerm, t2::SparseOpTerm)
    (t_out, phase) = mul(t1.ops, t2.ops)
    return OpTerm(t_out, phase * t1.coeff * t2.coeff)
end

"""
    mul(v1::SparseArraysN.SparseVector{T}, v2::SparseArraysN.SparseVector{T}) where {T <: AbstractOp}

Multiply `v1` and `v2`, factorwise, returning a tuple of a `SparseVector` and the phase accumulated
from each multiplication.
"""
function mul(v1::SparseArraysN.SparseVector{T}, v2::SparseArraysN.SparseVector{T}) where {T <: AbstractOp}
    vals = T[]
    inds = empty(v1.nzind)
    phase_data = AbstractOps.phase_data(T)
    i1 = 1  # index into terms in t1
    i2 = 1
    while i1 <= weight(v1) && i2 <= weight(v2)
        ind1 = v1.nzind[i1]
        ind2 = v2.nzind[i2]
        val1 = v1.nzval[i1]
        val2 = v2.nzval[i2]
        ## ind1 corresponds to non-structural identity in v2
        ## Multiplying v1[ind1] by identity gives v1[ind1], so we just push
        ## this value to the result.
        if ind1 < ind2
            push!(inds, ind1)
            push!(vals, val1)
            i1 += 1
        ## Symmetric to the case above.
        elseif ind1 > ind2
            push!(inds, ind2)
            push!(vals, val2)
            i2 += 1
        else  # tt1.op and tt2.op operate on the same index (degree of freedom)
            ## Only here do we multiply non-trivial factors.
            new_op = val1 * val2  # product up to a phase
            if iszero(new_op)  # if any factor vanishes, the term vanishes.
                _phase = 0
                return (SparseArraysN.SparseVector(0, empty!(inds), empty!(vals)), _phase)
            end
            i1 += 1
            i2 += 1
            if isone(new_op)  # Identity is not stored in sparse representation
                continue
            end
#            phase_data = AbstractOps.accumulate_phase(phase_data, val1, val2)
            AbstractOps.accumulate_phase!(phase_data, val1, val2)
            push!(inds, ind1)
            push!(vals, new_op)
        end
    end
    ## Include remaining factors from the longer string.
    if i1 <= weight(v1)
        for i in i1:weight(v1)
            push!(inds, v1.nzind[i])
            push!(vals, v1.nzval[i])
        end
    elseif i2 <= weight(v2)
        for i in i2:weight(v2)
            push!(inds, v2.nzind[i])
            push!(vals, v2.nzval[i])
        end
    end
    n = isempty(inds) ? 0 : inds[end]
    _phase = AbstractOps.compute_phase(phase_data, v1, v2)
    return (SparseArraysN.SparseVector(n, inds, vals), _phase)
end

function _show_sparse_term(io, term)
    ops = term.ops
    xnnz = length(ops.nzind)
    for i in eachindex(ops.nzind)
        _show_op_plain(io, ops.nzval[i])
        print(io, ops.nzind[i])
        if i < length(ops.nzind)
            print(io, " ")
        end
    end
    if xnnz != 0
        print(io, " * ")
    end
    if term.coeff isa Real
        print(io, term.coeff)
    else
        print(io, "(", term.coeff, ")")
    end
end

function Base.show(io::IO, term::OpTerm{T, <:SparseArraysN.SparseVector}) where {T}
    ops = term.ops
    xnnz = length(ops.nzind)
    print(io, length(term), "-element ", typeof(term), " with ", xnnz,
          " stored ", xnnz == 1 ? "entry" : "entries")
    print(io, ":\n")
    _show_sparse_term(io, term)
end

function Base.show(io::IO, opsum::OpSum{T, V}) where {T<:AbstractOp, V <: Vector{<:SparseArraysN.SparseVector}}
    (m, n) = size(opsum)
    print(io, m, "x", n, " ", typeof(opsum), ":\n")
    for i in eachindex(opsum)
        _show_sparse_term(io, opsum[i])
        if i != lastindex(opsum)
            print(io, "\n")
        end
    end
end

## TODO: should we remove this? I think the remaining method would then return a copy
sparse_op(term::OpTerm{T, V}) where {T<:AbstractOp, V<:SparseArraysN.SparseVector{T}} =
    term

"""
    sparse_op(x::OpTerm)
    sparse_op(x::OpSum)

Convert `x` to a sparse representation. Note that the function `sparse` instead may convert `x` to
a sparse matrix.
"""
function sparse_op(term::OpTerm)
    return OpTerm(SparseArraysN.sparse(term.ops), term.coeff)
end

"""
    dense_op(term::OpTerm)
    dense_op(_sum::OpSum)

Convert `term` or `_sum` to a dense representation.
"""
function dense_op(term::OpTerm)
    return OpTerm(Vector(term.ops), term.coeff)
end

function _convert_op(_sum::OpSum, convert_func)
    sum_out = OpSum([convert_func(first(_sum))]; already_sorted=true)
    @inbounds for i in 2:length(_sum)
        push!(sum_out, convert_func(_sum[i]))
    end
    return sum_out
end

dense_op(_sum::OpSum) = _convert_op(_sum, dense_op)
sparse_op(_sum::OpSum) = _convert_op(_sum, sparse_op)
